// 3D satellite visualisation.
//
// Frames: the engine streams TEME positions (Z = celestial north, X = vernal
// equinox). Three.js is Y-up, so every vector is remapped (x, y, z) -> (x, z, -y),
// a pure -90 deg rotation about X that preserves handedness. The Earth mesh is
// then spun by GMST so its texture lines up with the inertial satellite frame.

import * as THREE from '/vendor/three.module.js';
import { OrbitControls } from '/vendor/OrbitControls.js';

const EARTH_RADIUS_KM = 6378.135;      // WGS-72, matches the engine
const KM = 1 / EARTH_RADIUS_KM;        // scale: 1 scene unit == 1 Earth radius

// ---------------------------------------------------------------- scene setup
const canvas   = document.getElementById('scene');
const renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: false });
renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
renderer.setSize(innerWidth, innerHeight);
renderer.setClearColor(0x000000, 1);          // black space

const scene  = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(42, innerWidth / innerHeight, 0.01, 4000);
camera.position.set(0, 1.6, 4.2);

const controls = new OrbitControls(camera, canvas);
controls.enableDamping = true;
controls.dampingFactor = 0.06;
controls.rotateSpeed    = 0.45;
controls.zoomSpeed      = 0.8;
controls.minDistance    = 1.08;           // just above the surface
controls.maxDistance    = 120;
controls.enablePan      = false;

// ------------------------------------------------------------------- lighting
// The directional light *is* the Sun; its direction arrives with every frame,
// so the day/night terminator is physically placed.
const sunLight = new THREE.DirectionalLight(0xfff5e8, 2.6);
sunLight.position.set(5, 2, 5);
scene.add(sunLight);
scene.add(new THREE.AmbientLight(0x2a3a52, 0.55));   // keeps the night side readable

// ---------------------------------------------------------------------- Earth
const texLoader = new THREE.TextureLoader();
const loadTex = (url, srgb = false) => texLoader.load(
  url,
  (t) => { if (srgb) t.colorSpace = THREE.SRGBColorSpace; },
  undefined,
  () => console.warn('texture failed:', url)
);

const earthGeo = new THREE.SphereGeometry(1, 96, 64);
const earthMat = new THREE.MeshPhongMaterial({
  map:         loadTex('/assets/earth_atmos_2048.jpg', true),
  specularMap: loadTex('/assets/earth_specular_2048.jpg'),
  normalMap:   loadTex('/assets/earth_normal_2048.jpg'),
  normalScale: new THREE.Vector2(0.85, 0.85),
  specular:    new THREE.Color(0x2a4766),
  shininess:   18,
});
const earth = new THREE.Mesh(earthGeo, earthMat);
scene.add(earth);

// Clouds: a slightly larger translucent shell, drifting a touch faster than the
// surface so the planet doesn't look rigid.
const clouds = new THREE.Mesh(
  new THREE.SphereGeometry(1.006, 96, 64),
  new THREE.MeshLambertMaterial({
    map: loadTex('/assets/earth_clouds_1024.png', true),
    transparent: true, opacity: 0.42, depthWrite: false,
  })
);
scene.add(clouds);

// Atmosphere: rim glow rendered on the inside of a slightly larger sphere, so
// it only shows where the surface curves away from the camera.
const atmosphere = new THREE.Mesh(
  new THREE.SphereGeometry(1.055, 96, 64),
  new THREE.ShaderMaterial({
    uniforms: { uSun: { value: new THREE.Vector3(1, 0, 0) } },
    vertexShader: `
      varying vec3 vNormal; varying vec3 vWorld;
      void main() {
        vNormal = normalize(normalMatrix * normal);
        vWorld  = normalize(mat3(modelMatrix) * normal);
        gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
      }`,
    fragmentShader: `
      uniform vec3 uSun;
      varying vec3 vNormal; varying vec3 vWorld;
      void main() {
        // Strongest exactly at the limb.
        float rim = pow(1.0 - abs(dot(vNormal, vec3(0.0, 0.0, 1.0))), 3.2);
        // ...and only on the lit hemisphere, fading through the terminator.
        float lit = smoothstep(-0.35, 0.45, dot(vWorld, uSun));
        vec3 col = mix(vec3(0.16, 0.34, 0.72), vec3(0.42, 0.68, 1.0), lit);
        gl_FragColor = vec4(col, rim * (0.18 + 0.82 * lit));
      }`,
    blending: THREE.AdditiveBlending,
    side: THREE.BackSide,
    transparent: true,
    depthWrite: false,
  })
);
scene.add(atmosphere);

// ------------------------------------------------------------------ starfield
// Points on a large sphere, well outside any orbit we draw.
(function addStars() {
  const N = 6000, pos = new Float32Array(N * 3), col = new Float32Array(N * 3);
  for (let i = 0; i < N; i++) {
    // Uniform on the sphere: cos(theta) must be sampled uniformly, not theta.
    const u = Math.random() * 2 - 1, a = Math.random() * Math.PI * 2;
    const s = Math.sqrt(1 - u * u), r = 900;
    pos[i * 3] = r * s * Math.cos(a);
    pos[i * 3 + 1] = r * u;
    pos[i * 3 + 2] = r * s * Math.sin(a);
    const b = 0.55 + Math.random() * 0.45;
    const warm = Math.random() < 0.15;
    col[i * 3] = b; col[i * 3 + 1] = b * (warm ? 0.88 : 0.97); col[i * 3 + 2] = b * (warm ? 0.76 : 1.0);
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.BufferAttribute(pos, 3));
  g.setAttribute('color', new THREE.BufferAttribute(col, 3));
  scene.add(new THREE.Points(g, new THREE.PointsMaterial({
    size: 1.6, sizeAttenuation: false, vertexColors: true,
    transparent: true, opacity: 0.85, depthWrite: false,
  })));
})();

// ------------------------------------------------------------------ satellites
// One Points object for the whole catalogue: positions are rewritten in place
// from each stream frame, so N satellites cost one draw call.
const satGeo = new THREE.BufferGeometry();
const satMat = new THREE.ShaderMaterial({
  uniforms: { uScale: { value: 3.5 } },
  vertexShader: `
    attribute float aSize;
    attribute vec3  aColor;
    uniform   float uScale;
    varying   vec3  vColor;
    void main() {
      vColor = aColor;
      vec4 mv = modelViewMatrix * vec4(position, 1.0);
      // aSize is 0 for invalid objects, which zeroes the point out entirely.
      // Clamped so a close zoom can't turn 11k additive sprites into a whiteout.
      gl_PointSize = clamp(aSize * uScale * (3.0 / max(-mv.z, 0.05)), 0.0, 16.0);
      gl_Position  = projectionMatrix * mv;
    }`,
  fragmentShader: `
    varying vec3 vColor;
    void main() {
      vec2 c = gl_PointCoord - 0.5;
      float d = length(c);
      if (d > 0.5) discard;
      float a = smoothstep(0.5, 0.08, d);      // soft round dot
      gl_FragColor = vec4(vColor, a);
    }`,
  transparent: true,
  depthWrite: false,          // additive dots must not occlude each other
  depthTest: true,            // ...but the Earth must still occlude them
  blending: THREE.AdditiveBlending,
});
const satPoints = new THREE.Points(satGeo, satMat);
satPoints.frustumCulled = false;
scene.add(satPoints);

// Orbit path of the selected satellite.
const orbitLine = new THREE.Line(
  new THREE.BufferGeometry(),
  new THREE.LineBasicMaterial({ color: 0x4da3ff, transparent: true, opacity: 0.7 })
);
orbitLine.visible = false;
orbitLine.frustumCulled = false;
scene.add(orbitLine);

// Highlight ring around the selected satellite.
const marker = new THREE.Mesh(
  new THREE.RingGeometry(0.028, 0.038, 32),
  new THREE.MeshBasicMaterial({ color: 0x4da3ff, side: THREE.DoubleSide,
                                transparent: true, opacity: 0.95, depthTest: false })
);
marker.visible = false;
marker.renderOrder = 10;
scene.add(marker);

// --------------------------------------------------------------- catalog state
let catalog   = [];
let count     = 0;
let positions = null;      // Float32Array(count*3), scene coords
let sizes     = null;
let selected  = -1;

// Colour satellites by orbital regime -- the visual grouping people expect.
function regimeColor(periodMin) {
  if (periodMin <= 0)    return [0.55, 0.60, 0.70];
  if (periodMin < 128)   return [0.42, 0.85, 1.00];   // LEO  - cyan
  if (periodMin < 800)   return [1.00, 0.82, 0.35];   // MEO  - amber
  if (periodMin < 1600)  return [1.00, 0.52, 0.32];   // GEO  - orange
  return [0.78, 0.60, 1.00];                          // HEO / graveyard
}

function buildCatalog(sats) {
  catalog = sats;
  count = sats.length;
  positions = new Float32Array(count * 3);
  sizes = new Float32Array(count);
  const colors = new Float32Array(count * 3);
  for (let i = 0; i < count; i++) {
    const c = regimeColor(sats[i].period);
    colors[i * 3] = c[0]; colors[i * 3 + 1] = c[1]; colors[i * 3 + 2] = c[2];
    sizes[i] = sats[i].valid ? 1 : 0;
    positions[i * 3] = positions[i * 3 + 1] = positions[i * 3 + 2] = 0;
  }
  satGeo.setAttribute('position', new THREE.BufferAttribute(positions, 3));
  satGeo.setAttribute('aColor',   new THREE.BufferAttribute(colors, 3));
  satGeo.setAttribute('aSize',    new THREE.BufferAttribute(sizes, 1));
  document.getElementById('h-count').textContent = count.toLocaleString();
}

// ------------------------------------------------------------------- websocket
const WS_URL = `${location.protocol === 'https:' ? 'wss' : 'ws'}://${location.host}/ws`;
let socket = null, frames = 0, lastRateCheck = performance.now(), streamHz = 0;
let simDate = new Date();

function setLoader(msg, done = false) {
  const el = document.getElementById('loader');
  document.getElementById('loadmsg').textContent = msg;
  if (done) { el.classList.add('hidden'); setTimeout(() => (el.style.display = 'none'), 900); }
}

async function boot() {
  setLoader('fetching satellite catalog…');
  try {
    const r = await fetch('/api/catalog');
    if (!r.ok) throw new Error(`catalog HTTP ${r.status}`);
    const data = await r.json();
    buildCatalog(data.satellites);
    const m = data.meta;
    document.getElementById('h-engine').textContent =
      `C++ ${m.propagate_ms}ms${m.openmp ? ' · omp' : ''}`;
    document.querySelector('#title .sub').textContent =
      `${m.valid.toLocaleString()} objects · elements: ${m.source}`;
  } catch (e) {
    setLoader('engine unavailable — is the backend running?');
    console.error(e);
    return;
  }
  connect();
}

function connect() {
  setLoader('opening position stream…');
  socket = new WebSocket(WS_URL);
  socket.binaryType = 'arraybuffer';

  socket.onmessage = (ev) => {
    if (typeof ev.data === 'string') return;         // the JSON "hello"
    const buf = ev.data;
    const dv = new DataView(buf);
    const jd    = dv.getFloat64(0, true);
    const gmst  = dv.getFloat32(8, true);
    const sun   = [dv.getFloat32(12, true), dv.getFloat32(16, true), dv.getFloat32(20, true)];
    const n     = dv.getUint32(24, true);
    const src   = new Float32Array(buf, 28, n * 3);

    if (!positions || n !== count) return;           // catalog/stream mismatch

    // TEME (km) -> scene units, Z-up to Y-up.
    for (let i = 0; i < n; i++) {
      const x = src[i * 3], y = src[i * 3 + 1], z = src[i * 3 + 2];
      if (Number.isNaN(x)) { sizes[i] = 0; continue; }   // decayed / diverged
      positions[i * 3]     =  x * KM;
      positions[i * 3 + 1] =  z * KM;
      positions[i * 3 + 2] = -y * KM;
    }
    satGeo.attributes.position.needsUpdate = true;
    satGeo.attributes.aSize.needsUpdate = true;
    // Points.raycast() tests the ray against a cached bounding sphere. It is
    // computed lazily on first use, so without this a click landing before the
    // first frame would cache a zero-radius sphere and break picking for good;
    // it would also go stale as the satellites move. Recomputed on next pick.
    satGeo.boundingSphere = null;

    // Earth spins with sidereal time; the texture then matches the sky.
    earth.rotation.y  = gmst;
    clouds.rotation.y = gmst * 1.0008;

    // Sun direction, same frame remap.
    const sv = new THREE.Vector3(sun[0], sun[2], -sun[1]).normalize();
    sunLight.position.copy(sv).multiplyScalar(50);
    atmosphere.material.uniforms.uSun.value.copy(sv);

    // Julian date -> UTC calendar time for the HUD.
    simDate = new Date((jd - 2440587.5) * 86400000);
    frames++;

    if (loading) {
      loading = false;
      // Open the view on the daylit hemisphere: looking down the Sun vector
      // means the visible face is fully lit, swung round and raised slightly
      // so the shot reads as 3D rather than flat-on.
      const dist = camera.position.length();
      camera.position.copy(sv)
        .applyAxisAngle(new THREE.Vector3(0, 1, 0), -0.6)
        .setLength(dist);
      camera.position.y += dist * 0.28;
      camera.position.setLength(dist);
      controls.update();
      setLoader('', true);
    }
  };

  socket.onopen  = () => setLoader('receiving…');
  socket.onerror = () => setLoader('stream error');
  socket.onclose = () => {
    setLoader('stream closed — reconnecting…');
    setTimeout(connect, 2500);
  };
}
let loading = true;

// --------------------------------------------------------------------- picking
const raycaster = new THREE.Raycaster();
raycaster.params.Points.threshold = 0.02;
const pointer = new THREE.Vector2();
let dragged = false;

canvas.addEventListener('pointerdown', () => { dragged = false; });
canvas.addEventListener('pointermove', (e) => { if (e.buttons) dragged = true; });
canvas.addEventListener('pointerup', (e) => {
  if (dragged || !count) return;                 // an orbit drag, not a click
  pointer.x = (e.clientX / innerWidth) * 2 - 1;
  pointer.y = -(e.clientY / innerHeight) * 2 + 1;
  raycaster.setFromCamera(pointer, camera);
  // Scale the pick radius with zoom so distant dots stay clickable.
  raycaster.params.Points.threshold = 0.008 * camera.position.length();
  const hits = raycaster.intersectObject(satPoints);
  const hit = hits.find((h) => sizes[h.index] > 0);
  if (hit) select(hit.index); else deselect();
});

async function select(index) {
  selected = index;
  const s = catalog[index];
  const panel = document.getElementById('info');
  panel.style.display = 'block';
  document.getElementById('i-name').textContent = s.name;
  document.getElementById('i-id').textContent = s.id;
  document.getElementById('i-inc').textContent = `${s.incl.toFixed(2)}°`;
  document.getElementById('i-per').textContent = `${s.period.toFixed(1)} min`;
  marker.visible = true;
  if (showOrbits) loadOrbit(index);
  try {
    const r = await fetch(`/api/orbit/${index}?steps=2`);
    const d = await r.json();
    document.getElementById('i-alt').textContent = `${d.state.altitude_km.toFixed(0)} km`;
    document.getElementById('i-spd').textContent = `${d.state.speed_kms.toFixed(3)} km/s`;
  } catch { /* readout stays blank */ }
}

function deselect() {
  selected = -1;
  document.getElementById('info').style.display = 'none';
  marker.visible = false;
  orbitLine.visible = false;
}

async function loadOrbit(index) {
  try {
    const r = await fetch(`/api/orbit/${index}?steps=256`);
    const d = await r.json();
    const pts = new Float32Array(d.track.length * 3);
    d.track.forEach((p, i) => {
      pts[i * 3] = p[0] * KM; pts[i * 3 + 1] = p[2] * KM; pts[i * 3 + 2] = -p[1] * KM;
    });
    orbitLine.geometry.dispose();
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(pts, 3));
    orbitLine.geometry = g;
    orbitLine.visible = true;
  } catch { /* leave the previous path in place */ }
}

// -------------------------------------------------------------------- controls
// Rate slider is logarithmic: even steps map to human-meaningful speeds.
const RATES = [1, 5, 10, 30, 60, 120, 300, 600, 1800, 3600, 10800];
const rateEl = document.getElementById('rate');
const send = (o) => socket?.readyState === 1 && socket.send(JSON.stringify(o));

function applyRate() {
  const v = RATES[+rateEl.value];
  document.getElementById('rateval').textContent =
    v >= 3600 ? `${v / 3600}h/s` : v >= 60 ? `${v / 60}m/s` : `${v}×`;
  send({ cmd: 'rate', value: v });
}
rateEl.addEventListener('input', applyRate);

let paused = false;
document.getElementById('pause').addEventListener('click', (e) => {
  paused = !paused;
  e.target.textContent = paused ? 'resume' : 'pause';
  e.target.classList.toggle('on', paused);
  send({ cmd: paused ? 'pause' : 'resume' });
});
document.getElementById('now').addEventListener('click', () => send({ cmd: 'now' }));

let showOrbits = false;
document.getElementById('toggleOrbits').addEventListener('click', (e) => {
  showOrbits = !showOrbits;
  e.target.textContent = showOrbits ? 'on' : 'off';
  e.target.classList.toggle('on', showOrbits);
  if (showOrbits && selected >= 0) loadOrbit(selected);
  else orbitLine.visible = false;
});
document.getElementById('toggleClouds').addEventListener('click', (e) => {
  clouds.visible = !clouds.visible;
  e.target.textContent = clouds.visible ? 'on' : 'off';
  e.target.classList.toggle('on', clouds.visible);
});
document.getElementById('size').addEventListener('input', (e) => {
  satMat.uniforms.uScale.value = +e.target.value;
});

addEventListener('resize', () => {
  camera.aspect = innerWidth / innerHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(innerWidth, innerHeight);
});

// ----------------------------------------------------------------- render loop
let fpsFrames = 0, fpsLast = performance.now();

function tick() {
  requestAnimationFrame(tick);
  controls.update();

  if (selected >= 0 && positions) {
    marker.position.set(positions[selected * 3],
                        positions[selected * 3 + 1],
                        positions[selected * 3 + 2]);
    marker.quaternion.copy(camera.quaternion);        // always face the viewer
    const d = marker.position.distanceTo(camera.position);
    marker.scale.setScalar(0.8 * d);                  // constant on-screen size
  }

  renderer.render(scene, camera);

  // HUD, refreshed twice a second.
  fpsFrames++;
  const now = performance.now();
  if (now - fpsLast > 500) {
    const fps = (fpsFrames * 1000) / (now - fpsLast);
    document.getElementById('h-fps').textContent = `${fps.toFixed(0)} fps`;
    streamHz = (frames * 1000) / (now - lastRateCheck);
    document.getElementById('h-stream').textContent = `${streamHz.toFixed(0)} Hz`;
    document.getElementById('h-time').textContent =
      simDate.toISOString().replace('T', ' ').slice(0, 19);
    fpsFrames = 0; frames = 0; fpsLast = now; lastRateCheck = now;
  }
}

// Fade the hint out once the user has had a moment to read it.
setTimeout(() => document.getElementById('hint')?.classList.add('hidden'), 9000);

applyRate();
boot();
tick();
