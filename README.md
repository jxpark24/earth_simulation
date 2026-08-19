# Satellite Simulation

A real-time 3D visualisation of ~11,000 real satellites orbiting Earth, in the
browser. Orbital mechanics run in C++; rendering runs on the GPU.

![stack](https://img.shields.io/badge/engine-C%2B%2B%20SGP4%2FSDP4-blue) ![api](https://img.shields.io/badge/api-FastAPI%20%2B%20WebSocket-green) ![ui](https://img.shields.io/badge/ui-three.js-orange)

## Quick start

```bash
docker compose up --build
```

Then open **http://localhost:8080**.

First build takes a couple of minutes (it compiles the C++ extension and
vendors three.js + Earth textures into the images). Startup fetches current
orbital elements from Celestrak; a snapshot is baked into the image as a
fallback, so it also works with no internet.

## What it does

- **Earth in the middle, in 3D.** Drag to orbit the camera, scroll to zoom.
- **~11,000 real satellites**, positioned from live NORAD two-line elements —
  Starlink, the ISS, GPS, Galileo, weather and science satellites.
- Satellites are **coloured by orbital regime**: cyan LEO, amber MEO, orange
  GEO, violet highly-elliptical.
- **Click any satellite** to read its NORAD ID, altitude, speed, period and
  inclination, and to draw its orbit path.
- **Time control** from real-time up to 3 hours per second, pause, and reset.
- The Earth spins by true sidereal time and is lit by the **real Sun
  direction**, so the day/night terminator is physically placed.

## Architecture

```
browser ──HTTP──> nginx ──/api──> FastAPI ──> sat_engine (C++ .so)
   │                 │                            │
   │  three.js       └──/ws────> binary position stream
   └── WebGL: Earth, atmosphere, 11k-point cloud (1 draw call)
```

**Why this split:** propagating 11k satellites every frame is the expensive
part, so it lives in C++. Rendering must be WebGL to run in a browser at all.
The two talk over a compact binary WebSocket frame rather than JSON — 28 header
bytes plus 12 bytes per satellite, about 132 KB per frame at 20 Hz.

| Path | Contents |
|---|---|
| `backend/src/sgp4.{hpp,cpp}` | SGP4 near-Earth propagator, TLE parsing, Julian date, GMST |
| `backend/src/deepspace.{hpp,cpp}` | SDP4 deep-space: lunisolar periodics + geopotential resonance |
| `backend/src/engine.cpp` | pybind11 bindings; batched, OpenMP-parallel, GIL-released |
| `backend/app/main.py` | FastAPI, WebSocket stream, REST catalog/orbit endpoints |
| `backend/app/tle.py` | Celestrak fetch, disk cache, TLE parsing, de-duplication |
| `frontend/js/main.js` | three.js scene, stream decoding, picking, controls |

### Coordinate frames

The engine emits **TEME** positions (Z = celestial north, X = vernal equinox).
The frontend remaps `(x, y, z) → (x, z, −y)`, a −90° rotation about X, because
three.js is Y-up. The Earth mesh is then rotated by GMST, which is what makes
the texture line up with the inertial satellite frame.

## Accuracy

The propagator is validated against the reference Vallado/`python-sgp4`
implementation across the **entire 10,955-satellite catalog**:

| Population | Count | Max position error |
|---|---:|---:|
| Near-Earth (SGP4) | 10,858 | **0.45 m** |
| Deep-space (SDP4) | 97 | **7.3 m** |

Error-flag agreement with the reference is exact (0 disagreements). The residual
is float32 quantisation of the wire format, not model error — at 42,000 km one
float32 ulp is ≈2.5 m. Accuracy holds at ~3.7 m out to 7 days from epoch.

Spot checks against published values: ISS at 410.7 km / 7.668 km/s / 51.633°
inclination; METEOSAT-9 at 35,790 km with a 1436.1-minute period — exactly one
sidereal day, the defining property of a geostationary orbit.

## Performance

Full propagation of **10,955 satellites in ~0.5 ms** per frame (Apple Silicon,
OpenMP across all cores, mean of 200 warm iterations) — about 22 million
satellite-propagations per second. At the default 20 Hz stream that is roughly
1% of one core, so the cost is dominated by the network, not the physics.

The renderer draws all satellites as a single `THREE.Points` object with a
custom shader: one draw call regardless of catalog size.

## Configuration

Set in `docker-compose.yml`:

| Variable | Default | Meaning |
|---|---|---|
| `TLE_GROUPS` | `stations,starlink,gps-ops,galileo,weather,science` | Celestrak groups to load |
| `STREAM_HZ` | `20` | Position frames per second |
| `MAX_SATS` | `0` | Cap on satellite count (0 = unlimited) |
| `TLE_CACHE_TTL` | `21600` | Seconds before re-fetching elements |

To add more objects, append Celestrak group names — e.g. `oneweb`, `iridium-NEXT`,
`active` (~30,000 objects; the engine handles it, though bandwidth grows linearly).

## API

| Endpoint | Purpose |
|---|---|
| `GET /api/health` | Engine status, satellite count, timing |
| `GET /api/catalog` | Full catalog: names, NORAD IDs, inclination, period |
| `GET /api/orbit/{i}?steps=N` | One full orbital period, plus current state |
| `WS /ws` | Binary position stream |

WebSocket frame layout (little-endian): `float64 jd`, `float32 gmst`,
`float32 sun[3]`, `uint32 count`, then `float32 xyz[count]` in km.
Client commands are JSON: `{"cmd":"rate","value":60}`, `pause`, `resume`, `now`.

## Notes and limitations

- Deep-space objects use the full SDP4 model, but TLE accuracy itself degrades
  over days; refresh elements for best results (automatic every 6 hours).
- Satellite positions stream at 20 Hz and are not interpolated between frames.
  It is imperceptible at normal rates; raise `STREAM_HZ` if you drive time very
  fast and want smoother motion.
- Textures and three.js are vendored at image build time, so the running site
  makes no external requests.
