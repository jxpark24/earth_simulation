"""FastAPI service: propagates the satellite catalogue and streams it to the browser."""
from __future__ import annotations

import asyncio
import logging
import math
import os
import struct
import time
from contextlib import asynccontextmanager

import numpy as np
from fastapi import FastAPI, WebSocket, WebSocketDisconnect
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

import sat_engine          # the C++ extension
import tle

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)-7s %(name)s: %(message)s")
log = logging.getLogger("sim")

STREAM_HZ = float(os.environ.get("STREAM_HZ", 20))
MAX_SATS  = int(os.environ.get("MAX_SATS", 0))      # 0 == no limit
GROUPS    = os.environ.get("TLE_GROUPS", "").split(",") if os.environ.get("TLE_GROUPS") else None

STATE: dict = {"prop": None, "catalog": None, "meta": {}}


def unix_to_jd(unix_seconds: float) -> float:
    """Unix epoch seconds -> Julian date (UTC, treated as UT1)."""
    return unix_seconds / 86400.0 + 2440587.5


def sun_vector(jd: float) -> tuple[float, float, float]:
    """Low-precision unit vector to the Sun in the Earth-centred equatorial frame.

    Astronomical Almanac approximation, good to ~0.01 deg -- ample for placing
    the day/night terminator.
    """
    n = jd - 2451545.0
    L = math.radians((280.460 + 0.9856474 * n) % 360.0)      # mean longitude
    g = math.radians((357.528 + 0.9856003 * n) % 360.0)      # mean anomaly
    lam = L + math.radians(1.915) * math.sin(g) + math.radians(0.020) * math.sin(2 * g)
    eps = math.radians(23.439 - 4.0e-7 * n)                  # obliquity
    return (math.cos(lam),
            math.cos(eps) * math.sin(lam),
            math.sin(eps) * math.sin(lam))


def build_propagator() -> None:
    cat = tle.load(GROUPS)
    if MAX_SATS and len(cat.names) > MAX_SATS:
        log.info("limiting catalog %d -> %d", len(cat.names), MAX_SATS)
        cat.names, cat.line1, cat.line2 = (cat.names[:MAX_SATS],
                                           cat.line1[:MAX_SATS],
                                           cat.line2[:MAX_SATS])
    prop = sat_engine.Propagator()
    t0 = time.perf_counter()
    prop.load(cat.names, cat.line1, cat.line2)
    load_ms = (time.perf_counter() - t0) * 1000

    # Time one full batch so the startup log states real throughput.
    jd = unix_to_jd(time.time())
    t0 = time.perf_counter()
    prop.propagate_all(jd)
    prop_ms = (time.perf_counter() - t0) * 1000

    STATE["prop"] = prop
    STATE["catalog"] = prop.catalog()
    STATE["meta"] = {
        "count": prop.size(),
        "valid": prop.valid(),
        "source": cat.source,
        "openmp": bool(sat_engine.openmp),
        "load_ms": round(load_ms, 1),
        "propagate_ms": round(prop_ms, 3),
        "earth_radius_km": sat_engine.EARTH_RADIUS_KM,
    }
    log.info("loaded %d satellites (%d valid) in %.0f ms; "
             "full propagation %.2f ms/frame (openmp=%s)",
             prop.size(), prop.valid(), load_ms, prop_ms, bool(sat_engine.openmp))


@asynccontextmanager
async def lifespan(app: FastAPI):
    await asyncio.to_thread(build_propagator)
    yield


app = FastAPI(title="Satellite Simulation Engine", lifespan=lifespan)
app.add_middleware(CORSMiddleware, allow_origins=["*"],
                   allow_methods=["*"], allow_headers=["*"])


@app.get("/api/health")
def health():
    return {"status": "ok" if STATE["prop"] else "loading", **STATE["meta"]}


@app.get("/api/catalog")
def catalog():
    if not STATE["prop"]:
        return JSONResponse({"error": "still loading"}, status_code=503)
    return {"meta": STATE["meta"], "satellites": STATE["catalog"]}


@app.get("/api/orbit/{index}")
def orbit(index: int, steps: int = 180):
    """One full orbital period of a satellite, for drawing its path."""
    prop = STATE["prop"]
    if not prop:
        return JSONResponse({"error": "still loading"}, status_code=503)
    jd = unix_to_jd(time.time())
    track = prop.orbit_track(index, jd, steps)
    pts = [None if math.isnan(float(p[0])) else [float(p[0]), float(p[1]), float(p[2])]
           for p in track]
    return {"index": index, "state": prop.state_of(index, jd),
            "track": [p for p in pts if p]}


@app.websocket("/ws")
async def ws(sock: WebSocket):
    await sock.accept()
    prop = STATE["prop"]
    if not prop:
        await sock.close(code=1013, reason="engine still loading")
        return

    # Simulation clock: sim_time advances at `rate` x wall-clock.
    rate = 1.0
    paused = False
    sim_unix = time.time()
    last_wall = time.perf_counter()
    interval = 1.0 / STREAM_HZ

    await sock.send_json({"type": "hello", "meta": STATE["meta"], "hz": STREAM_HZ})

    async def read_commands():
        """Apply client control messages; runs until the socket closes."""
        nonlocal rate, paused, sim_unix
        try:
            while True:
                msg = await sock.receive_json()
                cmd = msg.get("cmd")
                if cmd == "rate":
                    rate = max(-100000.0, min(100000.0, float(msg.get("value", 1))))
                elif cmd == "pause":
                    paused = True
                elif cmd == "resume":
                    paused = False
                elif cmd == "now":
                    sim_unix = time.time()
                elif cmd == "seek":
                    sim_unix = float(msg.get("unix", time.time()))
        except Exception:
            return

    reader = asyncio.create_task(read_commands())
    frames = 0
    try:
        while True:
            now_wall = time.perf_counter()
            dt = now_wall - last_wall
            last_wall = now_wall
            if not paused:
                sim_unix += dt * rate

            jd = unix_to_jd(sim_unix)
            # Release the GIL inside C++; run off the event loop thread so the
            # socket stays responsive with large catalogues.
            pos = await asyncio.to_thread(prop.propagate_all, jd)
            sx, sy, sz = sun_vector(jd)

            # Header: jd f64 | gmst f32 | sun xyz f32 | count u32  (28 bytes)
            header = struct.pack("<dffffI", jd, sat_engine.gmst(jd),
                                 sx, sy, sz, pos.shape[0])
            await sock.send_bytes(header + np.ascontiguousarray(pos, dtype=np.float32).tobytes())

            frames += 1
            await asyncio.sleep(max(0.0, interval - (time.perf_counter() - now_wall)))
    except (WebSocketDisconnect, RuntimeError, ConnectionError):
        pass
    finally:
        reader.cancel()
        log.info("client disconnected after %d frames", frames)
