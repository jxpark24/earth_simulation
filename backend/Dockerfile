# ---------- stage 1: compile the C++ propagation engine ----------
FROM python:3.12-slim AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential cmake ca-certificates curl \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir pybind11==2.13.6 numpy==2.2.1

WORKDIR /build
COPY CMakeLists.txt .
COPY src/ src/

RUN cmake -B build -S . \
        -DCMAKE_BUILD_TYPE=Release \
        -Dpybind11_DIR="$(python -c 'import pybind11; print(pybind11.get_cmake_dir())')" \
    && cmake --build build -j "$(nproc)" \
    && ls -la build/*.so

# Bake a TLE snapshot into the image so the simulation still has real elements
# if the container later runs without internet. Never fails the build.
RUN mkdir -p /build/data && cd /build/data && \
    for g in stations starlink gps-ops galileo weather science; do \
        curl -sS -m 45 -A "satellite-sim/1.0" \
            -o "$g.tle" "https://celestrak.org/NORAD/elements/gp.php?GROUP=$g&FORMAT=tle" \
            && echo "bundled $g: $(wc -l < "$g.tle") lines" \
            || echo "skipped $g (offline)"; \
    done; rm -f $(find . -size -100c) 2>/dev/null; true

# ---------- stage 2: runtime ----------
FROM python:3.12-slim

# libgomp1 is the OpenMP runtime the engine links against.
RUN apt-get update && apt-get install -y --no-install-recommends libgomp1 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY --from=builder /build/build/*.so /app/
COPY --from=builder /build/data/ /app/data/
COPY app/ /app/

ENV PYTHONUNBUFFERED=1 \
    TLE_DATA_DIR=/app/data \
    STREAM_HZ=20

EXPOSE 8000
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000", "--ws", "websockets"]
