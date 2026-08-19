"""Fetching, caching and parsing of NORAD two-line element sets."""
from __future__ import annotations

import logging
import os
import time
from dataclasses import dataclass
from pathlib import Path

import httpx

log = logging.getLogger("tle")

CELESTRAK = "https://celestrak.org/NORAD/elements/gp.php?GROUP={group}&FORMAT=tle"

# Groups pulled by default. Chosen to give a visually varied sky: crewed
# stations, a big LEO constellation, MEO navigation and GEO weather birds.
DEFAULT_GROUPS = ["stations", "starlink", "gps-ops", "galileo", "weather", "science"]

DATA_DIR = Path(os.environ.get("TLE_DATA_DIR", Path(__file__).parent.parent / "data"))
CACHE_TTL = int(os.environ.get("TLE_CACHE_TTL", 6 * 3600))   # refresh every 6h


@dataclass
class Catalog:
    names: list[str]
    line1: list[str]
    line2: list[str]
    source: str          # "network", "cache" or "bundled"
    fetched_at: float

    def __len__(self) -> int:
        return len(self.names)


def parse_tle_text(text: str) -> tuple[list[str], list[str], list[str]]:
    """Parse 3-line (name + 2 element lines) TLE text, skipping malformed records."""
    names: list[str] = []
    l1: list[str] = []
    l2: list[str] = []
    lines = [ln.rstrip() for ln in text.splitlines() if ln.strip()]
    i = 0
    while i < len(lines) - 1:
        # A record is either "name / 1 ... / 2 ..." or a bare "1 ... / 2 ..." pair.
        if lines[i].startswith("1 ") and lines[i + 1].startswith("2 "):
            names.append(f"SAT {lines[i][2:7].strip()}")
            l1.append(lines[i])
            l2.append(lines[i + 1])
            i += 2
        elif (i + 2 < len(lines)
              and lines[i + 1].startswith("1 ") and lines[i + 2].startswith("2 ")):
            names.append(lines[i].strip())
            l1.append(lines[i + 1])
            l2.append(lines[i + 2])
            i += 3
        else:
            i += 1
    return names, l1, l2


def _cache_path(group: str) -> Path:
    return DATA_DIR / f"{group}.tle"


def _fetch_group(group: str, client: httpx.Client) -> str | None:
    url = CELESTRAK.format(group=group)
    try:
        r = client.get(url, timeout=25.0)
        r.raise_for_status()
        body = r.text
        # Celestrak answers rate-limit / unknown-group with a plain text notice.
        if "1 " not in body[:400] or len(body) < 100:
            log.warning("group %s returned no elements: %s", group, body[:120].strip())
            return None
        return body
    except Exception as exc:                       # network down, DNS, timeout...
        log.warning("fetch failed for %s: %s", group, exc)
        return None


def load(groups: list[str] | None = None, force_refresh: bool = False) -> Catalog:
    """Load TLEs, preferring fresh network data and falling back to disk cache."""
    groups = groups or DEFAULT_GROUPS
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    names: list[str] = []
    l1: list[str] = []
    l2: list[str] = []
    used_network = False
    used_cache = False

    with httpx.Client(follow_redirects=True,
                      headers={"User-Agent": "satellite-sim/1.0"}) as client:
        for group in groups:
            path = _cache_path(group)
            fresh = (path.exists()
                     and (time.time() - path.stat().st_mtime) < CACHE_TTL
                     and not force_refresh)

            text: str | None = None
            if fresh:
                text = path.read_text()
                used_cache = True
            else:
                text = _fetch_group(group, client)
                if text:
                    path.write_text(text)
                    used_network = True
                elif path.exists():
                    log.info("using stale cache for %s", group)
                    text = path.read_text()
                    used_cache = True

            if not text:
                log.error("no data at all for group %s", group)
                continue

            gn, g1, g2 = parse_tle_text(text)
            log.info("group %-10s -> %5d satellites", group, len(gn))
            names += gn
            l1 += g1
            l2 += g2

    # Drop duplicates (a satellite can appear in more than one group).
    seen: set[str] = set()
    dn: list[str] = []
    d1: list[str] = []
    d2: list[str] = []
    for n, a, b in zip(names, l1, l2):
        key = a[2:7]                                # NORAD catalogue number
        if key in seen:
            continue
        seen.add(key)
        dn.append(n)
        d1.append(a)
        d2.append(b)

    source = "network" if used_network else ("cache" if used_cache else "bundled")
    log.info("catalog: %d unique satellites (source=%s)", len(dn), source)
    return Catalog(dn, d1, d2, source, time.time())
