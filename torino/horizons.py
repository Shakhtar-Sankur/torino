"""State vectors from JPL Horizons, the reference that torino's propagator is checked against."""

from __future__ import annotations

import json
import urllib.parse

import numpy as np

from torino.sources import _get

HORIZONS_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"


def horizons_command(designation: str) -> str:
    """Horizons COMMAND for a small body: '99942;' for a number, 'DES=2008 JL3;' for a provisional designation."""
    designation = designation.strip()
    return f"{designation};" if designation.isdigit() else f"DES={designation};"


def state_vectors(designation: str, jd_tdb: list[float], center: str = "500@0") -> np.ndarray:
    """Barycentric (by default) ICRF state vectors [x, y, z, vx, vy, vz] in au and au/day at the given TDB Julian dates."""
    params = {
        "format": "json",
        "COMMAND": f"'{horizons_command(designation)}'",
        "EPHEM_TYPE": "VECTORS",
        "CENTER": f"'{center}'",
        "TLIST": " ".join(f"'{t:.9f}'" for t in jd_tdb),
        "TLIST_TYPE": "JD",
        "TIME_TYPE": "TDB",
        "VEC_TABLE": "2",
        "OUT_UNITS": "AU-D",
        "REF_PLANE": "FRAME",
        "REF_SYSTEM": "ICRF",
        "VEC_CORR": "NONE",
        "CSV_FORMAT": "YES",
        "OBJ_DATA": "NO",
    }
    result = json.loads(_get(f"{HORIZONS_URL}?{urllib.parse.urlencode(params)}"))["result"]
    return parse_vectors(result)


def parse_vectors(result: str) -> np.ndarray:
    """Parse the CSV block between $$SOE and $$EOE: JD, calendar date, x, y, z, vx, vy, vz."""
    start, end = result.find("$$SOE"), result.find("$$EOE")
    if start < 0 or end < 0:
        raise ValueError(f"no state vectors in the Horizons reply:\n{result[:1000]}")
    rows = []
    for line in result[start + 5:end].strip().splitlines():
        fields = [f.strip() for f in line.split(",")]
        rows.append([float(f) for f in fields[2:8]])
    return np.array(rows)

