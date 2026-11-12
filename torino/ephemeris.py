"""Positions of the Sun, the planets, the Moon and Pluto from JPL's DE440s planetary ephemeris.

DE440s is the 1849-2150 subset of DE440 (Park et al. 2021, AJ 161:105). The GM values are DE440's.
"""

from __future__ import annotations

import os
import urllib.request

import numpy as np

DE440S_URL = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp"
AU_KM = 149_597_870.7
DAY_S = 86_400.0

# GM in km^3/s^2 from the DE440 header, with the NAIF codes torino reads them under.
GM_KM3_S2 = {
    "Sun": 132_712_440_041.279419,
    "Mercury": 22_031.868551,
    "Venus": 324_858.592000,
    "Earth": 398_600.435507,
    "Moon": 4_902.800118,
    "Mars": 42_828.375816,
    "Jupiter": 126_712_764.100000,
    "Saturn": 37_940_584.841800,
    "Uranus": 5_794_556.400000,
    "Neptune": 6_836_527.100580,
    "Pluto": 975.500000,
}

# Segments of the SPK file that give each body's position relative to the solar-system barycentre (code 0).
# The planets other than the Earth are their system barycentres; the Earth and the Moon chain through the
# Earth-Moon barycentre (code 3).
_CHAINS = {
    "Sun": [(0, 10)],
    "Mercury": [(0, 1)],
    "Venus": [(0, 2)],
    "Earth": [(0, 3), (3, 399)],
    "Moon": [(0, 3), (3, 301)],
    "Mars": [(0, 4)],
    "Jupiter": [(0, 5)],
    "Saturn": [(0, 6)],
    "Uranus": [(0, 7)],
    "Neptune": [(0, 8)],
    "Pluto": [(0, 9)],
}


def gm_au3_day2(name: str) -> float:
    return GM_KM3_S2[name] * DAY_S**2 / AU_KM**3


def default_path() -> str:
    return os.environ.get("TORINO_EPHEMERIS",
                          os.path.join(os.path.expanduser("~"), ".cache", "torino", "de440s.bsp"))


def ensure_de440s(path: str | None = None) -> str:
    """Download DE440s once and return its path."""
    path = path or default_path()
    if not os.path.exists(path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        tmp = path + ".part"
        urllib.request.urlretrieve(DE440S_URL, tmp)
        os.replace(tmp, path)
    return path


class Ephemeris:
    """Barycentric ICRF positions in au, at TDB Julian dates."""

    def __init__(self, path: str | None = None, bodies: list[str] | None = None):
        from jplephem.spk import SPK

        self.kernel = SPK.open(ensure_de440s(path))
        self.bodies = bodies or list(_CHAINS)
        self.gm = np.array([gm_au3_day2(b) for b in self.bodies])

    def positions(self, jd_tdb: float) -> np.ndarray:
        """Array of shape (len(bodies), 3)."""
        out = np.empty((len(self.bodies), 3))
        for i, body in enumerate(self.bodies):
            km = sum(self.kernel[a, b].compute(jd_tdb) for a, b in _CHAINS[body])
            out[i] = km / AU_KM
        return out
