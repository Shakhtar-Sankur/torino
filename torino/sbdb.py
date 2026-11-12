"""Non-gravitational parameters of JPL's orbit solutions, from the JPL Small-Body Database (SBDB) API."""

from __future__ import annotations

import json
import urllib.parse

from torino.ephemeris import AU_KM, DAY_S
from torino.sources import _get

SBDB_URL = "https://ssd-api.jpl.nasa.gov/sbdb.api"
# Radiation pressure on a perfect absorber at 1 au: the nominal total solar irradiance (IAU 2015 B3) over c, in N/m^2.
SOLAR_PRESSURE_1AU = 1361.0 / 299_792_458.0


def nongrav_parameters(designation: str) -> tuple[float, float, float]:
    """(A1, A2, A3) in au/day^2 for the object's current JPL orbit solution; (0, 0, 0) if it has none."""
    data = json.loads(_get(f"{SBDB_URL}?{urllib.parse.urlencode({'sstr': designation})}"))
    return parse_model_pars(data["orbit"].get("model_pars") or [])


def parse_model_pars(model_pars: list[dict]) -> tuple[float, float, float]:
    """Turn SBDB's model_pars into (A1, A2, A3) for torino's (1 au / r)^2 model.

    The A1, A2 and A3 parameters are used as they are, if the solution uses the same law (NM = 2, NK = 0, ALN = 1,
    R0 = 1 au); other laws, such as the ones used for comets, raise ValueError. Radiation pressure from an
    area-to-mass ratio (AMRAT) is added to A1. A thermal Yarkovsky model given by a bulk density (RHO), as in
    Bennu's mission-based solution, is not modelled and is left out.
    """
    pars = {p["name"]: float(p["value"]) for p in model_pars}
    a123 = [pars.get(name, 0.0) for name in ("A1", "A2", "A3")]
    if any(a123):
        law = {name: pars.get(name, default) for name, default in (("NM", 2.0), ("NK", 0.0), ("ALN", 1.0), ("R0", 1.0))}
        if law != {"NM": 2.0, "NK": 0.0, "ALN": 1.0, "R0": 1.0}:
            raise ValueError(f"non-gravitational law {law} is not the (1 au / r)^2 law torino models")
    if "AMRAT" in pars:
        a123[0] += pars["AMRAT"] * SOLAR_PRESSURE_1AU * DAY_S**2 / (AU_KM * 1000.0)
    return tuple(a123)
