"""Propagate a small body's barycentric orbit under the gravity of the Sun, the planets, the Moon, Pluto and,
optionally, the largest main-belt asteroids (torino.perturbers), and optionally non-gravitational forces.

The small body is a test particle: it feels the bodies of the ephemeris but does not move them. The force model is
Newtonian point-mass gravity, plus, with relativity=True, the Sun's post-Newtonian (1PN) term in the PPN form with
beta = gamma = 1, as in the equations of motion JPL uses for small bodies (Moyer 2003, eq. 4-61, Sun only).

Non-gravitational forces use the model of JPL's asteroid orbits: constant radial, transverse and normal accelerations
A1, A2 and A3 (au/day^2) at 1 au from the Sun, scaled by (1 au / r)^2. A2 is the Yarkovsky effect; A1 is mostly
radiation pressure.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

C_AU_DAY = 299_792.458 * 86_400.0 / 149_597_870.7  # speed of light in au/day


def acceleration(r: np.ndarray, body_positions: np.ndarray, gm: np.ndarray) -> np.ndarray:
    """Point-mass gravity at position r (au) from bodies at body_positions (au), in au/day^2."""
    d = body_positions - r
    dist3 = np.einsum("ij,ij->i", d, d) ** 1.5
    return (gm[:, None] * d / dist3[:, None]).sum(axis=0)


def sun_relativity(r: np.ndarray, v: np.ndarray, gm_sun: float) -> np.ndarray:
    """1PN acceleration from the Sun (au/day^2); r and v are the body's position and velocity relative to the Sun."""
    rn = np.linalg.norm(r)
    return gm_sun / (C_AU_DAY**2 * rn**3) * ((4.0 * gm_sun / rn - v @ v) * r + 4.0 * (r @ v) * v)


def nongrav_acceleration(r: np.ndarray, v: np.ndarray, a123) -> np.ndarray:
    """Radial, transverse and normal acceleration (au/day^2) scaled by (1 au / r)^2; r and v relative to the Sun."""
    rn = np.linalg.norm(r)
    radial = r / rn
    normal = np.cross(r, v)
    normal /= np.linalg.norm(normal)
    transverse = np.cross(normal, radial)
    a1, a2, a3 = a123
    return (a1 * radial + a2 * transverse + a3 * normal) / rn**2


def _rhs(jd0: float, ephemeris, relativity: bool, nongrav=None):
    """Right-hand side of the equations of motion, with time in days from jd0."""
    sun = ephemeris.bodies.index("Sun") if relativity else None
    heliocentric = relativity or nongrav is not None

    def rhs(t, y):
        positions = ephemeris.positions(jd0 + t)
        a = acceleration(y[:3], positions, ephemeris.gm)
        if heliocentric:
            sun_r, sun_v = ephemeris.sun_state(jd0 + t)
            r, v = y[:3] - sun_r, y[3:] - sun_v
            if relativity:
                a = a + sun_relativity(r, v, ephemeris.gm[sun])
            if nongrav is not None:
                a = a + nongrav_acceleration(r, v, nongrav)
        return np.concatenate([y[3:], a])

    return rhs


def propagate(state0: np.ndarray, jd0: float, jd_out: list[float], ephemeris,
              rtol: float = 1e-12, atol: float = 1e-15, relativity: bool = False, nongrav=None) -> np.ndarray:
    """Integrate [x, y, z, vx, vy, vz] (au, au/day) from jd0 to each of jd_out (TDB Julian dates).

    ephemeris needs .positions(jd) -> (n, 3) array in au and .gm -> (n,) array in au^3/day^2; with relativity it also
    needs .bodies (names, including "Sun"), and with relativity or nongrav it needs .sun_state(jd) -> (position,
    velocity). nongrav is (A1, A2, A3) in au/day^2, or None.
    Time is measured in days from jd0 inside the integrator, so the step size is not limited by the size of the JD.
    Returns an array of shape (len(jd_out), 6).
    """
    rhs = _rhs(jd0, ephemeris, relativity, nongrav)
    t_out = np.asarray(jd_out, dtype=float) - jd0
    span = (0.0, float(t_out.max())) if t_out.max() > 0 else (0.0, float(t_out.min()))
    sol = solve_ivp(rhs, span, np.asarray(state0, dtype=float), method="DOP853",
                    t_eval=np.sort(t_out) if span[1] > 0 else np.sort(t_out)[::-1], rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    by_t = dict(zip(sol.t, sol.y.T))
    return np.array([by_t[t] for t in t_out])


def propagate_dense(state0: np.ndarray, jd0: float, jd1: float, ephemeris,
                    rtol: float = 1e-12, atol: float = 1e-15, relativity: bool = False):
    """Integrate from jd0 to jd1 and return a function giving the state at any TDB Julian date in between.

    It uses the integrator's own 7th-order interpolant, which is as accurate as the steps themselves.
    """
    sol = solve_ivp(_rhs(jd0, ephemeris, relativity), (0.0, jd1 - jd0), np.asarray(state0, dtype=float),
                    method="DOP853", dense_output=True, rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    lo, hi = sorted((jd0, jd1))

    def state(jd_tdb: float) -> np.ndarray:
        if not lo <= jd_tdb <= hi:
            raise ValueError(f"JD {jd_tdb} is outside the integrated interval {lo}-{hi}")
        return sol.sol(jd_tdb - jd0)

    return state
