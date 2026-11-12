"""Propagate a small body's barycentric orbit under the gravity of the Sun, the planets, the Moon and Pluto.

The small body is a test particle: it feels the bodies of the ephemeris but does not move them. The force model is
Newtonian point-mass gravity, plus, with relativity=True, the Sun's post-Newtonian (1PN) term in the PPN form with
beta = gamma = 1, as in the equations of motion JPL uses for small bodies (Moyer 2003, eq. 4-61, Sun only).
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


def propagate(state0: np.ndarray, jd0: float, jd_out: list[float], ephemeris,
              rtol: float = 1e-12, atol: float = 1e-15, relativity: bool = False) -> np.ndarray:
    """Integrate [x, y, z, vx, vy, vz] (au, au/day) from jd0 to each of jd_out (TDB Julian dates).

    ephemeris needs .positions(jd) -> (n, 3) array in au and .gm -> (n,) array in au^3/day^2; with relativity it also
    needs .bodies (names, including "Sun") and .sun_state(jd) -> (position, velocity).
    Time is measured in days from jd0 inside the integrator, so the step size is not limited by the size of the JD.
    Returns an array of shape (len(jd_out), 6).
    """
    sun = ephemeris.bodies.index("Sun") if relativity else None

    def rhs(t, y):
        positions = ephemeris.positions(jd0 + t)
        a = acceleration(y[:3], positions, ephemeris.gm)
        if relativity:
            sun_r, sun_v = ephemeris.sun_state(jd0 + t)
            a = a + sun_relativity(y[:3] - sun_r, y[3:] - sun_v, ephemeris.gm[sun])
        return np.concatenate([y[3:], a])

    t_out = np.asarray(jd_out, dtype=float) - jd0
    span = (0.0, float(t_out.max())) if t_out.max() > 0 else (0.0, float(t_out.min()))
    sol = solve_ivp(rhs, span, np.asarray(state0, dtype=float), method="DOP853",
                    t_eval=np.sort(t_out) if span[1] > 0 else np.sort(t_out)[::-1], rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    by_t = dict(zip(sol.t, sol.y.T))
    return np.array([by_t[t] for t in t_out])
