"""Offline checks of the integrator against the exact two-body solution (no ephemeris file or network needed)."""

import numpy as np
import pytest

from torino.ephemeris import gm_au3_day2
from torino.horizons import horizons_command, parse_vectors
from torino.propagate import propagate


class SunOnly:
    """An ephemeris with only the Sun, fixed at the origin."""

    gm = np.array([gm_au3_day2("Sun")])

    def positions(self, jd):
        return np.zeros((1, 3))


def kepler_orbit(a=1.3, e=0.4):
    """Perihelion state of an elliptical orbit around the fixed Sun, and its period in days."""
    mu = SunOnly.gm[0]
    rp = a * (1 - e)
    vp = np.sqrt(mu * (1 + e) / rp)
    return np.array([rp, 0.0, 0.0, 0.0, vp, 0.0]), 2 * np.pi * np.sqrt(a**3 / mu)


def test_returns_to_perihelion_after_one_period():
    state0, period = kepler_orbit()
    jd0 = 2461041.5
    end = propagate(state0, jd0, [jd0 + period], SunOnly())[0]
    assert np.linalg.norm(end[:3] - state0[:3]) < 1e-9  # au, about 150 m


def test_conserves_energy_and_angular_momentum():
    state0, period = kepler_orbit(a=2.1, e=0.6)
    mu = SunOnly.gm[0]
    jd0 = 2461041.5
    states = propagate(state0, jd0, list(jd0 + np.linspace(1, 3 * period, 7)), SunOnly())

    def energy(s):
        return 0.5 * s[3:] @ s[3:] - mu / np.linalg.norm(s[:3])

    for s in states:
        assert energy(s) == pytest.approx(energy(state0), rel=1e-10)
        assert np.cross(s[:3], s[3:]) == pytest.approx(np.cross(state0[:3], state0[3:]), rel=1e-10)


def test_backwards_and_unsorted_times():
    state0, period = kepler_orbit()
    jd0 = 2461041.5
    fwd = propagate(state0, jd0, [jd0 + 40.0, jd0 + 10.0], SunOnly())
    back = propagate(fwd[1], jd0 + 10.0, [jd0], SunOnly())[0]
    assert np.linalg.norm(back[:3] - state0[:3]) < 1e-11
    assert fwd.shape == (2, 6)


def test_horizons_helpers():
    assert horizons_command("99942") == "99942;"
    assert horizons_command("2008 JL3") == "DES=2008 JL3;"
    result = ("header\n$$SOE\n2461041.500000000, A.D. 2026-Jan-01 00:00:00.0000,  6.306953539299172E-01, "
              "-3.792320854325312E-01, -1.251234635815718E-01,  1.134397706979294E-02,  1.728653038821504E-02,  "
              "6.708869319172678E-03,\n$$EOE\n")
    v = parse_vectors(result)
    assert v.shape == (1, 6)
    assert v[0, 0] == pytest.approx(0.6306953539299172)


class SunOnlyWithVelocity(SunOnly):
    """The fixed Sun, with what the relativistic term needs."""

    bodies = ("Sun",)

    def sun_state(self, jd):
        return np.zeros(3), np.zeros(3)


def perihelion_angle(state, mu):
    """Direction of the Laplace-Runge-Lenz (eccentricity) vector in the orbital plane."""
    r, v = state[:3], state[3:]
    e_vec = np.cross(v, np.cross(r, v)) / mu - r / np.linalg.norm(r)
    return np.arctan2(e_vec[1], e_vec[0])


def test_relativistic_perihelion_precession():
    # A tight, eccentric orbit so the effect is large; compare with Einstein's 6 pi GM / (c^2 a (1 - e^2)) per orbit.
    from torino.propagate import C_AU_DAY

    a, e, orbits = 0.2, 0.8, 20
    state0, period = kepler_orbit(a=a, e=e)
    mu = SunOnly.gm[0]
    jd0 = 2461041.5
    end = propagate(state0, jd0, [jd0 + orbits * period], SunOnlyWithVelocity(), relativity=True)[0]
    expected = orbits * 6 * np.pi * mu / (C_AU_DAY**2 * a * (1 - e**2))
    assert perihelion_angle(end, mu) == pytest.approx(expected, rel=1e-3)
    newtonian = propagate(state0, jd0, [jd0 + orbits * period], SunOnlyWithVelocity())[0]
    assert abs(perihelion_angle(newtonian, mu)) < 1e-3 * expected


def test_asteroids_follow_their_orbits_and_add_their_pull():
    from torino.perturbers import PERTURBERS, Asteroids, WithAsteroids
    from torino.propagate import acceleration

    planets = SunOnlyWithVelocity()
    jd0 = 2461041.5
    ceres0, period = kepler_orbit(a=2.77, e=0.08)
    asteroids = Asteroids(planets, jd0, jd0 + period, names=["Ceres"], states=np.array([ceres0]), relativity=False)
    # The dense trajectory starts at the given state and, without relativity, comes back to it after one orbit.
    assert asteroids.positions(jd0)[0] == pytest.approx(ceres0[:3], abs=1e-12)
    assert np.linalg.norm(asteroids.positions(jd0 + period)[0] - ceres0[:3]) < 1e-9
    with pytest.raises(ValueError):
        asteroids.positions(jd0 + period + 1)

    both = WithAsteroids(planets, asteroids)
    assert both.bodies == ["Sun", "Ceres"]
    gm_ceres = PERTURBERS["Ceres"][1]
    r = np.array([2.0, 0.5, 0.0])
    extra = acceleration(r, both.positions(jd0), both.gm) - acceleration(r, planets.positions(jd0), planets.gm)
    d = ceres0[:3] - r
    assert extra == pytest.approx(gm_ceres * d / np.linalg.norm(d) ** 3, rel=1e-12)
