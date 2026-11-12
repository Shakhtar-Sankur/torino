"""The 16 most massive main-belt asteroids, which JPL includes as perturbers in small-body orbits.

The list and the GM values are JPL's SB441-N16 set: Table 1 of D. Farnocchia (2021), "Small-Body Perturber Files
SB441-N16 and SB441-N343", JPL IOM 392R-21-005, whose GM values are those of DE441 (Park et al. 2021).

JPL distributes the asteroids' trajectories as a 616 MB file. torino instead starts from each asteroid's state in
JPL Horizons and integrates its orbit under the Sun, the planets, the Moon, Pluto and relativity, which is enough for
their pull on another small body: the asteroids are far away, and an error of a few km in their positions changes
that pull negligibly. Their pull on each other is left out.
"""

from __future__ import annotations

import numpy as np

from torino.horizons import state_vectors
from torino.propagate import propagate_dense

# Name: (asteroid number, GM in au^3/day^2), in the order of Table 1 (by mass).
PERTURBERS = {
    "Ceres": ("1", 1.3964518123081070e-13),
    "Vesta": ("4", 3.8548000225257904e-14),
    "Pallas": ("2", 3.0471146330043200e-14),
    "Hygiea": ("10", 1.2542530761640810e-14),
    "Davida": ("511", 8.6836253492286545e-15),
    "Interamnia": ("704", 6.3110343420878887e-15),
    "Europa": ("52", 5.9824315264869841e-15),
    "Sylvia": ("87", 4.8345606546105521e-15),
    "Eunomia": ("15", 4.5107799051436795e-15),
    "Juno": ("3", 4.2823439677995011e-15),
    "Psyche": ("16", 3.5445002842488978e-15),
    "Camilla": ("107", 3.2191392075878588e-15),
    "Thisbe": ("88", 2.6529436610356353e-15),
    "Iris": ("7", 2.5416014973471498e-15),
    "Euphrosyne": ("31", 2.4067012218937576e-15),
    "Cybele": ("65", 2.0917175955133682e-15),
}


def horizons_states(jd: float, planets, names: list[str]) -> np.ndarray:
    """Barycentric states of the asteroids at jd: Horizons' heliocentric states plus the planetary ephemeris' Sun."""
    sun = np.concatenate(planets.sun_state(jd))
    return np.array([state_vectors(PERTURBERS[name][0], [jd], center="500@10")[0] + sun for name in names])


class Asteroids:
    """Positions of the asteroids between jd_start and jd_end, from orbits integrated by torino."""

    def __init__(self, planets, jd_start: float, jd_end: float, names: list[str] | None = None,
                 states: np.ndarray | None = None, relativity: bool = True):
        self.bodies = names or list(PERTURBERS)
        self.gm = np.array([PERTURBERS[name][1] for name in self.bodies])
        if states is None:
            states = horizons_states(jd_start, planets, self.bodies)
        self.orbits = [propagate_dense(s, jd_start, jd_end, planets, relativity=relativity) for s in states]

    def positions(self, jd_tdb: float) -> np.ndarray:
        return np.array([orbit(jd_tdb)[:3] for orbit in self.orbits])


class WithAsteroids:
    """A planetary ephemeris with the asteroids added as more gravitating bodies."""

    def __init__(self, planets, asteroids: Asteroids):
        self.planets, self.asteroids = planets, asteroids
        self.bodies = list(planets.bodies) + asteroids.bodies
        self.gm = np.concatenate([planets.gm, asteroids.gm])

    def positions(self, jd_tdb: float) -> np.ndarray:
        return np.vstack([self.planets.positions(jd_tdb), self.asteroids.positions(jd_tdb)])

    def sun_state(self, jd_tdb: float) -> tuple[np.ndarray, np.ndarray]:
        return self.planets.sun_state(jd_tdb)
