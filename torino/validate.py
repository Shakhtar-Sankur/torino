"""Check the propagator against JPL Horizons.

For each object: take Horizons' state at the start epoch, propagate it with torino, and compare the position with
Horizons' at later epochs. Horizons uses JPL's full force model, so the difference measures what torino's force model
leaves out (plus integration error).

Horizons' states are taken relative to the Sun and converted to barycentric with torino's own ephemeris, so that the
comparison does not depend on which planetary ephemeris a trajectory was computed with. This matters for Bennu, whose
Horizons trajectory comes from the OSIRIS-REx mission and was computed with DE424: its barycentric positions in
Horizons are offset by a nearly constant 116 km from its heliocentric positions plus DE440's Sun.

    python -m torino.validate --out validation/m1.md
"""

from __future__ import annotations

import argparse

import numpy as np

from torino.ephemeris import AU_KM, Ephemeris
from torino.horizons import state_vectors
from torino.perturbers import Asteroids, WithAsteroids
from torino.propagate import propagate
from torino.sbdb import nongrav_parameters

# Near-Earth asteroids from the risk lists and well-observed ones, with Horizons designations.
OBJECTS = {
    "99942": "99942 Apophis",
    "101955": "101955 Bennu",
    "29075": "29075 (1950 DA)",
    "433": "433 Eros",
    "2000 SG344": "2000 SG344",
    "2008 JL3": "2008 JL3",
}
EPOCH_JD = 2461041.5  # 2026-01-01 00:00 TDB
DAYS = [30, 365]

# Force models, in the order they were added: whether each uses relativity, the main-belt asteroids and the object's
# non-gravitational parameters from its JPL orbit.
MODELS = {
    "Newtonian": {"relativity": False, "asteroids": False, "nongrav": False},
    "+ relativity": {"relativity": True, "asteroids": False, "nongrav": False},
    "+ asteroids": {"relativity": True, "asteroids": True, "nongrav": False},
    "+ non-gravitational": {"relativity": True, "asteroids": True, "nongrav": True},
}
MODEL_NOTES = [
    "Newtonian: point masses of the Sun, the eight planets, the Moon and Pluto, from DE440s.",
    "+ relativity: adds the Sun's post-Newtonian (1PN) acceleration.",
    "+ asteroids: adds the 16 most massive main-belt asteroids (JPL's SB441-N16 set and GM values).",
    ("+ non-gravitational: adds the non-gravitational accelerations of the object's JPL orbit solution, from the JPL "
     "Small-Body Database: radial, transverse (Yarkovsky) and normal (A1, A2, A3), and radiation pressure from the "
     "area-to-mass ratio. For objects without them this column equals the previous one."),
]


def heliocentric_truth(designation: str, jd: list[float], ephemeris) -> np.ndarray:
    """Horizons' heliocentric states, made barycentric with the ephemeris torino integrates with."""
    helio = state_vectors(designation, jd, center="500@10")
    return helio + np.array([np.concatenate(ephemeris.sun_state(t)) for t in jd])


def validate(objects: dict[str, str], jd0: float, days: list[int], ephemeris,
             models: dict[str, dict] = MODELS) -> list[tuple[str, str, int, float]]:
    """Rows of (object, model, days, position difference in km)."""
    with_asteroids = None
    if any(m["asteroids"] for m in models.values()):
        with_asteroids = WithAsteroids(ephemeris, Asteroids(ephemeris, jd0, jd0 + max(days)))
    rows = []
    for designation, label in objects.items():
        jd = [jd0] + [jd0 + d for d in days]
        truth = heliocentric_truth(designation, jd, ephemeris)
        a123 = nongrav_parameters(designation) if any(m["nongrav"] for m in models.values()) else None
        for model, m in models.items():
            ours = propagate(truth[0], jd0, jd[1:], with_asteroids if m["asteroids"] else ephemeris,
                             relativity=m["relativity"], nongrav=a123 if m["nongrav"] and any(a123) else None)
            for d, mine, ref in zip(days, ours, truth[1:]):
                rows.append((label, model, d, float(np.linalg.norm(mine[:3] - ref[:3]) * AU_KM)))
    return rows


def markdown(rows: list[tuple[str, str, int, float]], days: list[int], notes: list[str] = MODEL_NOTES) -> str:
    by = {(label, model, d): err for label, model, d, err in rows}
    labels = list(dict.fromkeys(row[0] for row in rows))
    models = list(dict.fromkeys(row[1] for row in rows))
    columns = [(m, d) for m in models for d in days]
    lines = [("Position difference from JPL Horizons (km), starting from Horizons' heliocentric state at "
              "2026-01-01 TDB."), "",
             "| Object | " + " | ".join(f"{m}, {d} days" for m, d in columns) + " |",
             "| --- | " + " | ".join("---:" for _ in columns) + " |"]
    for label in labels:
        lines.append(f"| {label} | " + " | ".join(f"{by[(label, m, d)]:,.1f}" for m, d in columns) + " |")
    lines += [""] + [f"- {note}" for note in notes]
    return "\n".join(lines) + "\n"


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", help="write the markdown table here")
    args = parser.parse_args(argv)
    text = markdown(validate(OBJECTS, EPOCH_JD, DAYS, Ephemeris()), DAYS)
    print(text)
    if args.out:
        with open(args.out, "w") as f:
            f.write(text)


if __name__ == "__main__":
    main()
