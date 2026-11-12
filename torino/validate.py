"""Check the propagator against JPL Horizons.

For each object: take Horizons' barycentric state at the start epoch, propagate it with torino, and compare the
position with Horizons' at later epochs. Horizons uses JPL's full force model, so the difference measures what
torino's force model leaves out (plus integration error).

    python -m torino.validate --out validation/m1.md
"""

from __future__ import annotations

import argparse

import numpy as np

from torino.ephemeris import AU_KM, Ephemeris
from torino.horizons import state_vectors
from torino.propagate import propagate

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


def validate(objects: dict[str, str], jd0: float, days: list[int], ephemeris) -> list[tuple[str, int, float]]:
    rows = []
    for designation, label in objects.items():
        jd = [jd0] + [jd0 + d for d in days]
        truth = state_vectors(designation, jd)
        ours = propagate(truth[0], jd0, jd[1:], ephemeris)
        for d, mine, ref in zip(days, ours, truth[1:]):
            rows.append((label, d, float(np.linalg.norm(mine[:3] - ref[:3]) * AU_KM)))
    return rows


def markdown(rows: list[tuple[str, int, float]], days: list[int], model: str) -> str:
    by = {(label, d): err for label, d, err in rows}
    labels = list(dict.fromkeys(label for label, _, _ in rows))
    lines = [f"Position difference from JPL Horizons (km), starting from Horizons' state at 2026-01-01 TDB. "
             f"Force model: {model}.", "",
             "| Object | " + " | ".join(f"{d} days" for d in days) + " |",
             "| --- | " + " | ".join("---:" for _ in days) + " |"]
    for label in labels:
        lines.append(f"| {label} | " + " | ".join(f"{by[(label, d)]:,.1f}" for d in days) + " |")
    return "\n".join(lines) + "\n"


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", help="write the markdown table here")
    args = parser.parse_args(argv)
    model = "Newtonian point masses: Sun, eight planets, Moon, Pluto (DE440s)"
    text = markdown(validate(OBJECTS, EPOCH_JD, DAYS, Ephemeris()), DAYS, model)
    print(text)
    if args.out:
        with open(args.out, "w") as f:
            f.write(text)


if __name__ == "__main__":
    main()
