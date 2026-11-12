# torino

An independent, open-source check of the public asteroid impact-risk lists, updated every day.

Two centres publish impact probabilities for near-Earth asteroids: NASA/JPL's **Sentry** and ESA's **NEO Coordination Centre (NEOCC)**. Each computes its own orbits and its own list of possible future impacts, and the two lists do not always agree. torino compares them object by object every day, keeps the history, and is building the tools to recompute the probabilities from public observations.

Named after the Torino scale, the 0–10 scale used to communicate an asteroid's impact hazard.

## First result (12 Nov 2026)

| | |
|---|---|
| Objects on JPL Sentry's list | 2,246 |
| Objects on ESA NEOCC's list | 2,208 |
| On both lists | 1,978 |
| Only on JPL's / only on ESA's | 268 / 230 |
| On both, with impact probabilities 10× or more apart | 296 |

For example, JPL gives 2010 VQ an impact probability of 1.6 × 10⁻³ and ESA 3.3 × 10⁻⁸. A difference is not by itself an error: the centres fit orbits to observations up to different times and list potential impacts over different time spans. 101955 Bennu and 29075 (1950 DA), the two highest Palermo-scale objects on JPL's list, are not on ESA's list because their potential impacts (2178–2290 and 2880) fall outside ESA's time span, which ends in 2127.

The daily reports are in [`reports/`](reports/), the full table in [`data/compare.csv`](data/compare.csv), and the git history of `data/` is the archive.

## Status

| Milestone | State |
|---|---|
| M0: daily comparison of JPL Sentry and ESA NEOCC, object by object, with a dated report and the history in git | running daily |
| M1: n-body orbit propagator, validated against JPL Horizons | done: all six test objects within 1 km of JPL Horizons after a year (below) |
| M2: Monte Carlo impact probabilities from public orbit solutions, compared with Sentry and NEOCC for listed objects | planned |
| M3: public page with the daily comparison and its history | planned |

## M0: the daily comparison

`torino/sources.py` parses the two lists into the same record. Objects are matched by asteroid number when numbered and by provisional designation otherwise (JPL writes `2008 JL3`, ESA `2008JL3`). `torino/compare.py` joins them and compares:

- the cumulative impact probability (IP), as `log10(IP_JPL / IP_ESA)`;
- the cumulative and maximum Palermo scale;
- the maximum Torino scale;
- the years of the listed potential impacts;
- which objects are on one list only.

Two impact probabilities count as disagreeing when one is at least 10 times the other. `torino/daily.py` writes `data/compare.csv` (one row per object, sorted so that day-to-day diffs show only what changed), `data/summary.json`, and `reports/YYYY/YYYY-MM-DD.md` with the highest-risk objects, the largest disagreements, the objects on one list only, and what each centre added, removed or changed since the previous run.

```
python -m torino.daily                     # fetch both lists and write data/ and reports/
python -m pytest                           # tests, on real samples of both lists
```

It uses only the Python standard library. `.github/workflows/daily.yml` runs it every day after ESA's usual update and commits the results.

## M1: the propagator

`torino/propagate.py` integrates a small body's barycentric orbit with SciPy's DOP853 (8th-order Runge-Kutta; relative tolerance 10⁻¹², time counted from the start epoch). The force model is Newtonian point-mass gravity from the Sun, the eight planets, the Moon and Pluto, with positions and GM values from JPL's DE440s ephemeris (`torino/ephemeris.py`), plus, optionally, the Sun's general-relativistic (post-Newtonian) acceleration, the 16 most massive main-belt asteroids (`torino/perturbers.py`), and non-gravitational forces. The asteroids are JPL's SB441-N16 set with its GM values (Farnocchia 2021, JPL IOM 392R-21-005); instead of JPL's 616 MB trajectory file, torino integrates their orbits itself from their Horizons states. Non-gravitational forces follow JPL's model for asteroids: constant radial, transverse and normal accelerations A1, A2, A3 at 1 au, scaled by (1 au / r)². The transverse one, A2, is the Yarkovsky effect. torino takes them from the object's current JPL orbit solution through the [Small-Body Database API](https://ssd-api.jpl.nasa.gov/doc/sbdb.html) (`torino/sbdb.py`), and adds radiation pressure when the solution gives an area-to-mass ratio. The tests check it against exact results: in the two-body problem it returns to perihelion within 150 m after one orbit and conserves energy and angular momentum to 10⁻¹⁰, and with relativity the perihelion advances by Einstein's 6πGM/(c²a(1−e²)) per orbit to within 0.1%.

`torino/validate.py` starts from JPL Horizons' state of each object and compares torino's propagated position with Horizons'. It takes Horizons' states relative to the Sun and makes them barycentric with torino's own ephemeris, so the comparison does not depend on which planetary ephemeris a trajectory was computed with (see Bennu, below).

Position difference from JPL Horizons (km), starting from Horizons' heliocentric state at 2026-01-01 TDB:

| Object | Newtonian | + relativity | + asteroids | + non-gravitational |
| --- | ---: | ---: | ---: | ---: |
| 99942 Apophis | 132.9 | 3.1 | 3.4 | **0.1** |
| 101955 Bennu | 58.7 | 3.3 | 3.6 | **0.7** |
| 29075 (1950 DA) | 19.4 | 0.2 | 0.0 | **0.0** |
| 433 Eros | 49.1 | 0.3 | 0.0 | **0.0** |
| 2000 SG344 | 72.2 | 0.4 | 0.0 | **0.0** |
| 2008 JL3 | 5.7 | 3.6 | 0.0 | **0.0** |

After 365 days; each column adds one part of the force model, and 0.0 means under 50 m (the 30-day values are in [`validation/m1.md`](validation/m1.md)). Relativity is the largest thing Newtonian gravity leaves out. With it and the asteroids, four of the six objects agree with JPL to better than 50 m after a year. The other two, Apophis and Bennu, are the two whose JPL orbits include non-gravitational forces; with them, Apophis agrees to 0.1 km (57 m).

Bennu's JPL orbit models these forces physically rather than with A1 and A2: radiation pressure from an estimated area-to-mass ratio, and the Yarkovsky effect from a thermal model with an estimated bulk density (solution 118, Farnocchia et al. 2021). torino includes the radiation pressure, which brings Bennu from 3.6 to 0.7 km, but not the thermal model; the remaining 0.7 km is mainly that, together with the 343 perturbers and the DE424 ephemeris of Bennu's trajectory. The area-to-mass ratio is the SBDB value, 2.64 × 10⁻⁶ m²/kg; Bennu's Horizons record shows 0.00000026, ten times smaller, which does not fit its trajectory (3.2 km after a year).

**Bennu.** Compared in barycentric coordinates, as in the first version of this check, Bennu was 1,805 km from Horizons after a year even with relativity, while the other objects were within 4 km. The difference was radial at first and grew from day 0 as a constant acceleration, which points to a frame mismatch rather than a missing force. Bennu's Horizons trajectory is not a standard orbit solution: it is the reconstruction from the OSIRIS-REx mission (Farnocchia et al. 2021, Icarus 369:114594), computed with the DE424 planetary ephemeris. Its barycentric positions in Horizons are offset by a nearly constant 115.9 km from its heliocentric positions plus DE440's Sun; for Apophis the same offset is 0.0 km. In the barycentric frame torino therefore put the Sun 116 km from where Bennu's trajectory had it. Starting from and comparing with the heliocentric state removes the offset, and Bennu agrees with Horizons to 3.3 km after a year.

```
pip install -r requirements.txt
python -m torino.validate --out validation/m1.md   # downloads DE440s (33 MB) once
```

## Data

- NASA/JPL Center for Near-Earth Object Studies, Sentry: `https://ssd-api.jpl.nasa.gov/sentry.api`
- ESA NEO Coordination Centre, risk list: `https://neo.ssa.esa.int/PSDB-portlet/download?file=esa_risk_list`

torino republishes values from these lists, with credit, for comparison. The terms of reuse of ESA NEOCC's data are being confirmed; until then, cite the two centres as the source of the numbers, not torino.

## License

Apache-2.0. See [LICENSE](LICENSE).
