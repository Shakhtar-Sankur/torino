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
| M1: n-body orbit propagator, validated against JPL Horizons | next |
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

## Data

- NASA/JPL Center for Near-Earth Object Studies, Sentry: `https://ssd-api.jpl.nasa.gov/sentry.api`
- ESA NEO Coordination Centre, risk list: `https://neo.ssa.esa.int/PSDB-portlet/download?file=esa_risk_list`

torino republishes values from these lists, with credit, for comparison. The terms of reuse of ESA NEOCC's data are being confirmed; until then, cite the two centres as the source of the numbers, not torino.

## License

Apache-2.0. See [LICENSE](LICENSE).
