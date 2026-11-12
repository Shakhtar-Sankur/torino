"""Daily run: fetch both risk lists, compare them, and write the table, a summary and a dated report.

    python -m torino.daily --data data --reports reports

The previous data/compare.csv, if there is one, is read first so the report can list what changed since the last run.
The git history of data/ is the archive.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import os

from torino.compare import DISAGREEMENT_FACTOR, Row, changes, compare, summary
from torino.sources import ESA_URL, SENTRY_URL, fetch_esa, fetch_sentry, parse_esa, parse_sentry

COLUMNS = ["key", "designation", "listed_by",
           "ip_jpl", "ip_esa", "log10_ip_ratio",
           "ps_cum_jpl", "ps_cum_esa", "ts_jpl", "ts_esa",
           "years_jpl", "years_esa", "diameter_m_jpl", "diameter_m_esa"]


def _fmt(value, spec: str) -> str:
    return "" if value is None else format(value, spec)


def row_to_csv(r: Row) -> dict:
    j, e = r.jpl, r.esa
    return {
        "key": r.key,
        "designation": r.designation,
        "listed_by": r.listed_by,
        "ip_jpl": _fmt(j and j.ip_cum, ".4g"),
        "ip_esa": _fmt(e and e.ip_cum, ".4g"),
        "log10_ip_ratio": _fmt(r.log10_ip_ratio, ".3f"),
        "ps_cum_jpl": _fmt(j and j.ps_cum, ".2f"),
        "ps_cum_esa": _fmt(e and e.ps_cum, ".2f"),
        "ts_jpl": _fmt(j and j.ts_max, "d"),
        "ts_esa": _fmt(e and e.ts_max, "d"),
        "years_jpl": f"{j.first_year}-{j.last_year}" if j else "",
        "years_esa": f"{e.first_year}-{e.last_year}" if e else "",
        "diameter_m_jpl": _fmt(j and j.diameter_m, ".0f"),
        "diameter_m_esa": _fmt(e and e.diameter_m, ".0f"),
    }


def read_previous(path: str) -> dict[str, dict]:
    if not os.path.exists(path):
        return {}
    with open(path, newline="") as f:
        return {row["key"]: row for row in csv.DictReader(f)}


def write_csv(path: str, rows: list[Row]) -> None:
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS, lineterminator="\n")
        writer.writeheader()
        for r in rows:
            writer.writerow(row_to_csv(r))


def _ip(e) -> str:
    return "-" if e is None else f"{e.ip_cum:.2e}"


def _ps(e) -> str:
    return "-" if e is None else f"{e.ps_cum:.2f}"


def _ts(e) -> str:
    return "-" if e is None else ("n/a" if e.ts_max is None else str(e.ts_max))


def _years(e) -> str:
    return "-" if e is None else f"{e.first_year}-{e.last_year}"


def _table(rows: list[Row]) -> list[str]:
    out = ["| Object | Listed by | IP, JPL | IP, ESA | Palermo, JPL | Palermo, ESA | Torino, JPL | Torino, ESA "
           "| Years, JPL | Years, ESA |",
           "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for r in rows:
        out.append(f"| {r.designation} | {r.listed_by} | {_ip(r.jpl)} | {_ip(r.esa)} | {_ps(r.jpl)} | {_ps(r.esa)} "
                   f"| {_ts(r.jpl)} | {_ts(r.esa)} | {_years(r.jpl)} | {_years(r.esa)} |")
    return out


def report(day: str, rows: list[Row], stats: dict, diffs: list[str], esa_updated: str | None) -> str:
    both = [r for r in rows if r.listed_by == "both"]
    top = sorted(rows, key=lambda r: r.ps_cum_max, reverse=True)[:10]
    disagree = sorted((r for r in both if r.disagrees), key=lambda r: abs(r.log10_ip_ratio), reverse=True)[:15]
    jpl_only = sorted((r for r in rows if r.listed_by == "jpl"), key=lambda r: r.ps_cum_max, reverse=True)[:10]
    esa_only = sorted((r for r in rows if r.listed_by == "esa"), key=lambda r: r.ps_cum_max, reverse=True)[:10]
    lines = [
        f"# NASA vs ESA impact-risk lists, {day}",
        "",
        f"JPL Sentry lists {stats['jpl']} objects and ESA NEOCC lists {stats['esa']} "
        f"(ESA list updated {esa_updated or 'at an unknown time'} UTC). {stats['both']} are on both lists, "
        f"{stats['jpl_only']} only on JPL's and {stats['esa_only']} only on ESA's. For {stats['ip_disagree_10x']} of the "
        f"objects on both lists, one agency's impact probability is at least {DISAGREEMENT_FACTOR:g} times the other's, "
        f"and the Torino scale differs for {stats['torino_scale_differs']}.",
        "",
        "A difference is not by itself an error. The two centres compute their orbits independently, from "
        "observations up to different times, and list potential impacts over different time spans (the Years "
        "columns), so the same object can get different numbers, or appear on one list only.",
        "",
        "## Highest risk (by the higher cumulative Palermo scale)",
        "",
        *_table(top),
        "",
        f"## Largest disagreements in impact probability (factor {DISAGREEMENT_FACTOR:g} or more)",
        "",
        *(_table(disagree) if disagree else ["None today."]),
        "",
        "## Only on JPL's list (highest Palermo scale first)",
        "",
        *(_table(jpl_only) if jpl_only else ["None today."]),
        "",
        "## Only on ESA's list (highest Palermo scale first)",
        "",
        *(_table(esa_only) if esa_only else ["None today."]),
        "",
        "## Changes since the previous run",
        "",
        *([f"- {d}" for d in diffs] if diffs else ["None."]),
        "",
        "---",
        f"Data: NASA/JPL CNEOS Sentry ({SENTRY_URL}) and the ESA NEO Coordination Centre risk list ({ESA_URL}). "
        "IP is the cumulative impact probability; Palermo is the cumulative Palermo scale; Torino is the maximum "
        "Torino scale.",
        "",
    ]
    return "\n".join(lines)


def run(day: str, jpl, esa, esa_updated, data_dir: str, reports_dir: str) -> dict:
    rows = compare(jpl, esa)
    stats = summary(rows)
    csv_path = os.path.join(data_dir, "compare.csv")
    diffs = changes(read_previous(csv_path), rows) if os.path.exists(csv_path) else []

    os.makedirs(data_dir, exist_ok=True)
    write_csv(csv_path, rows)
    with open(os.path.join(data_dir, "summary.json"), "w") as f:
        json.dump({"date": day, "esa_list_updated_utc": esa_updated, **stats}, f, indent=1, sort_keys=True)
        f.write("\n")

    year_dir = os.path.join(reports_dir, day[:4])
    os.makedirs(year_dir, exist_ok=True)
    with open(os.path.join(year_dir, f"{day}.md"), "w") as f:
        f.write(report(day, rows, stats, diffs, esa_updated))
    return stats


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data", default="data")
    parser.add_argument("--reports", default="reports")
    parser.add_argument("--date", default=dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d"))
    parser.add_argument("--sentry-file", help="read Sentry JSON from a file instead of the API")
    parser.add_argument("--esa-file", help="read the ESA risk list from a file instead of the web")
    args = parser.parse_args(argv)

    if args.sentry_file:
        with open(args.sentry_file) as f:
            jpl = parse_sentry(json.load(f))
    else:
        jpl = fetch_sentry()
    if args.esa_file:
        with open(args.esa_file) as f:
            esa, esa_updated = parse_esa(f.read())
    else:
        esa, esa_updated = fetch_esa()

    # A short list means a partial download or a changed format; better to fail than to record it as real changes.
    if len(jpl) < 500 and not args.sentry_file or len(esa) < 500 and not args.esa_file:
        raise SystemExit(f"suspiciously short lists: JPL {len(jpl)}, ESA {len(esa)}")
    stats = run(args.date, jpl, esa, esa_updated, args.data, args.reports)
    print(json.dumps(stats))


if __name__ == "__main__":
    main()
