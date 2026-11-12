"""Compare NASA/JPL Sentry and ESA NEOCC risk lists object by object."""

from __future__ import annotations

import math
from dataclasses import dataclass

from torino.sources import RiskEntry

# Two impact probabilities "disagree" when one is at least this many times the other.
DISAGREEMENT_FACTOR = 10.0


@dataclass(frozen=True)
class Row:
    """One object, with what each agency says about it."""

    key: str
    jpl: RiskEntry | None
    esa: RiskEntry | None

    @property
    def listed_by(self) -> str:
        if self.jpl and self.esa:
            return "both"
        return "jpl" if self.jpl else "esa"

    @property
    def designation(self) -> str:
        return (self.jpl or self.esa).designation

    @property
    def log10_ip_ratio(self) -> float | None:
        """log10(JPL impact probability / ESA impact probability); positive when JPL gives the higher probability."""
        if not (self.jpl and self.esa) or self.jpl.ip_cum <= 0 or self.esa.ip_cum <= 0:
            return None
        return math.log10(self.jpl.ip_cum / self.esa.ip_cum)

    @property
    def ps_cum_max(self) -> float:
        """The higher of the two cumulative Palermo scales, for ranking objects by risk."""
        return max(e.ps_cum for e in (self.jpl, self.esa) if e)

    @property
    def disagrees(self) -> bool:
        ratio = self.log10_ip_ratio
        return ratio is not None and abs(ratio) >= math.log10(DISAGREEMENT_FACTOR)


def compare(jpl: list[RiskEntry], esa: list[RiskEntry]) -> list[Row]:
    """Join the two lists on the object key; rows are sorted by key so the output is stable from day to day."""
    by_jpl = {e.key: e for e in jpl}
    by_esa = {e.key: e for e in esa}
    return [Row(key, by_jpl.get(key), by_esa.get(key)) for key in sorted(by_jpl.keys() | by_esa.keys())]


def summary(rows: list[Row]) -> dict:
    both = [r for r in rows if r.listed_by == "both"]
    torino_differs = [r for r in both if r.jpl.ts_max is not None and r.esa.ts_max is not None
                      and r.jpl.ts_max != r.esa.ts_max]
    return {
        "objects": len(rows),
        "jpl": sum(1 for r in rows if r.jpl),
        "esa": sum(1 for r in rows if r.esa),
        "both": len(both),
        "jpl_only": sum(1 for r in rows if r.listed_by == "jpl"),
        "esa_only": sum(1 for r in rows if r.listed_by == "esa"),
        "ip_disagree_10x": sum(1 for r in both if r.disagrees),
        "torino_scale_differs": len(torino_differs),
    }


def changes(previous: dict[str, dict], rows: list[Row]) -> list[str]:
    """Describe what changed since the previous run, given the previous compare.csv rows keyed by object key."""
    out = []
    current = {r.key: r for r in rows}
    for agency in ("jpl", "esa"):
        name = "JPL Sentry" if agency == "jpl" else "ESA NEOCC"
        before = {k for k, p in previous.items() if p.get(f"ip_{agency}")}
        now = {k for k, r in current.items() if getattr(r, agency)}
        for key in sorted(now - before):
            out.append(f"{name} added {current[key].designation}")
        for key in sorted(before - now):
            out.append(f"{name} removed {previous[key]['designation']}")
        for key in sorted(now & before):
            entry = getattr(current[key], agency)
            old_ts = previous[key].get(f"ts_{agency}")
            new_ts = "" if entry.ts_max is None else str(entry.ts_max)
            if old_ts != new_ts:
                out.append(f"{name} changed the Torino scale of {entry.designation}: {old_ts or 'n/a'} -> {new_ts or 'n/a'}")
            old_ip = float(previous[key][f"ip_{agency}"])
            if old_ip > 0 and entry.ip_cum > 0 and abs(math.log10(entry.ip_cum / old_ip)) >= 1:
                out.append(f"{name} changed the impact probability of {entry.designation}: {old_ip:.3g} -> {entry.ip_cum:.3g}")
    return out
