"""Fetch and parse the two public impact-risk lists.

NASA/JPL's Sentry lists every near-Earth asteroid with a non-zero chance of hitting the Earth in the next century,
through its JSON API. ESA's NEO Coordination Centre (NEOCC) publishes its own risk list as a fixed-width text table.
Both are parsed into the same RiskEntry record, keyed so that the same object gets the same key in both lists.
"""

from __future__ import annotations

import json
import re
import time
import urllib.request
from dataclasses import dataclass

SENTRY_URL = "https://ssd-api.jpl.nasa.gov/sentry.api"
ESA_URL = "https://neo.ssa.esa.int/PSDB-portlet/download?file=esa_risk_list"
USER_AGENT = "torino (https://github.com/Shakhtar-Sankur/torino)"


@dataclass(frozen=True)
class RiskEntry:
    """One object on one agency's risk list."""

    key: str  # asteroid number if numbered, else the provisional designation without spaces
    designation: str  # as the agency writes it
    ip_cum: float  # cumulative impact probability over all listed potential impacts
    ps_cum: float  # cumulative Palermo scale
    ps_max: float  # Palermo scale of the most significant potential impact
    ts_max: int | None  # maximum Torino scale (None when the agency gives n/a)
    first_year: int  # first year of the listed potential impacts
    last_year: int  # last year of the listed potential impacts
    diameter_m: float | None  # estimated diameter in metres
    v_inf_kms: float | None  # velocity relative to the Earth (JPL: v_inf; ESA: impact velocity)


_PROVISIONAL_TAIL = re.compile(r"[A-Z]{1,2}\d*")


def object_key(designation: str) -> str:
    """Key for matching the same object across lists.

    '2008 JL3' (JPL) and '2008JL3' (ESA) -> '2008JL3'; '443104' (JPL), '443104 2013XK22' and '101955 Bennu' (ESA)
    -> the number. A provisional designation starts with a four-digit year, so a leading number only counts as an
    asteroid number when the next word is not the letters-and-digits tail of a provisional designation.
    """
    tokens = designation.split()
    if len(tokens) >= 2 and tokens[0].isdigit() and not _PROVISIONAL_TAIL.fullmatch(tokens[1].upper()):
        return tokens[0]
    return "".join(tokens).upper()


def _float(text: str | None) -> float | None:
    if text is None:
        return None
    text = str(text).strip()
    if text in ("", "n/a", "n/d", "-"):  # ESA writes n/a or n/d where it gives no Torino scale
        return None
    return float(text)


def _years(text: str) -> tuple[int, int]:
    first, _, last = text.strip().partition("-")
    return int(first), int(last or first)


def parse_sentry(data: dict) -> list[RiskEntry]:
    """Parse the JSON returned by the Sentry API's summary mode."""
    entries = []
    for row in data["data"]:
        first, last = _years(row["range"])
        diameter_km = _float(row.get("diameter"))
        ts = _float(row.get("ts_max"))
        entries.append(RiskEntry(
            key=object_key(row["des"]),
            designation=row["des"],
            ip_cum=float(row["ip"]),
            ps_cum=float(row["ps_cum"]),
            ps_max=float(row["ps_max"]),
            ts_max=None if ts is None else int(ts),
            first_year=first,
            last_year=last,
            diameter_m=None if diameter_km is None else diameter_km * 1000.0,
            v_inf_kms=_float(row.get("v_inf")),
        ))
    return entries


_ESA_UPDATED = re.compile(r"Last Update:\s*(\d{4}-\d{2}-\d{2} \d{2}:\d{2}) UTC")


def parse_esa(text: str) -> tuple[list[RiskEntry], str | None]:
    """Parse ESA NEOCC's risk list. Returns the entries and the list's 'Last Update' time (UTC), if given."""
    lines = text.splitlines()
    updated = None
    entries = []
    for line in lines:
        match = _ESA_UPDATED.match(line.strip())
        if match:
            updated = match.group(1)
            continue
        fields = [f.strip() for f in line.split("|")]
        # Data rows have 11 fields before the trailing '|'; header rows have text where the numbers go.
        if len(fields) < 11 or not re.fullmatch(r"\d{4}-\d{4}", fields[8]):
            continue
        first, last = _years(fields[8])
        ts = _float(fields[6])
        entries.append(RiskEntry(
            key=object_key(fields[0]),
            designation=" ".join(fields[0].split()),
            ip_cum=float(fields[9]),
            ps_cum=float(fields[10]),
            ps_max=float(fields[5]),
            ts_max=None if ts is None else int(ts),
            first_year=first,
            last_year=last,
            diameter_m=_float(fields[1]),
            v_inf_kms=_float(fields[7]),
        ))
    return entries, updated


def _get(url: str, attempts: int = 4) -> bytes:
    """GET with a few retries, since a daily job should not fail on one dropped connection."""
    for attempt in range(attempts):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(request, timeout=60) as response:
                return response.read()
        except OSError:
            if attempt == attempts - 1:
                raise
            time.sleep(2 ** (attempt + 1))
    raise AssertionError("unreachable")


def fetch_sentry() -> list[RiskEntry]:
    return parse_sentry(json.loads(_get(SENTRY_URL)))


def fetch_esa() -> tuple[list[RiskEntry], str | None]:
    return parse_esa(_get(ESA_URL).decode("utf-8", errors="replace"))
