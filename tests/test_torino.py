import csv
import json
import math
import os

import pytest

from torino.compare import changes, compare, summary
from torino.daily import read_previous, run
from torino.sources import object_key, parse_esa, parse_sentry

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures")


@pytest.fixture
def jpl():
    with open(os.path.join(FIXTURES, "sentry_sample.json")) as f:
        return parse_sentry(json.load(f))


@pytest.fixture
def esa():
    with open(os.path.join(FIXTURES, "esa_sample.txt")) as f:
        return parse_esa(f.read())


@pytest.mark.parametrize("designation, key", [
    ("2011 TO", "2011TO"),            # JPL provisional: a year, not an asteroid number
    ("2012 BA102", "2012BA102"),
    ("2008JL3", "2008JL3"),           # ESA provisional, no space
    ("443104", "443104"),             # JPL numbered
    ("443104               2013XK22", "443104"),  # ESA numbered, with its provisional designation
    ("101955 Bennu", "101955"),       # numbered with a name
])
def test_object_key(designation, key):
    assert object_key(designation) == key


def test_parse_sentry(jpl):
    by_key = {e.key: e for e in jpl}
    assert set(by_key) == {"2011TO", "2008JL3", "443104", "101955", "2023VD3", "29075"}
    bennu = by_key["101955"]
    assert bennu.ip_cum == pytest.approx(5.717e-4)
    assert bennu.ps_cum == -1.40
    assert by_key["2011TO"].first_year == 2064 and by_key["2011TO"].last_year == 2064
    assert by_key["2011TO"].diameter_m == pytest.approx(18.0)


def test_parse_esa(esa):
    entries, updated = esa
    assert updated == "2026-11-11 17:21"
    by_key = {e.key: e for e in entries}
    assert set(by_key) == {"2008JL3", "2023VD3", "2023TB2", "443104", "2003UQ25", "2024XC"}
    jl3 = by_key["2008JL3"]
    assert (jl3.ip_cum, jl3.ps_cum, jl3.ps_max, jl3.ts_max) == (1.84e-4, -2.10, -2.10, 0)
    assert (jl3.first_year, jl3.last_year, jl3.diameter_m) == (2027, 2126, 30.0)
    assert by_key["2023TB2"].ts_max is None  # 'n/a' in the list
    assert by_key["2024XC"].ts_max is None  # 'n/d' in the list
    assert by_key["443104"].designation == "443104 2013XK22"


def test_compare(jpl, esa):
    rows = {r.key: r for r in compare(jpl, esa[0])}
    assert rows["2008JL3"].listed_by == "both"
    assert rows["101955"].listed_by == "jpl"
    assert rows["2003UQ25"].listed_by == "esa"
    jl3 = rows["2008JL3"]
    assert jl3.log10_ip_ratio == pytest.approx(math.log10(jl3.jpl.ip_cum / 1.84e-4))
    stats = summary(list(rows.values()))
    assert (stats["jpl"], stats["esa"], stats["both"]) == (6, 6, 3)
    assert stats["jpl_only"] + stats["esa_only"] + stats["both"] == stats["objects"]


def test_run_writes_table_report_and_changes(tmp_path, jpl, esa):
    data, reports = tmp_path / "data", tmp_path / "reports"
    entries, updated = esa
    run("2026-11-11", jpl, entries, updated, str(data), str(reports))
    with open(data / "compare.csv", newline="") as f:
        table = list(csv.DictReader(f))
    assert [r["key"] for r in table] == sorted(r["key"] for r in table)
    assert (reports / "2026" / "2026-11-11.md").read_text().startswith("# NASA vs ESA impact-risk lists, 2026-11-11")

    # Next day ESA drops one object: the report and changes() must say so
    previous = read_previous(str(data / "compare.csv"))
    later = [e for e in entries if e.key != "2003UQ25"]
    assert changes(previous, compare(jpl, later)) == ["ESA NEOCC removed 2003UQ25"]
    run("2026-11-12", jpl, later, updated, str(data), str(reports))
    assert "- ESA NEOCC removed 2003UQ25" in (reports / "2026" / "2026-11-12.md").read_text()
