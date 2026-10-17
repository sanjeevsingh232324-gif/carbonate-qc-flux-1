"""Verifier for insight/ocean-carbonate-qc-flux. Ground truth lives only in this (separate) image."""
import json
import math
import os
import pathlib

import pytest

HERE = pathlib.Path(__file__).resolve().parent
REPORT = pathlib.Path(os.environ.get("REPORT_PATH", "/app/report.json"))
TRUTH = json.loads((HERE / "ground_truth.json").read_text())

LABELS = set(TRUTH["labels"])
FAULTS = TRUTH["faults"]
TOL = TRUTH["tolerances"]
TOP_KEYS = {"flagged", "diagnosis", "corrected", "monthly_flux"}


def _finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


@pytest.fixture(scope="module")
def report():
    assert REPORT.is_file(), f"{REPORT} was not created"
    try:
        data = json.loads(REPORT.read_text())
    except json.JSONDecodeError as exc:
        pytest.fail(f"report.json is not valid JSON: {exc}")
    assert isinstance(data, dict), "report.json must be a JSON object"
    return data


def test_top_level_keys(report):
    assert set(report) == TOP_KEYS, f"top-level keys must be exactly {sorted(TOP_KEYS)}, got {sorted(report)}"


def test_types_and_internal_consistency(report):
    flagged, diagnosis, corrected = report["flagged"], report["diagnosis"], report["corrected"]
    assert isinstance(flagged, list) and all(isinstance(s, str) for s in flagged), "flagged must be a list of strings"
    assert len(set(flagged)) == len(flagged), "flagged contains duplicates"
    assert isinstance(diagnosis, dict) and isinstance(corrected, dict)
    assert set(diagnosis) == set(flagged), "diagnosis keys must equal flagged"
    assert set(corrected) == set(flagged), "corrected keys must equal flagged"
    for sid, lab in diagnosis.items():
        assert lab in LABELS, f"{sid}: label {lab!r} not in the allowed list"
    for sid, obj in corrected.items():
        assert isinstance(obj, dict) and obj, f"{sid}: corrected entry must be a non-empty object"
        for var, val in obj.items():
            assert _finite(val), f"{sid}.{var}: value must be a finite number"


def test_flagged_set_exact(report):
    got, want = set(report["flagged"]), set(FAULTS)
    assert got == want, f"missing: {sorted(want - got)}; false positives: {sorted(got - want)}"


def test_diagnosis_labels(report):
    wrong = {sid: (report["diagnosis"].get(sid), f["label"]) for sid, f in FAULTS.items()
             if report["diagnosis"].get(sid) != f["label"]}
    assert not wrong, f"wrong diagnoses (got, expected): {wrong}"


def test_corrected_values(report):
    problems = []
    for sid, f in FAULTS.items():
        got = report["corrected"].get(sid, {})
        want = f["corrected"]
        if set(got) != set(want):
            problems.append(f"{sid}: variables {sorted(got)} != {sorted(want)}")
            continue
        for var, w in want.items():
            if abs(got[var] - w) > f["tolerance"]:
                problems.append(f"{sid}.{var}: {got[var]} vs {w} (tol {f['tolerance']})")
    assert not problems, "; ".join(problems)


def test_monthly_flux(report):
    flux, want = report["monthly_flux"], TRUTH["monthly_flux"]
    assert isinstance(flux, dict)
    assert set(flux) == set(want), f"month keys {sorted(flux)} != {sorted(want)}"
    problems = []
    for key, w in want.items():
        g = flux[key]
        if not _finite(g):
            problems.append(f"{key}: not a finite number")
            continue
        if abs(g - w) > max(TOL["flux_rel"] * abs(w), TOL["flux_abs"]):
            problems.append(f"{key}: {g:.5f} vs {w:.5f}")
    assert not problems, "; ".join(problems)
