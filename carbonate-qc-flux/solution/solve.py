#!/usr/bin/env python3
"""Reference solution for insight/ocean-carbonate-qc-flux.

Approach
--------
Each bottle has four carbonate parameters (TA, DIC, pH, fCO2), so the system is over-determined.
Using TA and fCO2 as the base pair, compute the expected pH (at the pH measurement temperature)
and the expected DIC, and compare with the measured values. Then test each fault hypothesis
(undo the fault, recompute the residuals) and keep the hypothesis with the lowest penalised
chi-square; "no fault" wins unless a fault explains the residuals clearly better.

Flux: Weiss (1974) K0, Weiss and Price (1980) pH2O, Weiss (1974) fugacity coefficient,
Wanninkhof (2014) k and Schmidt number, EOS-80 density.
"""
import calendar
import json
import math
import os
import sys

import numpy as np
import pandas as pd
import PyCO2SYS as pyco2

APP = os.environ.get("APP_DIR", "/app")
OUT = os.environ.get("REPORT_PATH", os.path.join(APP, "report.json"))
DEBUG = bool(os.environ.get("SOLVE_DEBUG"))

# VERIFY against the installed PyCO2SYS: 10 = Lueker 2000, 1 = Dickson 1990 (bisulfate),
# 2 = Lee 2010 (total boron), 1 = Dickson and Riley 1979 (fluoride), pH scale 1 = total.
OPTS = dict(opt_pH_scale=1, opt_k_carbonic=10, opt_k_bisulfate=1,
            opt_total_borate=2, opt_k_fluoride=1)

# Residual scatter of clean data. Set from the generator's [audit] output.
SIG_PH = 0.0017
SIG_DIC = 1.7
FAULT_PENALTY = 6.0      # extra chi-square cost for claiming a fault (prior against faults)

LABELS = ["ph_seawater_scale", "ph_temp_mislabeled", "dic_ta_swapped",
          "dic_volumetric_units", "ta_volumetric_units"]


def rho_sw(S, t):
    """EOS-80 density, kg m-3, 1 atm."""
    rw = (999.842594 + 6.793952e-2 * t - 9.095290e-3 * t ** 2 + 1.001685e-4 * t ** 3
          - 1.120083e-6 * t ** 4 + 6.536332e-9 * t ** 5)
    a = 0.824493 - 4.0899e-3 * t + 7.6438e-5 * t ** 2 - 8.2467e-7 * t ** 3 + 5.3875e-9 * t ** 4
    b = -5.72466e-3 + 1.0227e-4 * t - 1.6546e-6 * t ** 2
    return rw + a * S + b * S ** 1.5 + 4.8314e-4 * S ** 2


# ----------------------------------------------------------------------------- bottle QC
def residuals(ta, dic, ph, ph_t, fco2, S, T, si, po4):
    r = pyco2.sys(par1=ta, par2=fco2, par1_type=1, par2_type=5, salinity=S, temperature=T,
                  pressure=0.0, temperature_out=ph_t, pressure_out=0.0,
                  total_silicate=si, total_phosphate=po4, **OPTS)
    return (ph - r["pH_out"]) / SIG_PH, (dic - r["dic"]) / SIG_DIC


def sws_to_total(ph, ta, S, t, si, po4):
    r = pyco2.sys(par1=ph, par2=ta, par1_type=3, par2_type=1, salinity=S, temperature=t,
                  pressure=0.0, total_silicate=si, total_phosphate=po4,
                  **{**OPTS, "opt_pH_scale": 2})
    return np.asarray(r["pH_total"], float)


def diagnose(b):
    S, T = b["salinity"].values, b["temp_c"].values
    si, po4 = b["silicate_umolkg"].values, b["phosphate_umolkg"].values
    ta0, dic0 = b["ta_umolkg"].values, b["dic_umolkg"].values
    ph0, pht0, f = b["ph_meas"].values, b["ph_temp_c"].values, b["fco2_uatm"].values
    rho25 = rho_sw(S, 25.0) / 1000.0

    hyp = {
        "none": dict(ta=ta0, dic=dic0, ph=ph0, pht=pht0),
        "ph_seawater_scale": dict(ta=ta0, dic=dic0, ph=sws_to_total(ph0, ta0, S, pht0, si, po4), pht=pht0),
        "ph_temp_mislabeled": dict(ta=ta0, dic=dic0, ph=ph0, pht=np.full_like(pht0, 25.0)),
        "dic_ta_swapped": dict(ta=dic0, dic=ta0, ph=ph0, pht=pht0),
        "dic_volumetric_units": dict(ta=ta0, dic=dic0 / rho25, ph=ph0, pht=pht0),
        "ta_volumetric_units": dict(ta=ta0 / rho25, dic=dic0, ph=ph0, pht=pht0),
    }
    names = list(hyp)
    scores = []
    for name in names:
        h = hyp[name]
        z1, z2 = residuals(h["ta"], h["dic"], h["ph"], h["pht"], f, S, T, si, po4)
        s = z1 ** 2 + z2 ** 2 + (0.0 if name == "none" else FAULT_PENALTY)
        scores.append(s)
    scores = np.vstack(scores)                       # (n_hyp, n_samples)
    best = scores.argmin(axis=0)

    flagged, diagnosis, corrected = [], {}, {}
    for j, sid in enumerate(b["sample_id"].values):
        name = names[best[j]]
        if DEBUG:
            print(sid, name, np.round(scores[:, j], 1).tolist(), file=sys.stderr)
        if name == "none":
            continue
        flagged.append(str(sid))
        diagnosis[str(sid)] = name
        h = hyp[name]
        if name == "ph_seawater_scale":
            corrected[str(sid)] = {"ph_meas": round(float(h["ph"][j]), 4)}
        elif name == "ph_temp_mislabeled":
            corrected[str(sid)] = {"ph_temp_c": 25.0}
        elif name == "dic_ta_swapped":
            corrected[str(sid)] = {"dic_umolkg": round(float(ta0[j]), 2), "ta_umolkg": round(float(dic0[j]), 2)}
        elif name == "dic_volumetric_units":
            corrected[str(sid)] = {"dic_umolkg": round(float(h["dic"][j]), 2)}
        else:
            corrected[str(sid)] = {"ta_umolkg": round(float(h["ta"][j]), 2)}
    return flagged, diagnosis, corrected


# ----------------------------------------------------------------------------- flux
def k0_weiss(S, tc):
    t100 = (tc + 273.15) / 100.0
    return math.exp(-60.2409 + 93.4517 / t100 + 23.3585 * math.log(t100)
                    + S * (0.023517 - 0.023656 * t100 + 0.0047036 * t100 ** 2))


def flux_row(r):
    tc, S = r.sst_c, r.sss
    T = tc + 273.15
    p_atm = r.slp_hpa / 1013.25
    pw = math.exp(24.4543 - 67.4509 * (100.0 / T) - 4.8489 * math.log(T / 100.0) - 0.000544 * S)
    bvir = -1636.75 + 12.0408 * T - 0.0327957 * T ** 2 + 3.16528e-5 * T ** 3
    delta = 57.7 - 0.118 * T
    phi = math.exp((bvir + 2.0 * delta) * p_atm / (82.0578 * T))
    f_atm = r.xco2_atm_ppm * 1e-6 * (p_atm - pw) * phi
    sc = 2116.8 - 136.25 * tc + 4.7353 * tc ** 2 - 0.092307 * tc ** 3 + 0.0007555 * tc ** 4
    k = 0.251 * r.wind_10m_ms ** 2 * (sc / 660.0) ** -0.5      # cm h-1
    k0 = k0_weiss(S, tc) * rho_sw(S, tc)                         # mol m-3 atm-1
    return (k / 100.0) * k0 * (r.fco2_sw_uatm * 1e-6 - f_atm)   # mol m-2 h-1


def monthly_flux(u):
    ts = pd.to_datetime(u["timestamp"], utc=True)
    u = u.assign(year=ts.dt.year.values, month=ts.dt.month.values)
    u["flux_h"] = [flux_row(r) for r in u.itertuples()]
    multi_year = u["year"].nunique() > 1
    out = {}
    for (y, m), g in u.groupby(["year", "month"]):
        hours = 24 * calendar.monthrange(int(y), int(m))[1]
        key = f"{int(y):04d}-{int(m):02d}" if multi_year else str(int(m))
        out[key] = float(g["flux_h"].mean() * hours)
    return out


def main():
    bottle = pd.read_csv(os.path.join(APP, "bottle.csv"))
    underway = pd.read_csv(os.path.join(APP, "underway.csv"))
    flagged, diagnosis, corrected = diagnose(bottle)
    report = {"flagged": flagged, "diagnosis": diagnosis, "corrected": corrected,
              "monthly_flux": monthly_flux(underway)}
    with open(OUT, "w") as fh:
        json.dump(report, fh, indent=2)
    print(f"flagged {len(flagged)} samples; wrote {OUT}")


if __name__ == "__main__":
    main()
