"""
Why the optimiser gains +0.217 % and not 2-6 %: the diagnostics behind
docs/AEP_GAIN_AUDIT.md (2026-09-13).

Evaluation only -- no optimisation here (see reoptimise.py). Everything runs
through the repo's own objective chain (BladeParameterisation -> solve_rotor
-> aerodynamic_power) with the site resource and the 17 midpoint bins, so the
numbers are the objective's numbers, not a re-implementation of it.

Sections (each prints a table and lands in diagnostics.json):

  E1  bin-by-bin energy and gain attribution, x0 vs x* (adjoint optimum)
  E2  x* re-evaluated with the generator rating frozen at x0's 3822 W
  E8  Cp-lambda of x0 and x* at three wind speeds
  E4  x0 with the operating TSR swept (the MPPT-gain lever alone)
  E3  Schmitz designed at lambda_d and operated at lambda_d, lambda_d = 5..9
  E6  max-L/D point vs Re, and the energy-weighted 75 %-span Re
  S1  rotor-speed ceiling (Region 2.5), rating fixed: x0 and x* re-evaluated
  S2  same with the rating floating (the current model's P_aero(11 m/s))
  S3  what the literature's baselines score under this objective
  S4  manufacturability clipping of the Schmitz root
  SP  spanwise loading of x0 vs x* at 8.5 m/s
  X   the reoptimise.py optima cross-evaluated under every strategy

Under provisional bounds (chord_max_m = 0.45 m provisional). Run from the
repo root:

    python verification/aep_gain_audit/diagnostics.py      # ~3 min

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import glob
import json
import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "src"))
sys.path.insert(0, os.path.join(REPO, "tests"))

from bem.rotor import solve_rotor  # noqa: E402
from config import load_design_rotor, load_site  # noqa: E402
from design import BladeParameterisation, DesignBounds  # noqa: E402
from design.schmitz import (  # noqa: E402
    DEFAULT_DESIGN_REYNOLDS, max_lift_to_drag_point, schmitz_distribution)
from objective import WeibullResource  # noqa: E402
from objective.objective import HOURS_PER_YEAR  # noqa: E402
from objective.power import aerodynamic_power, wind_speed_bins  # noqa: E402
from polars.interpolant import PolarDomainError  # noqa: E402
from polars.polar import interpolant_for  # noqa: E402
from test_parameterisation import PROVISIONAL_BOUNDS  # noqa: E402

design = load_design_rotor()
site = load_site()
AIR = (site.air_density, site.kinematic_viscosity)
RES = WeibullResource.from_config()
PAR = BladeParameterisation()
R, LAM, VR = design.radius_m, design.design_tsr, design.rated_wind_speed_ms
EDGES, MIDS, _ = wind_speed_bins()
MASS = RES.probability_between(EDGES[:-1], EDGES[1:])
OUT = {}


def load_x0():
    path = os.path.join(REPO, "verification", "baseline", "x0.json")
    with open(path, encoding="utf-8") as handle:
        d = json.load(handle)
    return np.array(d["chord_control_points_m"] + d["twist_control_points_rad"], dtype=float)


def load_x_star():
    path = os.path.join(REPO, "verification", "adjoint_optimisation", "result.json")
    with open(path, encoding="utf-8") as handle:
        return np.array(json.load(handle)["x_star"], dtype=float)


def aep_sched(x, vtip_max=None, p_gen=None, lam=LAM):
    """
    AEP [MWh/yr] under: MPPT at `lam`; constant rpm above a tip-speed ceiling
    (Region 2.5); cap at `p_gen` (fixed rating) or at P_aero(V_rated) under
    the schedule (floating, the project's current model).
    """

    g = PAR.to_geometry(x)

    def tsr(v):
        return lam if vtip_max is None else min(lam, vtip_max / v)

    p_rated = p_gen if p_gen is not None else aerodynamic_power(g, VR, tsr(VR), *AIR)[0]
    rows, total = [], 0.0
    for v, m in zip(MIDS, MASS):
        p, r = aerodynamic_power(g, float(v), tsr(float(v)), *AIR)
        if not r["converged"]:
            raise RuntimeError(f"non-converged station at V = {v}")
        re = [s["reynolds"] for s in r["stations"]]
        rows.append({"v": float(v), "tsr": tsr(float(v)), "Cp": r["Cp"], "Ct": r["Ct"],
                     "P_w": p, "P_limited_w": min(p, p_rated), "capped": p > p_rated,
                     "re_min": min(re), "re_max": max(re)})
        total += min(p, p_rated) * m
    return total * HOURS_PER_YEAR / 1e6, p_rated, rows


def pct(a, b):
    return 100.0 * (a / b - 1.0)


def cp(geometry, tsr, v):
    return solve_rotor(geometry, tsr=tsr, v_inf=v, air_density=AIR[0],
                       kinematic_viscosity=AIR[1])["Cp"]


def main():
    x0, xs = load_x0(), load_x_star()
    A0, P0, rows0 = aep_sched(x0)
    As, Ps, rows_s = aep_sched(xs)
    OUT["site"] = {"k": RES.k, "c_ms": RES.c, "p_below_cut_in": RES.cdf(3.0),
                   "p_above_rated": 1 - RES.cdf(VR), "p_above_cut_out": 1 - RES.cdf(20.0)}
    print(f"site k={RES.k:.4f} c={RES.c:.4f}  P(V<3)={RES.cdf(3.0):.4f}  "
          f"P(V>{VR})={1 - RES.cdf(VR):.4f}  P(V>20)={1 - RES.cdf(20.0):.4f}")

    # -- E1 ------------------------------------------------------------------
    print("\n=== E1: bin-by-bin, x0 vs x* (fixed TSR 6.5, floating cap) ===")
    print(f"AEP x0={A0:.6f}  x*={As:.6f}  gain={pct(As, A0):+.4f} %   "
          f"P_rated {P0:.1f} -> {Ps:.1f} W ({pct(Ps, P0):+.3f} %)")
    print(f"{'V':>5} {'mass':>7} {'Cp0':>7} {'dCp%':>7} {'cap':>4} {'E0 MWh':>8} "
          f"{'share%':>7} {'dE/dAEP%':>9} {'Re0 range':>14}")
    below = above = d_below = d_above = 0.0
    bins = []
    for a, b, m in zip(rows0, rows_s, MASS):
        e0 = a["P_limited_w"] * m * HOURS_PER_YEAR / 1e6
        es = b["P_limited_w"] * m * HOURS_PER_YEAR / 1e6
        if a["capped"]:
            above += e0
            d_above += es - e0
        else:
            below += e0
            d_below += es - e0
        bins.append({"v": a["v"], "mass": float(m), "Cp_x0": a["Cp"], "Cp_xstar": b["Cp"],
                     "dCp_pct": pct(b["Cp"], a["Cp"]), "capped": a["capped"],
                     "energy_x0_mwh": e0, "share_pct": 100 * e0 / A0,
                     "share_of_gain_pct": 100 * (es - e0) / (As - A0),
                     "re_min_x0": a["re_min"], "re_max_x0": a["re_max"]})
        print(f"{a['v']:5.1f} {m:7.4f} {a['Cp']:7.4f} {pct(b['Cp'], a['Cp']):+7.3f} "
              f"{'cap' if a['capped'] else '':>4} {e0:8.4f} {100 * e0 / A0:7.2f} "
              f"{100 * (es - e0) / (As - A0):+9.2f} {a['re_min'] / 1e3:5.0f}k-{a['re_max'] / 1e3:4.0f}k")
    print(f"energy share: below rated {100 * below / A0:.1f} %, capped {100 * above / A0:.1f} %")
    print(f"gain from below-rated bins {100 * d_below / A0:+.4f} %, from capped bins {100 * d_above / A0:+.4f} %")
    OUT["E1"] = {"aep_x0": A0, "aep_xstar": As, "gain_pct": pct(As, A0), "p_rated_x0_w": P0,
                 "p_rated_xstar_w": Ps, "share_below_rated_pct": 100 * below / A0,
                 "share_capped_pct": 100 * above / A0, "gain_below_rated_pct": 100 * d_below / A0,
                 "gain_capped_pct": 100 * d_above / A0, "bins": bins}

    # -- E2 ------------------------------------------------------------------
    Asf = aep_sched(xs, p_gen=P0)[0]
    print(f"\n=== E2: x* with the rating frozen at {P0:.0f} W: AEP={Asf:.6f}, gain {pct(Asf, A0):+.4f} %")
    OUT["E2"] = {"aep_xstar_fixed_rating": Asf, "gain_pct": pct(Asf, A0)}

    # -- E8 ------------------------------------------------------------------
    print("\n=== E8: Cp-lambda, x0 / x* ===")
    g0, gs = PAR.to_geometry(x0), PAR.to_geometry(xs)
    OUT["E8"] = {}
    for v in (11.0, 6.5, 4.5):
        line = f" V={v:4.1f}:"
        OUT["E8"][str(v)] = {}
        for tsr in (5.5, 6.0, 6.25, 6.5, 6.75, 7.0, 7.5, 8.0):
            c0, cs = cp(g0, tsr, v), cp(gs, tsr, v)
            OUT["E8"][str(v)][str(tsr)] = [c0, cs]
            line += f"  {tsr}: {c0:.4f}/{cs:.4f}"
        print(line)

    # -- E4 ------------------------------------------------------------------
    print("\n=== E4: x0, operating TSR swept ===")
    OUT["E4"] = {}
    for lam in (5.5, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 8.0):
        a = aep_sched(x0, lam=lam)[0]
        OUT["E4"][str(lam)] = {"aep": a, "vs_x0_pct": pct(a, A0)}
        print(f" lambda_op={lam:5.2f}  AEP={a:.5f} ({pct(a, A0):+.3f} %)  tip speed at rated {lam * VR:.1f} m/s")

    # -- E3 ------------------------------------------------------------------
    interp = interpolant_for("sg6043")
    alpha_d, cl_d, ld = max_lift_to_drag_point(interp, DEFAULT_DESIGN_REYNOLDS)
    rR = PAR.radii / R
    print(f"\n=== E3: Schmitz at lambda_d, operated at lambda_d "
          f"(alpha={math.degrees(alpha_d):.2f} deg, Cl={cl_d:.4f}, L/D={ld:.1f}) ===")
    OUT["E3"] = {}
    for lam in (5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0):
        ch, tw = schmitz_distribution(rR, cl_d, alpha_d, tip_speed_ratio=lam)
        x, _, _ = PAR.fit(ch, tw)
        a, p_rated, rows = aep_sched(x, lam=lam)
        r105 = [r for r in rows if abs(r["v"] - 10.5) < 1e-6][0]
        re_min = min(r["re_min"] for r in rows)
        OUT["E3"][str(lam)] = {"aep": a, "vs_x0_pct": pct(a, A0), "root_chord_m": float(ch[0]),
                               "tip_chord_m": float(ch[-1]), "Cp_10p5": r105["Cp"],
                               "Ct_10p5": r105["Ct"], "p_rated_w": p_rated, "re_min": re_min}
        print(f" lambda_d={lam:3.1f}  root {ch[0] * 1000:4.0f} mm tip {ch[-1] * 1000:3.0f} mm  "
              f"Cp(10.5)={r105['Cp']:.4f} Ct={r105['Ct']:.3f}  Re_min={re_min / 1e3:3.0f}k  "
              f"P_rated={p_rated:.0f} W  AEP={a:.5f} ({pct(a, A0):+.3f} %)  tip {lam * VR:.0f} m/s")

    # -- E6 ------------------------------------------------------------------
    print("\n=== E6: max-L/D point vs Re ===")
    OUT["E6"] = {}
    for Re in (100e3, 150e3, 200e3, 250e3, 300e3, 400e3, 500e3):
        a_, c_, l_ = max_lift_to_drag_point(interp, Re)
        OUT["E6"][str(int(Re))] = {"alpha_deg": math.degrees(a_), "cl": c_, "ld": l_}
        print(f" Re={Re / 1e3:4.0f}k  alpha*={math.degrees(a_):.2f} deg  Cl={c_:.4f}  L/D={l_:.1f}")
    tot = wre = 0.0
    for row, m in zip(rows0, MASS):
        r = solve_rotor(g0, tsr=LAM, v_inf=row["v"], air_density=AIR[0], kinematic_viscosity=AIR[1])
        st = r["stations"]
        i = int(0.75 * len(st))
        e = row["P_limited_w"] * m
        tot += e
        wre += e * st[i]["reynolds"]
    OUT["E6"]["energy_weighted_re_75pct_span"] = wre / tot
    print(f" energy-weighted Re at ~75 % span (x0) = {wre / tot / 1e3:.0f}k")

    # -- S1 / S2 -------------------------------------------------------------
    print("\n=== S1: rotor-speed ceiling, rating FIXED at x0's P_rated ===")
    OUT["S1"] = {}
    for vt in (None, 75.0, 70.0, 65.0, 60.0, 55.0, 50.0):
        a0 = aep_sched(x0, vt, p_gen=P0)[0]
        a1 = aep_sched(xs, vt, p_gen=P0)[0]
        rpm = None if vt is None else vt / R * 60 / (2 * math.pi)
        lam_rated = LAM if vt is None else min(LAM, vt / VR)
        OUT["S1"][str(vt)] = {"aep_x0": a0, "vs_no_ceiling_pct": pct(a0, A0), "aep_xstar": a1,
                              "xstar_vs_x0_pct": pct(a1, a0), "omega_max_rpm": rpm,
                              "v_c_ms": None if vt is None else vt / LAM, "lambda_at_rated": lam_rated}
        print(f" Vtip={str(vt):>5}  Omega_max={'-' if rpm is None else f'{rpm:.0f}':>4} rpm  "
              f"AEP x0={a0:.5f} ({pct(a0, A0):+.3f} %)  x* re-evaluated: {pct(a1, a0):+.3f} % vs x0  "
              f"lambda@rated={lam_rated:.2f}")
    print("\n=== S2: rotor-speed ceiling, rating FLOATING ===")
    OUT["S2"] = {}
    for vt in (None, 65.0, 60.0, 55.0):
        a0, p0, _ = aep_sched(x0, vt)
        a1, _, _ = aep_sched(xs, vt)
        OUT["S2"][str(vt)] = {"aep_x0": a0, "vs_no_ceiling_pct": pct(a0, A0), "p_rated_w": p0,
                              "aep_xstar": a1, "xstar_vs_x0_pct": pct(a1, a0)}
        print(f" Vtip={str(vt):>5}  AEP x0={a0:.5f} ({pct(a0, A0):+.3f} %)  P_rated={p0:.0f} W  "
              f"x*: {pct(a1, a0):+.3f} % vs x0")

    # -- S3 ------------------------------------------------------------------
    print("\n=== S3: the literature's baselines under this objective ===")
    OUT["S3"] = {}

    def record(label, x):
        try:
            a = aep_sched(x)[0]
            OUT["S3"][label] = {"aep": a, "vs_x0_pct": pct(a, A0)}
            print(f" {label}: AEP={a:.5f} ({pct(a, A0):+.3f} %)")
        except PolarDomainError as error:
            OUT["S3"][label] = {"error": str(error)}
            print(f" {label}: PolarDomainError")

    lam_r = LAM * rR
    phi_b = np.arctan(2.0 / (3.0 * lam_r))
    c_b = 8 * math.pi * (rR * R) * np.sin(phi_b) / (3 * design.n_blades * cl_d * lam_r)
    record("Betz (no swirl), same polar point", PAR.fit(c_b, phi_b - alpha_d)[0])
    ch_s, tw_s = schmitz_distribution(rR, cl_d, alpha_d)
    record("linear chord+twist through Schmitz endpoints",
           PAR.fit(np.linspace(ch_s[0], ch_s[-1], len(rR)), np.linspace(tw_s[0], tw_s[-1], len(rR)))[0])
    record("generic linear 200->70 mm, 20->0 deg",
           PAR.fit(np.linspace(0.20, 0.07, len(rR)), np.radians(np.linspace(20.0, 0.0, len(rR))))[0])
    for ad, cld, lbl in ((8.0, 1.40, "Schmitz at alpha=8 deg, Cl=1.40 (datasheet-style)"),
                         (3.0, 1.10, "Schmitz at alpha=3 deg, Cl=1.10")):
        ch, tw = schmitz_distribution(rR, cld, math.radians(ad))
        record(lbl, PAR.fit(ch, tw)[0])
    for Re in (100e3, 400e3):
        a_, c_, _ = max_lift_to_drag_point(interp, Re)
        ch, tw = schmitz_distribution(rR, c_, a_)
        record(f"Schmitz at design Re={Re / 1e3:.0f}k (alpha={math.degrees(a_):.2f}, Cl={c_:.3f})",
               PAR.fit(ch, tw)[0])

    # -- S4 ------------------------------------------------------------------
    print("\n=== S4: manufacturability clipping of the Schmitz root ===")
    OUT["S4"] = {}
    for cmax in (0.25, 0.22, 0.20, 0.18, 0.15, 0.12):
        b = DesignBounds(n_chord=5, n_twist=5, chord_min_m=PROVISIONAL_BOUNDS["chord_min_m"],
                         chord_max_m=cmax, twist_min_rad=PROVISIONAL_BOUNDS["twist_min_rad"],
                         twist_max_rad=PROVISIONAL_BOUNDS["twist_max_rad"])
        xc, viol = b.clip_physical(x0)
        try:
            a = aep_sched(np.asarray(xc))[0]
            OUT["S4"][str(cmax)] = {"aep": a, "vs_x0_pct": pct(a, A0), "clipped_control_points": len(viol)}
            print(f" chord_max={cmax:.2f} m: {len(viol)} CP clipped, AEP={a:.5f} ({pct(a, A0):+.3f} %)")
        except PolarDomainError:
            OUT["S4"][str(cmax)] = {"error": "PolarDomainError"}
            print(f" chord_max={cmax:.2f} m: PolarDomainError")

    # -- SP ------------------------------------------------------------------
    print("\n=== SP: spanwise at 8.5 m/s, lambda 6.5 ===")
    OUT["SP"] = {}
    for lbl, g in (("x0", g0), ("x_star", gs)):
        r = solve_rotor(g, tsr=LAM, v_inf=8.5, air_density=AIR[0], kinematic_viscosity=AIR[1])
        OUT["SP"][lbl] = {"Cp": r["Cp"], "Ct": r["Ct"], "stations": [
            {"r_over_R": s["r"] / g.R, "chord_m": c, "twist_deg": math.degrees(t),
             "alpha_deg": math.degrees(s["alpha"]), "a": s["a"], "a_prime": s["a_prime"],
             "F": s["F"], "Cl": s["Cl"], "Cd": s["Cd"], "reynolds": s["reynolds"]}
            for s, c, t in zip(r["stations"], g.chord, g.twist)]}
        print(f" {lbl}: Cp={r['Cp']:.4f} Ct={r['Ct']:.4f}")

    # -- X: cross-evaluate reoptimise.py optima ---------------------------------
    print("\n=== X: reoptimise.py optima under every strategy (rating fixed at x0's) ===")
    OUT["X"] = {}
    ceilings = (None, 60.0, 55.0, 50.0)
    print(f"{'blade':>16} | " + " ".join(f"{str(c):>8}" for c in ceilings) + " | Ct(8.5, 6.5) | chords mm / twists deg")
    for path in [None] + sorted(glob.glob(os.path.join(_HERE, "opt_*.json"))):
        if path is None:
            name, x = "x0", x0
        else:
            name = os.path.basename(path)[4:-5]
            with open(path, encoding="utf-8") as handle:
                x = np.array(json.load(handle)["x_opt"], dtype=float)
        vals = [aep_sched(x, c, p_gen=P0)[0] for c in ceilings]
        g = PAR.to_geometry(x)
        ct = solve_rotor(g, tsr=LAM, v_inf=8.5, air_density=AIR[0], kinematic_viscosity=AIR[1])["Ct"]
        OUT["X"][name] = {"aep_by_ceiling": dict(zip(map(str, ceilings), vals)), "Ct_8p5_6p5": ct,
                          "cp_lambda_11": {str(t): cp(g, t, 11.0) for t in (4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0)},
                          "x": [float(v) for v in x]}
        print(f"{name:>16} | " + " ".join(f"{v:8.4f}" for v in vals) + f" | {ct:12.4f} | "
              + " ".join(f"{c * 1000:4.0f}" for c in x[:5]) + " / "
              + " ".join(f"{math.degrees(t):5.2f}" for t in x[5:]))

    with open(os.path.join(_HERE, "diagnostics.json"), "w", encoding="utf-8") as handle:
        json.dump(OUT, handle, indent=1,
                  default=lambda o: o.tolist() if isinstance(o, np.ndarray) else float(o))
    print("\nwrote diagnostics.json")


if __name__ == "__main__":
    main()
