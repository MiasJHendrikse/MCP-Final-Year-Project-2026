"""
Generate the committed residual histories for the Phase VI rotor.

Work order Task 5: "Solver converges across the full operating envelope for
baseline **and** perturbed geometries; residual histories committed." This is
the script that produces them, kept alongside its output as a versioned report
artefact rather than in `results/` (which is generated and partly gitignored).

Run from the repo root:

    python verification/phase_vi/generate_residual_histories.py

Redraw the committed figure from the committed JSON, without recomputing
anything:

    python verification/phase_vi/generate_residual_histories.py --replot

Writes `residual_histories.json` and `residual_convergence.png` next to this
file. The JSON holds each station's residual history -- the root-finder's own
sequence of `[phi, R]` evaluations -- and the initial residual each history is
normalised by, but not the convergence criterion: that is the solver's own
`RESIDUAL_RTOL`, imported from `bem.station` and labelled on the criterion
line.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import argparse
import json
import math
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "src"))

from bem.rotor import PHASE_VI_RATED_RPM, phase_vi_geometry  # noqa: E402
from bem.station import RESIDUAL_RTOL, StationParams, solve_station  # noqa: E402
from config import load_phase_vi_rotor  # noqa: E402
from polars.polar import CachedPolar, interpolant_for  # noqa: E402

PHASE_VI = load_phase_vi_rotor()
NU = PHASE_VI.kinematic_viscosity
OMEGA_RATED = PHASE_VI_RATED_RPM * 2.0 * math.pi / 60.0

#: Baseline operating points plus perturbed geometries -- the work order asks
#: for both, and the perturbed cases are the ones that matter, since an
#: optimiser line search visits geometries no fixed-rotor sweep ever does.
CASES = [
    ("v=5 m/s, rated rpm", 5.0, OMEGA_RATED, 1.0, 0.0),
    ("v=7 m/s, rated rpm", 7.0, OMEGA_RATED, 1.0, 0.0),
    ("v=10 m/s, rated rpm", 10.0, OMEGA_RATED, 1.0, 0.0),
    ("v=15 m/s, rated rpm", 15.0, OMEGA_RATED, 1.0, 0.0),
    ("v=20 m/s, rated rpm", 20.0, OMEGA_RATED, 1.0, 0.0),
    ("tsr=2.0, v=7", 7.0, 2.0 * 7.0 / 5.029, 1.0, 0.0),
    ("tsr=4.0, v=7", 7.0, 4.0 * 7.0 / 5.029, 1.0, 0.0),
    ("tsr=6.0, v=7", 7.0, 6.0 * 7.0 / 5.029, 1.0, 0.0),
    ("tsr=7.5, v=7", 7.0, 7.5 * 7.0 / 5.029, 1.0, 0.0),
    ("chord x0.7", 7.0, OMEGA_RATED, 0.7, 0.0),
    ("chord x1.3", 7.0, OMEGA_RATED, 1.3, 0.0),
    ("twist -6 deg", 7.0, OMEGA_RATED, 1.0, math.radians(-6.0)),
    ("twist +6 deg", 7.0, OMEGA_RATED, 1.0, math.radians(6.0)),
]


def stations(v_inf, omega, chord_scale, twist_delta):
    geometry = phase_vi_geometry()
    interpolant = interpolant_for(geometry.polar_cache)
    out = []
    for r, chord, twist in zip(geometry.r, geometry.chord, geometry.twist):
        chord = chord * chord_scale
        reynolds = math.hypot(v_inf, omega * r) * chord / NU
        out.append(StationParams(
            r=r, chord=chord, twist=twist + twist_delta,
            airfoil=CachedPolar(interpolant, reynolds),
            tsr=omega * r / v_inf, R=geometry.R,
            n_blades=geometry.n_blades, r_hub=geometry.r_hub,
        ))
    return out


def build():
    cases = []
    for name, v_inf, omega, chord_scale, twist_delta in CASES:
        entries = []
        for index, station in enumerate(stations(v_inf, omega, chord_scale, twist_delta)):
            result = solve_station(station, record_history=True)
            entries.append({
                "station": index,
                "r": station.r,
                "lambda_r": station.tsr,
                "converged": result["converged"],
                "failure": result["failure"],
                "phi": result["phi"],
                "iterations": result["iterations"],
                "function_calls": result["function_calls"],
                "residual_initial": result["residual_initial"],
                "residual_final": result["residual"],
                "residual_relative": (result["residual"] / result["residual_initial"]
                                      if result["residual_initial"] else None),
                # The root-finder's own trail: every (phi, R) it evaluated, in
                # order. This is what "residual histories" means -- not a
                # summary, the actual sequence, so a reader can see the
                # reduction rather than take the final number on trust.
                "history": [[phi, value] for phi, value in result["history"]],
            })
        cases.append({
            "name": name, "v_inf": v_inf, "omega": omega,
            "chord_scale": chord_scale, "twist_delta_rad": twist_delta,
            "stations": entries,
        })
    return cases


def summarise(cases):
    total = sum(len(c["stations"]) for c in cases)
    failed = [(c["name"], s["station"]) for c in cases for s in c["stations"]
              if not s["converged"]]
    relatives = [s["residual_relative"] for c in cases for s in c["stations"]
                 if s["residual_relative"] is not None]
    iterations = [s["iterations"] for c in cases for s in c["stations"]]
    return {
        "cases": len(cases),
        "stations_total": total,
        "stations_failed": len(failed),
        "failures": failed,
        "worst_relative_residual": max(relatives) if relatives else None,
        "max_iterations": max(iterations),
        "max_function_calls": max(s["function_calls"] for c in cases
                                  for s in c["stations"]),
    }


def require_histories(data):
    """
    Return the cases, or refuse the replot naming exactly what is missing.

    The figure needs each station's `history` and the `residual_initial` it
    is normalised by. The convergence criterion is not a JSON field: it is
    the solver's `RESIDUAL_RTOL`, imported from `bem.station`.
    """

    cases = data.get("cases")
    missing = []
    if not cases:
        missing.append("the top-level 'cases' array")
    else:
        for case in cases:
            for station in case.get("stations", []):
                where = (f"case {case.get('name', '?')!r}, "
                         f"station {station.get('station', '?')}")
                if not station.get("history"):
                    missing.append(f"{where}: 'history'")
                if not station.get("residual_initial"):
                    missing.append(f"{where}: 'residual_initial'")
    if missing:
        extra = (f" (and {len(missing) - 5} more)" if len(missing) > 5 else "")
        raise SystemExit("residual_histories.json cannot support the figure; "
                         "missing " + "; ".join(missing[:5]) + extra)
    return cases


def plot_histories(data, out_dir):
    """
    Draw `residual_convergence.png` from the per-station histories.

    Each history is divided by its own initial residual, so the figure shows
    the relative residual the solver drives to the criterion.
    """

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.ticker
        from plotting import figstyle
    except ImportError:
        print("matplotlib not available; skipping the figure")
        return

    figstyle.apply()
    fig, ax = plt.subplots(figsize=figstyle.SINGLE)

    for case in require_histories(data):
        for station in case["stations"]:
            initial = station["residual_initial"]
            best = math.inf
            running = []
            for _phi, value in station["history"]:
                best = min(best, abs(value) / initial)
                running.append(max(best, 1e-18))
            ax.semilogy(range(len(running)), running, linewidth=0.6,
                        alpha=0.35, color=figstyle.PALETTE[0])

    # `bem.station.RESIDUAL_RTOL`: the relative-residual criterion itself,
    # not a stored copy that could drift from the solver.
    ax.axhline(RESIDUAL_RTOL, color="0.3", linestyle="--", linewidth=1.0,
               label=r"$10^{-9} R_0$")
    ax.set_xlabel(figstyle.LABELS["evaluations"])
    # A count of evaluations: whole-number ticks only.
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    ax.set_ylabel(figstyle.LABELS["residual"])
    ax.legend()
    path = os.path.join(out_dir, "residual_convergence.png")
    figstyle.save(fig, path)
    plt.close(fig)
    print(f"wrote {path}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replot", action="store_true",
                        help="redraw residual_convergence.png from the "
                             "committed JSON, without recomputing anything")
    args = parser.parse_args(argv)

    if args.replot:
        with open(os.path.join(_HERE, "residual_histories.json"),
                  encoding="utf-8") as handle:
            data = json.load(handle)
        plot_histories(data, _HERE)
        return

    cases = build()
    summary = summarise(cases)

    out = os.path.join(_HERE, "residual_histories.json")
    with open(out, "w", encoding="utf-8") as handle:
        json.dump({"summary": summary, "cases": cases}, handle, indent=1)
    print(f"wrote {out}")

    plot_histories({"cases": cases}, _HERE)

    print()
    for key, value in summary.items():
        print(f"  {key}: {value}")


if __name__ == "__main__":
    main()
