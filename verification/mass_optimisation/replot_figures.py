"""
Redraw every figure of `verification/mass_optimisation` from the committed
JSON beside the scripts.

One command, and no optimisation, sweep or check: nothing here solves,
samples or evaluates a blade. Every number comes from the committed
artefacts; `_common.build_problem()` supplies the parameterisation the
figures place the curves with, and `_common.prepared_problem()` (config and
the polar cache only) the rendered blade surfaces. No JSON is written.

    python verification/mass_optimisation/replot_figures.py

Writes, next to this script:

    blade_delta0_geometry.png   result_delta0.json (x0 / x_c / x_m)
    blade_delta0_loads.png      result_delta0.json (x0 / x_c / x_m)
    reference_blades.png        reference_blades.json
    multistart_delta0.png       multistart_delta0.json
    pareto_front.png            pareto.json (+ verification/baseline/x0.json)
    pareto_ablation.png         ablation.json
    blades_rendered_iso.png     reference_blades.json
    blades_rendered_plan.png    reference_blades.json
    blades_rendered_edge.png    reference_blades.json

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import _common as C  # noqa: E402
import run_mass_checks  # noqa: E402
import run_multistart  # noqa: E402
import run_sweep  # noqa: E402

BLADE_KEYS = ("x0", "x_c", "x_m")


def main():
    written = []

    # One parameterisation serves the planform, the loads and the front.
    problem = C.build_problem()

    # The planform and the rated-point loads of the three committed blades.
    result = C.load_json(run_mass_checks.RESULT_PATH)
    written += C.plot_blade(problem,
                            [(key, result["blades"][key]) for key in BLADE_KEYS],
                            _HERE, stem="blade_delta0")

    # The reference-blades table.
    table = C.load_json(run_mass_checks.TABLE_PATH)["blades"]
    written.append(run_mass_checks.plot_reference_blades(
        table, run_mass_checks.FIGURE_PATH))

    # The rendered blades (cosmetic root cylinder and tip rounding; the
    # surfaces need the parameterisation and the root fraction only).
    render_problem, _u0 = C.prepared_problem()
    by_label = {b["label"]: b for b in table}
    written += C.render_blades(render_problem,
                               [(key, by_label[key]) for key in BLADE_KEYS],
                               _HERE)

    # The multi-start figure.
    multistart = C.load_json(run_multistart.OUT_PATH)
    written.append(run_multistart.plot_multistart(
        multistart["runs"], multistart["multistart_u_spread_inf"],
        run_multistart.FIGURE_PATH))

    # The front and the ablation.
    pareto = C.load_json(run_sweep.PARETO_PATH)
    ablation = C.load_json(run_sweep.ABLATION_PATH)
    written.append(run_sweep.plot_pareto_front(
        problem.parameterisation, C.load_x0(), pareto["front"], pareto["steps"],
        run_sweep.FRONT_FIGURE_PATH))
    written.append(run_sweep.plot_pareto_ablation(
        ablation["cases"], run_sweep.ABLATION_FIGURE_PATH))

    for path in written:
        print(f"wrote {path}")
    return written


if __name__ == "__main__":
    main()
