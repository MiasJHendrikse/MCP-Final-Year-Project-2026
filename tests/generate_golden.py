"""
Regenerate the golden files in tests/golden/.

    python tests/generate_golden.py

Read this before running it
----------------------------
The golden files are the regression net under the whole remediation work
order. Regenerating them is not a routine step -- it is how a *deliberate*
behaviour change gets recorded, and the procedure is:

    A task is not complete while the golden regression is failing, unless the
    change is intended, in which case the golden file is regenerated in its
    OWN commit, separate from the code change, with the reason in the commit
    message.

So: commit the code change first, run the test, and only if it fails for a
reason you can state does the regeneration happen -- as a second commit whose
message says which task changed what and why. Never fold a regeneration into
a code commit; that is exactly the silent absorption the net exists to stop.

Float precision
----------------
json.dump writes Python's repr for floats, which round-trips IEEE-754 doubles
exactly. The golden files therefore carry full precision and diff cleanly.

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import json
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from golden_reference import GOLDEN_DIR, build_golden, golden_path  # noqa: E402


def main():
    os.makedirs(GOLDEN_DIR, exist_ok=True)

    print("Building golden snapshot...")
    t0 = time.perf_counter()
    snapshot = build_golden()
    print(f"Built in {time.perf_counter() - t0:.1f} s")

    for key, payload in snapshot.items():
        path = golden_path(key)
        with open(path, "w", newline="\n") as f:
            json.dump(payload, f, indent=2, sort_keys=True)
            f.write("\n")
        print(f"Wrote {os.path.relpath(path, os.path.join(GOLDEN_DIR, '..', '..'))}")

    cross = snapshot["cross_tool_summary"]["mean_abs_cp_deviation_pct"]
    print("\nCross-tool mean |Cp| deviation against the STORED external values:")
    print(f"  vs stored CCBlade : {cross['vs_stored_ccblade']:.2f} %")
    print(f"  vs stored pyBEMT  : {cross['vs_stored_pybemt']:.2f} %")
    print("  (tracking figures, not validation -- see tests/golden/README.md)")


if __name__ == "__main__":
    main()
