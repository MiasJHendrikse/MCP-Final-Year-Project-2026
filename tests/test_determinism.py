"""
Bitwise reproducibility of the solver (work order Task 7).

The brief's determinism requirement is stated on AEP: the same design vector
must give a bitwise identical answer, including "via a different code path
(fresh process)". `src/objective/` does not exist yet -- it lands with plan
step 1.5 -- so the AEP-level test lands with it.

What is testable now is the layer everything above it inherits: `solve_rotor`,
and specifically that its answer does not depend on call order, on what ran
before it in the process, or on process-local state. That is the property
Task 4 bought by deleting `xfoil.polar_lookup`'s `ACTIVE_AIRFOIL` global and
`bem.airfoil.S809Polar`'s module-level lookup, and it is worth pinning at this
level rather than waiting: if it fails here, the AEP test would fail for a
reason that has nothing to do with AEP.

Both halves of the brief's requirement are covered:

  * same process, repeated and interleaved with other work
  * fresh subprocess, compared byte-for-byte against this one

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import json
import math
import os
import subprocess
import sys

import pytest

from bem.rotor import PHASE_VI_RATED_RPM, phase_vi_geometry, solve_rotor
from config import load_phase_vi_rotor
from polars.polar import interpolant_for

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(_HERE, ".."))

_AIR = load_phase_vi_rotor()
AIR = {"air_density": _AIR.air_density,
       "kinematic_viscosity": _AIR.kinematic_viscosity}

#: The operating point every check below uses. Fixed so the subprocess and the
#: in-process run are comparing the same thing.
V_INF = 7.0
TSR = PHASE_VI_RATED_RPM * 2.0 * math.pi / 60.0 * 5.029 / V_INF

#: Per-station fields compared bitwise. Cp and Ct alone would not catch a
#: spanwise redistribution that leaves the integral unchanged.
FIELDS = ("phi", "a", "a_prime", "alpha", "Cl", "Cd", "Cn", "Ct", "F",
          "reynolds", "w")


def _signature(result):
    """A flat, JSON-round-trippable summary of everything worth comparing."""

    return {
        "Cp": result["Cp"],
        "Ct": result["Ct"],
        "stations": [[station[field] for field in FIELDS]
                     for station in result["stations"]],
    }


def _solve():
    return _signature(solve_rotor(phase_vi_geometry(), tsr=TSR, v_inf=V_INF, **AIR))


def test_repeated_calls_are_bitwise_identical():
    """The same call, five times, in one process."""

    first = _solve()
    for _ in range(4):
        assert _solve() == first


def test_result_does_not_depend_on_call_order():
    """
    Interleaving other work must not move the answer.

    This is the direct guard on what `ACTIVE_AIRFOIL` cost. Under a
    module-level "current airfoil", solving a different rotor in between --
    which from Phase 1.4 means a different *cache* -- could change what the
    next caller got. Both caches are touched between the two solves here.
    """

    first = _solve()

    interpolant_for("sg6043")
    interpolant_for("s809")
    solve_rotor(phase_vi_geometry(), tsr=4.0, v_inf=9.0, **AIR)
    interpolant_for("sg6043")

    assert _solve() == first


#: Run in a *fresh* interpreter, printing the signature as JSON. Written as a
#: string rather than a committed script because it must stay in lockstep with
#: the constants above -- a separate file could drift and the test would then
#: be comparing two different operating points and calling them identical.
_SUBPROCESS_SOURCE = """
import json, math, sys
sys.path.insert(0, {src!r})
from bem.rotor import phase_vi_geometry, solve_rotor
from config import load_phase_vi_rotor

air_config = load_phase_vi_rotor()
result = solve_rotor(
    phase_vi_geometry(), tsr={tsr!r}, v_inf={v_inf!r},
    air_density=air_config.air_density,
    kinematic_viscosity=air_config.kinematic_viscosity,
)
fields = {fields!r}
print(json.dumps({{
    "Cp": result["Cp"],
    "Ct": result["Ct"],
    "stations": [[s[f] for f in fields] for s in result["stations"]],
}}))
"""


def test_fresh_subprocess_reproduces_the_same_bits():
    """
    The brief's "same vector via a different code path (fresh process)".

    A fresh interpreter shares no memoised interpolant, no imported module
    state and no accumulated call history with this one, so agreement here
    rules out the whole class of process-local state that a repeated in-process
    call cannot.

    Compared through JSON with `repr`-exact floats: `json.dumps` uses
    `float.__repr__`, which round-trips IEEE doubles exactly, so this is a
    bitwise comparison and not a 17-significant-digit one.
    """

    source = _SUBPROCESS_SOURCE.format(
        src=os.path.join(REPO_ROOT, "src"), tsr=TSR, v_inf=V_INF, fields=FIELDS)

    completed = subprocess.run(
        [sys.executable, "-c", source],
        capture_output=True, text=True, cwd=REPO_ROOT, timeout=300,
    )
    assert completed.returncode == 0, completed.stderr

    assert json.loads(completed.stdout) == _solve()


def test_json_round_trip_is_exact():
    """
    The comparison above is only bitwise if the transport is lossless.

    Asserted rather than assumed, because if `json.dumps` ever rounded, the
    subprocess test would quietly weaken from "bitwise" to "close" and would
    still pass.
    """

    signature = _solve()
    assert json.loads(json.dumps(signature)) == signature
