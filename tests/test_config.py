"""
The Task 1 acceptance criteria, as executable checks.

Task 1 of the remediation work order removed two hard-coded air properties that
produced *silently wrong answers rather than errors*: `air_density=1.225` as a
default argument on four solver signatures, and a sea-level
`AIR_KINEMATIC_VISCOSITY` module constant. An AEP call that omitted the density
came back 25 % high and still landed inside the 4-6 MWh/yr plausibility band;
the viscosity would have put every design-rotor Reynolds number 23 % high.

Its "done when" list is four items. Three of them are properties of the code
rather than of one run, so they are pinned here instead of being checked once
by hand and then quietly regressing:

  * no air-property literal survives anywhere under `src/`      -> test_no_air_property_literals_in_src
  * omitting rho or nu is a TypeError, not a wrong number       -> test_*_requires_air_properties
  * an unresolved TODO raises rather than being substituted     -> test_unresolved_*

The fourth -- that the Phase VI cross-check figures are unchanged after the
rewiring -- is what tests/test_golden_regression.py already covers, at a
tolerance, across the whole Cp(lambda) curve and the spanwise state. It is not
duplicated here.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""

import os
import re

import pytest

import config
from bem.powercurve import cp_lambda_curve, operating_point, power_curve
from bem.rotor import phase_vi_geometry, solve_rotor

_HERE = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.abspath(os.path.join(_HERE, "..", "src"))

#: The two literals Task 1 removed. Written as regexes rather than substrings
#: so `1.2255` or a Reynolds number ending in the same digits cannot match by
#: accident, and so both spellings of the exponent are caught.
AIR_PROPERTY_LITERALS = (
    r"(?<![\d.])1\.225(?![\d])",          # sea-level density, kg/m^3
    r"(?<![\d.])1\.5e-0?5(?![\d])",       # sea-level kinematic viscosity, m^2/s
)


def _python_files(root):
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d != "__pycache__"]
        for name in filenames:
            if name.endswith(".py"):
                yield os.path.join(dirpath, name)


def test_no_air_property_literals_in_src():
    """
    No sea-level air property appears as a literal anywhere under `src/`.

    This is the work order's own acceptance check ("grep is the check"), and it
    covers comments and docstrings as well as code deliberately: a comment
    carrying the number is how it gets copied back into code later. The values
    live in `config/`, which is input data, not source.
    """

    offenders = []
    for path in _python_files(SRC_DIR):
        with open(path, encoding="utf-8") as f:
            for lineno, line in enumerate(f, start=1):
                for pattern in AIR_PROPERTY_LITERALS:
                    if re.search(pattern, line):
                        rel = os.path.relpath(path, os.path.dirname(SRC_DIR))
                        offenders.append(f"{rel}:{lineno}: {line.rstrip()}")

    assert not offenders, (
        "hard-coded sea-level air properties found under src/. They belong in "
        "config/ and must be threaded through as required arguments:\n  "
        + "\n  ".join(offenders)
    )


@pytest.mark.parametrize("call", [
    pytest.param(lambda g: solve_rotor(g, tsr=6.0), id="solve_rotor"),
    pytest.param(lambda g: operating_point(g, v_inf=7.0, tsr=6.0), id="operating_point"),
    pytest.param(lambda g: power_curve(g, [7.0], tsr=6.0), id="power_curve"),
    pytest.param(lambda g: cp_lambda_curve(g, [6.0]), id="cp_lambda_curve"),
])
def test_solver_entry_points_require_air_properties(call):
    """
    Every solver entry point refuses to run without rho and nu.

    A `TypeError` here is the entire point of Task 1: the previous behaviour
    was a plausible-looking number computed at the wrong site.
    """

    geometry = phase_vi_geometry()
    with pytest.raises(TypeError, match="air_density|kinematic_viscosity"):
        call(geometry)


def test_solve_rotor_rejects_air_properties_positionally():
    """
    rho and nu are keyword-only, so they cannot be passed by position.

    Positional passing is how the two would get swapped -- and a density used
    as a viscosity is off by five orders of magnitude, which at least fails
    loudly; the reverse quietly does not.
    """

    geometry = phase_vi_geometry()
    with pytest.raises(TypeError):
        solve_rotor(geometry, 6.0, 7.0, 1.0, 1e-5)


def test_phase_vi_config_supplies_the_validation_condition():
    """The validation rotor's sea-level condition is config, not a default."""

    phase_vi = config.load_phase_vi_rotor()

    assert phase_vi.air_density == pytest.approx(1.225)
    assert phase_vi.kinematic_viscosity == pytest.approx(1.5e-5)
    assert phase_vi.rated_rpm == pytest.approx(71.63)
    assert phase_vi.tip_pitch_deg == pytest.approx(3.0)
    # mu is derived from the pair above rather than recorded, so it cannot
    # drift from them. compare_ccblade.py passes it to CCBlade as `mu=`.
    assert phase_vi.dynamic_viscosity == (
        phase_vi.air_density * phase_vi.kinematic_viscosity
    )


def test_site_atmosphere_is_internally_consistent():
    """
    The site's derived atmosphere agrees with the barometric relation.

    `load_site` recomputes p, rho and nu from elevation, mean temperature, mu
    and the ISA constants and rejects the file on disagreement, so simply
    loading it is the check. Re-asserted here against the plan's exit band so
    a change to the tolerance cannot quietly widen it.
    """

    site = config.load_site()

    assert 0.95 < site.air_density < 1.03, (
        "plan item 1.1's exit criterion: site density between 0.95 and 1.03 kg/m^3"
    )
    assert site.kinematic_viscosity == pytest.approx(
        site.dynamic_viscosity_pa_s / site.air_density, rel=config.ATMOSPHERE_TOLERANCE
    )
    # The site is markedly less dense than sea level -- the point of plan 1.2.
    assert site.air_density < 1.0


def test_unresolved_wind_resource_raises_rather_than_substituting():
    """
    The Global Wind Atlas fields are TODO and behave like it.

    Ground rule 3, and plan 1.3: a fabricated wind resource propagates silently
    into every AEP figure downstream. So these must not be numbers, must not be
    `None` (which fails late, far from the config, complaining about
    `NoneType`), and must not be quietly falsy.
    """

    site = config.load_site()

    for name in ("weibull_k", "weibull_c_ms", "mean_wind_speed_ms"):
        value = getattr(site, name)
        assert not config.is_resolved(value), f"{name} should still be TODO"

        with pytest.raises(config.UnresolvedConfigError, match=name):
            float(value)
        with pytest.raises(config.UnresolvedConfigError, match=name):
            2.0 * value
        with pytest.raises(config.UnresolvedConfigError, match=name):
            # `if cfg.weibull_k:` is a silent fallback waiting to happen.
            bool(value)


def test_unresolved_design_bounds_raise():
    """The design-variable bounds are TODO too (plan 7.1), and behave the same."""

    bounds = config.load_design_rotor().parameterisation

    for name in ("chord_min_m", "chord_max_m", "twist_min_deg", "twist_max_deg"):
        with pytest.raises(config.UnresolvedConfigError, match=name):
            float(getattr(bounds, name))


def test_resolved_design_values_are_present():
    """Everything that *is* known is loaded, so the TODOs are not load failures."""

    design = config.load_design_rotor()

    assert design.radius_m == pytest.approx(2.0)   # fixed, not a design variable
    assert design.n_blades == 3
    assert design.design_tsr == pytest.approx(6.5)
    assert design.cut_in_wind_speed_ms == pytest.approx(3.0)
    assert design.cut_out_wind_speed_ms == pytest.approx(20.0)
    # Strip count and design-variable count are independent quantities (plan
    # 4.1); conflating them is the plan's named most-likely misreading.
    assert design.parameterisation.n_bem_strips == 25
    assert design.parameterisation.n_design_variables == 10


def test_polar_cache_metadata_matches_the_committed_cache():
    """
    `polars_s809.yaml` describes the cache that is actually on disk.

    Descriptive metadata that has drifted from the files it describes is worse
    than none, so the Reynolds list is checked against the directory rather
    than taken on trust. Ncrit is per-cache by design -- SG6043 will differ --
    so there is no global default to compare against.
    """

    cache = config.load_polar_cache("s809")
    assert cache.ncrit == pytest.approx(5.0)
    assert cache.n_panel == 240

    cache_dir = os.path.join(os.path.dirname(SRC_DIR), *cache.directory.split("/"))
    on_disk = {
        int(m.group(1))
        for name in os.listdir(cache_dir)
        if (m := re.fullmatch(r"S809_Re(\d+)\.csv", name))
    }
    assert set(cache.reynolds) == on_disk


def test_config_dir_is_the_only_reader_of_the_yaml_files():
    """
    Nothing outside `src/config/` opens `config/`.

    This is what makes "nothing downstream hard-codes a site value" hold going
    forward rather than only today: a module that reads the YAML itself would
    reintroduce a second, unchecked path for a site value to enter the code.
    """

    offenders = []
    for path in _python_files(SRC_DIR):
        if os.path.dirname(path) == os.path.join(SRC_DIR, "config"):
            continue
        with open(path, encoding="utf-8") as f:
            source = f.read()
        if "yaml" in source.lower() and re.search(r"\byaml\.(safe_)?load\b", source):
            offenders.append(os.path.relpath(path, os.path.dirname(SRC_DIR)))

    assert not offenders, (
        "these modules parse YAML directly; config/ must be read only through "
        f"src/config: {offenders}"
    )
