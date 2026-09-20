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

import contextlib
import math
import os
import re
import shutil
import tempfile

import pytest
import yaml

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


# ---------------------------------------------------------------------------
# Driving the loader's rejection paths
# ---------------------------------------------------------------------------
# The loader rejects an internally inconsistent site.yaml. Testing that means
# handing it a deliberately inconsistent one, which means a throwaway config
# directory -- the real `config/` is versioned input data and is never written
# to by a test.
#
# `CONFIG_DIR` is read into a module constant at import, so the env var alone
# is too late; `load_site` is `lru_cache`d, so the cache has to be cleared on
# the way in AND on the way out, or a poisoned entry leaks into later tests.


def _load_raw(filename):
    """The config file as plain nested dicts, for mutating."""

    with open(os.path.join(config.CONFIG_DIR, filename), encoding="utf-8") as f:
        return yaml.safe_load(f)


@contextlib.contextmanager
def _temporary_config(overrides):
    """Point the loader at a copy of `config/` with `overrides` applied."""

    from config import loader

    original = loader.CONFIG_DIR
    scratch = tempfile.mkdtemp(prefix="blade-config-")
    try:
        shutil.copytree(original, scratch, dirs_exist_ok=True)
        for filename, data in overrides.items():
            with open(os.path.join(scratch, filename), "w",
                      encoding="utf-8") as f:
                yaml.safe_dump(data, f)

        loader.CONFIG_DIR = scratch
        _clear_config_caches()
        yield
    finally:
        loader.CONFIG_DIR = original
        _clear_config_caches()
        shutil.rmtree(scratch, ignore_errors=True)


def _clear_config_caches():
    """Every `lru_cache`d loader, so an override is actually read."""

    from config import loader

    for load in (loader.load_site, loader.load_design_rotor,
                 loader.load_phase_vi_rotor, loader.load_polar_cache):
        load.cache_clear()


def test_an_inconsistent_atmosphere_is_rejected():
    """
    The atmosphere check has guarded site.yaml since Task 1 but its rejection
    path was never actually exercised. Driving it here, alongside the new wind
    resource check, so both are known to fire rather than assumed to.
    """

    site_yaml = _load_raw("site.yaml")
    site_yaml["atmosphere"]["air_density"] *= 1.05

    with _temporary_config({"site.yaml": site_yaml}):
        with pytest.raises(config.ConfigError, match="air_density"):
            config.load_site()


def test_resolved_wind_resource_is_a_usable_number():
    """
    Replaces `test_unresolved_wind_resource_raises_rather_than_substituting`,
    whose job ended when the resource landed on 2026-09-13.

    Ground rule 3 has not gone away -- it has been satisfied. What is checked
    now is that the three fields are genuinely resolved, positive, and physical,
    so a `TODO` reintroduced by a bad merge fails here rather than three modules
    downstream.
    """

    site = config.load_site()

    for name in ("weibull_k", "weibull_c_ms", "mean_wind_speed_ms"):
        value = getattr(site, name)
        assert config.is_resolved(value), f"{name} should be resolved"
        assert float(value) > 0.0


def test_gwa_area_is_still_unresolved_and_behaves_like_it():
    """
    `gwa_area` is deliberately still TODO, and this is the guard on that.

    The 2026-09-13 extraction was a GASP *point*, not a Global Wind Atlas
    *area*, so there is no area selection to record. Writing the point's
    coordinates there would assert an extraction that was never performed. This
    keeps the sentinel's behaviour honest in the meantime: not a number, not
    `None`, and not quietly falsy.
    """

    value = config.load_site().gwa_area
    assert not config.is_resolved(value)

    with pytest.raises(config.UnresolvedConfigError, match="gwa_area"):
        float(value)
    with pytest.raises(config.UnresolvedConfigError, match="gwa_area"):
        # `if site.gwa_area:` is a silent fallback waiting to happen.
        bool(value)


def test_coordinates_are_the_khomas_hochland_site():
    """
    Latitude and longitude landed with the same extraction.

    The longitude SIGN is the assertion that matters. GASP labelled the point
    "W16.558"; west 16.558 is open ocean several hundred km off Namibia, where
    the panel's terrain-derived statistics could not have come from. East is
    what the pinned satellite views show and east is what is recorded.
    """

    site = config.load_site()

    assert site.latitude_deg == pytest.approx(-22.427972, abs=1e-5)
    assert site.longitude_deg == pytest.approx(16.557833, abs=1e-5)
    assert site.longitude_deg > 0.0, "the site is in EASTERN Namibia-longitude"


def test_mean_wind_speed_must_agree_with_k_and_c():
    """
    The loader rejects a wind resource whose derived mean has drifted.

    This is the check the 2026-09-10 resumption checklist assumed already
    existed. It did not -- only the atmosphere was checked -- so it was added
    with the resource. Editing `k` and forgetting the mean is otherwise a
    silent few-percent error in every AEP figure, which is exactly the failure
    mode the atmosphere check exists to prevent.
    """

    site_yaml = _load_raw("site.yaml")
    site_yaml["wind_resource"]["weibull_k"] *= 1.10

    with _temporary_config({"site.yaml": site_yaml}):
        with pytest.raises(config.ConfigError, match="mean_wind_speed_ms"):
            config.load_site()


def test_a_partially_resolved_wind_resource_is_rejected():
    """
    All three fields, or none. Two out of three is an error.

    Filling `k` and `c` but leaving the mean as TODO would otherwise load
    cleanly and hand a sentinel to whatever read the mean -- a failure far from
    its cause.
    """

    site_yaml = _load_raw("site.yaml")
    site_yaml["wind_resource"]["mean_wind_speed_ms"] = "TODO: forgotten"

    with _temporary_config({"site.yaml": site_yaml}):
        with pytest.raises(config.ConfigError, match="partially resolved"):
            config.load_site()


def test_design_bounds_are_resolved_with_their_recorded_values():
    """
    The design-variable bounds were TODO (plan 7.1) until 2026-09-19; they
    are resolved now, with the basis recorded in the YAML. Pinned to the
    decided numbers so an edit has to restate its basis, and so a `TODO`
    creeping back in would fail here rather than three modules downstream.
    """

    bounds = config.load_design_rotor().parameterisation

    assert float(bounds.chord_min_m) == pytest.approx(0.045)
    assert float(bounds.chord_max_m) == pytest.approx(0.30)      # O4, 2026-09-19
    assert float(bounds.twist_min_deg) == pytest.approx(-2.0)
    assert float(bounds.twist_max_deg) == pytest.approx(35.0)


def test_a_todo_bound_still_raises_when_used():
    """The unresolved mechanism is kept: a TODO in the YAML raises on use."""

    rotor_yaml = _load_raw("rotor_design.yaml")
    rotor_yaml["parameterisation"]["bounds"]["chord_max_m"] = "TODO: reopened"

    with _temporary_config({"rotor_design.yaml": rotor_yaml}):
        bounds = config.load_design_rotor().parameterisation
        with pytest.raises(config.UnresolvedConfigError, match="chord_max_m"):
            float(bounds.chord_max_m)


def test_the_machine_facts_are_recorded():
    """
    B1 (rotor-speed ceiling) and the solidity cap, resolved provisionally on
    2026-09-19. `max_tip_speed_ms` is derived from the rpm and the radius so
    the two cannot drift apart.
    """

    design = config.load_design_rotor()

    assert design.max_rotor_speed_rpm == pytest.approx(300.0)
    assert design.max_tip_speed_ms == pytest.approx(300.0 * 2.0 * math.pi / 60.0 * 2.0)
    assert design.max_local_solidity == pytest.approx(0.5)


def test_no_ceiling_is_spelled_null():
    """`max_rotor_speed_rpm: null` means no ceiling -- the control case."""

    rotor_yaml = _load_raw("rotor_design.yaml")
    rotor_yaml["operating"]["max_rotor_speed_rpm"] = None

    with _temporary_config({"rotor_design.yaml": rotor_yaml}):
        design = config.load_design_rotor()
        assert design.max_rotor_speed_rpm is None
        assert design.max_tip_speed_ms is None


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


def test_the_mass_problem_is_recorded():
    """
    The 2026-09-20 re-pitch (docs/PLAN-mass-objective-2026-09-20.md): the
    objective is the shell material and both manufacturability rows are on.
    """

    design = config.load_design_rotor()

    assert design.mass_model == "shell"
    assert design.monotone_chord is True
    assert design.monotone_twist is True
    assert design.min_chord_m == pytest.approx(0.060)      # the buildable-tip floor, 2026-09-20
    assert design.min_chord_m > design.parameterisation.chord_min_m   # a row, not the box


def test_the_laminate_is_recorded_with_its_source():
    """
    The structural inputs resolved 2026-09-20 (evening): one laminate, the
    E-LT-5500/EP-3 row of Griffith & Ashwill (2011) Table 19, a 2 mm skin
    tied to the 60 mm tip, and the GL hand-lay-up safety factor
    (docs/MATERIALS-STRUCTURAL-INPUTS.md). Pinned so an edit has to restate
    its basis, as for the bounds.
    """

    design = config.load_design_rotor()

    assert design.laminate_density_kg_m3 == pytest.approx(1920.0)
    assert design.shell_thickness_m == pytest.approx(0.002)
    assert design.youngs_modulus_pa == pytest.approx(41.8e9)
    assert design.allowable_stress_pa == pytest.approx(702.0e6)
    # gamma_f * gamma_M0 * C1a * C2a * C3a(hand lay-up) * C4a(non post-cured)
    assert design.safety_factor == pytest.approx(1.35 * 1.35 * 1.35 * 1.1 * 1.2 * 1.1, abs=5e-5)
    assert design.design_allowable_stress_pa == pytest.approx(702.0e6 / 3.5725, rel=1e-12)
    # A 2 mm skin on the 60 mm floor: two skins fit inside the 10 % section.
    assert 2.0 * design.shell_thickness_m < 0.10 * design.min_chord_m
    for name in ("laminate_density_kg_m3", "shell_thickness_m", "youngs_modulus_pa",
                 "allowable_stress_pa", "safety_factor"):
        assert isinstance(getattr(design, name), float), name


def test_the_tip_clearance_stays_todo_and_behaves_like_it():
    """
    `tip_clearance_m` is machine geometry with no defensible number yet: not
    a number, not `None`, not falsy, and the derived allowable it would gate
    is untouched by it. The `Unresolved` mechanism is what the absolute
    deflection row relies on to refuse assembly.
    """

    design = config.load_design_rotor()

    value = design.tip_clearance_m
    assert not config.is_resolved(value)
    with pytest.raises(config.UnresolvedConfigError, match="tip_clearance_m"):
        float(value)
    with pytest.raises(config.UnresolvedConfigError, match="tip_clearance_m"):
        bool(value)
    assert config.is_resolved(design.design_allowable_stress_pa)


def test_a_structural_input_that_is_neither_number_nor_todo_is_rejected():
    """`41.8e9` without an exponent sign is a YAML string, not a float; caught at load."""

    rotor_yaml = _load_raw("rotor_design.yaml")
    rotor_yaml["structure"]["youngs_modulus_pa"] = "41.8e9"

    with _temporary_config({"rotor_design.yaml": rotor_yaml}):
        with pytest.raises(config.ConfigError, match="exponent"):
            config.load_design_rotor()


def test_a_reopened_structural_input_comes_back_unresolved():
    rotor_yaml = _load_raw("rotor_design.yaml")
    rotor_yaml["structure"]["laminate_density_kg_m3"] = "TODO: laminate re-opened"

    with _temporary_config({"rotor_design.yaml": rotor_yaml}):
        design = config.load_design_rotor()
        assert not config.is_resolved(design.laminate_density_kg_m3)
        with pytest.raises(config.UnresolvedConfigError, match="re-opened"):
            float(design.laminate_density_kg_m3)


def test_an_unknown_mass_model_is_rejected():
    rotor_yaml = _load_raw("rotor_design.yaml")
    rotor_yaml["objective"]["mass_model"] = "hollow"

    with _temporary_config({"rotor_design.yaml": rotor_yaml}):
        with pytest.raises(config.ConfigError, match="mass_model"):
            config.load_design_rotor()


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
