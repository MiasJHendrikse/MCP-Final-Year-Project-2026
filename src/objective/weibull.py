"""
The Weibull wind-speed distribution.

    f(V) = (k/c) * (V/c)^(k-1) * exp(-(V/c)^k)

`k` and `c` come from the Global Wind Atlas extraction at the site's 20 m hub
height. They are still `TODO` in `config/site.yaml`, so `from_config()` raises
rather than substituting anything -- see `docs/OUTSTANDING-INPUTS.md` section 1.

This is the mechanism, complete and tested, waiting on two numbers. Everything
downstream of it that does not need those numbers (the bin scheme, the
operating strategy, per-bin power, and the smoothness gate's surrogate) is
built and works today.

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import math
from dataclasses import dataclass

import numpy as np

from config import is_resolved, load_site
from config.unresolved import UnresolvedConfigError


@dataclass(frozen=True)
class WeibullResource:
    """
    A two-parameter Weibull wind-speed distribution.

    Parameters
    ----------
    k : float
        Shape parameter, dimensionless. Higher k means a narrower spread of
        wind speeds about the mean.
    c : float
        Scale parameter, m/s.
    """

    k: float
    c: float

    def __post_init__(self):
        if self.k <= 0.0:
            raise ValueError(f"Weibull k must be positive, got {self.k}")
        if self.c <= 0.0:
            raise ValueError(f"Weibull c must be positive, got {self.c}")

    @classmethod
    def from_config(cls):
        """
        Build from `config/site.yaml`.

        Raises
        ------
        config.unresolved.UnresolvedConfigError
            While `weibull_k` and `weibull_c_ms` are `TODO`, which is the
            current state. Named fields in the message, so a caller that
            reaches here is told exactly what is missing.
        """

        site = load_site()
        missing = [name for name, value in (("weibull_k", site.weibull_k),
                                            ("weibull_c_ms", site.weibull_c_ms))
                   if not is_resolved(value)]
        if missing:
            raise UnresolvedConfigError(
                "the site wind resource is still TODO in config/site.yaml: "
                f"{', '.join(missing)}. These come from the Global Wind Atlas "
                "extraction at 20 m hub height; see "
                "docs/OUTSTANDING-INPUTS.md section 1. Construct "
                "WeibullResource(k, c) explicitly for a study that states its "
                "own provisional values."
            )

        return cls(k=float(site.weibull_k), c=float(site.weibull_c_ms))

    def pdf(self, v):
        """Probability density at wind speed `v`, per (m/s)."""

        v = np.asarray(v, dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            density = ((self.k / self.c) * (v / self.c) ** (self.k - 1.0)
                       * np.exp(-((v / self.c) ** self.k)))
        return np.where(v > 0.0, density, 0.0)

    def cdf(self, v):
        """Cumulative probability up to `v`."""

        v = np.asarray(v, dtype=float)
        return np.where(v > 0.0, 1.0 - np.exp(-((v / self.c) ** self.k)), 0.0)

    def probability_between(self, lower, upper):
        """
        Exact probability mass in [lower, upper], from the CDF.

        Bin masses come from the CDF rather than from the PDF times the bin
        width: the CDF is exact for any bin width, so the bin scheme's accuracy
        is decided by how well the *power* curve is resolved, not by how well
        the distribution is. Those are separate concerns and mixing them makes
        a bin-width study uninterpretable.
        """

        return self.cdf(upper) - self.cdf(lower)

    @property
    def mean_speed(self):
        """V-bar = c * Gamma(1 + 1/k), m/s."""

        return self.c * math.gamma(1.0 + 1.0 / self.k)
