"""
Central finite differences and the step-size sweep.

Central only, never forward
----------------------------
A forward difference has truncation error O(h); central has O(h^2). With an
objective that costs ~0.2 s and a noise floor set by the polar interpolant's
C2 breaks rather than by round-off, the extra evaluation per variable buys a
much wider usable band of `h`, and it is the central estimate at a chosen
`h*` that the discrete adjoint is measured against (Tier 3). A
forward difference is never good enough to be that reference.

The sweep, and why it saves everything
---------------------------------------
`step_size_sweep` evaluates the central gradient at every step in `steps`
and returns all of them. The choice of `h*` is made *afterwards* by the
caller (see `verification/fd_step_size/run_sweep.py`), from the recorded
values, so the criterion can be changed or the curves re-plotted against a
later reference without re-running the ~300 evaluations.

If a stencil point leaves the polar cache, `CachedPolar` raises
`PolarDomainError`. That is recorded -- the affected component is `nan` and
the `(h, j)` pair is listed -- and the sweep continues. It is never replaced
with a forward difference or a clamped value: an out-of-cache step is a fact
about that `h`, not a number to be patched over.

Both functions are pure: no config access, no rotor, no knowledge of what
`fun` is. `fun` is any `ndarray -> float`.

Author: MJ Hendrikse
Project: MCP820S -- Gradient-Based Aerodynamic Optimisation of a Small Wind Turbine Blade
"""

import numpy as np

from polars.interpolant import PolarDomainError

#: The label recorded in a sweep for a step whose stencil left the polar
#: cache. A string, deliberately: it is written into the sweep JSON as-is.
OUT_OF_CACHE = "out of cache"


def _step_vector(h, n):
    """`h` as an (n,) array of strictly positive steps, scalar or per-variable."""

    steps = np.broadcast_to(np.asarray(h, dtype=float), (n,)).copy()
    if not np.all(steps > 0.0):
        raise ValueError(f"finite-difference steps must be positive, got {h}")
    return steps


def _central_component(fun, u, j, h):
    """(f(u + h e_j) - f(u - h e_j)) / (2h) for one variable `j`."""

    forward = np.array(u, dtype=float)
    backward = np.array(u, dtype=float)
    forward[j] += h
    backward[j] -= h
    return (fun(forward) - fun(backward)) / (2.0 * h)


def central_difference(fun, u, h):
    """
    Central-difference gradient of `fun` at `u`.

    Parameters
    ----------
    fun : callable
        `ndarray (n,) -> float`.
    u : array_like
        The point, shape (n,).
    h : float or array_like
        Step size: one scalar for every variable, or an (n,) array with one
        step per variable (the `h*_j` from the step-size study).

    Returns
    -------
    (grad, n_evals)
        `grad` shape (n,); `n_evals` is the number of `fun` calls made,
        always `2n`.

    Raises
    ------
    polars.interpolant.PolarDomainError
        Propagated unchanged if any stencil point leaves the polar cache.
        This function makes no attempt to recover; `step_size_sweep` is the
        place where that is recorded.
    """

    u = np.asarray(u, dtype=float)
    n = u.shape[0]
    steps = _step_vector(h, n)

    grad = np.empty(n)
    for j in range(n):
        grad[j] = _central_component(fun, u, j, steps[j])

    return grad, 2 * n


def step_size_sweep(fun, u, steps):
    """
    The central gradient at every step in `steps`, all values kept.

    Parameters
    ----------
    fun, u
        As for `central_difference`.
    steps : iterable of float
        Scalar steps to try, e.g. `np.logspace(-2, -9, 15)`. Each is applied
        to every variable.

    Returns
    -------
    dict
        `{float(h): grad}` in the order given, `grad` shape (n,). A component
        whose stencil raised `PolarDomainError` is `nan`; the list of such
        `(h, j)` pairs is attached under the key `OUT_OF_CACHE` so the record
        survives serialisation. Nothing is substituted.
    """

    u = np.asarray(u, dtype=float)
    n = u.shape[0]

    results = {}
    out_of_cache = []
    for h in steps:
        h = float(h)
        grad = np.empty(n)
        for j in range(n):
            try:
                grad[j] = _central_component(fun, u, j, h)
            except PolarDomainError:
                grad[j] = np.nan
                out_of_cache.append((h, j))
        results[h] = grad

    results[OUT_OF_CACHE] = out_of_cache
    return results
