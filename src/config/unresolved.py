"""
The `TODO` mechanism: how an input that has not arrived yet is represented,
and why it is an object rather than a number.

Plan section 1.3 and the remediation work order's ground rule 3 both say the
same thing: a fabricated wind resource propagates silently into every AEP
figure downstream and invalidates the whole results chapter. So a config field
that is still `TODO` must not be a plausible-looking float, and must not be
`None` either -- `None` fails late, somewhere far from the config, with a
message about `NoneType`.

`Unresolved` fails at the point of first use, with the field name and the note
from the YAML file in the message. Every arithmetic and conversion protocol is
wired to raise, so there is no path by which one of these reaches a result.

Author: MJ Hendrikse
Project: DSP810S -- Inverse Design of Small Wind Turbine Blades
"""


class UnresolvedConfigError(RuntimeError):
    """A config field still marked `TODO` was used as if it had a value."""


class Unresolved:
    """
    Placeholder for a config field whose YAML value is a `TODO:` string.

    Carries the field's dotted path and the note from the file, and raises
    `UnresolvedConfigError` naming both on any attempt to use it as a value.

    Truthiness raises too, deliberately: `if cfg.weibull_k:` is a silent
    fallback waiting to happen. Use `config.is_resolved(...)` to test.
    """

    __slots__ = ("field", "note")

    def __init__(self, field, note):
        self.field = field
        self.note = note

    def __repr__(self):
        return f"<Unresolved {self.field}: {self.note}>"

    def _fail(self, *_args, **_kwargs):
        raise UnresolvedConfigError(
            f"config field '{self.field}' is still TODO ({self.note}) and has "
            f"no value. It must be supplied from real data before anything "
            f"that depends on it is run -- it is never defaulted or invented."
        )

    # Every route by which this could be mistaken for a number.
    __float__ = __int__ = __index__ = __bool__ = _fail
    __add__ = __radd__ = __sub__ = __rsub__ = _fail
    __mul__ = __rmul__ = __truediv__ = __rtruediv__ = _fail
    __pow__ = __rpow__ = __neg__ = __abs__ = _fail
    __lt__ = __le__ = __gt__ = __ge__ = _fail
    __format__ = _fail


def is_resolved(value):
    """True if `value` is a real value rather than an outstanding `TODO`."""

    return not isinstance(value, Unresolved)
