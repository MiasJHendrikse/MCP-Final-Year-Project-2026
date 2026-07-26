Selig-format coordinate files for airfoils with no NACA parametric equivalent, loaded
into XFOIL via `run_xfoil_polar(airfoil_cmd="LOAD <path>")` (see `src/xfoil/xfoil_runner.py`).

- `s809.dat` — NREL S809, fetched from the airfoiltools.com mirror of the UIUC
  Applied Aerodynamics coordinate database (the canonical
  m-selig.ae.illinois.edu path 404'd at fetch time). 66 points, unit chord.
  Geometry re-verified 2026-07-26 against an independently downloaded
  percent-chord copy of the same section — see the 2026-07-26 journal entry for
  the numbers.

- `s809_qblade.dat` — the same section, formatted for import into QBlade
  (Airfoil Design > Import Airfoil). QBlade reads the same Selig layout XFOIL
  does — a name line followed by `x y` pairs on unit chord, running from the
  trailing edge forward over the upper surface to the leading edge and back
  along the lower surface — so this file differs from `s809.dat` only in the
  name line (`S809`, with the apostrophe of "NREL's S809 Airfoil" removed since
  it is carried straight through into QBlade's object name) and in using a
  fixed-width 7-decimal column format. Coordinates are identical to
  `s809.dat` to machine precision; keep the two in step if either is ever
  regenerated.
