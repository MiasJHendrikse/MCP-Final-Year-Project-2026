Selig-format coordinate files for airfoils with no NACA parametric equivalent, loaded
into XFOIL via `run_xfoil_polar(airfoil_cmd="LOAD <path>")` (see `src/xfoil/xfoil_runner.py`).

- `s809.dat` — NREL S809, fetched from the airfoiltools.com mirror of the UIUC
  Applied Aerodynamics coordinate database (the canonical
  m-selig.ae.illinois.edu path 404'd at fetch time). 66 points, unit chord.
  Geometry re-verified 2026-07-26 against an independently downloaded
  percent-chord copy of the same section — see the 2026-07-26 journal entry for
  the numbers.

- `sg6043.dat` — SG6043, the design rotor's section (plan §3.1). Fetched
  2026-08-29 from the canonical UIUC Applied Aerodynamics coordinate database,
  `https://m-selig.ae.illinois.edu/ads/coord/sg6043.dat`, and kept **verbatim**
  — the file is byte-for-byte what UIUC serves, including its leading-dot
  number format and the `.999999` closing x, rather than reformatted to match
  `s809.dat`. 81 points, unit chord, TE → upper → LE → lower → TE.

  Verified before use, not assumed: the coordinates give **max thickness
  10.015 % at x/c 0.321** and **max camber 5.500 % at x/c 0.497**, against
  published figures of 10.02 % at 32.10 % and 5.50 % at 49.70 % for this
  section. Agreement to three decimals in both magnitudes and both locations
  is what identifies the file as genuinely SG6043; the contour is also closed,
  x-monotonic on each surface, and has zero trailing-edge gap.

  This closes one of the two inputs plan §1.2 was waiting on. The other — the
  UIUC *Summary of Low-Speed Airfoil Data* Vol. 3 experimental tables at
  Re = 60k/100k/200k/300k — is still outstanding, and without it the `n_crit`
  sensitivity study can sweep values but cannot *select* one against measured
  data the way S809's `Ncrit = 5` was selected. See
  `config/polars_sg6043.yaml`.

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
