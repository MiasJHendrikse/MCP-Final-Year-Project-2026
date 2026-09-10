# `verification/`

Versioned report figures and data — the artefacts that back a claim in the
writeup. Work order Task 8, and the plan's working conventions.

Distinct from `results/`, which is generated output and partly gitignored.
The rule is: if regenerating it would change what the dissertation says, it
lives here and is committed. If it is scratch, a sweep you ran once, or an
intermediate you would not cite, it belongs in `results/`.

Each subdirectory carries its own README stating what the artefact shows, how
it was produced, and what would make it wrong.

| Directory | Status | Backs |
|---|---|---|
| `polar_interpolant/` | populated (Task 3) | that the C¹ interpolant removed the bilinear derivative staircase |
| `phase_vi/` | populated (Task 5) | that the solver converges across the envelope with checked residuals |
| `smoothness_gate/` | awaiting plan step 1.8 | that the objective is smooth enough to differentiate |
| `representation_study/` | awaiting plan step 1.6 | the choice of blade parameterisation and control-point count |

## Regenerating

Every subdirectory has a generator script next to its output, and the output is
committed alongside it. Regeneration is deliberate, not routine — the same rule
`tests/golden/README.md` states for the golden files. A figure that changes
without an explanation in the commit message is a defect, not an update.
