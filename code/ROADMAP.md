# girt branch — roadmap

Goal of this branch: **redo and clean the whole Rt project from the
generation-interval (renewal) perspective**, replacing the mechanistic (SEIR
compartmental) framing as the primary one. `mechrt` becomes a *special case*
(the GI implied by `pi_lat`/`pi_IR`), not the foundation. We proved the two are
identical when `g = zeta^EI/mu^IR` (retro + real-time, to ~1e-14), so this is a
reframing + cleanup, not a re-derivation — except on **real data**, where a
free/literature GI is a genuinely different (and more defensible) estimator.

Status legend: [ ] todo  [~] partial  [x] done

---

## Why GI-first (what it buys us)

1. **One epi input, not three.** The method needs a generation interval `g` and
   a reporting delay `pi_EY`. Gone: latent `pi_EI`, infectious `pi_IR`, recovery
   rate, mean-infectious-period, the `sum(P_still_I)` discretization gotcha,
   susceptible tracking `S/N`, prevalence `P_hat`/`K_EP`.
2. **Cleaner code.** mechrt's `design.R` is ~650 lines carrying beta-mode, kink
   regressors, regression-spline + KKT tail/left constraints, S/I imputation
   variants, `REMOVED_OPTIONS.md`. The GI design is ~180 lines, natural-spline
   only. Big surface reduction.
3. **Cleaner exposition.** The renewal equation is the standard object; the
   compartmental connection is one proposition (appendix3 `compartmental-generation`).
4. **Real-data novelty.** Only under GI misspecification (real data) do girt and
   mechrt differ. Decoupling `g` from SEIR delays (literature/estimated GI) is
   the actual scientific contribution this branch can make.

---

## Guiding principles

- **Parallel, then replace.** New girt pipelines write to `results|figures/girt/...`
  trees and are validated against cached mechrt outputs before anything mechrt is
  retired. Nothing in `code/mechrt/` is edited on this branch until girt is the
  proven canonical.
- **Golden-master tested.** Lock the equivalence (cached mechrt numbers) as
  regression tests so cleanup/refactors can't silently drift.
- **No dead options.** Every knob girt keeps must be used by a real pipeline.

---

## Phase 0 — Package hardening  [~]

- [x] Self-contained `code/` (loader, R/, `examples/demo.R`, README).
- [x] Validated retro + real-time equivalence vs mechrt to ~1e-14 (the
      comparison scripts live on `main`, which retains the mechanistic codebase).
- [ ] **Settle the name** (girt is provisional) — sweep folder, `girt.R`, `gi_*`,
      `GIRT_DIR`.
- [ ] **Unify the solver.** With natural-spline (empty tail-constraint A), the
      KKT path (`gi_solve`) and null-space taper path (`.gi_solve_taper`) collapse
      to the same solve; merge into one `gi_solve(..., lam_taper, P_taper)`.
- [ ] **Decide DoW scope** (keep for real data; verify the profiling path on a
      DoW dataset).
- [ ] Trim copied-engine cruft that the GI method never hits (regression-spline
      branches in extract/taper comments, `param_mode="beta"` guards, kink).
- [ ] Optional: convert to a real R package (DESCRIPTION/NAMESPACE, `tests/testthat`).

## Phase 1 — Simulation parity  [ ]

Port the retrospective + real-time + conformal benchmarks to girt, reproduce the
canonical numbers, write to `results|figures/girt/sim/...`.

- [ ] flu3 (canonical wiggly+smooth): retro_benchmark, daily real-time vintages,
      conformal, wiggly/smooth tables. (retro+RT spot-checks already match.)
- [ ] flu4 (overdispersed NB1): same.
- [ ] rtestim (4 scenarios): same, incl. the per-scenario EpiNow2 baselines.
- [ ] flu2 / flu_based: port or mark for deprecation (flu_based already deprecating).
- [ ] One girt-native sim config object (`g`, `pi_EY`, `severity`) replacing the
      mechanistic `(pi_lat, pi_IR, pi_EH, mip)` configs.

## Phase 2 — Real data  [ ]

This is where GI-first is a *different* model, not a rename.

- [ ] Real flu (HHS hosp): port the retro + real-time + conformal pipeline; first
      with the SEIR-implied `g` (parity check), then with a literature/estimated
      GI (the new experiment).
- [ ] covidestim comparison: same; compare girt(literature GI) vs
      girt(SEIR-implied) vs covidestim point estimate.
- [ ] **Experiment:** does decoupling `g` from the fitted delays improve
      calibration / point accuracy on real data? (the headline real-data result.)

## Phase 3 — Paper  [ ]

- [ ] Reframe `paper3.tex` GI-first: renewal equation + `g` as the primary model;
      compartmental version as the `compartmental-generation` special case.
- [ ] Regenerate all figures/tables from girt outputs.
- [ ] Update notation table (drop mechanistic-only symbols; `g`, `pi_EY` central).

## Phase 4 — Cleanup / retire mechanistic  [ ]

- [ ] Once girt reproduces every result, make girt canonical: repoint pipelines,
      retire `code/mechrt/` (or keep a thin `gi_from_compartmental` shim).
- [ ] Delete `deprecated/`, `REMOVED_OPTIONS.md`, dead branches.
- [ ] Single deconvolution (identity-weighted exposures); remove log/plain legacy.
- [ ] Reconcile `code/mechrt/renewal/` (the `_gi_obs_lib.R` prototype + the
      latent-exposure `_renewal_lib.R`) into girt or `examples/`.

---

## Open decisions (need user input)

1. **Package name** (girt → ?).
2. **Data/results/figures layout**: parallel `…/girt/` trees during migration,
   then replace in place? Or girt canonical from the start?
3. **mechrt's fate**: delete after parity, or keep as the compartmental shim?
4. **Real R package** (DESCRIPTION + testthat) vs source-loader as now?
5. **Scope of "everything"**: include flu2/flu_based, or let them stay deprecated?

## Immediate next steps (when we resume)

1. Lock golden-master tests: girt vs cached mechrt on flu3 retro + RT (already
   ~1e-14) so the cleanup is safe.
2. Unify the solver + trim engine cruft (Phase 0).
3. Port flu3 retro_benchmark to girt end-to-end (Phase 1, canonical first).
