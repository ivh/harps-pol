# CLAUDE.md

Guidance for Claude Code working in this repository. Things that cost time to
find out, and that reading the code alone will not tell you.

## Overview

`espdr_demod_pol` is a pyesorex recipe that demodulates HARPSpol sequences from
the S2D products of `espdr_sci_red`. Layout follows `~/pycr2res.git`: thin CPL
interface in `pyrecipes/`, real work in `pyespdr/`, EDPS extension in
`workflows/harpspol/`.

- Run anything through `uv`. On macOS prefix with `env -u DYLD_LIBRARY_PATH`, or
  an installed ESO pipeline's CPL shadows the one PyCPL bundles and imports die
  with `Symbol not found: _cpl_wcs_duplicate`.
- `uv run pytest` — 20 tests. 17 need the reference data next door
  (`../112.25MG.001/reduc`) and skip without it; the 3 workflow tests need
  `~/pipes/harps-3.6.0/workflows`.

## HARPSpol header facts

Verified against programme 112.25MG.001, 2024-01-02.

- **`ESO INS RET25 POS`** is the quarter-wave plate angle: 45/135/225/315 for a
  4-exposure cycle, 45/135 for 2. RET25 is the only moving element recorded —
  HARPS has two fibres fed by a beam-splitter, so there is no nodding or any
  other degree of freedom. Don't go looking for one.
- **`ESO DPR TECH = ECHELLE,CIRPOL`** exists in the *raw* frames only. It is not
  propagated to S2D.
- **`ESO INS MODE`** is `HARPSPOL` in raw but **`HARPS`** in the S2D. Never use
  it to detect polarimetry in a reduced product.
- **`ESO TPL NAME = HARPS_pol_obs_cir`** *does* survive into S2D, and is the
  reliable discriminator. Note `TPL ID` is the generic `HARPS_pol_obs_all`;
  `TPL NAME` carries the `_cir` / `_lin` distinction. A linear template still
  cannot tell you Q from U — that is in the half-wave plate angles.
- `TPL EXPNO` / `TPL NEXP` / `TPL START` are all present in S2D, so template
  grouping needs no filename parsing.
- **Fibre A has 71 orders, fibre B has 70.** The one B lacks is A index 44,
  starting at 5246.8 Å. Matched orders agree to 0.009 Å in start wavelength
  against a ~40 Å order spacing, so matching by wavelength is unambiguous. Why
  the order is missing is still unknown.
- S2D extensions: `SCIDATA`, `ERRDATA`, `QUALDATA`, `WAVEDATA_VAC_BARY`,
  `WAVEDATA_AIR_BARY`, `DLLDATA_*`. Data is float32, so "bit-identical"
  comparisons against stored products bottom out at ~1e-7 relative.

## astropy header traps

- `header["PRO_CATG"] = x` silently writes a literal 8-character `PRO_CATG`
  keyword. It is *not* `ESO PRO CATG`, and the inherited one from the template
  survives untouched — the product then still claims to be whatever it was
  copied from. This was a real bug in the original script.
- Assigning a new long keyword emits `VerifyWarning` even though the HIERARCH
  card is created correctly. Write `header["HIERARCH ESO PRO ..."] = x` to keep
  the log clean. Updating an *existing* long keyword does not warn.

## Demodulation maths

For `X = (R**p - 1) / (R**p + 1)`, `dX/dR = 2p * R**(p-1) / (R**p + 1)**2`.
p = 1/2 for a 2-exposure cycle, 1/4 for 4. The original `error_helper` was a
factor 2 too large at p = 1/4; 2-exposure errors were right. Confirmed against
the stored products: the ratio is exactly 0.5000 for the 4-exposure target and
1.0000 for all six 2-exposure ones.

Every spectrum enters the ratio R exactly once, as numerator or denominator, so
its relative error contributes the same either way and they add in quadrature.
The original code had `Rae`/`Rbe` built from the wrong fibre's errors, which
cancelled in the total and was therefore harmless but unreadable.

**The blaze cancels in the Stokes ratios but not in the intensity.** So
`S2D_POL_STOKES` is identical whether you feed it `S2D_*` or `S2D_BLAZE_*`, and
only `S2D_POL_I` tells you which was used. That is how the flavour each
reference product was made with was recovered — `demod_all.sh` uses plain S2D
for the first four targets and blaze for the last three, and `SEQUENCES` in
`tests/test_demod.py` records it.

## Pipelines on this machine

- `~/pipes/espdr-3.3.0` — ESPRESSO. Has an EDPS workflow, gated on
  `INSTRUME == "ESPRESSO"` throughout. No polarimetric mode.
- `~/pipes/harps-3.3.0` — HARPS, **no EDPS workflow**, Reflex OCA only.
- `~/pipes/harps-3.6.0` — HARPS, **has** the EDPS workflow. Check the newest
  tarball before concluding a pipeline lacks something.

The HARPS 3.6.0 workflow already handles HARPSpol:

- `harps_rules.is_obj_cirpol` / `is_obj_linpol` classify pol science, both onto
  the ordinary `OBJ_SKY` tag.
- `harps_task_functions.set_sky_subtraction` turns off sky correction when
  `INS MODE == HARPSPOL`.
- `should_combine` explicitly excludes HARPSPOL from `combine_science` — that
  gap is exactly where `demod_pol` belongs.
- The `object` task carries `.with_min_group_size(2)` and
  `.with_grouping_keywords([kwd.tpl_start])`, so per-template grouping comes for
  free downstream.
- It does **not** classify the science `S2D_A` / `S2D_B` products; only the
  calibration S2Ds and `S1D_FINAL_A`. `harpspol_wkf.py` adds those rules.

## EDPS (1.7.1)

- A workflow module is extended, not forked, by importing another one:
  `create_workflow` recurses into any module member where
  `isinstance(obj, ModuleType) and 'wkf' in obj.__name__` and unions its tasks,
  data sources and classification rules. `edps/workflow/meta_wkf.py` uses the
  same trick, so it is intended.
- Discovery (`edps/scripts/server.py`): walks each `workflow_dir` looking for a
  file matching `{package}.*_wkf.py$` where `package` is the containing
  directory's name, registers it as `{package}.{stem}`, and puts the
  **parent** of that directory on `sys.path`. So the layout must be
  `workflows/<name>/<name>_wkf.py`, and the pipeline's own workflow must be
  importable as `harps` — a symlink named `harps` is enough.
- `workflow_dir` is a comma-separated list, so ours sits alongside the
  pipeline's with no installation into its tree.
- To build a workflow without a server:
  `WorkflowManager(config=None).create_workflow(module)`.
- `Task` attributes are not what you would guess: `.command` (not `.recipe`),
  `.main_input.name`, `.input_filter` is a set of **strings**,
  `.grouping_keywords`, `.min_group_size`. Classification rules have no `.name`
  on every subclass.

## pyesorex / PyCPL

- pyesorex loads C recipes as well as Python ones (`cpl.ui.CRecipe`; it even
  handles both existing under one name). But **CPL versions cannot be mixed in
  one process** — pointing it at plugins built against CPL 7.3.2 while PyCPL
  bundles 7.4 aborts with
  `CplCore-ERROR: CPL 7.3.2 memory management subsystem is not initialized`.
  This matters because EDPS has a single global `esorex_path`, so running a
  Python recipe inside a workflow means running *everything* under pyesorex.
- PyCPL here comes from `https://ivh.github.io/pycpl/simple/`, the user's own
  repackaging with the C libraries bundled (`~/pycpl.git`).

## Still open

See the end of README.md. The short version: the `S2D_POL_*` PRO.CATG values
are invented and need registering with ESO, product provenance is hand-written
rather than from `cpl.dfs`, and the linear-polarimetry template name and its
Q/U convention are unknown for want of data.
