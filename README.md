# pyespdr

Polarimetric demodulation for the [ESPRESSO
pipeline](https://www.eso.org/sci/software/pipelines/espresso/) (`espdr`), which
also reduces HARPS and HARPSpol data. Built on ESO's
[pyesorex](https://www.eso.org/sci/software/pyesorex/) framework and
[PyCPL](https://ivh.github.io/pycpl/), so it needs no changes to, and no build
of, the C pipeline.

## espdr_demod_pol

A HARPSpol observation is a cycle of 2 or 4 exposures at different retarder
angles, optionally repeated for signal to noise. The science recipe extracts the
two polarimetric beams into separate products, fibre A and fibre B, so a cycle
arrives as 4 or 8 S2D files. This recipe combines them with the ratio method.

**Input**
- SOF with either `S2D_A` / `S2D_B` or `S2D_BLAZE_A` / `S2D_BLAZE_B`, all
  exposures of one observing template. Mixing the two flavours is rejected.
- Parameters
  - `--null` (default: true): compute the null spectrum; ignored for
    2-exposure cycles
  - `--stokes` (default: AUTO): Stokes parameter, or derive it from the
    observing template

**Output**, sharing the prefix and timestamp of the first exposure:
- `_S2D_POL_I.fits` — intensity, the mean of all beams
- `_S2D_POL_STOKES.fits` — the Stokes parameter
- `_S2D_POL_NULL.fits` — null spectrum, 4-exposure cycles only

**Sequence ordering.** `espdr_sci_red` does not classify its products by
retarder angle, and the workflow is free to hand the frames over in any order,
so the recipe takes each exposure's place in the sequence from `ESO INS RET<n>
POS` in its own header. Sorting by angle gives the 90-degree pairs — (45, 135)
and (225, 315) — that the ratio method needs.

**Repeated cycles.** A longer template is the same cycle observed again, so the
exposures are walked in `ESO TPL EXPNO` order and a new cycle is started
whenever an angle comes round again. Each cycle is demodulated on its own — so
one taken in worse conditions cannot quietly drag the others, and each gets its
own null — and the Stokes and null spectra are then combined by inverse-variance
weighting. The intensity is the plain mean, as it already is within a cycle. The
number of cycles lands in `ESO QC POL NCYCLE`. All input frames must come from
one template (`ESO TPL ID` + `ESO TPL START`).

**Stokes parameter.** Taken from the observing template, `ESO TPL NAME`:
`HARPS_pol_obs_cir` means Stokes V. A linear template does not say whether it is
Q or U — that depends on the half-wave plate angles — so that needs `--stokes`.
If the template name is missing the recipe falls back on the retarder unit,
where RET25 is the quarter-wave plate and hence V.

**Order matching.** Fibre B is extracted with one order fewer than fibre A. The
orders are matched by their first wavelength rather than by a hard-coded index,
and the products are written on fibre B's order set with its wavelength
extensions passed through unchanged.

## Usage

```bash
uv sync
export PYESOREX_PLUGIN_DIR="$(pwd)/pyrecipes"
uv run pyesorex --recipes
uv run pyesorex --man-page espdr_demod_pol
uv run pyesorex espdr_demod_pol pol.sof
```

On macOS, if you have ESO pipelines installed, run under `env -u
DYLD_LIBRARY_PATH` so the bundled CPL is not shadowed by the system one.

There is also `demod.py` for use without pyesorex; it takes the S2D files
directly and works out the fibre from the filename.

```bash
uv run demod.py reduc/r.HARPS.2024-01-02T0[23]*_S2D_?.fits
```

## Tests

`tests/test_demod.py` reproduces the reference products in
`../112.25MG.001/reduc`, made by the pre-recipe version of this code. Point
`$HARPSPOL_TEST_DATA` elsewhere if the data lives somewhere else; the tests skip
if it is missing.

The 8-exposure path is tested against a synthetic template: the 29 CMa cycle
re-observed as exposures 5-8. Two identical cycles must give back the
single-cycle spectrum with the error smaller by sqrt(2), which is what the test
asserts. No real repeated-cycle data was available.

`tests/test_workflow.py` builds the merged EDPS workflow the way EDPS would and
checks the task wiring. It needs the pipeline's workflow directory; point
`$HARPS_WORKFLOW_DIR` at it if it is not in `~/pipes/harps-3.6.0/workflows`.

One deliberate difference: the old code's `error_helper` used a `dX/dR` twice
too large, so the Stokes and null **errors** of 4-exposure sequences were
overestimated by a factor 2. The tests assert the new errors are half the
reference ones there, and equal everywhere else. Fluxes are unchanged.

## EDPS workflow

`workflows/harpspol/harpspol_wkf.py` adds a `demod_pol` task to the HARPS
pipeline's own EDPS workflow (HARPS 3.6.0 and later, which is the first release
to ship one). It does not replace or fork it: importing `harps_wkf` makes EDPS
merge every task and data source of the pipeline workflow into ours, so only the
extra task and the product classifications it needs are defined here.

Add this repository's `workflows` directory to EDPS's `workflow_dir`, alongside
the pipeline's own:

```ini
workflow_dir = /path/to/share/esopipes/workflows, /path/to/harpspol.git/workflows
```

The task hangs off the pipeline's `object` task filtered to `S2D_A` / `S2D_B`,
exactly as the pipeline's own `combine_science` hangs off it filtered to
`S1D_FINAL_A`. `object` is already grouped on `TPL START` with a minimum group
size of 2, so one `demod_pol` job is one polarimetric template. The task is
conditional on `INS MODE == HARPSPOL`, which is the same switch
`harps_task_functions.should_combine` uses to *skip* HARPSpol in
`combine_science`.

Since EDPS runs one executable for the whole workflow, `esorex_path` has to
point at `pyesorex` for this to work. pyesorex does load C recipes, but they
must be built against the same CPL that PyCPL bundles — mixing CPL versions in
one process aborts with `CPL 7.x memory management subsystem is not
initialized`.

## Structure

- `pyrecipes/` — recipe files for pyesorex discovery, CPL interface only
- `pyespdr/` — the package with the actual work
- `workflows/harpspol/` — the EDPS workflow extension
- `old/` — the plugin for the old python-based HARPS reduction pipeline

## Still to settle with ESO

- `S2D_POL_I`, `S2D_POL_STOKES` and `S2D_POL_NULL` are made-up PRO.CATG values;
  they need registering for EDPS association and archive ingestion.
- The HARPS workflow does not classify the science `S2D_A` / `S2D_B` products,
  since nothing consumed them until now; `harpspol_wkf.py` adds those rules. They
  would sit more naturally in `harps_classification.py`.
- The task reads the un-corrected `S2D_A` / `S2D_B`. Whether the blaze-corrected
  pair is the better input is still open — the Stokes spectrum is unaffected,
  the intensity is not.
- Product provenance is written by hand in `_set_product_header`. Full DFS
  headers would come from `cpl.dfs` once the recipe is built against the
  pipeline proper.
- Why fibre B is extracted with one order fewer, and whether it is always the
  same order.
- `TEMPLATE_STOKES` only knows `HARPS_pol_obs_cir`. The linear template's name
  and its Q/U angle convention need filling in once there is linear data.
