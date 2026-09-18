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

**Resampling.** The two fibres are sampled about 0.31 pixels apart (0.0044 A
against a 0.0134 A pixel), so fibre B has to be put on fibre A's grid. The
Stokes ratio is interpolated as a ratio, which has the continuum divided out
and so interpolates better than either flux alone; the intensity needs the
fluxes themselves, and those go through a flux-conserving rebin -- differencing
a monotonic cubic interpolant of the cumulative flux, which is `espdr_rebin`'s
algorithm with `gsl_interp_cspline` swapped for PCHIP so it cannot overshoot
into negative flux. Both use PCHIP rather than linear interpolation: at HARPS's
3.18 px/FWHM sampling and this grid offset, linear broadens a line by 6.3% and
PCHIP by 2.4%. Pixels of fibre A's grid that fall outside fibre B's are masked
rather than extrapolated.

The ratio's *error* also comes from the aligned spectra, so it lands on the
same grid as the value it belongs to; every product and its error share a mask.

Exposure to exposure within one fibre the grids drift by only ~0.03 px (the
barycentric velocity moves ~30 m/s across a sequence), which is not worth a
further interpolation. Variance goes through the same operator as the flux,
which ignores the correlation resampling introduces between neighbouring
pixels -- as HDRL and `espdr_rebin` also do.

## Trying it out

Nothing here needs the C pipeline built or installed, so the recipe can be run
on existing S2D products on any machine.

**1. Get the code.** Needs Python 3.13 and [uv](https://docs.astral.sh/uv/).

```bash
git clone https://github.com/ivh/harps-pol
cd harps-pol
uv sync
```

**2. Get some HARPSpol S2D products.** Any 2- or 4-exposure polarimetric
template reduced with `espdr_sci_red` will do — what matters is the `S2D_A` and
`S2D_B` pair per exposure. The dataset this was developed against is ESO
programme **112.25MG.001**, night of **2024-01-02**, seven targets in
`INS MODE = HARPSPOL` (`DPR TYPE = STAR,CIRPOL,...`, template
`HARPS_pol_obs_all`): 29 CMa in a 4-exposure cycle, the rest in 2-exposure
cycles. Raw frames come from the ESO archive; run them through `espdr_sci_red`
as usual.

**3. Write a SOF** listing both fibres of every exposure in the template:

```
/data/r.HARPS.2024-01-02T00:36:59.492_S2D_A.fits  S2D_A
/data/r.HARPS.2024-01-02T00:36:59.492_S2D_B.fits  S2D_B
/data/r.HARPS.2024-01-02T00:47:32.940_S2D_A.fits  S2D_A
/data/r.HARPS.2024-01-02T00:47:32.940_S2D_B.fits  S2D_B
/data/r.HARPS.2024-01-02T00:58:05.829_S2D_A.fits  S2D_A
/data/r.HARPS.2024-01-02T00:58:05.829_S2D_B.fits  S2D_B
/data/r.HARPS.2024-01-02T01:08:38.837_S2D_A.fits  S2D_A
/data/r.HARPS.2024-01-02T01:08:38.837_S2D_B.fits  S2D_B
```

The order of the lines does not matter; the recipe sorts by retarder angle.

**4. Run it.**

```bash
uv run pyesorex --recipe-dir=pyrecipes --recipes
uv run pyesorex --recipe-dir=pyrecipes --man-page espdr_demod_pol
uv run pyesorex --recipe-dir=pyrecipes espdr_demod_pol pol.sof
```

`--recipe-dir` can be replaced by `export PYESOREX_PLUGIN_DIR="$(pwd)/pyrecipes"`.

On macOS, if ESO pipelines are installed, prefix with `env -u
DYLD_LIBRARY_PATH` so the bundled CPL is not shadowed by the system one.

### Using ESO's own PyCPL

`pyproject.toml` pulls PyCPL from <https://ivh.github.io/pycpl/simple/>, an
unofficial repackaging that bundles the C libraries so nothing has to be
installed first. The recipe itself uses only `cpl.core` and `cpl.ui`, so ESO's
own PyCPL should serve equally well — drop the `[tool.uv.sources]` block and the
`pycpl` index from `pyproject.toml` and let it resolve from the `eso` index that
is already listed. That path has not been tested here, since it needs a matching
CPL installation.

### Without pyesorex

`demod.py` does the same thing straight from the command line, working out the
fibre from the filename:

```bash
uv run demod.py /data/r.HARPS.2024-01-02T0[23]*_S2D_?.fits
```

## Tests

```bash
uv run pytest
```

`tests/test_demod.py` compares against `*_S2D_POL_*_LINEAR.fits`, the frozen
output of the pre-recipe version of this code, for all seven targets of
programme 112.25MG.001. Keeping those separate from the current
`*_S2D_POL_*.fits` is what stops the test comparing the code against its own
output. The data
is not in the repository: the tests look for a directory of `espdr_sci_red`
output next to it (`../112.25MG.001/reduc`) and skip if it is missing. Point
`$HARPSPOL_TEST_DATA` at your own copy — the timestamps the tests expect are
listed in `SEQUENCES` at the top of the file, and the reference
`_S2D_POL_*.fits` products have to sit alongside the `S2D_*` inputs. Without
them the remaining tests still run.

The 8-exposure path is tested against a synthetic template: the 29 CMa cycle
re-observed as exposures 5-8. Two identical cycles must give back the
single-cycle spectrum with the error smaller by sqrt(2), which is what the test
asserts. No real repeated-cycle data was available.

`tests/test_workflow.py` builds the merged EDPS workflow the way EDPS would and
checks the task wiring. It needs the pipeline's workflow directory; point
`$HARPS_WORKFLOW_DIR` at it if it is not in `~/pipes/harps-3.6.0/workflows`.

Three deliberate differences from the reference products:

- The old `error_helper` used a `dX/dR` twice too large, so the Stokes and null
  **errors** of 4-exposure sequences were overestimated by a factor 2. The
  tests assert the new errors are half the reference ones there, and equal
  everywhere else.
- The ratio is now interpolated with PCHIP rather than linearly.
- The intensity now co-adds the two fibres after aligning them.

The last two move the products by 2.4e-3 of peak rms at worst, so the reference
comparison is an rms tolerance rather than a bit-level check. The numerics
themselves are tested directly: flux conservation of the rebin, that only the
extrapolated order edges get masked, and that PCHIP does not broaden a shifted
line more than linear does.

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
