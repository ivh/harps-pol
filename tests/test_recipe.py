"""The recipe: its parameter declaration, which EDPS depends on, and runs of it."""

import importlib.util
import os

import numpy as np
import pytest
from astropy.io import fits

cpl = pytest.importorskip("cpl")
import cpl.ui  # noqa: E402

from pyespdr.demod import demodulate_cycles, read_s2d, split_cycles  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RECIPE = os.path.join(os.path.dirname(HERE), "pyrecipes", "espdr_demod_pol.py")


def _module():
    spec = importlib.util.spec_from_file_location("espdr_demod_pol", RECIPE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_parameters_follow_the_espdr_naming_convention():
    """EDPS writes a recipe config keyed by the full CPL name.

    Every other espdr parameter is ``espdr.<recipe>.<param>`` -- see
    ``harps_parameters.yaml`` -- and pyesorex keys its ``settings`` dict the
    same way, so a bare name is unreachable from the workflow.
    """
    module = _module()
    names = {p.name for p in module.DemodPol().parameters}
    assert names == {"espdr.espdr_demod_pol.null"}


def test_the_command_line_keeps_the_short_names():
    module = _module()
    aliases = {p.cli_alias for p in module.DemodPol().parameters}
    assert aliases == {"null"}



DATA = os.environ.get(
    "HARPSPOL_TEST_DATA",
    os.path.join(HERE, "..", "..", "112.25MG.001", "reduc"),
)
needs_data = pytest.mark.skipif(
    not os.path.isdir(DATA), reason=f"reference data not found at {DATA}")

STAMPS_4 = ["00:36:59.492", "00:47:32.940", "00:58:05.829", "01:08:38.837"]
STAMPS_2 = ["01:20:47.880", "01:56:20.986"]


def _s2d(stamp, flavour, fiber, root=DATA):
    return os.path.join(root, f"r.HARPS.2024-01-02T{stamp}_S2D_{flavour}{fiber}.fits")


def _frameset(stamps, flavour="", root=DATA):
    frames = cpl.ui.FrameSet()
    for stamp in stamps:
        for fiber in "AB":
            frames.append(cpl.ui.Frame(_s2d(stamp, flavour, fiber, root),
                                       tag=f"S2D_{flavour}{fiber}"))
    return frames


def _run(frames, settings=None):
    """Run the recipe; products land in the working directory."""
    out = _module().DemodPol().run(frames, settings or {})
    return {(fits.getval(f.file, "ESO QC POL STOKES"), f.tag): f.file
            for f in out}


def _library(stamps, flavour="", root=DATA):
    specs_a = [read_s2d(_s2d(t, flavour, "A", root), "A") for t in stamps]
    specs_b = [read_s2d(_s2d(t, flavour, "B", root), "B") for t in stamps]
    return demodulate_cycles(split_cycles(specs_a), split_cycles(specs_b),
                             null=len(stamps) == 4)


def _assert_matches(product, expected, key):
    data = fits.getdata(product, "SCIDATA")
    err = fits.getdata(product, "ERRDATA")
    assert np.allclose(data, expected[key].filled(np.nan), equal_nan=True)
    assert np.allclose(err, expected[f"{key}_ERR"].filled(np.nan),
                       equal_nan=True)


@needs_data
def test_four_exposure_circular_run_gives_v_with_null(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    products = _run(_frameset(STAMPS_4))
    assert set(products) == {("V", "S2D_POL_I"), ("V", "S2D_POL_STOKES"),
                             ("V", "S2D_POL_NULL")}

    expected = _library(STAMPS_4)
    for key, catg in (("I", "S2D_POL_I"), ("STOKES", "S2D_POL_STOKES"),
                      ("NULL", "S2D_POL_NULL")):
        _assert_matches(products["V", catg], expected, key)

    header = fits.getheader(products["V", "S2D_POL_STOKES"])
    assert header["ESO PRO CATG"] == "S2D_POL_STOKES"
    assert header["ESO QC POL ANGLES"] == "45,135,225,315"
    assert header["ESO QC POL NCYCLE"] == 1


@needs_data
def test_two_exposure_run_has_no_null(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    products = _run(_frameset(STAMPS_2))
    assert set(products) == {("V", "S2D_POL_I"), ("V", "S2D_POL_STOKES")}
    _assert_matches(products["V", "S2D_POL_STOKES"], _library(STAMPS_2),
                    "STOKES")


@needs_data
def test_null_can_be_switched_off(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    products = _run(_frameset(STAMPS_4),
                    {"espdr.espdr_demod_pol.null": False})
    assert ("V", "S2D_POL_NULL") not in products


@needs_data
def test_blaze_input_gives_the_same_stokes(tmp_path, monkeypatch):
    """The blaze cancels in the ratio, so only the intensity may differ.

    Up to float32 rounding of the inputs: measured 1.0e-7 at worst, 1.6e-5 of
    the error.  Absolute, since Stokes crosses zero.
    """
    monkeypatch.chdir(tmp_path)
    plain = fits.getdata(_run(_frameset(STAMPS_2))["V", "S2D_POL_STOKES"],
                         "SCIDATA")
    blaze = _run(_frameset(STAMPS_2, "BLAZE_"))
    assert blaze["V", "S2D_POL_STOKES"].endswith("_S2D_POL_STOKES.fits")
    assert np.allclose(fits.getdata(blaze["V", "S2D_POL_STOKES"], "SCIDATA"),
                       plain, rtol=0, atol=1e-6, equal_nan=True)


@needs_data
def test_mixed_s2d_flavours_are_rejected(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    frames = _frameset(STAMPS_2)
    for frame in _frameset(STAMPS_2, "BLAZE_"):
        frames.append(frame)
    with pytest.raises(ValueError, match="mixes blaze-corrected"):
        _run(frames)


def _fake_linear(root, angles):
    """Copy the 4-exposure circular sequence, relabelled as half-wave plate data."""
    for stamp, angle in zip(STAMPS_4, angles):
        for fiber in "AB":
            with fits.open(_s2d(stamp, "", fiber)) as hdul:
                header = hdul[0].header
                del header["ESO INS RET25 POS"]
                header["HIERARCH ESO INS RET50 POS"] = angle
                header["ESO TPL NAME"] = "HARPS_pol_obs_lin"
                hdul.writeto(_s2d(stamp, "", fiber, root))


@needs_data
def test_interleaved_linear_template_gives_q_and_u(tmp_path, monkeypatch):
    raw = tmp_path / "in"
    raw.mkdir()
    _fake_linear(raw, [0.0, 22.5, 45.0, 67.5])
    monkeypatch.chdir(tmp_path)

    products = _run(_frameset(STAMPS_4, root=raw))
    assert set(products) == {(s, c) for s in "QU"
                             for c in ("S2D_POL_I", "S2D_POL_STOKES")}
    # Each set is named after its own first exposure, so they cannot collide.
    assert len(set(products.values())) == 4

    for stokes, stamps, angles in (("Q", STAMPS_4[0::2], "0,45"),
                                   ("U", STAMPS_4[1::2], "22.5,67.5")):
        product = products[stokes, "S2D_POL_STOKES"]
        assert os.path.basename(product).startswith(
            f"r.HARPS.2024-01-02T{stamps[0]}")
        assert fits.getval(product, "ESO QC POL ANGLES") == angles
        _assert_matches(product, _library(stamps, root=raw), "STOKES")


@needs_data
def test_half_wave_plate_at_circular_angles_is_rejected(tmp_path, monkeypatch):
    """45/135/225/315 on the half-wave plate never swaps the beams."""
    raw = tmp_path / "in"
    raw.mkdir()
    _fake_linear(raw, [45.0, 135.0, 225.0, 315.0])
    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError, match="as many of each"):
        _run(_frameset(STAMPS_4, root=raw))
