"""Regression tests against the products of the original demod.py script.

The reference products in ``112.25MG.001/reduc`` were made by the pre-recipe
version of this code, which resampled linearly and co-added the two fibres
without aligning them.  Both were changed deliberately, so the comparison is a
"same spectrum, known numerics change" check rather than a bit-level one; see
``test_agrees_with_reference_products`` for the tolerances and why.
"""

import dataclasses
import os

import numpy as np
import pytest
from astropy.io import fits

from pyespdr.demod import (
    _bin_edges,
    _rebin_conservative,
    _resample_to,
    _select_orders,
    demodulate,
    demodulate_cycles,
    match_orders,
    order_sequence,
    read_s2d,
    resample_spectrum,
    retarder_angle,
    split_cycles,
    stokes_parameter,
)

DATA = os.environ.get(
    "HARPSPOL_TEST_DATA",
    os.path.join(os.path.dirname(__file__), "..", "..", "112.25MG.001", "reduc"),
)
PREFIX = "r.HARPS.2024-01-02T"

# target -> (timestamps in acquisition order, S2D flavour used for the reference)
SEQUENCES = {
    "29 CMa": (["00:36:59.492", "00:47:32.940", "00:58:05.829", "01:08:38.837"], ""),
    "HD 54879": (["01:20:47.880", "01:56:20.986"], ""),
    "HD 61954": (["02:36:31.739", "03:22:04.861"], ""),
    "HD 92206": (["04:13:18.603", "04:58:51.855"], ""),
    "HD 62000": (["05:49:41.087", "06:35:13.848"], "BLAZE_"),
    "HD 101545": (["07:26:58.372", "07:47:31.865"], "BLAZE_"),
    "del Cir": (["08:12:11.782", "08:22:44.850"], "BLAZE_"),
}

pytestmark = pytest.mark.skipif(
    not os.path.isdir(DATA), reason=f"reference data not found at {DATA}")


def _path(stamp, suffix):
    return os.path.join(DATA, f"{PREFIX}{stamp}_{suffix}.fits")


def _load(stamps, flavour):
    seq_a = order_sequence(
        [read_s2d(_path(t, f"S2D_{flavour}A"), "A") for t in stamps])
    seq_b = order_sequence(
        [read_s2d(_path(t, f"S2D_{flavour}B"), "B") for t in stamps])
    return seq_a, seq_b


#: The deliberate numerics changes move the products by this much, as a
#: fraction of each spectrum's peak: PCHIP instead of linear interpolation of
#: the ratio, and fibre B resampled onto fibre A before co-adding the
#: intensity.  Measured worst case over the seven targets is 2.4e-3.
NUMERICS_RMS = 5e-3


@pytest.mark.parametrize("target", sorted(SEQUENCES))
def test_agrees_with_reference_products(target):
    stamps, flavour = SEQUENCES[target]
    seq_a, seq_b = _load(stamps, flavour)
    products = demodulate(seq_a, seq_b, null=len(stamps) == 4)

    for key, catg in (("I", "S2D_POL_I"), ("STOKES", "S2D_POL_STOKES"),
                      ("NULL", "S2D_POL_NULL")):
        if key not in products:
            continue
        with fits.open(_path(stamps[0], catg)) as hdul:
            ref_flux = hdul["SCIDATA"].data
            ref_err = hdul["ERRDATA"].data

        new_flux = products[key].filled(np.nan)
        rms = np.sqrt(np.nanmean((new_flux - ref_flux) ** 2))
        assert rms < NUMERICS_RMS * np.nanmax(np.abs(ref_flux)), \
            f"{target} {catg} flux"

        # The old 4-exposure code used dX/dR twice too large; everything else
        # must still come out at the same scale.
        expected = 0.5 if len(stamps) == 4 and key != "I" else 1.0
        new_err = products[key + "_ERR"].filled(np.nan)
        ratio = np.nanmedian(new_err / ref_err)
        assert ratio == pytest.approx(expected, rel=1e-3), \
            f"{target} {catg} error scale"


def test_conservative_rebin_conserves_flux():
    """Onto the same range, differencing the cumulative flux telescopes."""
    stamps, flavour = SEQUENCES["HD 54879"]
    _, seq_b = _load(stamps, flavour)
    spec = seq_b[0]
    for order in (0, 35, 69):
        edges = _bin_edges(spec.wave[order], spec.dll[order])
        flux = np.ma.getdata(spec.flux)[order]
        total = _rebin_conservative(flux, edges, np.array([edges[0], edges[-1]]))
        assert total.sum() == pytest.approx(flux.sum(), rel=1e-5)


def test_resampling_masks_the_extrapolated_edges():
    """Fibre A's grid runs past fibre B's, and those pixels are not invented."""
    stamps, flavour = SEQUENCES["HD 54879"]
    seq_a, seq_b = _load(stamps, flavour)
    keep = match_orders(seq_a[0].wave, seq_b[0].wave)
    target = _select_orders(seq_a[0], keep)

    before = np.ma.getmaskarray(seq_b[0].flux).sum()
    after = resample_spectrum(seq_b[0], target.wave, target.dll)
    mask = np.ma.getmaskarray(after.flux)
    assert mask.sum() > before
    # Only order edges, a handful of pixels out of 4096.
    assert mask.sum() < 3 * mask.shape[0]
    assert mask[:, 5:-5].sum() == 0


def test_pchip_broadens_a_shifted_line_less_than_linear():
    """The reason for the interpolant, at HARPS sampling and grid offset."""
    sampling, shift = 3.18, 0.31  # px per FWHM, px between the fibres' grids
    x = np.arange(81.0)
    sigma = sampling / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    line = np.exp(-0.5 * ((x - 40.0) / sigma) ** 2)

    def fwhm(y):
        half = y.max() / 2.0
        above = np.where(y > half)[0]
        return above[-1] - above[0]

    wave = np.array([x])
    shifted = np.array([x + shift])
    values = np.ma.array(np.array([line]), mask=np.zeros((1, 81), bool))
    pchip = _resample_to(values, wave, shifted)[0]
    linear = np.interp(x + shift, x, line)

    assert fwhm(np.asarray(pchip)) <= fwhm(linear)
    # and it stays positive, which a natural cubic spline need not
    assert np.all(np.asarray(pchip) > -1e-12)


def test_sequence_order_comes_from_the_angle():
    """Shuffling the input frames must not change the result."""
    stamps, flavour = SEQUENCES["29 CMa"]
    seq_a, seq_b = _load(stamps, flavour)
    assert [s.angle for s in seq_a] == [45.0, 135.0, 225.0, 315.0]

    shuffled_a = order_sequence([seq_a[i] for i in (2, 0, 3, 1)])
    shuffled_b = order_sequence([seq_b[i] for i in (3, 1, 2, 0)])
    assert [s.filename for s in shuffled_a] == [s.filename for s in seq_a]
    assert [s.filename for s in shuffled_b] == [s.filename for s in seq_b]


def test_rejects_a_sequence_that_is_not_ninety_degree_pairs():
    stamps, flavour = SEQUENCES["29 CMa"]
    seq_a, _ = _load(stamps, flavour)
    seq_a[1].angle = 200.0
    with pytest.raises(ValueError, match="90-degree pair"):
        order_sequence(seq_a)


def test_orders_are_matched_by_wavelength_not_by_index():
    stamps, flavour = SEQUENCES["HD 54879"]
    seq_a, seq_b = _load(stamps, flavour)
    keep = match_orders(seq_a[0].wave, seq_b[0].wave)
    assert len(keep) == seq_b[0].norders == 70
    assert seq_a[0].norders == 71
    # The order fibre B is missing is the one at ~5246.8 A, index 44 in A.
    assert set(range(71)) - set(keep) == {44}


def test_stokes_parameter_from_the_retarder_unit():
    stamps, flavour = SEQUENCES["HD 54879"]
    seq_a, _ = _load(stamps, flavour)
    assert retarder_angle(seq_a[0].header) == (25, 45.0)
    assert stokes_parameter(seq_a) == "V"


def _repeat_as_second_cycle(cycle):
    """Fake an 8-exposure template by re-observing the same cycle."""
    return list(cycle) + [
        dataclasses.replace(s, expno=s.expno + len(cycle)) for s in cycle]


def test_eight_exposure_template_splits_into_two_cycles():
    stamps, flavour = SEQUENCES["29 CMa"]
    seq_a, _ = _load(stamps, flavour)
    cycles = split_cycles(_repeat_as_second_cycle(seq_a))
    assert len(cycles) == 2
    for cycle in cycles:
        assert [s.angle for s in cycle] == [45.0, 135.0, 225.0, 315.0]


def test_repeated_cycle_averages_down_the_error():
    """Two identical cycles must give the same spectrum, sqrt(2) better."""
    stamps, flavour = SEQUENCES["29 CMa"]
    seq_a, seq_b = _load(stamps, flavour)
    single = demodulate_cycles([seq_a], [seq_b], null=True)
    doubled = demodulate_cycles(split_cycles(_repeat_as_second_cycle(seq_a)),
                                split_cycles(_repeat_as_second_cycle(seq_b)),
                                null=True)

    assert single["NCYCLE"] == 1
    assert doubled["NCYCLE"] == 2
    assert doubled["NEXP"] == 8

    for key in ("I", "STOKES", "NULL"):
        one = single[key].filled(np.nan)
        two = doubled[key].filled(np.nan)
        assert np.allclose(one, two, rtol=1e-10, equal_nan=True), key

        one_err = single[key + "_ERR"].filled(np.nan)
        two_err = doubled[key + "_ERR"].filled(np.nan)
        assert np.allclose(two_err, one_err / np.sqrt(2.0),
                           rtol=1e-10, equal_nan=True), key + "_ERR"


def test_split_cycles_rejects_two_templates():
    seq_a, _ = _load(SEQUENCES["29 CMa"][0], "")
    other, _ = _load(SEQUENCES["HD 54879"][0], "")
    with pytest.raises(ValueError, match="different\n?\\s*observing templates"):
        split_cycles(seq_a + other)


def test_stokes_comes_from_the_template_name():
    seq_a, _ = _load(SEQUENCES["HD 54879"][0], "")
    assert seq_a[0].tpl_name == "HARPS_pol_obs_cir"
    assert stokes_parameter(seq_a) == "V"


def test_linear_template_demands_an_explicit_stokes_parameter():
    seq_a, _ = _load(SEQUENCES["HD 54879"][0], "")
    linear = [dataclasses.replace(s, tpl_name="HARPS_pol_obs_lin")
              for s in seq_a]
    with pytest.raises(ValueError, match="Q or U"):
        stokes_parameter(linear)


def test_stokes_falls_back_to_the_retarder_unit():
    seq_a, _ = _load(SEQUENCES["HD 54879"][0], "")
    no_template = [dataclasses.replace(s, tpl_name="") for s in seq_a]
    assert stokes_parameter(no_template) == "V"


def test_errors_share_the_grid_and_mask_of_their_values():
    """The ratio's error comes from the aligned spectra, so it lines up."""
    stamps, flavour = SEQUENCES["29 CMa"]
    seq_a, seq_b = _load(stamps, flavour)
    products = demodulate(seq_a, seq_b, null=True)

    for key in ("I", "STOKES", "NULL"):
        value_mask = np.ma.getmaskarray(products[key])
        error_mask = np.ma.getmaskarray(products[key + "_ERR"])
        assert np.array_equal(value_mask, error_mask), key
        # Order edges fibre B does not reach, plus whatever the inputs
        # already had masked; either way a small fraction of the frame.
        assert value_mask.sum() < 0.01 * value_mask.size
