"""Polarimetric demodulation of S2D spectra produced by the espdr pipeline.

The espdr science recipe extracts the two polarimetric beams into separate
products (fibre A and fibre B).  A HARPSpol observation is a cycle of 2 or 4
exposures taken at different retarder angles, possibly repeated for signal to
noise; combining them with the ratio method gives the Stokes parameter, the
intensity, and -- for 4-exposure cycles -- a null spectrum.

The sequence position of each exposure is taken from the retarder angle in the
header, not from the filename or the exposure time, because the science recipe
does not classify its products by angle and the workflow is free to hand us the
frames in any order.
"""

from __future__ import annotations

import os
import re
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np
from astropy.io import fits
from scipy.interpolate import interp1d

PRO_CATG_I = "S2D_POL_I"
PRO_CATG_STOKES = "S2D_POL_STOKES"
PRO_CATG_NULL = "S2D_POL_NULL"

TAGS_A = ("S2D_A", "S2D_BLAZE_A")
TAGS_B = ("S2D_B", "S2D_BLAZE_B")

RECIPE_NAME = "espdr_demod_pol"

_RET_KEY = re.compile(r"^ESO INS RET(\d+) POS$")

#: Observing template suffix -> Stokes parameter.  The template says what the
#: observer set out to measure, which beats inferring it from the hardware.  A
#: linear template does not say whether it is Q or U -- that depends on the
#: half-wave plate angles -- so it has to be given explicitly.
TEMPLATE_STOKES = {"cir": "V"}

#: Fallback when the template name is missing: retarder unit -> Stokes
#: parameter.  RET25 is the quarter-wave plate (DPR.TECH = ECHELLE,CIRPOL in
#: the raw frames).
RETARDER_STOKES = {25: "V"}

#: Two orders are the same order if their first wavelength agrees to this many
#: Angstrom.  Adjacent echelle orders start ~40 A apart, so this is very loose.
ORDER_MATCH_TOL = 1.0

#: Retarder angles are compared modulo 360 with this tolerance, in degrees.
ANGLE_TOL = 1.0


@dataclass
class Spectrum:
    """One S2D file: the beam of one fibre at one retarder angle."""

    filename: str
    fiber: str
    angle: float
    retarder: int
    mjd: float
    wave: np.ndarray
    flux: np.ma.MaskedArray
    err: np.ma.MaskedArray
    header: fits.Header = field(repr=False)
    tpl_id: str = ""
    tpl_name: str = ""
    tpl_start: str = ""
    expno: int = 0

    @property
    def norders(self) -> int:
        return self.flux.shape[0]

    @property
    def template(self) -> tuple[str, str]:
        """Identity of the observing template this exposure belongs to."""
        return (self.tpl_id, self.tpl_start)


def retarder_angle(header: fits.Header) -> tuple[int, float]:
    """Return (retarder unit number, angle in degrees) from an S2D header."""
    found = [(int(m.group(1)), header[key]) for key in header
             if (m := _RET_KEY.match(str(key)))]
    if not found:
        raise ValueError("No 'ESO INS RET<n> POS' keyword found; "
                         "is this a polarimetric observation?")
    if len(found) > 1:
        raise ValueError(f"More than one retarder unit in header: "
                         f"{[n for n, _ in found]}")
    unit, angle = found[0]
    return unit, float(angle) % 360.0


def read_s2d(filename: str, fiber: str) -> Spectrum:
    """Load flux, error and barycentric vacuum wavelengths from an S2D file."""
    with fits.open(filename) as hdul:
        header = hdul[0].header.copy()
        flux = hdul["SCIDATA"].data
        wave = hdul["WAVEDATA_VAC_BARY"].data
        err = hdul["ERRDATA"].data

    mask = (flux == 0) | np.isnan(flux)
    retarder, angle = retarder_angle(header)
    return Spectrum(
        filename=filename,
        fiber=fiber,
        angle=angle,
        retarder=retarder,
        mjd=float(header["MJD-OBS"]),
        wave=wave,
        flux=np.ma.array(flux, mask=mask),
        err=np.ma.array(err, mask=mask),
        header=header,
        tpl_id=str(header.get("ESO TPL ID", "")),
        tpl_name=str(header.get("ESO TPL NAME", "")),
        tpl_start=str(header.get("ESO TPL START", "")),
        expno=int(header.get("ESO TPL EXPNO", 0)),
    )


def match_orders(wave_a: np.ndarray, wave_b: np.ndarray) -> np.ndarray:
    """Index array selecting the orders of A that also exist in B.

    Fibre B is extracted with one order fewer than fibre A.  Rather than
    hard-coding which one is missing, match the orders by their first
    wavelength; the wavelength scales are not rewritten by the demodulation, so
    they have to correspond anyway.
    """
    dist = np.abs(wave_b[:, 0][:, None] - wave_a[:, 0][None, :])
    idx = dist.argmin(axis=1)
    worst = dist[np.arange(len(idx)), idx].max()
    if worst > ORDER_MATCH_TOL:
        raise ValueError(f"Could not match orders between fibres: worst "
                         f"start-wavelength mismatch is {worst:.3f} A")
    if len(set(idx)) != len(idx):
        raise ValueError("Order matching between fibres is not one-to-one")
    return idx


def order_sequence(specs: Sequence[Spectrum]) -> list[Spectrum]:
    """Put the exposures of one fibre into demodulation order, by angle.

    The ratio method pairs consecutive exposures whose retarder angles differ by
    90 degrees: (45, 135) and, for a 4-exposure sequence, (225, 315).  Sorting
    by angle produces exactly that pairing.
    """
    n = len(specs)
    if n not in (2, 4):
        raise ValueError(f"Need 2 or 4 exposures per fibre, got {n}")

    units = {s.retarder for s in specs}
    if len(units) != 1:
        raise ValueError(f"Exposures use different retarder units: {units}")

    ordered = sorted(specs, key=lambda s: (s.angle, s.mjd))
    angles = [s.angle for s in ordered]
    if len(set(angles)) != n:
        raise ValueError(f"Repeated retarder angles in sequence: {angles}. "
                         f"Split the template into sub-sequences first.")

    for first, second in zip(ordered[0::2], ordered[1::2]):
        sep = (second.angle - first.angle) % 360.0
        if abs(sep - 90.0) > ANGLE_TOL:
            raise ValueError(
                f"Exposures at {first.angle:g} and {second.angle:g} deg are "
                f"not a 90-degree pair ({sep:g} deg apart)")
    return ordered


def stokes_parameter(specs: Sequence[Spectrum]) -> str:
    """Which Stokes parameter this sequence measures.

    Taken from the observing template (``HARPS_pol_obs_cir`` and friends), which
    records what the observer meant to do, and falls back on the retarder unit
    when the template name is missing.
    """
    name = specs[0].tpl_name
    suffix = name.rsplit("_", 1)[-1].lower() if name else ""
    if suffix in TEMPLATE_STOKES:
        return TEMPLATE_STOKES[suffix]
    if suffix == "lin":
        raise ValueError(
            f"Template {name!r} is linear polarimetry, which does not say "
            f"whether this is Stokes Q or U; give the Stokes parameter "
            f"explicitly.")

    unit = specs[0].retarder
    try:
        return RETARDER_STOKES[unit]
    except KeyError:
        raise ValueError(
            f"Template {name!r} is not a known polarimetric template and "
            f"retarder unit RET{unit} is unknown; give the Stokes parameter "
            f"explicitly.") from None


def split_cycles(specs: Sequence[Spectrum]) -> list[list[Spectrum]]:
    """Split one fibre's exposures into demodulation cycles.

    A longer template is just the 2- or 4-exposure cycle repeated for signal to
    noise, so walk the exposures in template order and start a new cycle
    whenever a retarder angle comes round again.  Each cycle is then put in
    angle order by :func:`order_sequence`.
    """
    if not specs:
        raise ValueError("No exposures to split into cycles")

    templates = {s.template for s in specs}
    if len(templates) > 1:
        raise ValueError(f"Exposures come from {len(templates)} different "
                         f"observing templates: {sorted(templates)}")

    cycles: list[list[Spectrum]] = []
    current: list[Spectrum] = []
    seen: set[float] = set()
    for spec in sorted(specs, key=lambda s: (s.expno, s.mjd)):
        if spec.angle in seen:
            cycles.append(current)
            current, seen = [], set()
        current.append(spec)
        seen.add(spec.angle)
    cycles.append(current)

    lengths = {len(c) for c in cycles}
    if len(lengths) > 1:
        raise ValueError(f"Template splits into cycles of unequal length: "
                         f"{[len(c) for c in cycles]}")
    return [order_sequence(c) for c in cycles]


def _ratio(seq: Sequence[Spectrum], invert: bool,
           swap: bool) -> np.ma.MaskedArray:
    """Product of the intra-pair flux ratios of one fibre.

    Fibre B is built the other way up than fibre A (``invert``) because the
    retarder swaps the two beams between the fibres.  ``swap`` inverts every
    pair but the first, which is what turns the Stokes demodulation into the
    null demodulation.
    """
    out = None
    for k, (first, second) in enumerate(zip(seq[0::2], seq[1::2])):
        flip = invert ^ (swap and k > 0)
        ratio = first.flux / second.flux if flip else second.flux / first.flux
        out = ratio if out is None else out * ratio
    return out


def _relative_error(specs: Sequence[Spectrum]) -> np.ma.MaskedArray:
    """Relative error of the ratio R.

    Every spectrum enters R exactly once, as a numerator or a denominator, so
    the relative errors simply add in quadrature regardless of which.
    """
    return np.sqrt(sum((s.err / s.flux) ** 2 for s in specs))


def demodulate(seq_a: Sequence[Spectrum], seq_b: Sequence[Spectrum],
               null: bool = True) -> dict[str, np.ndarray]:
    """Ratio-method demodulation of a 2- or 4-exposure polarimetric sequence.

    ``seq_a`` and ``seq_b`` must already be in angle order and matched to each
    other exposure by exposure.  Returns intensity, Stokes and (for a
    4-exposure sequence) null spectra with their errors, all on the wavelength
    scale of fibre B's orders.
    """
    npairs = len(seq_a) // 2
    if len(seq_b) != len(seq_a):
        raise ValueError("Fibre A and B sequences have different lengths")
    if null and npairs < 2:
        raise ValueError("A null spectrum needs 4 exposures")

    # Restrict fibre A to the orders that fibre B also has.
    keep = match_orders(seq_a[0].wave, seq_b[0].wave)
    seq_a = [_select_orders(s, keep) for s in seq_a]
    _check_wavelengths([s.wave for s in seq_a] + [s.wave for s in seq_b])

    wave_a = seq_a[0].wave
    wave_b = seq_b[0].wave
    allspec = [s for pair in zip(seq_a, seq_b) for s in pair]

    nspec = len(allspec)
    intensity = sum(s.flux for s in allspec) / nspec
    intensity_err = np.sqrt(sum(s.err**2 for s in allspec)) / nspec

    rel_err = _relative_error(allspec)
    out: dict[str, np.ndarray] = {
        "I": intensity,
        "I_ERR": intensity_err,
        "NEXP": len(seq_a),
    }

    exponent = 1.0 / (2 * npairs)
    for name, swap in (("STOKES", False), ("NULL", True)):
        if name == "NULL" and not null:
            continue
        ratio_a = _ratio(seq_a, invert=False, swap=swap)
        ratio_b = _resample_to(_ratio(seq_b, invert=True, swap=swap),
                               wave_b, wave_a)
        ratio = ratio_a * ratio_b

        root = ratio**exponent
        value = (root - 1.0) / (root + 1.0)
        # dX/dR for X = (R**p - 1)/(R**p + 1)
        dvalue = 2.0 * exponent * ratio ** (exponent - 1.0) / (root + 1.0) ** 2
        out[name] = value
        out[f"{name}_ERR"] = np.abs(dvalue) * ratio * rel_err

    return out


def demodulate_cycles(cycles_a: Sequence[Sequence[Spectrum]],
                      cycles_b: Sequence[Sequence[Spectrum]],
                      null: bool = True) -> dict[str, np.ndarray]:
    """Demodulate each cycle of a template separately, then combine them.

    Each cycle is demodulated on its own rather than folded into one long
    ratio, so that a cycle taken in worse conditions cannot quietly drag the
    others and so that every cycle gets its own null.  The Stokes and null
    spectra are then combined by inverse-variance weighting; the intensity is
    the plain mean, as it already is within a cycle.
    """
    if len(cycles_a) != len(cycles_b):
        raise ValueError(f"{len(cycles_a)} cycles in fibre A but "
                         f"{len(cycles_b)} in fibre B")

    per_cycle = [demodulate(a, b, null=null)
                 for a, b in zip(cycles_a, cycles_b)]
    if len(per_cycle) == 1:
        out = dict(per_cycle[0])
        out["NCYCLE"] = 1
        return out

    _check_wavelengths([c[0].wave for c in cycles_b])

    n = len(per_cycle)
    out: dict[str, np.ndarray] = {
        "I": sum(p["I"] for p in per_cycle) / n,
        "I_ERR": np.sqrt(sum(p["I_ERR"] ** 2 for p in per_cycle)) / n,
        "NEXP": sum(p["NEXP"] for p in per_cycle),
        "NCYCLE": n,
    }
    for key in ("STOKES", "NULL"):
        if key not in per_cycle[0]:
            continue
        weights = [1.0 / p[f"{key}_ERR"] ** 2 for p in per_cycle]
        total = sum(weights)
        out[key] = sum(w * p[key] for w, p in zip(weights, per_cycle)) / total
        out[f"{key}_ERR"] = 1.0 / np.sqrt(total)
    return out


def _select_orders(spec: Spectrum, keep: np.ndarray) -> Spectrum:
    return Spectrum(
        filename=spec.filename,
        fiber=spec.fiber,
        angle=spec.angle,
        retarder=spec.retarder,
        mjd=spec.mjd,
        wave=spec.wave[keep],
        flux=spec.flux[keep],
        err=spec.err[keep],
        header=spec.header,
        tpl_id=spec.tpl_id,
        tpl_name=spec.tpl_name,
        tpl_start=spec.tpl_start,
        expno=spec.expno,
    )


def _resample_to(values: np.ma.MaskedArray, wave_from: np.ndarray,
                 wave_to: np.ndarray) -> np.ma.MaskedArray:
    """Put a per-order quantity from one wavelength scale onto another."""
    out = values.copy()
    for i in range(out.shape[0]):
        out[i] = interp1d(wave_from[i], out[i].filled(np.nan),
                          fill_value="extrapolate", bounds_error=False)(wave_to[i])
    return out


def _check_wavelengths(waves: Sequence[np.ndarray], delta: float = 0.1) -> None:
    """All frames and fibres must share the same order start wavelengths."""
    reference = waves[0]
    for order in range(reference.shape[0]):
        starts = [w[order][0] for w in waves]
        if not np.allclose(starts, reference[order][0], atol=delta):
            raise ValueError(f"Order {order} wavelengths disagree between "
                             f"frames: {starts}")


def write_products(template: Spectrum, products: dict[str, np.ndarray],
                   inputs: Sequence[Spectrum], stokes: str,
                   outdir: str = ".", version: str = "0.1") -> list[tuple[str, str]]:
    """Write the demodulated spectra as S2D files, one per product.

    The first exposure of fibre B is used as the template, because fibre B
    carries the order set the products are on and the wavelength extensions are
    passed through unchanged.
    """
    base = os.path.basename(template.filename)
    for tag in TAGS_B:
        if base.endswith(f"_{tag}.fits"):
            suffix = tag
            break
    else:
        raise ValueError(f"Template {base!r} is not a fibre B S2D product")

    written: list[tuple[str, str]] = []
    todo = [(PRO_CATG_I, "I", "I_ERR"), (PRO_CATG_STOKES, "STOKES", "STOKES_ERR")]
    if "NULL" in products:
        todo.append((PRO_CATG_NULL, "NULL", "NULL_ERR"))

    for catg, flux_key, err_key in todo:
        outname = os.path.join(outdir, base.replace(f"_{suffix}.fits", f"_{catg}.fits"))
        with fits.open(template.filename) as hdul:
            hdul["SCIDATA"].data = _filled(products[flux_key])
            hdul["ERRDATA"].data = _filled(products[err_key])
            _set_product_header(hdul[0].header, catg, inputs, stokes,
                                version, int(products.get("NCYCLE", 1)))
            hdul.writeto(outname, overwrite=True)
        written.append((outname, catg))
    return written


def _filled(a) -> np.ndarray:
    return a.filled(np.nan) if isinstance(a, np.ma.MaskedArray) else a


def _next_rec_index(header: fits.Header) -> int:
    used = {int(m.group(1)) for key in header
            if (m := re.match(r"^ESO PRO REC(\d+) ID$", str(key)))}
    return max(used, default=0) + 1


def _set_product_header(header: fits.Header, catg: str,
                        inputs: Sequence[Spectrum], stokes: str,
                        version: str, ncycles: int = 1) -> None:
    """Stamp DFS product keywords on a header inherited from the template.

    Without this the products keep the template's ``ESO PRO CATG`` and claim to
    be S2D_B.  Full provenance would come from ``cpl.dfs`` once the recipe is
    built against the pipeline proper; what is set here is the subset that
    matters for classification and for tracing the inputs.
    """
    # An earlier version of this code wrote a literal 'PRO_CATG' keyword, which
    # is not the hierarchical DFS one and left the inherited S2D_B in place.
    if "PRO_CATG" in header:
        del header["PRO_CATG"]
    header["ESO PRO CATG"] = catg
    header["ESO PRO TECH"] = "ECHELLE,POLARIMETRY"
    header["ESO PRO TYPE"] = "REDUCED"

    rec = _next_rec_index(header)
    new = {
        f"ESO PRO REC{rec} ID": RECIPE_NAME,
        f"ESO PRO REC{rec} PIPE ID": f"pyespdr/{version}",
        f"ESO PRO REC{rec} PARAM1 NAME": "stokes",
        f"ESO PRO REC{rec} PARAM1 VALUE": stokes,
        f"ESO PRO REC{rec} PARAM2 NAME": "nexp",
        f"ESO PRO REC{rec} PARAM2 VALUE": str(len(inputs) // 2),
        "ESO QC POL NCYCLE": ncycles,
        "ESO QC POL ANGLES": ",".join(
            f"{s.angle:g}" for s in inputs if s.fiber == "A"),
        "ESO QC POL STOKES": stokes,
    }
    for i, spec in enumerate(inputs, start=1):
        new[f"ESO PRO REC{rec} RAW{i} NAME"] = os.path.basename(spec.filename)
        new[f"ESO PRO REC{rec} RAW{i} CATG"] = spec.header.get("ESO PRO CATG", "")

    for key, value in new.items():
        # The explicit prefix keeps astropy from warning about long keywords.
        header[f"HIERARCH {key}"] = value
