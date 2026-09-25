"""PyEsoRex recipe: polarimetric demodulation of espdr S2D products."""

from typing import Any

import cpl.core
import cpl.ui

from pyespdr.demod import (
    TAGS_A,
    TAGS_B,
    demodulate_cycles,
    read_s2d,
    split_cycles,
    stokes_parameter,
    write_products,
)

VERSION = "0.1"
RECIPE = "espdr_demod_pol"

#: espdr names its parameters "espdr.<recipe>.<param>" -- see the
#: espdr.espdr_mflat.* entries in the pipeline's harps_parameters.yaml.  EDPS
#: sets them by that full name, so ours have to match, and pyesorex keys
#: ``settings`` by it too; the short form survives as the command-line alias.
PREFIX = f"espdr.{RECIPE}."

STOKES_CHOICES = ("AUTO", "I", "Q", "U", "V")


def _aliased(parameter, alias: str):
    parameter.cli_alias = alias
    parameter.cfg_alias = alias
    return parameter


def validate_stokes(value: str) -> str:
    """pyesorex does not enforce a ParameterEnum's alternatives itself."""
    if value not in STOKES_CHOICES:
        raise ValueError(
            f"--stokes={value!r} is not one of "
            f"{', '.join(STOKES_CHOICES)}")
    return value


class DemodPol(cpl.ui.PyRecipe):
    _name = RECIPE
    _version = VERSION
    _author = "Thomas Marquart"
    _email = "thomas.marquart@physics.uu.se"
    _copyright = "GPL-3.0-or-later"
    _synopsis = "Demodulate a HARPSpol sequence into Stokes, intensity and null"
    _description = (
        "Combines the fibre A and fibre B S2D spectra of a polarimetric\n"
        "sequence with the ratio method.  A template that repeats the 2- or\n"
        "4-exposure cycle is split on ESO TPL EXPNO, demodulated cycle by\n"
        "cycle, and combined by inverse-variance weighting.\n\n"
        "INPUT FRAMES\n"
        "  S2D_A / S2D_B              all exposures of one template, or\n"
        "  S2D_BLAZE_A / S2D_BLAZE_B  the blaze-corrected equivalent\n\n"
        "OUTPUT FRAMES\n"
        "  S2D_POL_I                  intensity\n"
        "  S2D_POL_STOKES             Stokes parameter\n"
        "  S2D_POL_NULL               null spectrum (4-exposure cycles only)\n\n"
        "The position of each exposure within the sequence is taken from the\n"
        "retarder angle (ESO INS RET<n> POS), not from the order in which the\n"
        "frames arrive, since the science recipe does not distinguish them."
    )

    def __init__(self):
        self.parameters = cpl.ui.ParameterList(
            [
                _aliased(cpl.ui.ParameterValue(
                    name=f"{PREFIX}null",
                    context=RECIPE,
                    description="Compute the null spectrum (4-exposure cycles only)",
                    default=True,
                ), "null"),
                _aliased(cpl.ui.ParameterEnum(
                    name=f"{PREFIX}stokes",
                    context=RECIPE,
                    description=(
                        "Stokes parameter, or AUTO to derive it from the "
                        "observing template"
                    ),
                    default="AUTO",
                    alternatives=list(STOKES_CHOICES),
                ), "stokes"),
            ]
        )

    def run(self, frameset: cpl.ui.FrameSet,
            settings: dict[str, Any]) -> cpl.ui.FrameSet:
        do_null = settings.get(f"{PREFIX}null", True)
        stokes = validate_stokes(settings.get(f"{PREFIX}stokes", "AUTO"))

        specs_a, specs_b, blaze = _sort_frames(frameset)
        cpl.core.Msg.info(
            self._name,
            f"{len(specs_a)} exposures, "
            f"{'blaze-corrected' if blaze else 'un-corrected'} S2D",
        )

        cycles_a = split_cycles(specs_a)
        cycles_b = split_cycles(specs_b)
        if len(cycles_a) != len(cycles_b):
            raise ValueError(f"{len(cycles_a)} cycles in fibre A but "
                             f"{len(cycles_b)} in fibre B")
        for cycle_a, cycle_b in zip(cycles_a, cycles_b):
            _check_pairing(cycle_a, cycle_b)

        template = cycles_a[0][0].tpl_name or cycles_a[0][0].tpl_id or "unknown"
        cpl.core.Msg.info(self._name, f"template {template}")
        for i, cycle in enumerate(cycles_a, start=1):
            cpl.core.Msg.info(
                self._name,
                f"cycle {i}/{len(cycles_a)}, retarder angles: "
                + ", ".join(f"{s.angle:g}" for s in cycle),
            )

        if stokes == "AUTO":
            stokes = stokes_parameter(cycles_a[0])
        cpl.core.Msg.info(self._name, f"Stokes parameter: {stokes}")

        nexp = len(cycles_a[0])
        want_null = do_null and nexp == 4
        if do_null and not want_null:
            cpl.core.Msg.warning(
                self._name, "No null spectrum from a 2-exposure cycle")

        products = demodulate_cycles(cycles_a, cycles_b, null=want_null)

        inputs = [s for cycle in zip(cycles_a, cycles_b)
                  for pair in zip(*cycle) for s in pair]
        # pyesorex collects products from the working directory.
        written = write_products(cycles_b[0][0], products, inputs, stokes,
                                 outdir=".", version=VERSION)

        out = cpl.ui.FrameSet()
        for filename, catg in written:
            cpl.core.Msg.info(self._name, f"wrote {filename}")
            out.append(cpl.ui.Frame(file=filename, tag=catg,
                                    group=cpl.ui.Frame.FrameGroup.PRODUCT))
        return out


def _sort_frames(frameset: cpl.ui.FrameSet):
    """Split the input frames by fibre, and check they are all one S2D flavour."""
    by_tag: dict[str, list] = {}
    for frame in frameset:
        by_tag.setdefault(frame.tag, []).append(frame.file)

    unknown = set(by_tag) - set(TAGS_A) - set(TAGS_B)
    if unknown:
        raise ValueError(f"Unexpected frame tags in input: {sorted(unknown)}")

    blaze = {"S2D_BLAZE_A", "S2D_BLAZE_B"} & set(by_tag)
    plain = {"S2D_A", "S2D_B"} & set(by_tag)
    if blaze and plain:
        raise ValueError("Input mixes blaze-corrected and un-corrected S2D "
                         "frames; give one flavour only")
    if not blaze and not plain:
        raise ValueError("No S2D frames in input")

    tag_a, tag_b = ("S2D_BLAZE_A", "S2D_BLAZE_B") if blaze else ("S2D_A", "S2D_B")
    files_a, files_b = by_tag.get(tag_a, []), by_tag.get(tag_b, [])
    if len(files_a) != len(files_b):
        raise ValueError(f"{len(files_a)} frames tagged {tag_a} but "
                         f"{len(files_b)} tagged {tag_b}")

    return ([read_s2d(f, "A") for f in files_a],
            [read_s2d(f, "B") for f in files_b],
            bool(blaze))


def _check_pairing(seq_a, seq_b) -> None:
    """Each fibre A exposure must have a fibre B exposure from the same frame."""
    for a, b in zip(seq_a, seq_b):
        if abs(a.angle - b.angle) > 1e-6 or abs(a.mjd - b.mjd) > 1e-9:
            raise ValueError(
                f"Fibre A and B frames do not pair up: "
                f"{a.filename} (angle {a.angle:g}, MJD {a.mjd}) vs "
                f"{b.filename} (angle {b.angle:g}, MJD {b.mjd})")
