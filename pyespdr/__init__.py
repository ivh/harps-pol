"""Python recipes and support code for the ESPRESSO/HARPS pipeline (espdr)."""

from pyespdr.demod import (
    PRO_CATG_I,
    PRO_CATG_NULL,
    PRO_CATG_STOKES,
    Spectrum,
    demodulate,
    demodulate_cycles,
    order_sequence,
    read_s2d,
    resample_spectrum,
    split_cycles,
    stokes_parameter,
    write_products,
)

__all__ = [
    "PRO_CATG_I",
    "PRO_CATG_NULL",
    "PRO_CATG_STOKES",
    "Spectrum",
    "demodulate",
    "demodulate_cycles",
    "order_sequence",
    "read_s2d",
    "resample_spectrum",
    "split_cycles",
    "stokes_parameter",
    "write_products",
]
