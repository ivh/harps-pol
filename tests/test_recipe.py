"""The recipe's parameter declaration, which EDPS depends on."""

import importlib.util
import os

import pytest

pytest.importorskip("cpl.ui")

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

