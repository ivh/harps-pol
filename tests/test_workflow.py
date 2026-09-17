"""The EDPS workflow extension builds and carries the pipeline's own tasks.

Needs the HARPS pipeline's workflow directory, which is not part of this
repository; point ``$HARPS_WORKFLOW_DIR`` at it, otherwise these tests skip.
"""

import importlib
import os
import sys

import pytest

HARPS_WORKFLOWS = os.environ.get(
    "HARPS_WORKFLOW_DIR", os.path.expanduser("~/pipes/harps-3.6.0/workflows"))

pytestmark = pytest.mark.skipif(
    not os.path.isfile(os.path.join(HARPS_WORKFLOWS, "harps_wkf.py")),
    reason=f"HARPS workflow not found at {HARPS_WORKFLOWS}")

HERE = os.path.dirname(os.path.abspath(__file__))
OURS = os.path.join(os.path.dirname(HERE), "workflows")


@pytest.fixture
def workflow(tmp_path):
    """Build the merged workflow the way EDPS would.

    EDPS discovers ``<package>/<package>*_wkf.py`` and puts each package's
    parent directory on the path, so the pipeline's workflow has to be
    importable as ``harps``.
    """
    edps = pytest.importorskip("edps.generator.workflow_manager")
    (tmp_path / "harps").symlink_to(HARPS_WORKFLOWS)

    sys.path[:0] = [str(tmp_path), OURS]
    try:
        module = importlib.import_module("harpspol.harpspol_wkf")
        yield edps.WorkflowManager(config=None).create_workflow(module)
    finally:
        for path in (str(tmp_path), OURS):
            sys.path.remove(path)
        for name in [n for n in sys.modules
                     if n.startswith(("harps.", "harpspol.", "harps_wkf"))
                     or n in ("harps", "harpspol")]:
            del sys.modules[name]


def test_the_pipeline_tasks_come_along(workflow):
    names = {t.name for t in workflow.tasks}
    assert "demod_pol" in names
    # A sample of the HARPS pipeline's own tasks, inherited by importing its
    # workflow module rather than copied.
    assert {"object", "combine_science", "bias", "wave_thar_fp"} <= names


def test_demod_runs_on_the_science_products(workflow):
    demod = next(t for t in workflow.tasks if t.name == "demod_pol")
    assert demod.command == "espdr_demod_pol"
    assert demod.main_input.name == "object"
    assert sorted(demod.input_filter) == ["S2D_A", "S2D_B"]


def test_science_is_already_grouped_by_template(workflow):
    """Our per-template grouping is inherited, not re-declared."""
    science = next(t for t in workflow.tasks if t.name == "object")
    assert science.grouping_keywords == ["tpl.start"]
    assert science.min_group_size == 2
