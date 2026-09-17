"""EDPS workflow adding polarimetric demodulation to the HARPS pipeline.

This does not replace the HARPS workflow, it extends it.  Importing
``harps_wkf`` makes EDPS merge every task and data source of the pipeline's own
workflow into this one (see ``edps.generator.workflow_manager.create_workflow``,
which recurses into any imported module whose name contains ``wkf``), so all
that is defined here is the one extra task and the product classifications it
needs.

To use it, add this directory's parent to EDPS's ``workflow_dir``, alongside the
pipeline's own:

    workflow_dir = /path/to/esopipes/workflows, /path/to/harpspol.git/workflows
"""

from edps import SCIENCE, JobParameters, classification_rule, task
from harps import harps_keywords as kwd
from harps import harps_wkf
from harps.harps_classification import harps
from harps.harps_task_functions import which_ins_mode

__title__ = "HARPSpol demodulation workflow"

# The HARPS workflow classifies the calibration S2Ds but not the science ones,
# since nothing downstream consumed them until now.
S2D_A = classification_rule("S2D_A", {**harps, kwd.pro_catg: "S2D_A"})
S2D_B = classification_rule("S2D_B", {**harps, kwd.pro_catg: "S2D_B"})
S2D_BLAZE_A = classification_rule("S2D_BLAZE_A", {**harps, kwd.pro_catg: "S2D_BLAZE_A"})
S2D_BLAZE_B = classification_rule("S2D_BLAZE_B", {**harps, kwd.pro_catg: "S2D_BLAZE_B"})

# Products of the demodulation.  These PRO.CATG values are not registered with
# ESO yet.
S2D_POL_I = classification_rule("S2D_POL_I", {**harps, kwd.pro_catg: "S2D_POL_I"})
S2D_POL_STOKES = classification_rule("S2D_POL_STOKES", {**harps, kwd.pro_catg: "S2D_POL_STOKES"})
S2D_POL_NULL = classification_rule("S2D_POL_NULL", {**harps, kwd.pro_catg: "S2D_POL_NULL"})


def is_harpspol(params: JobParameters) -> bool:
    return params.get_workflow_param("which_ins_mode") == "HARPSPOL"


# The counterpart of the pipeline's own combine_science, which explicitly skips
# HARPSPOL (see harps_task_functions.should_combine).  The science task is
# already grouped on TPL START with a minimum group size of 2, so one job here
# is one polarimetric template.
demod_pol = (task("demod_pol")
             .with_recipe("espdr_demod_pol")
             .with_condition(is_harpspol)
             .with_main_input(harps_wkf.science)
             .with_dynamic_parameter("which_ins_mode", which_ins_mode)
             .with_input_filter(S2D_A, S2D_B)
             .with_meta_targets([SCIENCE])
             .build())
