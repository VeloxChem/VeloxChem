import pytest

from veloxchem.xtbdriver import XtbDriver

__all__ = ["using_cppe", "using_xtb", "using_mdanalysis"]


def pymodule_available(module):
    from importlib import util

    try:
        module_spec = util.find_spec(module, package=None)
    except ModuleNotFoundError:
        module_spec = None

    if module_spec is None:
        return False
    else:
        return True


using_cppe = pytest.mark.skipif(
    not pymodule_available("cppe"),
    reason="Not detecting module cppe. Install package if necessary",
)

using_xtb = pytest.mark.skipif(
    not XtbDriver.is_available(),
    reason="Not detecting xTB. Install package if necessary",
)

using_mdanalysis = pytest.mark.skipif(
    not pymodule_available("MDAnalysis"),
    reason="Not detecting MDAnalysis. Install package if necessary",
)
