import tempfile
import geometric

from .molecule import Molecule
from .scfrestdriver import ScfRestrictedDriver
from .gradientdriver import GradientDriver
from .optimizationengine import OptimizationEngine


class OptimizationDriver:
    """
    Implements optimization driver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - scf_drv: The SCF driver.
        - grad_drv: The gradient driver.
    """

    def __init__(self, comm, ostream):
        """
        Initializes optimization driver.
        """

        self.comm = comm
        self.ostream = ostream

        self.scf_drv = ScfRestrictedDriver(self.comm, self.ostream)
        self.grad_drv = GradientDriver(self.comm, self.ostream)

    def update_settings(self, scf_dict, method_dict=None):
        """
        Updates settings in optimization driver.

        :param scf_dict:
            The input dictionary of scf group.
        :param method_dict:
            The input dicitonary of method settings group.
        """

        self.scf_drv.update_settings(scf_dict, method_dict)
        self.grad_drv.update_settings(scf_dict, method_dict)

    def compute(self, molecule, ao_basis, min_basis=None):
        """
        Performs geometry optimization.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param min_basis:
            The minimal AO basis set.

        :return:
            The molecule with final geometry.
        """

        opt_engine = OptimizationEngine(molecule, ao_basis, self.scf_drv,
                                        self.grad_drv)

        tmpf = tempfile.mktemp()
        m = geometric.optimize.run_optimizer(customengine=opt_engine,
                                             check=1,
                                             input=tmpf)

        coords = m.xyzs[-1] / geometric.nifty.bohr2ang
        labels = molecule.get_labels()

        final_mol = Molecule(labels, coords.reshape(-1, 3), units='au')
        return final_mol
