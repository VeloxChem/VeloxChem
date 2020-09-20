import os
import tempfile
import geometric

from .veloxchemlib import mpi_master
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
        - rank: The rank of MPI process.
        - coordsys: The coordinate system.
        - constraints: The constraints.
        - check_interval: The interval (number of steps) for checking
          coordinate system.
        - transition: The flag for transition state searching.
        - hessian: The flag for computing Hessian.
        - scf_drv: The SCF driver.
        - grad_drv: The gradient driver.
    """

    def __init__(self, comm, ostream):
        """
        Initializes optimization driver.
        """

        self.comm = comm
        self.rank = comm.Get_rank()
        self.ostream = ostream

        self.coordsys = 'tric'
        self.check_interval = 0
        self.constraints = None

        self.transition = False
        self.hessian = 'never'

        self.scf_drv = ScfRestrictedDriver(self.comm, self.ostream)
        self.grad_drv = GradientDriver(self.comm, self.ostream)

    def update_settings(self, opt_dict, scf_dict, method_dict=None):
        """
        Updates settings in optimization driver.

        :param opt_dict:
            The input dictionary of optimize group.
        :param scf_dict:
            The input dictionary of scf group.
        :param method_dict:
            The input dicitonary of method settings group.
        """

        if 'coordsys' in opt_dict:
            self.coordsys = opt_dict['coordsys'].lower()
        if 'check_interval' in opt_dict:
            self.check_interval = int(opt_dict['check_interval'])
        if 'constraints' in opt_dict:
            self.constraints = opt_dict['constraints']

        if 'transition' in opt_dict:
            key = opt_dict['transition'].lower()
            self.transition = True if key == 'yes' else False

        if 'hessian' in opt_dict:
            self.hessian = opt_dict['hessian'].lower()
        elif self.transition:
            self.hessian = 'first'

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

        opt_engine = OptimizationEngine(molecule, ao_basis, min_basis,
                                        self.scf_drv, self.grad_drv)

        if self.rank == mpi_master():
            suffix = '.scf.h5'
            if self.scf_drv.checkpoint_file[-len(suffix):] == suffix:
                temp_f = self.scf_drv.checkpoint_file[:-len(suffix)]
            else:
                temp_f = self.scf_drv.checkpoint_file

        with tempfile.TemporaryDirectory() as temp_d:

            if self.rank != mpi_master():
                temp_f = self.scf_drv.checkpoint_file
                temp_f = os.path.join(temp_d, 'tmp_{:d}'.format(self.rank))

            m = geometric.optimize.run_optimizer(customengine=opt_engine,
                                                 coordsys=self.coordsys,
                                                 check=self.check_interval,
                                                 constraints=self.constraints,
                                                 transition=self.transition,
                                                 hessian=self.hessian,
                                                 input=temp_f)

        coords = m.xyzs[-1] / geometric.nifty.bohr2ang
        labels = molecule.get_labels()

        return Molecule(labels, coords.reshape(-1, 3), units='au')
