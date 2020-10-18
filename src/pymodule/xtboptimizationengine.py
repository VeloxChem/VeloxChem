import geometric

from .veloxchemlib import mpi_master
from .molecule import Molecule
from .xtbdriver import XTBDriver
from .xtbgradientdriver import XTBGradientDriver

class XTBOptimizationEngine(geometric.engine.Engine):
    """
    Implements optimization engine for geomeTRIC.

    :param molecule:
        The molecule. 

    Instance variables
        - vlx_molecule: The molecule.
        - vlx_ao_basis: The AO basis set.
        - vlx_min_basis: The minimal AO basis set.
        - scf_drv: The SCF driver.
        - grad_drv: The gradient driver.
        - comm: The MPI communicator.
        - rank: The rank of MPI process.
    """

    def __init__(self, comm, molecule, scf_dict, method_dict, ostream):
        """
        Initializes XTB optimization engine for geomeTRIC.
        """

        g_molecule = geometric.molecule.Molecule()
        g_molecule.elem = molecule.get_labels()
        g_molecule.xyzs = [
            molecule.get_coordinates() * geometric.nifty.bohr2ang
        ]

        super().__init__(g_molecule)
         
        self.comm = comm
        self.rank = comm.Get_rank()
        self.ostream = ostream

        self.vlx_molecule = molecule

        self.scf_dict = scf_dict 
        self.method_dict = method_dict

    def calc_new(self, coords, dirname):
        """
        Implements calc_new method for the engine.

        :param coords:
            The coordinates.
        :param dirname:
            The relative path.

        :return:
            A dictionary containing energy and gradient.
        """

        labels = self.vlx_molecule.get_labels()

        if self.rank == mpi_master():
            new_mol = Molecule(labels, coords.reshape(-1, 3), units='au')
        else:
            new_mol = Molecule()
        new_mol.broadcast(self.rank, self.comm)

        xtb_drv = XTBDriver(self.comm)
        xtb_drv.compute_energy(new_mol, self.scf_dict, self.method_dict,
                               self.ostream)
        energy = xtb_drv.get_energy()

        grad_drv = XTBGradientDriver(self.comm, xtb_drv, self.ostream)
        grad_drv.compute(new_mol)
        gradient = grad_drv.get_gradient()

        return {
            'energy': energy,
            'gradient': gradient.flatten(),
        }

    def copy_scratch(self, src, dest):
        """
        Implements copy_scratch method for the engine.

        :param src:
            The source.
        :param dest:
            The destination.
        """

        return
