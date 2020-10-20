import geometric

from .veloxchemlib import mpi_master
from .molecule import Molecule


class OptimizationEngine(geometric.engine.Engine):
    """
    Implements optimization engine for geomeTRIC.

    :param molecule:
        The molecule.
    :param ao_basis:
        The AO basis set.
    :param min_basis:
        The minimal AO basis set.
    :param scf_drv:
        The SCF driver.
    :param grad_drv:
        The gradient driver.

    Instance variables
        - vlx_molecule: The molecule.
        - vlx_ao_basis: The AO basis set.
        - vlx_min_basis: The minimal AO basis set.
        - ener_drv: The energy driver.
        - grad_drv: The gradient driver.
        - flag: The type of the optimization engine.
        - comm: The MPI communicator.
        - rank: The rank of MPI process.
    """

    def __init__(self, molecule, ao_basis, min_basis, ener_drv, grad_drv, flag):
        """
        Initializes optimization engine for geomeTRIC.
        """

        g_molecule = geometric.molecule.Molecule()
        g_molecule.elem = molecule.get_labels()
        g_molecule.xyzs = [
            molecule.get_coordinates() * geometric.nifty.bohr2ang
        ]

        super().__init__(g_molecule)

        self.vlx_molecule = molecule
        self.vlx_ao_basis = ao_basis
        self.vlx_min_basis = min_basis

        self.ener_drv = ener_drv
        self.grad_drv = grad_drv
        self.flag = flag

        self.comm = grad_drv.comm
        self.rank = grad_drv.comm.Get_rank()

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
            new_mol.set_charge(self.vlx_molecule.get_charge())
            new_mol.set_multiplicity(self.vlx_molecule.get_multiplicity())
        else:
            new_mol = Molecule()
        new_mol.broadcast(self.rank, self.comm)

        if self.flag.upper() == 'XTB':
            self.ener_drv.compute(new_mol)
            energy = self.ener_drv.get_energy()

            self.grad_drv.compute(new_mol)
            gradient = self.grad_drv.get_gradient()

        elif self.flag.upper() == 'SCF':
            self.ener_drv.compute(new_mol, self.vlx_ao_basis,
                                  self.vlx_min_basis)
            energy = self.ener_drv.get_scf_energy()

            self.grad_drv.compute(new_mol, self.vlx_ao_basis,
                                  self.vlx_min_basis)
            gradient = self.grad_drv.get_gradient()

        energy = self.comm.bcast(energy, root=mpi_master())
        gradient = self.comm.bcast(gradient, root=mpi_master())

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
