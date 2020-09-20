import geometric

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
        - scf_drv: The SCF driver.
        - grad_drv: The gradient driver.
    """

    def __init__(self, molecule, ao_basis, min_basis, scf_drv, grad_drv):
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

        self.scf_drv = scf_drv
        self.grad_drv = grad_drv

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
        new_mol = Molecule(labels, coords.reshape(-1, 3), units='au')

        self.scf_drv.compute(new_mol, self.vlx_ao_basis, self.vlx_min_basis)
        energy = self.scf_drv.get_scf_energy()

        self.grad_drv.compute(new_mol, self.vlx_ao_basis, self.vlx_min_basis)
        gradient = self.grad_drv.get_gradient()

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
