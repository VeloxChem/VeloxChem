import geometric

from .molecule import Molecule


class OptimizationEngine(geometric.engine.Engine):

    def __init__(self, molecule, basis, scf_drv, grad_drv):

        g_molecule = geometric.molecule.Molecule()
        g_molecule.elem = molecule.get_labels()
        g_molecule.xyzs = [
            molecule.get_coordinates() * geometric.nifty.bohr2ang
        ]

        super().__init__(g_molecule)

        self.vlx_molecule = molecule
        self.vlx_basis = basis

        self.scf_drv = scf_drv
        self.grad_drv = grad_drv

    def calc_new(self, coords, dirname):

        labels = self.vlx_molecule.get_labels()
        new_mol = Molecule(labels, coords.reshape(-1, 3), units='au')

        self.scf_drv.compute(new_mol, self.vlx_basis)
        energy = self.scf_drv.get_scf_energy()

        self.grad_drv.compute(new_mol, self.vlx_basis)
        gradient = self.grad_drv.get_gradient()

        return {
            'energy': energy,
            'gradient': gradient.flatten(),
        }
