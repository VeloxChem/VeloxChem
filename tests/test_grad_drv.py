import numpy as np
import pytest

from veloxchem.veloxchemlib import FockGeom1000Driver, T4CScreener
from veloxchem.veloxchemlib import make_matrix
from veloxchem.veloxchemlib import mpi_master, mat_t
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


@pytest.mark.solvers
class TestGradientDriver:

    def run_grad_drv(self, molecule, xcfun_label, basis_label):

        basis = MolecularBasis.read(molecule, basis_label)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_results = scf_drv.compute(molecule, basis)

        if scf_drv.rank == mpi_master():
            D = scf_results['D_alpha']
        else:
            D = None
        D = scf_drv.comm.bcast(D, root=mpi_master())

        den_mat_sym = make_matrix(basis, mat_t.symmetric)
        den_mat_sym.set_values(D)

        den_mat_gen = make_matrix(basis, mat_t.general)
        den_mat_gen.set_values(D)

        fock_grad_drv = FockGeom1000Driver()

        screener = T4CScreener()
        screener.partition(basis, molecule, 'eri')

        natoms = molecule.number_of_atoms()

        local_atoms = list(range(natoms))[scf_drv.rank::scf_drv.nodes]

        for iatom in local_atoms:

            screener_atom = T4CScreener()
            screener_atom.partition_atom(basis, molecule, 'eri', iatom)

            for fock_type in ['2jk', '2jkx', 'j', 'k', 'kx']:

                if fock_type in ['2jkx', 'kx']:
                    exchange_scaling_factor = 0.2
                else:
                    exchange_scaling_factor = 0.0

                thresh_int = 12

                gmats_sym = fock_grad_drv.compute(basis, screener_atom,
                                                  screener, den_mat_sym, iatom,
                                                  fock_type,
                                                  exchange_scaling_factor, 0.0,
                                                  thresh_int)

                gmats_gen = fock_grad_drv.compute(basis, screener_atom,
                                                  screener, den_mat_gen, iatom,
                                                  fock_type,
                                                  exchange_scaling_factor, 0.0,
                                                  thresh_int)

                for i, label in enumerate(['X', 'Y', 'Z']):
                    gmat_sym = gmats_sym.matrix_to_numpy(label)
                    gmat_gen = gmats_gen.matrix_to_numpy(label)

                    assert np.max(np.abs(gmat_sym - gmat_gen)) < 1.0e-10

    def test_nh3_svp(self):

        molstr = """
        N         -1.96309        1.59755       -0.01963
        H         -1.95876        2.61528        0.03109
        H         -2.48929        1.27814        0.79244
        H         -2.52930        1.35928       -0.83265
        """
        mol = Molecule.read_molecule_string(molstr, units='bohr')

        self.run_grad_drv(mol, 'hf', 'def2-svp')
        self.run_grad_drv(mol, 'slda', 'def2-svp')
        self.run_grad_drv(mol, 'b3lyp', 'def2-svp')
