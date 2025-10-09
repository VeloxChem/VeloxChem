import numpy as np

from veloxchem.veloxchemlib import T4CScreener
from veloxchem.veloxchemlib import make_matrix, mat_t
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.rijkfockdriver import RIJKFockDriver
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.fockdriver import FockDriver


class TestRIMODriver:

    def get_data_h2o(self):

        h2ostr = """
            O    0.000000000000        0.000000000000        0.000000000000
            H    0.000000000000        0.740848095288        0.582094932012
            H    0.000000000000       -0.740848095288        0.582094932012
        """

        mol = Molecule.read_str(h2ostr, 'angstrom')
        bas = MolecularBasis.read(mol, 'sto-3g')
        aux_bas = MolecularBasis.read(mol, 'def2-universal-jkfit')

        return mol, bas, aux_bas

    def test_h2o_compute_fock(self):

        mol_h2o, bas_sto3g, bas_aux = self.get_data_h2o()

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol_h2o, bas_sto3g)

        dmat = make_matrix(bas_sto3g, mat_t.symmetric)
        dmat.set_values(scf_results['D_alpha'])

        molorbs = scf_drv.molecular_orbitals

        ri_fock_drv = RIJKFockDriver()
        ri_fock_drv.compute_metric(mol_h2o, bas_aux, verbose=False)
        ri_fock_drv.compute_bq_vectors(mol_h2o,
                                       bas_sto3g,
                                       bas_aux,
                                       verbose=False)

        jmat = ri_fock_drv.compute_j_fock(dmat, 'j', verbose=False)
        kmat = ri_fock_drv.compute_k_fock(dmat, molorbs, verbose=False)

        # screen basis function pairs
        t4c_drv = T4CScreener()
        t4c_drv.partition(bas_sto3g, mol_h2o, "eri")

        # compute Fock matrix
        fock_drv = FockDriver()
        jmat_ref = fock_drv._compute_fock_omp(t4c_drv, dmat, "j", 0.0, 0.0, 15)
        kmat_ref = fock_drv._compute_fock_omp(t4c_drv, dmat, "k", 0.0, 0.0, 15)

        assert np.max(np.abs(kmat_ref.to_numpy() - kmat.to_numpy())) < 1e-3
        assert np.max(np.abs(jmat_ref.to_numpy() - jmat.to_numpy())) < 1e-3
