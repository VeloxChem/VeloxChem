import numpy as np
import math

from veloxchem import MolecularBasis
from veloxchem import Molecule
from veloxchem import FockGeom1000Driver
from veloxchem import RIFockDriver
from veloxchem import TwoCenterElectronRepulsionDriver
from veloxchem import T4CScreener
from veloxchem import SubMatrix
from veloxchem import make_matrix
from veloxchem import mat_t
from veloxchem.rigradientdriver import RIFockGradDriver
from veloxchem.scfrestdriver import ScfRestrictedDriver


class TestRIJFockGeomExGradDriver:

    def get_data_h2o(self):

        h2ostr = """
           O   1.206898105612   1.227024188072  -0.113458607566
           H   0.490554304233   1.138628097468   0.508036765414
           H   1.992568601678   1.153431301285   0.406410355109
        """

        mol = Molecule.read_str(h2ostr, 'angstrom')
        bas = MolecularBasis.read(mol, 'def2-svp')
        aux_bas = MolecularBasis.read(mol, 'def2-universal-jfit')

        return mol, bas, aux_bas

    def test_h2o_ri_j_grad_screened(self):

        mol, bas, bas_aux = self.get_data_h2o()

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_res = scf_drv.compute(mol, bas)
        density = scf_res['D_alpha']

        mo = scf_res['C_alpha']
        nocc = mol.number_of_alpha_electrons()
        mo_occ = mo[:, :nocc]
        mo_vir = mo[:, nocc:]
        nvir = mo_vir.shape[1]

        exc_mat = np.zeros((nocc, nvir))
        exc_mat[nocc - 1, 0] = 1.0
        exc_density = np.linalg.multi_dot([mo_occ, exc_mat, mo_vir.T])
        sym_exc_density = 0.5 * (exc_density + exc_density.T)

        gs_den_mat = make_matrix(bas, mat_t.symmetric)
        gs_den_mat.set_values(density)

        rw_den_sym_mat = make_matrix(bas, mat_t.symmetric)
        rw_den_sym_mat.set_values(sym_exc_density)

        rw_den_gen_mat = make_matrix(bas, mat_t.general)
        rw_den_gen_mat.set_values(exc_density)

        fock_grad_drv = FockGeom1000Driver()

        ket_t4c = T4CScreener()
        ket_t4c.partition(bas, mol, 'eri')

        # compute J metric
        t2c_drv = TwoCenterElectronRepulsionDriver()
        matj = t2c_drv.compute(mol, bas_aux)
        rmatj = np.linalg.inv(matj.full_matrix().to_numpy())
        invmatj = SubMatrix([0, 0, rmatj.shape[0], rmatj.shape[0]])
        invmatj.set_values(rmatj)

        ri_fock_drv = RIFockDriver(invmatj)
        ri_fock_drv.prepare_buffers(mol, bas, bas_aux)

        gs_bq = ri_fock_drv.compute_bq_vector(gs_den_mat)
        rw_bq = ri_fock_drv.compute_bq_vector(rw_den_sym_mat)

        ri_grad_drv = RIFockGradDriver()

        for iatom in range(mol.number_of_atoms()):

            bra_t4c = T4CScreener()
            bra_t4c.partition_atom(bas, mol, 'eri', iatom)

            rgrad = fock_grad_drv.compute(bas, bra_t4c, ket_t4c, gs_den_mat,
                                          rw_den_gen_mat, iatom, 'j', 0.0, 0.0,
                                          12)

            cgrad = ri_grad_drv.direct_compute(ket_t4c, bas, bas_aux, mol,
                                               rw_bq, gs_bq, rw_den_sym_mat,
                                               gs_den_mat, iatom, 12)
            cg_xyz = cgrad.coordinates()

            for i in range(3):
                assert math.isclose(rgrad[i],
                                    cg_xyz[i],
                                    rel_tol=1.0e-4,
                                    abs_tol=1.0e-4)
