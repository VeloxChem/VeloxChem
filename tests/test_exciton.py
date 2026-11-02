from mpi4py import MPI
from pathlib import Path
import numpy as np
import pytest

from veloxchem.veloxchemlib import OverlapDriver
from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import hartree_in_ev
from veloxchem.veloxchemlib import rotatory_strength_in_cgs
from veloxchem.veloxchemlib import get_dimer_ao_indices
from veloxchem.mpitask import MpiTask
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.excitondriver import ExcitonModelDriver
from veloxchem.errorhandler import VeloxChemError


@pytest.mark.solvers
class TestExciton:

    @staticmethod
    def assemble_matrices(ao_inds_1, ao_inds_2, s11, s12, s21, s22):

        n1 = len(ao_inds_1)
        n2 = len(ao_inds_2)
        smat = np.zeros((n1 + n2, n1 + n2))

        for row in range(n1):
            for col in range(n1):
                smat[ao_inds_1[row], ao_inds_1[col]] = s11[row, col]

        for row in range(n1):
            for col in range(n2):
                smat[ao_inds_1[row], ao_inds_2[col]] = s12[row, col]

        for row in range(n2):
            for col in range(n1):
                smat[ao_inds_2[row], ao_inds_1[col]] = s21[row, col]

        for row in range(n2):
            for col in range(n2):
                smat[ao_inds_2[row], ao_inds_2[col]] = s22[row, col]

        return smat

    @staticmethod
    def get_plot_results():

        return {
            'adiabatic_eigenvalues': np.array([0.20, 0.35]),
            'oscillator_strengths': np.array([0.10, 0.25]),
            'rotatory_strengths': np.array([0.01, -0.02]),
        }

    def test_assemble_matrices(self):

        here = Path(__file__).parent
        inpfile = here / 'data' / 'dimer.inp'
        outfile = inpfile.with_suffix('.out')

        task = MpiTask([str(inpfile), str(outfile)])

        molecule = task.molecule
        basis = task.ao_basis

        # build sub molecules

        atomlist_1 = list(range(0, 4))
        atomlist_2 = list(range(4, 9))

        mol_1 = molecule.slice(atomlist_1)
        mol_2 = molecule.slice(atomlist_2)

        basis_1 = basis.slice(atomlist_1)
        basis_2 = basis.slice(atomlist_2)

        # get indices of AOs from sub molecules

        ao_inds_1, ao_inds_2 = get_dimer_ao_indices(mol_1, mol_2, basis_1,
                                                    basis_2)

        # compute overlap

        if task.mpi_rank == mpi_master():

            ovl_drv = OverlapDriver()

            Smat = ovl_drv.compute(molecule, basis)

            S11mat = ovl_drv.compute(mol_1, basis_1)
            S12mat = ovl_drv.compute(mol_1, mol_2, basis_1, basis_2)
            S21mat = ovl_drv.compute(mol_2, mol_1, basis_2, basis_1)
            S22mat = ovl_drv.compute(mol_2, basis_2)

            S = Smat.full_matrix().to_numpy()

            S11 = S11mat.full_matrix().to_numpy()
            S12 = S12mat.full_matrix().to_numpy()
            S21 = S21mat.full_matrix().to_numpy()
            S22 = S22mat.full_matrix().to_numpy()

            S_exmod = self.assemble_matrices(ao_inds_1, ao_inds_2, S11, S12,
                                             S21, S22)

            assert np.max(np.abs(S - S_exmod)) < 1.0e-10

    @staticmethod
    def record_monomer_scf(exciton_drv, seen):

        # NOTE: the exciton model driver keeps no SCF results; those of the
        # monomers are recorded here so that the RI case can check them.
        monomer_scf = exciton_drv.monomer_scf

        def recording_monomer_scf(*args, **kwargs):
            scf_results = monomer_scf(*args, **kwargs)
            seen.append(scf_results)
            return scf_results

        exciton_drv.monomer_scf = recording_monomer_scf

    @staticmethod
    def check_ri_used(seen):

        assert seen
        for scf_results in seen:
            assert scf_results['ri_jk'] or scf_results['ri_coulomb']
            assert not scf_results['ri_fallback']

    def run_exciton_model(self, method_dict, ref_H, threshold, ri):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'exciton.inp')

        task = MpiTask([inpfile, None])
        exciton_dict = task.input_dict['exciton']
        # filename is necessary for checkpoint file
        exciton_dict['filename'] = task.input_dict['filename']

        exciton_drv = ExcitonModelDriver(task.mpi_comm, task.ostream)

        # fragments: 3
        # atoms_per_fragment: 3
        # nstates: 2
        # ct_nocc: 1
        # ct_nvir: 1

        # set attributes directly
        exciton_drv.fragments = '3'
        exciton_drv.atoms_per_fragment = '3'
        exciton_drv.nstates = 2
        exciton_drv.ct_nocc = 1
        exciton_drv.ct_nvir = 1
        exciton_drv.ri = ri
        exciton_drv.filename = exciton_dict['filename']

        monomer_scf_results = []
        self.record_monomer_scf(exciton_drv, monomer_scf_results)
        if 'xcfun' in method_dict:
            exciton_drv.xcfun = method_dict['xcfun']

        first_result = exciton_drv.compute(task.molecule, task.ao_basis)

        if task.mpi_rank == mpi_master():
            if ri:
                self.check_ri_used(monomer_scf_results)

            diag_diff = np.max(np.abs(np.diag(exciton_drv.H) - np.diag(ref_H)))
            abs_diff = np.max(np.abs(np.abs(exciton_drv.H) - np.abs(ref_H)))

            assert diag_diff < threshold
            assert abs_diff < threshold

            ref_eigvals, ref_eigvecs = np.linalg.eigh(ref_H)
            eigvals, eigvecs = np.linalg.eigh(exciton_drv.H)

            eigval_diff = np.max(np.abs(eigvals - ref_eigvals))
            assert eigval_diff < threshold

            backup_H = np.array(exciton_drv.H)

            np.testing.assert_allclose(first_result['hamiltonian'],
                                       exciton_drv.H)
            assert first_result['num_states'] == exciton_drv.H.shape[0]

        task.mpi_comm.barrier()

        # reset attributes and test update_settings
        exciton_drv.fragments = None
        exciton_drv.atoms_per_fragment = None
        exciton_drv.nstates = None
        exciton_drv.ct_nocc = None
        exciton_drv.ct_nvir = None
        exciton_drv.filename = None
        exciton_drv.xcfun = None
        exciton_drv.update_settings(exciton_dict, method_dict)

        exciton_drv.restart = True
        second_result = exciton_drv.compute(task.molecule, task.ao_basis)

        if task.mpi_rank == mpi_master():
            assert np.max(np.abs(backup_H - exciton_drv.H)) < 1.0e-10

            np.testing.assert_allclose(second_result['hamiltonian'],
                                       exciton_drv.H)
            assert second_result['num_states'] == first_result['num_states']

            exciton_h5 = Path(exciton_drv.checkpoint_file)

            for ind in range(len(exciton_drv.monomers)):

                final_h5 = Path(
                    str(exciton_h5).replace('_exciton.h5',
                                            f'_exciton_monomer_{ind + 1}.h5'))
                if final_h5.is_file():
                    final_h5.unlink()

                scf_h5 = Path(
                    str(exciton_h5).replace(
                        '_exciton.h5', f'_exciton_monomer_{ind + 1}_scf.h5'))
                if scf_h5.is_file():
                    scf_h5.unlink()

                rsp_h5 = Path(
                    str(exciton_h5).replace(
                        '_exciton.h5', f'_exciton_monomer_{ind + 1}_rsp.h5'))
                if rsp_h5.is_file():
                    rsp_h5.unlink()

            if exciton_h5.is_file():
                exciton_h5.unlink()

    @pytest.mark.parametrize('ri', [False, True])
    def test_exciton_model_rhf(self, ri):

        # vlxtag: RHF, Exciton_Model

        ref_H = {
            False: np.array([
                [
                    0.34427707, 0.00000000, 0.01212927, 0.00440677, 0.00045881,
                    0.00001154, -0.05006625, -0.14129907, -0.01428315,
                    0.00120021, 0.00000000, 0.00000000
                ],
                [
                    0.00000000, 0.41055831, -0.00268369, 0.00946559, 0.00000658,
                    -0.00002995, 0.00799698, 0.00104151, -0.00053603,
                    -0.00000206, 0.00000000, 0.00000000
                ],
                [
                    0.01212927, -0.00268369, 0.31644641, 0.00000000,
                    -0.01350857, -0.00365320, -0.13733722, -0.04998405,
                    0.00000000, 0.00000000, -0.05364967, 0.13888201
                ],
                [
                    0.00440677, 0.00946559, 0.00000000, 0.40065260, -0.00690413,
                    0.01135558, -0.02749723, -0.01394045, 0.00000000,
                    0.00000000, -0.01687523, 0.02684388
                ],
                [
                    0.00045881, 0.00000658, -0.01350857, -0.00690413,
                    0.34427707, 0.00000000, 0.00000000, 0.00000000, 0.00118706,
                    -0.01435363, 0.14447657, -0.05768827
                ],
                [
                    0.00001154, -0.00002995, -0.00365320, 0.01135558,
                    0.00000000, 0.41055831, 0.00000000, 0.00000000, 0.00000341,
                    -0.00065488, 0.00039640, -0.00841610
                ],
                [
                    -0.05006625, 0.00799698, -0.13733722, -0.02749723,
                    0.00000000, 0.00000000, 0.42953658, 0.00351461, 0.03426403,
                    0.00000000, 0.00000000, 0.00122229
                ],
                [
                    -0.14129907, 0.00104151, -0.04998405, -0.01394045,
                    0.00000000, 0.00000000, 0.00351461, 0.43767789, 0.00000000,
                    0.18627958, -0.00769894, 0.00000000
                ],
                [
                    -0.01428315, -0.00053603, 0.00000000, 0.00000000,
                    0.00118706, 0.00000341, 0.03426403, 0.00000000, 0.55734486,
                    0.00002348, -0.18322905, 0.00000000
                ],
                [
                    0.00120021, -0.00000206, 0.00000000, 0.00000000,
                    -0.01435363, -0.00065488, 0.00000000, 0.18627958,
                    0.00002348, 0.55770178, 0.00000000, 0.03043711
                ],
                [
                    0.00000000, 0.00000000, -0.05364967, -0.01687523,
                    0.14447657, 0.00039640, 0.00000000, -0.00769894,
                    -0.18322905, 0.00000000, 0.42274769, -0.00501119
                ],
                [
                    0.00000000, 0.00000000, 0.13888201, 0.02684388, -0.05768827,
                    -0.00841610, 0.00122229, 0.00000000, 0.00000000, 0.03043711,
                    -0.00501119, 0.41495367
                ],
            ]),
            True: np.array(
                [[
                    0.34427656, 0.00000000, 0.01212991, 0.00440704, -0.00045875,
                    -0.00001154, 0.05006975, 0.14130055, 0.01428380, 0.00120029,
                    0.00000000, 0.00000000
                ],
                 [
                     0.00000000, 0.41054536, -0.00268370, 0.00946663,
                     -0.00000658, 0.00002995, -0.00799449, -0.00104165,
                     0.00053606, -0.00000206, 0.00000000, 0.00000000
                 ],
                 [
                     0.01212991, -0.00268370, 0.31645038, 0.00000000,
                     0.01350912, 0.00365325, 0.13733578, 0.04998735, 0.00000000,
                     0.00000000, -0.05365336, -0.13887987
                 ],
                 [
                     0.00440704, 0.00946663, 0.00000000, 0.40064202, 0.00690441,
                     -0.01135659, 0.02751176, 0.01394163, 0.00000000,
                     0.00000000, -0.01687692, -0.02685831
                 ],
                 [
                     -0.00045875, -0.00000658, 0.01350912, 0.00690441,
                     0.34427656, 0.00000000, 0.00000000, 0.00000000, 0.00118714,
                     0.01435430, -0.14447751, -0.05769319
                 ],
                 [
                     -0.00001154, 0.00002995, 0.00365325, -0.01135659,
                     0.00000000, 0.41054536, 0.00000000, 0.00000000, 0.00000341,
                     0.00065493, -0.00039649, -0.00841326
                 ],
                 [
                     0.05006975, -0.00799449, 0.13733578, 0.02751176,
                     0.00000000, 0.00000000, 0.42953018, 0.00351597, 0.03426394,
                     0.00000000, 0.00000000, 0.00122239
                 ],
                 [
                     0.14130055, -0.00104165, 0.04998735, 0.01394163,
                     0.00000000, 0.00000000, 0.00351597, 0.43767284, 0.00000000,
                     -0.18628157, 0.00769924, 0.00000000
                 ],
                 [
                     0.01428380, 0.00053606, 0.00000000, 0.00000000, 0.00118714,
                     0.00000341, 0.03426394, 0.00000000, 0.55734066,
                     -0.00002348, 0.18323197, 0.00000000
                 ],
                 [
                     0.00120029, -0.00000206, 0.00000000, 0.00000000,
                     0.01435430, 0.00065493, 0.00000000, -0.18628157,
                     -0.00002348, 0.55769787, 0.00000000, -0.03043677
                 ],
                 [
                     0.00000000, 0.00000000, -0.05365336, -0.01687692,
                     -0.14447751, -0.00039649, 0.00000000, 0.00769924,
                     0.18323197, 0.00000000, 0.42274030, 0.00501271
                 ],
                 [
                     0.00000000, 0.00000000, -0.13887987, -0.02685831,
                     -0.05769319, -0.00841326, 0.00122239, 0.00000000,
                     0.00000000, -0.03043677, 0.00501271, 0.41494496
                 ]])
        }[ri]

        self.run_exciton_model({}, ref_H, 1.0e-6, ri)

    @pytest.mark.parametrize('ri', [False, True])
    def test_exciton_model_blyp(self, ri):

        # vlxtag: RKS, Exciton_Model

        ref_H = {
            False: np.array([
                [
                    0.26550384, 0.00000000, -0.00054035, -0.00007100,
                    0.00051242, 0.00001167, 0.02127844, -0.08555976, 0.00973521,
                    -0.00104548, 0.00000000, 0.00000000
                ],
                [
                    0.00000000, 0.33377631, 0.00003374, 0.00011106, -0.00001378,
                    0.00004700, 0.00471365, 0.00007118, -0.00037418, 0.00000237,
                    0.00000000, 0.00000000
                ],
                [
                    -0.00054035, 0.00003374, 0.24687169, 0.00000000, 0.00065576,
                    0.00001511, -0.08456819, 0.02188022, 0.00000000, 0.00000000,
                    -0.02468705, -0.08723805
                ],
                [
                    -0.00007100, 0.00011106, 0.00000000, 0.32511120, 0.00025808,
                    -0.00016134, -0.00174433, 0.00025917, 0.00000000,
                    0.00000000, -0.00011455, -0.00208295
                ],
                [
                    0.00051242, -0.00001378, 0.00065576, 0.00025808, 0.26550384,
                    0.00000000, 0.00000000, 0.00000000, -0.00103147, 0.00972945,
                    -0.08737873, -0.02491830
                ],
                [
                    0.00001167, 0.00004700, 0.00001511, -0.00016134, 0.00000000,
                    0.33377631, 0.00000000, 0.00000000, -0.00000088, 0.00041671,
                    -0.00001771, -0.00476247
                ],
                [
                    0.02127844, 0.00471365, -0.08456819, -0.00174433,
                    0.00000000, 0.00000000, 0.20903741, 0.00876112, -0.02327806,
                    0.00000000, 0.00000000, -0.00093677
                ],
                [
                    -0.08555976, 0.00007118, 0.02188022, 0.00025917, 0.00000000,
                    0.00000000, 0.00876112, 0.23442793, 0.00000000, -0.08476348,
                    0.00959230, 0.00000000
                ],
                [
                    0.00973521, -0.00037418, 0.00000000, 0.00000000,
                    -0.00103147, -0.00000088, -0.02327806, 0.00000000,
                    0.25032577, 0.00001953, -0.08250407, 0.00000000
                ],
                [
                    -0.00104548, 0.00000237, 0.00000000, 0.00000000, 0.00972945,
                    0.00041671, 0.00000000, -0.08476348, 0.00001953, 0.25234794,
                    0.00000000, 0.02018038
                ],
                [
                    0.00000000, 0.00000000, -0.02468705, -0.00011455,
                    -0.08737873, -0.00001771, 0.00000000, 0.00959230,
                    -0.08250407, 0.00000000, 0.21467743, -0.00882228
                ],
                [
                    0.00000000, 0.00000000, -0.08723805, -0.00208295,
                    -0.02491830, -0.00476247, -0.00093677, 0.00000000,
                    0.00000000, 0.02018038, -0.00882228, 0.21791698
                ],
            ]),
            True: np.array(
                [[
                    0.26543009, 0.00000000, -0.00054087, -0.00007121,
                    0.00051247, 0.00001171, 0.02127045, -0.08555833,
                    -0.00973726, 0.00104568, 0.00000000, 0.00000000
                ],
                 [
                     0.00000000, 0.33352648, -0.00003344, -0.00011186,
                     0.00001378, -0.00004697, -0.00471296, -0.00007439,
                     -0.00037431, 0.00000230, 0.00000000, 0.00000000
                 ],
                 [
                     -0.00054087, -0.00003344, 0.24684996, 0.00000000,
                     0.00065626, 0.00001487, -0.08456646, 0.02187197,
                     0.00000000, 0.00000000, 0.02467770, 0.08723521
                 ],
                 [
                     -0.00007121, -0.00011186, 0.00000000, 0.32481636,
                     0.00025849, -0.00016255, -0.00175019, 0.00027184,
                     0.00000000, 0.00000000, 0.00012471, 0.00208866
                 ],
                 [
                     0.00051247, 0.00001378, 0.00065626, 0.00025849, 0.26543009,
                     0.00000000, 0.00000000, 0.00000000, 0.00103167,
                     -0.00973150, 0.08737672, 0.02490917
                 ],
                 [
                     0.00001171, -0.00004697, 0.00001487, -0.00016255,
                     0.00000000, 0.33352648, 0.00000000, 0.00000000, 0.00000084,
                     -0.00041682, 0.00002737, 0.00475739
                 ],
                 [
                     0.02127045, -0.00471296, -0.08456646, -0.00175019,
                     0.00000000, 0.00000000, 0.20902250, 0.00876248, 0.02327158,
                     0.00000000, 0.00000000, 0.00093696
                 ],
                 [
                     -0.08555833, -0.00007439, 0.02187197, 0.00027184,
                     0.00000000, 0.00000000, 0.00876248, 0.23440618, 0.00000000,
                     0.08476145, -0.00959362, 0.00000000
                 ],
                 [
                     -0.00973726, -0.00037431, 0.00000000, 0.00000000,
                     0.00103167, 0.00000084, 0.02327158, 0.00000000, 0.25031052,
                     0.00001954, -0.08250274, 0.00000000
                 ],
                 [
                     0.00104568, 0.00000230, 0.00000000, 0.00000000,
                     -0.00973150, -0.00041682, 0.00000000, 0.08476145,
                     0.00001954, 0.25233280, 0.00000000, 0.02017458
                 ],
                 [
                     0.00000000, 0.00000000, 0.02467770, 0.00012471, 0.08737672,
                     0.00002737, 0.00000000, -0.00959362, -0.08250274,
                     0.00000000, 0.21465466, -0.00882329
                 ],
                 [
                     0.00000000, 0.00000000, 0.08723521, 0.00208866, 0.02490917,
                     0.00475739, 0.00093696, 0.00000000, 0.00000000, 0.02017458,
                     -0.00882329, 0.21790385
                 ]])
        }[ri]

        self.run_exciton_model({'xcfun': 'blyp'}, ref_H, 1.0e-5, ri)

    @pytest.mark.parametrize('ri', [False, True])
    def test_exciton_model_ecd(self, ri):

        mol_xyz = """12
        c2h4-dimer
        C         -1.37731        1.01769       -0.71611
        C         -0.04211        1.07142       -0.72602
        H         -1.96225        1.74636       -0.16458
        H         -1.90859        0.23094       -1.24174
        H          0.49049        1.84498       -0.18262
        H          0.54315        0.32947       -1.25941
        C         -1.17537       -1.48468        2.37427
        C          0.06813       -1.06658        2.62697
        H         -1.35657       -2.40378        1.82687
        H          0.92893       -1.63558        2.29127
        H         -2.03527       -0.90348        2.69157
        H          0.24803       -0.13578        3.15527
        """
        molecule = Molecule.read_xyz_string(mol_xyz)
        basis = MolecularBasis.read(molecule, 'def2-svp')

        exmod_drv = ExcitonModelDriver()

        exmod_drv.fragments = '2'
        exmod_drv.atoms_per_fragment = '6'
        exmod_drv.charges = '0'
        exmod_drv.nstates = 5
        exmod_drv.ct_nocc = 1
        exmod_drv.ct_nvir = 1
        exmod_drv.ri = ri

        monomer_scf_results = []
        self.record_monomer_scf(exmod_drv, monomer_scf_results)

        exmod_drv.ostream.mute()
        exmod_results = exmod_drv.compute(molecule, basis)

        if exmod_drv.rank == mpi_master():
            if ri:
                self.check_ri_used(monomer_scf_results)

            eigvals, eigvecs = np.linalg.eigh(exmod_drv.H)

            elec_trans_dipoles = np.matmul(eigvecs.T,
                                           exmod_drv.elec_trans_dipoles)
            velo_trans_dipoles = np.matmul(eigvecs.T,
                                           exmod_drv.velo_trans_dipoles)
            magn_trans_dipoles = np.matmul(eigvecs.T,
                                           exmod_drv.magn_trans_dipoles)

            excitation_energies = []
            oscillator_strengths = []
            rotatory_strengths = []

            for s in range(eigvals.size):
                ene = eigvals[s]
                dip_strength = np.sum(elec_trans_dipoles[s, :]**2)
                f = (2.0 / 3.0) * dip_strength * ene

                velo_trans_dipoles[s, :] /= (-ene)
                R = np.vdot(velo_trans_dipoles[s, :], magn_trans_dipoles[s, :])

                excitation_energies.append(ene)
                oscillator_strengths.append(f)
                rotatory_strengths.append(R)

            excitation_energies = np.array(excitation_energies)
            oscillator_strengths = np.array(oscillator_strengths)
            rotatory_strengths = np.array(rotatory_strengths)

            ref_ene = {
                False: np.array([
                    0.28932589, 0.31669640, 0.34117880, 0.34154634, 0.34448158,
                    0.34474205, 0.35545359, 0.35550180, 0.37509079, 0.37587732,
                    0.40570065, 0.41752625
                ]),
                True: np.array([
                    0.28934303, 0.31671347, 0.34118154, 0.34154941, 0.34449241,
                    0.34475279, 0.35546895, 0.35551709, 0.37511313, 0.37589993,
                    0.40570240, 0.41752876
                ])
            }[ri]
            ref_osc = {
                False: np.array([
                    0.031064, 1.241580, 0.000108, 0.013926, 0.000704, 0.000126,
                    0.000001, 0.000025, 0.000134, 0.001633, 0.024283, 0.005320
                ]),
                True: np.array([
                    0.0310628, 1.2415517, 0.0001073, 0.0138891, 0.0007042,
                    0.0001265, 1.099318e-06, 2.527011e-05, 0.0001339, 0.0016322,
                    0.0242953, 0.0053207
                ])
            }[ri]
            ref_rot = {
                False: np.array([
                    0.149686, -0.157138, -0.000080, 0.000033, 0.021327,
                    -0.001754, -0.000358, -0.002835, -0.010791, -0.007212,
                    0.001252, 0.016260
                ]),
                True: np.array([
                    0.1496247, -0.1570825, -8.001903e-05, 3.289724e-05,
                    0.0213306, -0.0017562, -0.0003583, -0.0028375, -0.0107963,
                    -0.0072060, 0.0012529, 0.0162616
                ])
            }[ri]

            assert np.max(np.abs(excitation_energies - ref_ene)) < 1.0e-6
            assert np.max(np.abs(oscillator_strengths - ref_osc)) < 1.0e-5
            assert np.max(np.abs(rotatory_strengths - ref_rot)) < 1.0e-5

            assert np.max(
                np.abs(exmod_results['excitation_energies_in_ev'] /
                       hartree_in_ev() - ref_ene)) < 1.0e-6
            assert np.max(
                np.abs(exmod_results['oscillator_strengths'] -
                       ref_osc)) < 1.0e-5
            assert np.max(
                np.abs(exmod_results['rotatory_strengths'] /
                       rotatory_strength_in_cgs() - ref_rot)) < 1.0e-5

    def test_get_spectrum_plot_dict(self):

        exciton_drv = ExcitonModelDriver()
        exciton_results = self.get_plot_results()

        plot_dict = exciton_drv._get_spectrum_plot_dict(exciton_results)

        np.testing.assert_allclose(plot_dict['eigenvalues'],
                                   exciton_results['adiabatic_eigenvalues'])
        np.testing.assert_allclose(plot_dict['oscillator_strengths'],
                                   exciton_results['oscillator_strengths'])
        np.testing.assert_allclose(plot_dict['rotatory_strengths'],
                                   exciton_results['rotatory_strengths'])

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_get_spectrum_plot_dict_missing_key(self):

        exciton_drv = ExcitonModelDriver()
        exciton_results = self.get_plot_results()
        exciton_results.pop('rotatory_strengths')

        with pytest.raises(VeloxChemError,
                           match='Missing key "rotatory_strengths"'):
            exciton_drv._get_spectrum_plot_dict(exciton_results)

    def test_plot_uv_vis(self):

        exciton_drv = ExcitonModelDriver()
        exciton_results = self.get_plot_results()
        plt = pytest.importorskip('matplotlib.pyplot')

        fig, ax = plt.subplots()
        try:
            exciton_drv.plot_uv_vis(exciton_results,
                                    broadening_type='gaussian',
                                    broadening_value=0.15,
                                    ax=ax)
        finally:
            plt.close(fig)

    def test_plot_ecd(self):

        exciton_drv = ExcitonModelDriver()
        exciton_results = self.get_plot_results()
        plt = pytest.importorskip('matplotlib.pyplot')

        fig, ax = plt.subplots()
        try:
            exciton_drv.plot_ecd(exciton_results,
                                 broadening_type='lorentzian',
                                 broadening_value=0.20,
                                 ax=ax)
        finally:
            plt.close(fig)
