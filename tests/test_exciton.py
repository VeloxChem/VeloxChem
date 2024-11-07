from pathlib import Path
import numpy as np
import pytest

from veloxchem.veloxchemlib import OverlapDriver
from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import get_dimer_ao_indices
from veloxchem.mpitask import MpiTask
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.excitondriver import ExcitonModelDriver


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

    def run_exciton_model(self, method_dict, ref_H, threshold):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'exciton.inp')

        task = MpiTask([inpfile, None])
        exciton_dict = task.input_dict['exciton']
        # filename is necessary for checkpoint file
        exciton_dict['filename'] = task.input_dict['filename']

        exciton_drv = ExcitonModelDriver(task.mpi_comm, task.ostream)
        exciton_drv.update_settings(exciton_dict, method_dict)
        exciton_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        if task.mpi_rank == mpi_master():

            diag_diff = np.max(np.abs(np.diag(exciton_drv.H) - np.diag(ref_H)))
            abs_diff = np.max(np.abs(np.abs(exciton_drv.H) - np.abs(ref_H)))

            assert diag_diff < threshold
            assert abs_diff < threshold

            ref_eigvals, ref_eigvecs = np.linalg.eigh(ref_H)
            eigvals, eigvecs = np.linalg.eigh(exciton_drv.H)

            eigval_diff = np.max(np.abs(eigvals - ref_eigvals))
            assert eigval_diff < threshold

            backup_H = np.array(exciton_drv.H)

        task.mpi_comm.barrier()

        exciton_drv.restart = True
        exciton_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        if task.mpi_rank == mpi_master():
            assert np.max(np.abs(backup_H - exciton_drv.H)) < 1.0e-10

            exciton_h5 = Path(exciton_drv.checkpoint_file)

            for ind in range(len(exciton_drv.monomers)):
                scf_h5 = exciton_h5.with_suffix(f'.monomer_{ind + 1}.scf.h5')
                if scf_h5.is_file():
                    scf_h5.unlink()
                scf_final_h5 = scf_h5.with_suffix('.results.h5')
                if scf_final_h5.is_file():
                    scf_final_h5.unlink()
                rsp_h5 = exciton_h5.with_suffix(f'.monomer_{ind + 1}.rsp.h5')
                if rsp_h5.is_file():
                    rsp_h5.unlink()
                rsp_final_h5 = rsp_h5.with_suffix('.solutions.h5')
                if rsp_final_h5.is_file():
                    rsp_final_h5.unlink()

            if exciton_h5.is_file():
                exciton_h5.unlink()

    def test_exciton_model_rhf(self):

        # vlxtag: RHF, Exciton_Model

        ref_H = np.array([
            [
                0.34427707, 0.00000000, 0.01212927, 0.00440677, 0.00045881,
                0.00001154, -0.05006625, -0.14129907, -0.01428315, 0.00120021,
                0.00000000, 0.00000000
            ],
            [
                0.00000000, 0.41055831, -0.00268369, 0.00946559, 0.00000658,
                -0.00002995, 0.00799698, 0.00104151, -0.00053603, -0.00000206,
                0.00000000, 0.00000000
            ],
            [
                0.01212927, -0.00268369, 0.31644641, 0.00000000, -0.01350857,
                -0.00365320, -0.13733722, -0.04998405, 0.00000000, 0.00000000,
                -0.05364967, 0.13888201
            ],
            [
                0.00440677, 0.00946559, 0.00000000, 0.40065260, -0.00690413,
                0.01135558, -0.02749723, -0.01394045, 0.00000000, 0.00000000,
                -0.01687523, 0.02684388
            ],
            [
                0.00045881, 0.00000658, -0.01350857, -0.00690413, 0.34427707,
                0.00000000, 0.00000000, 0.00000000, 0.00118706, -0.01435363,
                0.14447657, -0.05768827
            ],
            [
                0.00001154, -0.00002995, -0.00365320, 0.01135558, 0.00000000,
                0.41055831, 0.00000000, 0.00000000, 0.00000341, -0.00065488,
                0.00039640, -0.00841610
            ],
            [
                -0.05006625, 0.00799698, -0.13733722, -0.02749723, 0.00000000,
                0.00000000, 0.42953658, 0.00351461, 0.03426403, 0.00000000,
                0.00000000, 0.00122229
            ],
            [
                -0.14129907, 0.00104151, -0.04998405, -0.01394045, 0.00000000,
                0.00000000, 0.00351461, 0.43767789, 0.00000000, 0.18627958,
                -0.00769894, 0.00000000
            ],
            [
                -0.01428315, -0.00053603, 0.00000000, 0.00000000, 0.00118706,
                0.00000341, 0.03426403, 0.00000000, 0.55734486, 0.00002348,
                -0.18322905, 0.00000000
            ],
            [
                0.00120021, -0.00000206, 0.00000000, 0.00000000, -0.01435363,
                -0.00065488, 0.00000000, 0.18627958, 0.00002348, 0.55770178,
                0.00000000, 0.03043711
            ],
            [
                0.00000000, 0.00000000, -0.05364967, -0.01687523, 0.14447657,
                0.00039640, 0.00000000, -0.00769894, -0.18322905, 0.00000000,
                0.42274769, -0.00501119
            ],
            [
                0.00000000, 0.00000000, 0.13888201, 0.02684388, -0.05768827,
                -0.00841610, 0.00122229, 0.00000000, 0.00000000, 0.03043711,
                -0.00501119, 0.41495367
            ],
        ])

        self.run_exciton_model({}, ref_H, 1.0e-6)

    def test_exciton_model_blyp(self):

        # vlxtag: RKS, Exciton_Model

        ref_H = np.array([
            [
                0.26550384, 0.00000000, -0.00054035, -0.00007100, 0.00051242,
                0.00001167, 0.02127844, -0.08555976, 0.00973521, -0.00104548,
                0.00000000, 0.00000000
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
                -0.00016134, -0.00174433, 0.00025917, 0.00000000, 0.00000000,
                -0.00011455, -0.00208295
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
                0.02127844, 0.00471365, -0.08456819, -0.00174433, 0.00000000,
                0.00000000, 0.20903741, 0.00876112, -0.02327806, 0.00000000,
                0.00000000, -0.00093677
            ],
            [
                -0.08555976, 0.00007118, 0.02188022, 0.00025917, 0.00000000,
                0.00000000, 0.00876112, 0.23442793, 0.00000000, -0.08476348,
                0.00959230, 0.00000000
            ],
            [
                0.00973521, -0.00037418, 0.00000000, 0.00000000, -0.00103147,
                -0.00000088, -0.02327806, 0.00000000, 0.25032577, 0.00001953,
                -0.08250407, 0.00000000
            ],
            [
                -0.00104548, 0.00000237, 0.00000000, 0.00000000, 0.00972945,
                0.00041671, 0.00000000, -0.08476348, 0.00001953, 0.25234794,
                0.00000000, 0.02018038
            ],
            [
                0.00000000, 0.00000000, -0.02468705, -0.00011455, -0.08737873,
                -0.00001771, 0.00000000, 0.00959230, -0.08250407, 0.00000000,
                0.21467743, -0.00882228
            ],
            [
                0.00000000, 0.00000000, -0.08723805, -0.00208295, -0.02491830,
                -0.00476247, -0.00093677, 0.00000000, 0.00000000, 0.02018038,
                -0.00882228, 0.21791698
            ],
        ])

        self.run_exciton_model({'xcfun': 'blyp'}, ref_H, 1.0e-5)

    def test_exciton_model_ecd(self):

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

        exmod_settings = {
            'fragments': '2',
            'atoms_per_fragment': '6',
            'charges': '0',
            'nstates': '5',
            'ct_nocc': '1',
            'ct_nvir': '1',
        }

        method_settings = {}
        exmod_drv = ExcitonModelDriver()
        exmod_drv.update_settings(exmod_settings, method_settings)
        exmod_drv.ostream.mute()
        exmod_drv.compute(molecule, basis)

        if exmod_drv.rank == mpi_master():
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

            # TODO: get e, f and R from exmod_results

            excitation_energies = np.array(excitation_energies)
            oscillator_strengths = np.array(oscillator_strengths)
            rotatory_strengths = np.array(rotatory_strengths)

            ref_ene = np.array([
                0.28932589, 0.31669640, 0.34117880, 0.34154634, 0.34448158,
                0.34474205, 0.35545359, 0.35550180, 0.37509079, 0.37587732,
                0.40570065, 0.41752625
            ])
            ref_osc = np.array([
                0.031064, 1.241580, 0.000108, 0.013926, 0.000704, 0.000126,
                0.000001, 0.000025, 0.000134, 0.001633, 0.024283, 0.005320
            ])
            ref_rot = np.array([
                0.149686, -0.157138, -0.000080, 0.000033, 0.021327, -0.001754,
                -0.000358, -0.002835, -0.010791, -0.007212, 0.001252, 0.016260
            ])

            assert np.max(np.abs(excitation_energies - ref_ene)) < 1.0e-6
            assert np.max(np.abs(oscillator_strengths - ref_osc)) < 1.0e-5
            assert np.max(np.abs(rotatory_strengths - ref_rot)) < 1.0e-5
