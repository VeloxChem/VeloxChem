from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import DenseMatrix
from veloxchem.veloxchemlib import OverlapMatrix
from veloxchem.veloxchemlib import OverlapIntegralsDriver
from veloxchem.veloxchemlib import KineticEnergyMatrix
from veloxchem.veloxchemlib import KineticEnergyIntegralsDriver
from veloxchem.veloxchemlib import NuclearPotentialMatrix
from veloxchem.veloxchemlib import NuclearPotentialIntegralsDriver
from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.veloxchemlib import get_dimer_ao_indices
from veloxchem.mpitask import MpiTask
from veloxchem.excitondriver import ExcitonModelDriver


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
        inpfile = here / 'inputs' / 'dimer.inp'
        outfile = inpfile.with_suffix('.out')

        task = MpiTask([str(inpfile), str(outfile)])

        molecule = task.molecule
        basis = task.ao_basis

        # build sub molecules

        mol_1 = molecule.get_sub_molecule(0, 4)
        mol_2 = molecule.get_sub_molecule(4, 5)

        # get indices of AOs from sub molecules

        ao_inds_1, ao_inds_2 = get_dimer_ao_indices(mol_1, mol_2, basis, basis)

        # compute overlap

        ovldrv = OverlapIntegralsDriver(task.mpi_comm)
        S = ovldrv.compute(molecule, basis)
        S11 = ovldrv.compute(mol_1, basis)
        S22 = ovldrv.compute(mol_2, basis)
        S12 = ovldrv.compute(mol_1, mol_2, basis)
        S21 = ovldrv.compute(mol_2, mol_1, basis)

        if is_mpi_master(task.mpi_comm):

            smat = self.assemble_matrices(ao_inds_1, ao_inds_2, S11.to_numpy(),
                                          S12.to_numpy(), S21.to_numpy(),
                                          S22.to_numpy())
            S_exmod = OverlapMatrix(DenseMatrix(smat))
            assert S == S_exmod

        # compute kinetic energy

        kindrv = KineticEnergyIntegralsDriver(task.mpi_comm)
        T = kindrv.compute(molecule, basis)
        T11 = kindrv.compute(mol_1, basis)
        T22 = kindrv.compute(mol_2, basis)
        T12 = kindrv.compute(mol_1, mol_2, basis)
        T21 = kindrv.compute(mol_2, mol_1, basis)

        if is_mpi_master(task.mpi_comm):

            tmat = self.assemble_matrices(ao_inds_1, ao_inds_2, T11.to_numpy(),
                                          T12.to_numpy(), T21.to_numpy(),
                                          T22.to_numpy())
            T_exmod = KineticEnergyMatrix(DenseMatrix(tmat))
            assert T == T_exmod

        # compute nuclear potential

        npotdrv = NuclearPotentialIntegralsDriver(task.mpi_comm)
        V = npotdrv.compute(molecule, basis)
        V11 = npotdrv.compute(mol_1, basis, molecule)
        V22 = npotdrv.compute(mol_2, basis, molecule)
        V12 = npotdrv.compute(mol_1, mol_2, basis, molecule)
        V21 = npotdrv.compute(mol_2, mol_1, basis, molecule)

        if is_mpi_master(task.mpi_comm):

            vmat = self.assemble_matrices(ao_inds_1, ao_inds_2, V11.to_numpy(),
                                          V12.to_numpy(), V21.to_numpy(),
                                          V22.to_numpy())
            V_exmod = NuclearPotentialMatrix(DenseMatrix(vmat))
            assert V == V_exmod

    def run_exciton_model(self, method_dict, ref_H, threshold):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'exciton.inp')

        task = MpiTask([inpfile, None])
        exciton_dict = task.input_dict['exciton']
        # filename is necessary for checkpoint file
        exciton_dict['filename'] = task.input_dict['filename']

        exciton_drv = ExcitonModelDriver(task.mpi_comm, task.ostream)
        exciton_drv.update_settings(exciton_dict, method_dict)
        exciton_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        if is_mpi_master(task.mpi_comm):

            diag_diff = np.max(np.abs(np.diag(exciton_drv.H) - np.diag(ref_H)))
            abs_diff = np.max(np.abs(np.abs(exciton_drv.H) - np.abs(ref_H)))

            assert diag_diff < threshold
            assert abs_diff < threshold

            ref_eigvals, ref_eigvecs = np.linalg.eigh(ref_H)
            eigvals, eigvecs = np.linalg.eigh(exciton_drv.H)

            eigval_diff = np.max(np.abs(eigvals - ref_eigvals))
            assert eigval_diff < threshold

            backup_H = np.array(exciton_drv.H)

        exciton_drv.restart = True
        exciton_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        if is_mpi_master(task.mpi_comm):
            assert np.max(np.abs(backup_H - exciton_drv.H)) < 1.0e-10

            exciton_h5 = Path(exciton_drv.checkpoint_file)

            for ind in range(len(exciton_drv.monomers)):
                scf_h5 = exciton_h5.with_suffix(f'.monomer_{ind + 1}.scf.h5')
                if scf_h5.is_file():
                    scf_h5.unlink()
                scf_final_h5 = scf_h5.with_suffix('.tensors.h5')
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
