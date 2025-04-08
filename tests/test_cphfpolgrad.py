from pathlib import Path
import numpy as np
import pytest
import h5py
import sys

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lrsolver import LinearResponseSolver
from veloxchem.cppsolver import ComplexResponse
from veloxchem.polorbitalresponse import PolOrbitalResponse


@pytest.mark.solvers
class TestCphfPolgrad:

    def run_cphfpolgrad_cg_real(self, molecule, basis, xcfun=None, label1=None, label2=None):
        # conjugate gradient solver
        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = xcfun
        scf_drv.ostream.mute()
        scf_tensors = scf_drv.compute(molecule, basis)

        # linear response
        rsp_settings = {'conv_thresh': 1.0e-5, 'frequencies': (0.0, 0.4)}
        lr_drv = LinearResponseSolver()
        lr_drv.a_operator = "electric dipole"
        lr_drv.b_operator = "electric dipole"
        lr_drv.update_settings(rsp_settings)
        lr_drv.ostream.mute()
        lr_results = lr_drv.compute(molecule, basis, scf_tensors)

        # test real analytical gradient
        cphfpolgrad_solver = PolOrbitalResponse()
        cphfpolgrad_settings = {'conv_thresh':2e-7, 'frequencies': (0.0, 0.4),
                                'use_subspace_solver': 'no'}
        cphfpolgrad_solver.update_settings(cphfpolgrad_settings)
        cphfpolgrad_solver.ostream.mute()
        cphfpolgrad_solver.compute(molecule, basis, scf_tensors, lr_results)
        cphfpolgrad_solver.compute_omega(molecule, basis, scf_tensors, lr_results)

        # get CPHF omega coefficients from dist. array
        dist_omega_coefficients = cphfpolgrad_solver.cphf_results['dist_omega_ao']
        dof = len(dist_omega_coefficients)

        cphfpolgrad_omega_coefficients = []
        for x in range(dof):
            solution_vec = dist_omega_coefficients[x].get_full_vector(0)
            if scf_drv.rank == mpi_master():
                cphfpolgrad_omega_coefficients.append(solution_vec)

        if scf_drv.rank == mpi_master():
            cphfpolgrad_results = cphfpolgrad_solver.cphf_results
            cphfpolgrad_coefficients = cphfpolgrad_results['cphf_ov']
            #cphfpolgrad_omega_coefficients = np.array([cphfpolgrad_results[0.0]['omega_ao'], 
            #                                          cphfpolgrad_results[0.4]['omega_ao']])
            np.set_printoptions(suppress=True, precision=10)
            here = Path(__file__).parent
            hf_file_name = str(here /'data'/'cphfpolgrad_coefficients.h5')
            hf = h5py.File(hf_file_name, 'r')
            cphfpolgrad_reference = np.array(hf.get(label1))
            cphfpolgrad_omega_reference = np.array(hf.get(label2))
            hf.close()

            # Here we are comparing the CPHF coefficients in MO basis, so
            # there might be sign differences; we compare absolute values instead.
            assert np.max(np.abs(cphfpolgrad_coefficients) - np.abs(cphfpolgrad_reference)) < 1.0e-6
            assert np.max(np.abs(cphfpolgrad_omega_coefficients)
                                   - np.abs(cphfpolgrad_omega_reference)) < 1.0e-6 

    def run_cphfpolgrad_subspace_real(self, molecule, basis, xcfun=None, label1=None, label2=None):
        # subspace solver
        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = xcfun
        scf_drv.ostream.mute()
        scf_tensors = scf_drv.compute(molecule, basis)

        # linear response
        rsp_settings = {'conv_thresh': 1.0e-5, 'frequencies': (0.0, 0.4)}
        lr_drv = LinearResponseSolver()
        lr_drv.a_operator = "electric dipole"
        lr_drv.b_operator = "electric dipole"
        lr_drv.update_settings(rsp_settings)
        lr_drv.ostream.mute()
        lr_results = lr_drv.compute(molecule, basis, scf_tensors)

        # test real analytical gradient
        cphfpolgrad_solver = PolOrbitalResponse()
        cphfpolgrad_settings = {'conv_thresh':2e-7, 'frequencies': (0.0, 0.4),
                                'use_subspace_solver': 'yes'}
        cphfpolgrad_solver.update_settings(cphfpolgrad_settings)
        cphfpolgrad_solver.ostream.mute()
        cphfpolgrad_solver.compute(molecule, basis, scf_tensors, lr_results)
        cphfpolgrad_solver.compute_omega(molecule, basis, scf_tensors, lr_results)

        # get CPHF lambda coefficients from dist array
        dist_cphf_coefficients = cphfpolgrad_solver.cphf_results['dist_cphf_ov']
        dof = len(dist_cphf_coefficients)

        cphfpolgrad_coefficients = []
        for x in range(dof):
            solution_vec = dist_cphf_coefficients[x].get_full_vector(0)
            if scf_drv.rank == mpi_master():
                cphfpolgrad_coefficients.append(solution_vec)

        # get CPHF omega coefficients from dist. array
        dist_omega_coefficients = cphfpolgrad_solver.cphf_results['dist_omega_ao']
        dof = len(dist_omega_coefficients)

        cphfpolgrad_omega_coefficients = []
        for x in range(dof):
            solution_vec = dist_omega_coefficients[x].get_full_vector(0)
            if scf_drv.rank == mpi_master():
                cphfpolgrad_omega_coefficients.append(solution_vec)

        if scf_drv.rank == mpi_master():
            cphfpolgrad_results = cphfpolgrad_solver.cphf_results
            #cphfpolgrad_omega_coefficients = np.array([cphfpolgrad_results[0.0]['omega_ao'], 
            #                                          cphfpolgrad_results[0.4]['omega_ao']])
            np.set_printoptions(suppress=True, precision=10)
            here = Path(__file__).parent
            hf_file_name = str(here /'data'/'cphfpolgrad_coefficients.h5')
            hf = h5py.File(hf_file_name, 'r')
            cphfpolgrad_reference = np.array(hf.get(label1))
            cphfpolgrad_omega_reference = np.array(hf.get(label2))
            hf.close()

            cphfpolgrad_coefficients = np.array(cphfpolgrad_coefficients)
            cphfpolgrad_coefficients = cphfpolgrad_coefficients.reshape(cphfpolgrad_reference.shape)

            # Here we are comparing the CPHF coefficients in MO basis, so
            # there might be sign differences; we compare absolute values instead.
            assert np.max(np.abs(cphfpolgrad_coefficients) - np.abs(cphfpolgrad_reference)) < 1.0e-6
            assert np.max(np.abs(cphfpolgrad_omega_coefficients)
                                   - np.abs(cphfpolgrad_omega_reference)) < 1.0e-6 

    def run_cphfpolgrad_cg_complex(self, molecule, basis, xcfun=None, label1=None, label2=None):
        # conjugate gradient solver
        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = xcfun
        scf_drv.ostream.mute()
        scf_tensors = scf_drv.compute(molecule, basis)

        # linear response
        rsp_settings = {'conv_thresh': 1.0e-5, 'frequencies': (0.0, 0.4),
                        'damping': 0.5}
        lr_drv = ComplexResponse()
        lr_drv.a_operator = "electric dipole"
        lr_drv.b_operator = "electric dipole"
        lr_drv.update_settings(rsp_settings)
        lr_drv.ostream.mute()
        lr_results = lr_drv.compute(molecule, basis, scf_tensors)

        # test complex analytical gradient
        cphfpolgrad_solver = PolOrbitalResponse()
        cphfpolgrad_settings = {'conv_thresh':2e-7, 'frequencies': (0.0, 0.4),
                                'is_complex': 'yes', 'damping': 0.5,
                                'use_subspace_solver': 'no'}
        cphfpolgrad_solver.update_settings(cphfpolgrad_settings)
        cphfpolgrad_solver.ostream.mute()
        cphfpolgrad_solver.compute(molecule, basis, scf_tensors, lr_results)
        cphfpolgrad_solver.compute_omega(molecule, basis, scf_tensors, lr_results)

        # get CPHF omega coefficients from dist. array
        dist_omega_coefficients = cphfpolgrad_solver.cphf_results['dist_omega_ao']
        dof = len(dist_omega_coefficients)

        cphfpolgrad_omega_coefficients = []
        for x in range(dof):
            solution_vec = dist_omega_coefficients[x].get_full_vector(0)
            if scf_drv.rank == mpi_master():
                cphfpolgrad_omega_coefficients.append(solution_vec)

        if scf_drv.rank == mpi_master():
            cphfpolgrad_results = cphfpolgrad_solver.cphf_results
            cphfpolgrad_coefficients = cphfpolgrad_results['cphf_ov']
            #cphfpolgrad_omega_coefficients = np.array([cphfpolgrad_results[0.0]['omega_ao'], 
            #                                          cphfpolgrad_results[0.4]['omega_ao']])
            np.set_printoptions(suppress=True, precision=10)
            here = Path(__file__).parent
            hf_file_name = str(here /'data'/'cphfpolgrad_coefficients.h5')
            hf = h5py.File(hf_file_name, 'r')
            cphfpolgrad_reference = np.array(hf.get(label1))[12:]
            cphfpolgrad_omega_reference = np.array(hf.get(label2))[12:]
            hf.close()

            # Here we ar comparing the CPHF coefficients in MO basis, so
            # there might be sign differences; we compare absolute values instead.
            assert np.max(np.abs(cphfpolgrad_coefficients) - np.abs(cphfpolgrad_reference)) < 1.0e-6
            assert np.max(np.abs(cphfpolgrad_omega_coefficients)
                                   - np.abs(cphfpolgrad_omega_reference)) < 1.0e-6 

    def run_cphfpolgrad_subspace_complex(self, molecule, basis, xcfun=None, label1=None, label2=None):
        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = xcfun
        scf_drv.ostream.mute()
        scf_tensors = scf_drv.compute(molecule, basis)

        # linear response
        rsp_settings = {'conv_thresh': 1.0e-5, 'frequencies': (0.0, 0.4),
                        'damping': 0.5}
        lr_drv = ComplexResponse()
        lr_drv.a_operator = "electric dipole"
        lr_drv.b_operator = "electric dipole"
        lr_drv.update_settings(rsp_settings)
        lr_drv.ostream.mute()
        lr_results = lr_drv.compute(molecule, basis, scf_tensors)

        # test complex analytical gradient
        cphfpolgrad_solver = PolOrbitalResponse()
        cphfpolgrad_settings = {'conv_thresh':2e-7, 'frequencies': (0.0, 0.4),
                                'is_complex': 'yes', 'damping': 0.5,
                                'use_subspace_solver': 'yes'}
        cphfpolgrad_solver.update_settings(cphfpolgrad_settings)
        cphfpolgrad_solver.ostream.mute()
        cphfpolgrad_solver.compute(molecule, basis, scf_tensors, lr_results)
        cphfpolgrad_solver.compute_omega(molecule, basis, scf_tensors, lr_results)

        # get CPHF lambda coefficients from dist. array
        dist_cphf_coefficients = cphfpolgrad_solver.cphf_results['dist_cphf_ov']
        dof = len(dist_cphf_coefficients)

        cphfpolgrad_coefficients = []
        for x in range(dof):
            solution_vec = dist_cphf_coefficients[x].get_full_vector(0)
            if scf_drv.rank == mpi_master():
                cphfpolgrad_coefficients.append(solution_vec)

        # get CPHF omega coefficients from dist. array
        dist_omega_coefficients = cphfpolgrad_solver.cphf_results['dist_omega_ao']
        dof = len(dist_omega_coefficients)

        cphfpolgrad_omega_coefficients = []
        for x in range(dof):
            solution_vec = dist_omega_coefficients[x].get_full_vector(0)
            if scf_drv.rank == mpi_master():
                cphfpolgrad_omega_coefficients.append(solution_vec)

        if scf_drv.rank == mpi_master():
            cphfpolgrad_results = cphfpolgrad_solver.cphf_results
            #cphfpolgrad_omega_coefficients = np.array([cphfpolgrad_results[0.0]['omega_ao'], 
            #                                          cphfpolgrad_results[0.4]['omega_ao']])
            np.set_printoptions(suppress=True, precision=10)
            here = Path(__file__).parent
            hf_file_name = str(here /'data'/'cphfpolgrad_coefficients.h5')
            hf = h5py.File(hf_file_name, 'r')
            cphfpolgrad_reference = np.array(hf.get(label1))[12:]
            cphfpolgrad_omega_reference = np.array(hf.get(label2))[12:]
            hf.close()

            cphfpolgrad_coefficients = np.array(cphfpolgrad_coefficients)
            cphfpolgrad_coefficients = cphfpolgrad_coefficients.reshape(cphfpolgrad_reference.shape)

            # Here we ar comparing the CPHF coefficients in MO basis, so
            # there might be sign differences; we compare absolute values instead.
            assert np.max(np.abs(cphfpolgrad_coefficients) - np.abs(cphfpolgrad_reference)) < 1.0e-6
            assert np.max(np.abs(cphfpolgrad_omega_coefficients)
                                   - np.abs(cphfpolgrad_omega_reference)) < 1.0e-6 

    def test_cphfpolgrad_coefficients_real(self):
        h2o_xyz = """3

        O     0.000000    0.000000    0.000000
        H     0.000000    0.504284    0.758602
        H     0.000000   -0.504284    0.758602 
        """
        basis_set_label = "sto-3g"

        molecule = Molecule.from_xyz_string(h2o_xyz)
        basis = MolecularBasis.read(molecule, basis_set_label)

        self.run_cphfpolgrad_cg_real(molecule, basis, None, "cphfpolgrad_coefficients_red_real",
                                  "cphfpolgrad_omega_coefficients_red_real")
        self.run_cphfpolgrad_subspace_real(molecule, basis, None, "cphfpolgrad_coefficients_red_real",
                                  "cphfpolgrad_omega_coefficients_red_real")

    def test_cphfpolgrad_coefficients_complex(self):
        h2o_xyz = """3

        O     0.000000    0.000000    0.000000
        H     0.000000    0.504284    0.758602
        H     0.000000   -0.504284    0.758602 
        """
        basis_set_label = "sto-3g"

        molecule = Molecule.from_xyz_string(h2o_xyz)
        basis = MolecularBasis.read(molecule, basis_set_label)

        self.run_cphfpolgrad_cg_complex(molecule, basis, None, "cphfpolgrad_coefficients_red_complex",
                                     "cphfpolgrad_omega_coefficients_red_complex")
        self.run_cphfpolgrad_subspace_complex(molecule, basis, None, "cphfpolgrad_coefficients_red_complex",
                                     "cphfpolgrad_omega_coefficients_red_complex")

    def test_cpkspolgrad_coefficients_real(self):
        h2o_xyz = """3

        O     0.000000    0.000000    0.000000
        H     0.000000    0.504284    0.758602
        H     0.000000   -0.504284    0.758602 
        """
        basis_set_label = "sto-3g"

        molecule = Molecule.from_xyz_string(h2o_xyz)
        basis = MolecularBasis.read(molecule, basis_set_label)


        self.run_cphfpolgrad_cg_real(molecule, basis, "b3lyp", "cpkspolgrad_coefficients_red_b3lyp_real",
                                  "cpkspolgrad_omega_coefficients_red_b3lyp_real")
        self.run_cphfpolgrad_subspace_real(molecule, basis, "b3lyp", "cpkspolgrad_coefficients_red_b3lyp_real",
                                  "cpkspolgrad_omega_coefficients_red_b3lyp_real")

    def test_cpkspolgrad_coefficients_complex(self):
        h2o_xyz = """3

        O     0.000000    0.000000    0.000000
        H     0.000000    0.504284    0.758602
        H     0.000000   -0.504284    0.758602 
        """
        basis_set_label = "sto-3g"

        molecule = Molecule.from_xyz_string(h2o_xyz)
        basis = MolecularBasis.read(molecule, basis_set_label)

        self.run_cphfpolgrad_cg_complex(molecule, basis, "b3lyp", "cpkspolgrad_coefficients_red_b3lyp_complex",
                                     "cpkspolgrad_omega_coefficients_red_b3lyp_complex")
        self.run_cphfpolgrad_subspace_complex(molecule, basis, "b3lyp", "cpkspolgrad_coefficients_red_b3lyp_complex",
                                     "cpkspolgrad_omega_coefficients_red_b3lyp_complex")
