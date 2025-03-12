from pathlib import Path
import numpy as np
import h5py
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.hessianorbitalresponse import HessianOrbitalResponse


@pytest.mark.solvers
class TestCphfSolver:

    def run_cphfsolver(self,
                       molecule,
                       basis,
                       xcfun=None,
                       label=None,
                       use_subcomms=False):
        scf_drv = ScfRestrictedDriver()
        method_settings = {}
        if xcfun is not None:
            scf_drv._dft = True
            scf_drv.xcfun = xcfun
            method_settings = {'xcfun': xcfun}

        scf_drv.ostream.mute()
        scf_tensors = scf_drv.compute(molecule, basis)

        hess_orbrsp_drv = HessianOrbitalResponse()
        orbrsp_settings = {'conv_thresh': 2e-7}
        hess_orbrsp_drv.update_settings(orbrsp_settings, method_settings)
        hess_orbrsp_drv.ostream.mute()

        hess_orbrsp_drv.use_subcomms = use_subcomms

        # TODO: hess_orbrsp_drv should return cphf_results
        hess_orbrsp_drv.compute(molecule, basis, scf_tensors)

        dist_cphf_coefficients = hess_orbrsp_drv.cphf_results['dist_cphf_ov']
        dof = len(dist_cphf_coefficients)

        cphf_coefficients = []
        for x in range(dof):
            solution_vec = dist_cphf_coefficients[x].get_full_vector(0)
            if scf_drv.rank == mpi_master():
                cphf_coefficients.append(solution_vec)

        if scf_drv.rank == mpi_master():
            here = Path(__file__).parent
            hf_file_name = str(here / 'data' / 'cphf_coefficients.h5')
            hf = h5py.File(hf_file_name, 'r')
            cphf_reference = np.array(hf.get(label))
            hf.close()

            cphf_coefficients = np.array(cphf_coefficients)
            cphf_coefficients = cphf_coefficients.reshape(cphf_reference.shape)

            # Here we are comparing the CPHF coefficients in MO basis, so
            # there might be sign differences; we compare absolute values instead.
            assert np.max(np.abs(cphf_coefficients) -
                          np.abs(cphf_reference)) < 1.0e-6

    def test_cphf_coefficients(self):
        nh3_xyz = """4

        N     0.000000000     0.000000000     0.000000000
        H    -0.653401663     0.309213352     0.817609879
        H     0.695693936     0.071283622    -0.702632331
        H     0.330018952    -0.826347607     0.052270023
        """
        basis_set_label = "sto-3g"

        molecule = Molecule.from_xyz_string(nh3_xyz)
        basis = MolecularBasis.read(molecule, basis_set_label)

        self.run_cphfsolver(molecule, basis, "hf", "cphf_coefficients")

        self.run_cphfsolver(molecule,
                            basis,
                            "hf",
                            "cphf_coefficients",
                            use_subcomms=True)

    def test_cpks_coefficients(self):
        nh3_xyz = """4

        N     0.000000000     0.000000000     0.000000000
        H    -0.653401663     0.309213352     0.817609879
        H     0.695693936     0.071283622    -0.702632331
        H     0.330018952    -0.826347607     0.052270023
        """
        basis_set_label = "sto-3g"

        molecule = Molecule.from_xyz_string(nh3_xyz)
        basis = MolecularBasis.read(molecule, basis_set_label)

        self.run_cphfsolver(molecule, basis, "b3lyp", "cpks_coefficients_b3lyp")

        self.run_cphfsolver(molecule,
                            basis,
                            "b3lyp",
                            "cpks_coefficients_b3lyp",
                            use_subcomms=True)
