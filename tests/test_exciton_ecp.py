from mpi4py import MPI
import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.outputstream import OutputStream
from veloxchem.excitondriver import ExcitonModelDriver


@pytest.mark.solvers
@pytest.mark.timeconsuming
class TestExcitonECP:

    @staticmethod
    def get_h2te_h2o_dimer():

        xyz_string = """6
        H2Te_H2O
        Te   0.000   0.000   0.000
        H    1.520   0.000   0.600
        H   -0.400   1.350   0.600
        O    0.700   0.300   3.600
        H    1.500   0.900   3.900
        H    0.100   1.100   4.300
        """

        molecule = Molecule.read_xyz_string(xyz_string)
        molecule.set_charge(0)
        molecule.set_multiplicity(1)

        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)

        return molecule, basis

    @staticmethod
    def get_reference_data():

        return {
            'exciton_energies': np.array([
                0.18067744, 0.24082588, 0.24689075, 0.32415115, 0.34479102,
                0.48593087
            ]),
            'exciton_osc': np.array(
                [0.0000, 0.0025, 0.0089, 0.0055, 0.0034, 0.0003]),
            'hamiltonian': np.array([
                [
                    1.81784041e-01, 0.00000000e+00, -2.70605074e-05,
                    1.45806827e-04, -3.04521417e-03, 1.77324361e-02
                ],
                [
                    0.00000000e+00, 2.42049633e-01, 3.34217740e-04,
                    -2.91999608e-04, 8.52026398e-03, -1.54251080e-03
                ],
                [
                    -2.70605074e-05, 3.34217740e-04, 2.49864837e-01,
                    0.00000000e+00, 1.55947521e-02, -3.72929334e-03
                ],
                [
                    1.45806827e-04, -2.91999608e-04, 0.00000000e+00,
                    3.44368426e-01, -3.65899970e-03, 5.49083164e-03
                ],
                [
                    -3.04521417e-03, 8.52026398e-03, 1.55947521e-02,
                    -3.65899970e-03, 3.20588375e-01, -3.99689206e-05
                ],
                [
                    1.77324361e-02, -1.54251080e-03, -3.72929334e-03,
                    5.49083164e-03, -3.99689206e-05, 4.84611797e-01
                ],
            ]),
        }

    def run_exciton_ecp(self):

        comm = MPI.COMM_WORLD
        ostream = OutputStream(None)

        molecule, basis = self.get_h2te_h2o_dimer()

        exciton_drv = ExcitonModelDriver(comm, ostream)
        exciton_drv.fragments = '2'
        exciton_drv.atoms_per_fragment = '3'
        exciton_drv.charges = '0'
        exciton_drv.nstates = 2
        exciton_drv.ct_nocc = 1
        exciton_drv.ct_nvir = 1
        exciton_drv.scf_conv_thresh = 1.0e-7
        exciton_drv.tda_conv_thresh = 1.0e-5
        exciton_drv.restart = False

        exciton_results = exciton_drv.compute(molecule, basis)

        return comm, exciton_results

    def test_exciton_ecp_h2te_h2o(self):
        comm, exciton_results = self.run_exciton_ecp()
        ref_data = self.get_reference_data()

        if comm.Get_rank() == mpi_master():
            exciton_energies = exciton_results['adiabatic_eigenvalues']
            exciton_osc = exciton_results['oscillator_strengths']

            assert len(exciton_energies) == 6
            assert np.all(np.isfinite(exciton_energies))
            assert np.all(np.isfinite(exciton_osc))

            np.testing.assert_allclose(exciton_energies,
                                       ref_data['exciton_energies'],
                                       atol=1.0e-6)
            np.testing.assert_allclose(exciton_osc,
                                       ref_data['exciton_osc'],
                                       atol=1.0e-4)

            hamiltonian = exciton_results['hamiltonian']
            ref_hamiltonian = ref_data['hamiltonian']

            np.testing.assert_allclose(np.diag(hamiltonian),
                                       np.diag(ref_hamiltonian),
                                       atol=5.0e-5)
            np.testing.assert_allclose(np.abs(hamiltonian),
                                       np.abs(ref_hamiltonian),
                                       atol=5.0e-5)
