from mpi4py import MPI
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.quadraticresponsedriver import QuadraticResponseDriver
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


@pytest.mark.solvers
class TestQrf:

    def run_scf(self, basis_set_label, *, ri):

        molecule_string = """
            O  0.0           0.0  0.0
            H   .7586020000  0.0  -.5042840000
            H   .7586020000  0.0   .5042840000
        """

        molecule = Molecule.read_molecule_string(molecule_string,
                                                 units='angstrom')
        basis = MolecularBasis.read(molecule, basis_set_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        return molecule, basis, scf_results

    def run_qrf(self, basis_set_label, ref_result, *, ri):

        molecule, basis, scf_results = self.run_scf(basis_set_label, ri=ri)

        qrf_drv = QuadraticResponseDriver()
        qrf_drv.ostream.mute()

        rsp_settings = {
            'conv_thresh': 1.0e-4,
            'b_frequencies': [0.2],
            'c_frequencies': [0.2],
            'damping': 0.1,
            'a_component': 'x',
            'b_component': 'x',
            'c_component': 'x'
        }
        qrf_drv.update_settings(rsp_settings)
        qrf_result_xxx = qrf_drv.compute(molecule, basis, scf_results)

        rsp_settings = {
            'conv_thresh': 1.0e-4,
            'b_frequencies': [0.2],
            'c_frequencies': [0.2],
            'damping': 0.1,
            'a_component': 'z',
            'b_component': 'z',
            'c_component': 'x'
        }
        qrf_drv.update_settings(rsp_settings)
        qrf_result_zzx = qrf_drv.compute(molecule, basis, scf_results)

        rsp_settings = {
            'conv_thresh': 1.0e-4,
            'b_frequencies': [0.2],
            'c_frequencies': [0.2],
            'damping': 0.1,
            'a_component': 'y',
            'b_component': 'y',
            'c_component': 'x'
        }
        qrf_drv.update_settings(rsp_settings)
        qrf_result_yyx = qrf_drv.compute(molecule, basis, scf_results)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            thresh = 1.0e-5

            # x-component
            assert abs(qrf_result_xxx[('qrf', 0.2, 0.2)].real -
                       ref_result['xxx'].real) < thresh
            assert abs(qrf_result_xxx[('qrf', 0.2, 0.2)].imag -
                       ref_result['xxx'].imag) < thresh

            # y-component
            assert abs(qrf_result_yyx[('qrf', 0.2, 0.2)].real -
                       ref_result['yyx'].real) < thresh
            assert abs(qrf_result_yyx[('qrf', 0.2, 0.2)].imag -
                       ref_result['yyx'].imag) < thresh

            # z-component
            assert abs(qrf_result_zzx[('qrf', 0.2, 0.2)].real -
                       ref_result['zzx'].real) < thresh
            assert abs(qrf_result_zzx[('qrf', 0.2, 0.2)].imag -
                       ref_result['zzx'].imag) < thresh

    @pytest.mark.parametrize('ri', [False, True])
    def test_svp_qrf(self, ri):

        ref_result = {
            False: {
                'xxx': -19.93302793247569 - 20.644727090814307j,
                'zzx': -17.755885250507482 - 28.454822496543663j,
                'yyx': 7.431953463379896 - 2.5194096972044324j,
            },
            True: {
                'xxx': -19.93617607 - 20.65138401j,
                'zzx': -17.76129431 - 28.46074480j,
                'yyx': 7.43141364 - 2.51908745j
            }
        }[ri]

        self.run_qrf('def2-svp', ref_result, ri=ri)

    @pytest.mark.parametrize('ri', [False, True])
    def test_tzvp_qrf(self, ri):

        ref_result = {
            False: {
                'xxx': -16.15272019367957 - 24.535238439985143j,
                'zzx': -12.50114052353982 - 30.600755240976852j,
                'yyx': 13.085839529577688 + 0.1218738912151196j,
            },
            True: {
                'xxx': -16.15539689 - 24.53581832j,
                'zzx': -12.50887849 - 30.60553853j,
                'yyx': 13.09102179 + 0.12223456j
            }
        }[ri]

        self.run_qrf('def2-tzvp', ref_result, ri=ri)
