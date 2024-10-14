from mpi4py import MPI
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.cubicresponsedriver import CubicResponseDriver
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


@pytest.mark.solvers
class TestCrf:

    def run_scf(self, xcfun_label, basis_set_label):

        molecule_string = """
            O  0.0           0.0  0.0
            H   .7586020000  0.0  -.5042840000
            H   .7586020000  0.0   .5042840000
        """
        molecule = Molecule.read_molecule_string(molecule_string,
                                                 units='angstrom')
        molecule.set_charge(0)
        molecule.set_multiplicity(1)

        basis = MolecularBasis.read(molecule, basis_set_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.conv_thresh = 1.0e-8
        scf_drv.xcfun = xcfun_label
        scf_results = scf_drv.compute(molecule, basis)

        return molecule, basis, scf_results

    def run_crf(self, xcfun_label, basis_set_label, ref_result):

        molecule, ao_basis, scf_results = self.run_scf(xcfun_label,
                                                       basis_set_label)

        wb = 0.11
        wc = -0.2
        wd = 0.06

        rsp_settings = {
            'conv_thresh': 1.0e-8,
            'b_frequencies': [wb],
            'c_frequencies': [wc],
            'd_frequencies': [wd],
            'a_component': 'z',
            'b_component': 'z',
            'c_component': 'z',
            'd_component': 'z',
            'damping': 0.1,
        }

        crf_prop = CubicResponseDriver()
        crf_prop.ostream.mute()
        crf_prop.update_settings(rsp_settings, {'xcfun': xcfun_label})
        crf_result = crf_prop.compute(molecule, ao_basis, scf_results)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            for key in ref_result:
                assert abs(crf_result[(key, wb, wc, wd)].real -
                           ref_result[key].real) < 1.0e-6
                assert abs(crf_result[(key, wb, wc, wd)].imag -
                           ref_result[key].imag) < 1.0e-6

    def test_lda_crf(self):

        ref_result = {'crf': -219.79050033492334 + 25.111107933316774j}

        self.run_crf('slda', 'def2-svp', ref_result)

    def test_gga_hyb_crf(self):

        ref_result = {'crf': -652.4805189222037+60.48380758714163j}

        self.run_crf('b3lyp', 'def2-svpd', ref_result)
