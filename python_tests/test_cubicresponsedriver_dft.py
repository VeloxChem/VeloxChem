import pytest

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.cubicresponsedriver import CubicResponseDriver
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


@pytest.mark.solvers
class TestCrf:

    def run_scf(self):

        molecule_string = """
            O  0.0           0.0  0.0
            H   .7586020000  0.0  -.5042840000
            H   .7586020000  0.0   .5042840000
        """

        basis_set_label = 'def2-svpd'

        scf_settings = {'conv_thresh': 1.0e-8}

        molecule = Molecule.read_str(molecule_string, units='ang')
        molecule.set_charge(0)
        molecule.set_multiplicity(1)

        basis = MolecularBasis.read(molecule, basis_set_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.state = False
        method_settings = {'xcfun': 'SLDA', 'grid_level': 4}
        scf_drv.update_settings(scf_settings,method_settings)
        scf_drv.compute(molecule, basis)

        return scf_drv.scf_tensors, molecule, basis,method_settings

    def run_crf(self, ref_result):

        scf_tensors, molecule, ao_basis,method_settings = self.run_scf()

        wb = 0
        wc = 0
        wd = 0

        rsp_settings = {
            'conv_thresh': 1.0e-8,
            'b_frequencies': [wb],
            'c_frequencies': [wc],
            'd_frequencies': [wd],
            'a_components': 'z',
            'b_components': 'z',
            'c_components': 'z',
            'd_components': 'z',
            'damping': 0
        }

        crf_prop = CubicResponseDriver()
        crf_prop.ostream.state = False
        crf_prop.update_settings(rsp_settings,method_settings)
        crf_result_yyzz = crf_prop.compute(molecule, ao_basis, scf_tensors)

        if is_mpi_master():

            assert abs(crf_result_yyzz[('T4', wb, wc, wd)].real -
                       ref_result['T4'].real) < 1.0e-6

            assert abs(crf_result_yyzz[('T4', wb, wc, wd)].imag -
                       ref_result['T4'].imag) < 1.0e-6

            assert abs(crf_result_yyzz[('X3', wb, wc, wd)].real -
                       ref_result['X3'].real) < 1.0e-6

            assert abs(crf_result_yyzz[('X3', wb, wc, wd)].imag -
                       ref_result['X3'].imag) < 1.0e-6

            assert abs(crf_result_yyzz[('A3', wb, wc, wd)].real -
                       ref_result['A3'].real) < 1.0e-6

            assert abs(crf_result_yyzz[('A3', wb, wc, wd)].imag -
                       ref_result['A3'].imag) < 1.0e-6

    def test_crf(self):

        ref_result = {
            'E3':   103.05627199,
            'T4':   58.87127176 ,
            'X2': -505.15579639 ,
            'X3':  133.84179900 ,
            'A2': -505.15579639 ,
            'A3':   44.61393300 ,
        }

        self.run_crf(ref_result)
