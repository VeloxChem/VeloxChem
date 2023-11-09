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

        molecule = Molecule.read_molecule_string(molecule_string,
                                                 units='angstrom')
        molecule.set_charge(0)
        molecule.set_multiplicity(1)

        basis = MolecularBasis.read(molecule, basis_set_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        method_settings = {'xcfun': 'SLDA', 'grid_level': 4}
        scf_drv.update_settings(scf_settings, method_settings)
        scf_results = scf_drv.compute(molecule, basis)

        return scf_results, molecule, basis, method_settings

    def run_crf(self, ref_result):

        scf_results, molecule, ao_basis, method_settings = self.run_scf()

        wb = 0
        wc = 0
        wd = 0

        rsp_settings = {
            'conv_thresh': 1.0e-8,
            'b_frequencies': [wb],
            'c_frequencies': [wc],
            'd_frequencies': [wd],
            'a_component': 'z',
            'b_component': 'z',
            'c_component': 'z',
            'd_component': 'z',
            'damping': 0
        }

        crf_prop = CubicResponseDriver()
        crf_prop.ostream.mute()
        crf_prop.update_settings(rsp_settings, method_settings)
        crf_result = crf_prop.compute(molecule, ao_basis, scf_results)

        if is_mpi_master():
            for key in ref_result:
                assert abs(crf_result[(key, wb, wc, wd)].real -
                           ref_result[key].real) < 1.0e-6
                assert abs(crf_result[(key, wb, wc, wd)].imag -
                           ref_result[key].imag) < 1.0e-6

    def test_crf(self):

        ref_result = {
            'E3': 103.05625002,
            'T4': 58.87171313,
            'X2': -505.15455288,
            'X3': 133.84194721,
            'A2': -505.15455288,
            'A3': 44.61398240,
        }

        self.run_crf(ref_result)
