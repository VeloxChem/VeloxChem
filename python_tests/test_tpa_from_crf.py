import pytest

from veloxchem.veloxchemlib import is_single_node, is_mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.cubicresponsedriver import CubicResponseDriver
from veloxchem.tpafulldriver import TpaFullDriver


class TestTpaFromCrf:

    def run_tpa_from_crf(self, xcfun_label):

        molecule_string = """
            O   0.0   0.0   0.0
            H   0.0   1.4   1.1
            H   0.0  -1.4   1.1
        """
        basis_set_label = 'def2-svp'

        molecule = Molecule.read_molecule_string(molecule_string, units='au')
        basis = MolecularBasis.read(molecule, basis_set_label, ostream=None)

        scf_settings = {}
        method_settings = {'xcfun': xcfun_label, 'grid_level': 1}

        scfdrv = ScfRestrictedDriver()
        scfdrv.ostream.mute()
        scfdrv.update_settings(scf_settings, method_settings)
        scf_results = scfdrv.compute(molecule, basis)

        w1 = 0.05
        w2 = -w1
        w3 = w1
        damping = 0.1
        components = 'xyz'

        ref_tpa_results = {
            'gamma': 0.0,
        }

        crf = CubicResponseDriver()
        crf.ostream.mute()
        rsp_settings = {}

        for a in components:
            for b in components:
                rsp_settings.update({
                    'b_frequencies': [w1],
                    'c_frequencies': [w2],
                    'd_frequencies': [w3],
                    'a_component': a,
                    'b_component': a,
                    'c_component': b,
                    'd_component': b,
                    'damping': damping,
                })
                crf.update_settings(rsp_settings, method_settings)
                crf_results = crf.compute(molecule, basis, scf_results)
                if is_mpi_master():
                    ref_tpa_results['gamma'] += crf_results[('crf', w1, w2, w3)]

        for a in components:
            for b in components:
                rsp_settings.update({
                    'b_frequencies': [w1],
                    'c_frequencies': [w2],
                    'd_frequencies': [w3],
                    'a_component': a,
                    'b_component': b,
                    'c_component': a,
                    'd_component': b,
                    'damping': damping,
                })
                crf.update_settings(rsp_settings, method_settings)
                crf_results = crf.compute(molecule, basis, scf_results)
                if is_mpi_master():
                    ref_tpa_results['gamma'] += crf_results[('crf', w1, w2, w3)]

        for a in components:
            for b in components:
                rsp_settings.update({
                    'b_frequencies': [w1],
                    'c_frequencies': [w2],
                    'd_frequencies': [w3],
                    'a_component': a,
                    'b_component': b,
                    'c_component': b,
                    'd_component': a,
                    'damping': damping,
                })
                crf.update_settings(rsp_settings, method_settings)
                crf_results = crf.compute(molecule, basis, scf_results)
                if is_mpi_master():
                    ref_tpa_results['gamma'] += crf_results[('crf', w1, w2, w3)]

        rsp_settings = {
            'frequencies': [w1],
            'damping': damping,
        }

        tpa = TpaFullDriver()
        tpa.ostream.mute()
        tpa.update_settings(rsp_settings, method_settings)
        tpa_results = tpa.compute(molecule, basis, scf_results)

        if is_mpi_master():
            for key in ref_tpa_results:
                ref_tpa_results[key] /= 15.0
                ref_tpa_results[key] *= -1.0  # rsp func. -> gamma

            tol = 1.0e-5

            for key in ref_tpa_results:
                ref_val = ref_tpa_results[key]
                calc_val = tpa_results[key][(w1, w2, w3)]
                assert abs(abs(calc_val.real / ref_val.real) - 1.0) < tol
                assert abs(abs(calc_val.imag / ref_val.imag) - 1.0) < tol

    @pytest.mark.skipif(is_single_node(), reason="multi-node only")
    def test_tpa_from_crf_lda(self):

        self.run_tpa_from_crf('slda')

    @pytest.mark.skipif(is_single_node(), reason="multi-node only")
    def test_tpa_from_crf_gga(self):

        self.run_tpa_from_crf('pbe0')

    @pytest.mark.skipif(is_single_node(), reason="multi-node only")
    def test_tpa_from_crf_mgga(self):

        self.run_tpa_from_crf('tpssh')
