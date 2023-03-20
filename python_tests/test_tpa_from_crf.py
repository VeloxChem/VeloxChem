import pytest

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.cubicresponsedriver import CubicResponseDriver
from veloxchem.tpafulldriver import TpaFullDriver


@pytest.mark.solvers
class TestTpaFromCrf:

    def test_tpa_from_crf(self):

        molecule_string = """
            O   0.0   0.0   0.0
            H   0.0   1.4   1.1
            H   0.0  -1.4   1.1
        """
        basis_set_label = 'def2-svp'

        molecule = Molecule.read_str(molecule_string, units='au')
        basis = MolecularBasis.read(molecule, basis_set_label, ostream=None)

        scf_settings = {}
        method_settings = {'xcfun': 'pbe0', 'grid_level': 1}

        scfdrv = ScfRestrictedDriver()
        scfdrv.ostream.mute()
        scfdrv.update_settings(scf_settings, method_settings)
        scfdrv.compute(molecule, basis)

        w1 = 0.05
        w2 = -w1
        w3 = w1
        damping = 0.1
        components = 'xyz'

        ref_tpa_results = {
            'E3': 0.0,
            'T4': 0.0,
            'A3': 0.0,
            'A2': 0.0,
            'X3': 0.0,
            'X2': 0.0,
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
                    'a_components': a,
                    'b_components': a,
                    'c_components': b,
                    'd_components': b,
                    'damping': damping,
                })
                crf.update_settings(rsp_settings, method_settings)
                crf_results = crf.compute(molecule, basis, scfdrv.scf_tensors)
                if is_mpi_master():
                    ref_tpa_results['gamma'] += crf_results[('gamma', w1, w2,
                                                             w3)]
                    ref_tpa_results['E3'] += crf_results[('E3', w1, w2, w3)]
                    ref_tpa_results['T4'] += crf_results[('T4', w1, w2, w3)]
                    ref_tpa_results['A3'] += crf_results[('A3', w1, w2, w3)]
                    ref_tpa_results['A2'] += crf_results[('A2', w1, w2, w3)]
                    ref_tpa_results['X3'] += crf_results[('X3', w1, w2, w3)]
                    ref_tpa_results['X2'] += crf_results[('X2', w1, w2, w3)]

        for a in components:
            for b in components:
                rsp_settings.update({
                    'b_frequencies': [w1],
                    'c_frequencies': [w2],
                    'd_frequencies': [w3],
                    'a_components': a,
                    'b_components': b,
                    'c_components': a,
                    'd_components': b,
                    'damping': damping,
                })
                crf.update_settings(rsp_settings, method_settings)
                crf_results = crf.compute(molecule, basis, scfdrv.scf_tensors)
                if is_mpi_master():
                    ref_tpa_results['gamma'] += crf_results[('gamma', w1, w2,
                                                             w3)]
                    ref_tpa_results['E3'] += crf_results[('E3', w1, w2, w3)]
                    ref_tpa_results['T4'] += crf_results[('T4', w1, w2, w3)]
                    ref_tpa_results['A3'] += crf_results[('A3', w1, w2, w3)]
                    ref_tpa_results['A2'] += crf_results[('A2', w1, w2, w3)]
                    ref_tpa_results['X3'] += crf_results[('X3', w1, w2, w3)]
                    ref_tpa_results['X2'] += crf_results[('X2', w1, w2, w3)]

        for a in components:
            for b in components:
                rsp_settings.update({
                    'b_frequencies': [w1],
                    'c_frequencies': [w2],
                    'd_frequencies': [w3],
                    'a_components': a,
                    'b_components': b,
                    'c_components': b,
                    'd_components': a,
                    'damping': damping,
                })
                crf.update_settings(rsp_settings, method_settings)
                crf_results = crf.compute(molecule, basis, scfdrv.scf_tensors)
                if is_mpi_master():
                    ref_tpa_results['gamma'] += crf_results[('gamma', w1, w2,
                                                             w3)]
                    ref_tpa_results['E3'] += crf_results[('E3', w1, w2, w3)]
                    ref_tpa_results['T4'] += crf_results[('T4', w1, w2, w3)]
                    ref_tpa_results['A3'] += crf_results[('A3', w1, w2, w3)]
                    ref_tpa_results['A2'] += crf_results[('A2', w1, w2, w3)]
                    ref_tpa_results['X3'] += crf_results[('X3', w1, w2, w3)]
                    ref_tpa_results['X2'] += crf_results[('X2', w1, w2, w3)]

        rsp_settings = {
            'frequencies': [w1],
            'damping': damping,
        }

        tpa = TpaFullDriver()
        tpa.ostream.mute()
        tpa.update_settings(rsp_settings, method_settings)
        tpa_results = tpa.compute(molecule, basis, scfdrv.scf_tensors)

        if is_mpi_master():
            for key in ref_tpa_results:
                ref_tpa_results[key] /= 15.0
                ref_tpa_results[key] *= -1.0  # rsp func. -> gamma

            key_aliases = {
                'E3': 't3_dict',
                'T4': 't4_dict',
                'A3': 'NaA3NxNy',
                'A2': 'NxA2Nyz',
                'X3': 'NaX3NyNz',
                'X2': 'NaX2Nyz',
                'gamma': 'gamma',
            }
            tol = 1.0e-5

            for key in ref_tpa_results:
                ref_val = ref_tpa_results[key]
                calc_val = tpa_results[key_aliases[key]][(w1, w2, w3)]
                assert abs(abs(calc_val.real / ref_val.real) - 1.0) < tol
                assert abs(abs(calc_val.imag / ref_val.imag) - 1.0) < tol
