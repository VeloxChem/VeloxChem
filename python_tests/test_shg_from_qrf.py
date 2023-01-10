import numpy as np
import pytest

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.firstorderprop import FirstOrderProperties
from veloxchem.quadraticresponsedriver import QuadraticResponseDriver
from veloxchem.shgdriver import ShgDriver


@pytest.mark.solvers
class TestShgFromQrf:

    def test_shg_from_qrf(self):

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
        scfdrv.ostream.state = False
        scfdrv.update_settings(scf_settings, method_settings)
        scfdrv.compute(molecule, basis)

        scf_prop = FirstOrderProperties()
        scf_prop.compute_scf_prop(molecule, basis, scfdrv.scf_tensors)

        w1 = 0.05
        w2 = w1
        damping = 0.1
        components = 'xyz'

        ref_shg_results = np.zeros(3, dtype='complex128')

        qrf = QuadraticResponseDriver()
        qrf.ostream.state = False
        rsp_settings = {}

        for ind, a in enumerate(components):

            for b in components:
                rsp_settings.update({
                    'b_frequencies': [w1],
                    'c_frequencies': [w2],
                    'a_components': a,
                    'b_components': b,
                    'c_components': b,
                    'damping': damping,
                })
                qrf.update_settings(rsp_settings, method_settings)
                qrf_results = qrf.compute(molecule, basis, scfdrv.scf_tensors)
                if is_mpi_master():
                    ref_shg_results[ind] += qrf_results[(w1, w2)]

            for b in components:
                rsp_settings.update({
                    'b_frequencies': [w1],
                    'c_frequencies': [w2],
                    'a_components': b,
                    'b_components': a,
                    'c_components': b,
                    'damping': damping,
                })
                qrf.update_settings(rsp_settings, method_settings)
                qrf_results = qrf.compute(molecule, basis, scfdrv.scf_tensors)
                if is_mpi_master():
                    ref_shg_results[ind] += qrf_results[(w1, w2)]

            for b in components:
                rsp_settings.update({
                    'b_frequencies': [w1],
                    'c_frequencies': [w2],
                    'a_components': b,
                    'b_components': b,
                    'c_components': a,
                    'damping': damping,
                })
                qrf.update_settings(rsp_settings, method_settings)
                qrf_results = qrf.compute(molecule, basis, scfdrv.scf_tensors)
                if is_mpi_master():
                    ref_shg_results[ind] += qrf_results[(w1, w2)]

        rsp_settings = {
            'frequencies': [w1],
            'damping': damping,
        }

        shg = ShgDriver()
        shg.ostream.state = False
        shg.update_settings(rsp_settings, method_settings)
        shg_results = shg.compute(molecule, basis, scfdrv.scf_tensors)

        if is_mpi_master():
            for ind in range(len(ref_shg_results)):
                ref_shg_results[ind] /= 5.0
                ref_shg_results[ind] *= -1.0  # rsp func. -> beta

            tol = 1.0e-5

            for ind in range(len(ref_shg_results)):
                ref_val = ref_shg_results[ind]
                calc_val = shg_results['beta'][w1][ind]
                if abs(ref_val) > tol:
                    assert abs(abs(calc_val.real / ref_val.real) - 1.0) < tol
                    assert abs(abs(calc_val.imag / ref_val.imag) - 1.0) < tol
                else:
                    assert abs(calc_val.real - ref_val.real) < tol
                    assert abs(calc_val.imag - ref_val.imag) < tol

            dipole = scf_prop.get_property('dipole moment')
            dipole_norm = np.linalg.norm(dipole)

            ref_beta_bar = np.sum(dipole * ref_shg_results) / dipole_norm
            calc_beta_bar = shg_results['beta_bar'][w1]

            assert abs(abs(calc_beta_bar.real / ref_beta_bar.real) - 1.0) < tol
            assert abs(abs(calc_beta_bar.imag / ref_beta_bar.imag) - 1.0) < tol
