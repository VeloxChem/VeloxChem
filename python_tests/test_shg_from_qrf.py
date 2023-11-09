import numpy as np

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.firstorderprop import FirstOrderProperties
from veloxchem.quadraticresponsedriver import QuadraticResponseDriver
from veloxchem.shgdriver import ShgDriver


class TestShgFromQrf:

    def run_shg_from_qrf(self, xcfun_label):

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

        scf_prop = FirstOrderProperties()
        scf_prop.compute_scf_prop(molecule, basis, scf_results)

        w1 = 0.05
        w2 = w1
        damping = 0.1
        components = 'xyz'

        ref_shg_results = np.zeros(3, dtype='complex128')

        qrf = QuadraticResponseDriver()
        qrf.ostream.mute()
        rsp_settings = {}

        for ind, a in enumerate(components):

            for b in components:
                rsp_settings.update({
                    'b_frequencies': [w1],
                    'c_frequencies': [w2],
                    'a_component': a,
                    'b_component': b,
                    'c_component': b,
                    'damping': damping,
                })
                qrf.update_settings(rsp_settings, method_settings)
                qrf_results = qrf.compute(molecule, basis, scf_results)
                if is_mpi_master():
                    ref_shg_results[ind] += qrf_results[(w1, w2)]

            for b in components:
                rsp_settings.update({
                    'b_frequencies': [w1],
                    'c_frequencies': [w2],
                    'a_component': b,
                    'b_component': a,
                    'c_component': b,
                    'damping': damping,
                })
                qrf.update_settings(rsp_settings, method_settings)
                qrf_results = qrf.compute(molecule, basis, scf_results)
                if is_mpi_master():
                    ref_shg_results[ind] += qrf_results[(w1, w2)]

            for b in components:
                rsp_settings.update({
                    'b_frequencies': [w1],
                    'c_frequencies': [w2],
                    'a_component': b,
                    'b_component': b,
                    'c_component': a,
                    'damping': damping,
                })
                qrf.update_settings(rsp_settings, method_settings)
                qrf_results = qrf.compute(molecule, basis, scf_results)
                if is_mpi_master():
                    ref_shg_results[ind] += qrf_results[(w1, w2)]

        rsp_settings = {
            'frequencies': [w1],
            'damping': damping,
        }

        shg = ShgDriver()
        shg.ostream.mute()
        shg.update_settings(rsp_settings, method_settings)
        shg_results = shg.compute(molecule, basis, scf_results)

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

    def test_shg_from_qrf_lda(self):

        self.run_shg_from_qrf('slda')

    def test_shg_from_qrf_gga(self):

        self.run_shg_from_qrf('pbe0')

    def test_shg_from_qrf_mgga(self):

        self.run_shg_from_qrf('tpssh')
