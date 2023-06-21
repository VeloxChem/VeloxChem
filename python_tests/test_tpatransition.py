import pytest

from veloxchem.veloxchemlib import is_single_node, is_mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tpatransitiondriver import TpaTransitionDriver


class TestTpaTransition:

    def run_scf(self, xcfun_label):

        molecule_string = """
            O  0.0           0.0  0.0
            H   .7586020000  0.0  -.5042840000
            H   .7586020000  0.0   .5042840000
        """
        basis_set_label = 'def2-svpd'
        scf_conv_thresh = 1.0e-8

        molecule = Molecule.read_molecule_string(molecule_string)
        basis = MolecularBasis.read(molecule, basis_set_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = xcfun_label
        scf_drv.conv_thresh = scf_conv_thresh
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        return scf_results, molecule, basis

    def run_tpatransition(self, xcfun_label, ref_results):

        tpa_nstates = 2
        tpa_conv_thresh = 1.0e-7

        scf_results, molecule, ao_basis = self.run_scf(xcfun_label)

        tpa_drv = TpaTransitionDriver()
        tpa_drv.xcfun = xcfun_label
        tpa_drv.nstates = tpa_nstates
        tpa_drv.conv_thresh = tpa_conv_thresh
        tpa_drv.ostream.mute()
        tpa_results = tpa_drv.compute(molecule, ao_basis, scf_results)

        if is_mpi_master():
            tpa_str = tpa_results['tpa_strengths']['linear']
            for (freq, val), (ref_freq, ref_val) in zip(tpa_str.items(),
                                                        ref_results.items()):
                assert abs(freq / ref_freq - 1.0) < 1.0e-6
                assert abs(val / ref_val - 1.0) < 1.0e-6

    @pytest.mark.solvers
    def test_tpatransition_hf(self):

        ref_result = {
            -0.1656922537003149: 4.130284073574981,
            -0.21190963394768278: 15.83816511497292,
        }
        self.run_tpatransition('hf', ref_result)

    def test_tpatransition_lda(self):

        ref_result = {
            -0.13269116242012727: 12.922686324743287,
            -0.18051942347838879: 55.35603386346646,
        }
        self.run_tpatransition('slda', ref_result)

    @pytest.mark.solvers
    def test_tpatransition_gga(self):

        ref_result = {
            -0.13304132289794054: 13.068064328595725,
            -0.179389130413857: 57.062788455439296,
        }
        self.run_tpatransition('bp86', ref_result)

    @pytest.mark.skipif(is_single_node(), reason="multi-node only")
    def test_tpatransition_mgga(self):

        ref_result = {
            -0.1381414072229231: 9.579957827267322,
            -0.18265658100098003: 41.94340887233736,
        }
        self.run_tpatransition('tpssh', ref_result)
