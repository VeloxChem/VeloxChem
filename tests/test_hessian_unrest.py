import numpy as np
import pytest
from pathlib import Path

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.scfhessiandriver import ScfHessianDriver
from veloxchem.vibrationalanalysis import VibrationalAnalysis


class TestScfUnrestrictedHessian:

    @staticmethod
    def get_embedded_water_molecule_and_basis():

        mol = Molecule.read_xyz_string("""3
        xyz
        O    1.2361419   1.0137761  -0.0612424
        H    0.5104418   0.8944555   0.5514190
        H    1.9926927   1.1973129   0.4956931
        """)
        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)
        return mol, bas

    def run_hessian_unrest(self, xcfun_label, ref_vib_freqs):

        mol = Molecule.read_xyz_string("""4
        xyz
        C       -0.8246083505   -2.6364102075   -0.0255277471
        H       -0.4137775504   -3.6364807662   -0.1702755824
        H       -1.8972335969   -2.4672317696   -0.1290827158
        H       -0.1628531453   -1.8055908756    0.2228478391
        """)
        mol.set_multiplicity(2)

        bas = MolecularBasis.read(mol, 'def2-svp')

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_drv.compute(mol, bas)

        vibanalysis_drv = VibrationalAnalysis(scf_drv)
        vibanalysis_drv.ostream.mute()
        vibanalysis_drv.compute(mol, bas)

        atom_pairs = [(0, 3), (2, 3)]
        hess_drv = ScfHessianDriver(scf_drv)
        hess_drv.ostream.mute()
        hess_drv.atom_pairs = atom_pairs
        hess_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            # compare calculated frequencies with reference
            calc_vib_freqs = vibanalysis_drv.vib_frequencies
            assert np.max(np.abs(calc_vib_freqs - ref_vib_freqs)) < 0.1

            # compare full Hessian with atom pair Hessian
            for a, b in atom_pairs:
                start_a, end_a = a * 3, (a + 1) * 3
                start_b, end_b = b * 3, (b + 1) * 3
                diff_aa = np.max(
                    np.abs(vibanalysis_drv.hessian[start_a:end_a,
                                                   start_a:end_a] -
                           hess_drv.hessian[start_a:end_a, start_a:end_a]))
                diff_ab = np.max(
                    np.abs(vibanalysis_drv.hessian[start_a:end_a,
                                                   start_b:end_b] -
                           hess_drv.hessian[start_a:end_a, start_b:end_b]))
                diff_bb = np.max(
                    np.abs(vibanalysis_drv.hessian[start_b:end_b,
                                                   start_b:end_b] -
                           hess_drv.hessian[start_b:end_b, start_b:end_b]))
                assert diff_aa < 1e-8
                assert diff_ab < 1e-8
                assert diff_bb < 1e-8

    @pytest.mark.solvers
    def test_hessian_unrest_pbe(self):

        xcfun_label = 'pbe'

        ref_vib_freqs = np.array(
            [393.84, 1309.27, 1309.32, 3105.28, 3311.26, 3311.30])

        self.run_hessian_unrest(xcfun_label, ref_vib_freqs)

    @pytest.mark.solvers
    def test_hessian_unrest_b3lyp(self):

        xcfun_label = 'b3lyp'

        ref_vib_freqs = np.array(
            [470.01, 1373.86, 1373.90, 3096.68, 3296.69, 3296.72])

        self.run_hessian_unrest(xcfun_label, ref_vib_freqs)

    def run_hessian_unrest_with_ecp(self, ref_vib_freqs, ref_ir_intens):

        mol = Molecule.read_xyz_string("""2
        xyz
        Au 0 0 0
        H  0 0 1.55
        """)
        mol.set_charge(1)
        mol.set_multiplicity(2)

        bas = MolecularBasis.read(mol, 'def2-svp')

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.compute(mol, bas)

        vib_drv = VibrationalAnalysis(scf_drv)
        vib_drv.ostream.mute()
        vib_results = vib_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            assert np.max(np.abs(vib_results['vib_frequencies'] -
                                 ref_vib_freqs)) < 0.1
            assert np.max(np.abs(vib_results['ir_intensities'] -
                                 ref_ir_intens)) < 0.1

    @pytest.mark.solvers
    def test_hessian_unrest_with_ecp(self):

        ref_vib_freqs = np.array([1949.51])
        ref_ir_intens = np.array([278.8583])
        self.run_hessian_unrest_with_ecp(ref_vib_freqs, ref_ir_intens)

    @pytest.mark.solvers
    def test_hessian_unrest_with_point_charges_and_vdw(self):

        mol, bas = self.get_embedded_water_molecule_and_basis()

        here = Path(__file__).parent
        potfile = str(here / 'data' / 'pe_water.pot')
        vdwfile = str(here / 'data' / 'pe_water.qm_vdw_params.txt')
        reffile = str(here / 'data' /
                      'water_analytical_hessian_pointcharges_unrest_b3lyp.txt')

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = 'b3lyp'
        scf_drv.point_charges = potfile
        scf_drv.qm_vdw_params = vdwfile
        scf_drv.compute(mol, bas)

        hess_drv = ScfHessianDriver(scf_drv)
        hess_drv.ostream.mute()
        hess_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            ref_hessian = np.loadtxt(reffile)
            np.testing.assert_allclose(hess_drv.hessian,
                                       ref_hessian,
                                       rtol=1.0e-8,
                                       atol=1.0e-10)
