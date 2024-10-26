import numpy as np

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver


class TestNoradrenaline:

    def run_scf_and_rpa(self, mol, bas, min_bas, xcfun_label, nstates,
                        ref_scf_energy, ref_eigvals, ref_osc_str, ref_rot_str,
                        tol_ene, tol_osc_str, tol_rot_str):

        scf_drv = ScfRestrictedDriver()
        scf_drv.filename = None
        scf_drv.checkpoint_file = None
        scf_drv.xcfun = xcfun_label
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas, min_bas)

        lr_drv = LinearResponseEigenSolver()
        lr_drv.filename = None
        lr_drv.checkpoint_file = None
        lr_drv.xcfun = xcfun_label
        lr_drv.nstates = 5
        lr_drv.nto = False
        lr_drv.ostream.mute()
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == 0:

            scf_energy = scf_results['scf_energy']
            eigvals = lr_results['eigenvalues']
            osc_str = lr_results['oscillator_strengths']
            rot_str = lr_results['rotatory_strengths']

            assert abs(scf_energy - ref_scf_energy) < tol_ene
            assert np.max(np.abs(eigvals - ref_eigvals)) < tol_ene
            assert np.max(np.abs(osc_str - ref_osc_str)) < tol_osc_str
            assert np.max(np.abs(rot_str - ref_rot_str)) < tol_rot_str

    def test_hf(self):

        here = Path(__file__).parent
        xyz_file = str(here / 'data' / 'noradrenaline.xyz')
        mol = Molecule.read_xyz_file(xyz_file)

        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)
        min_bas = MolecularBasis.read(mol, 'ao-start-guess', ostream=None)

        xcfun_label = 'hf'
        nstates = 5
        ref_scf_energy = -587.9305634780
        ref_eigvals = np.array(
            [0.21220288, 0.21905503, 0.27927333, 0.28147272, 0.28292503])
        ref_osc_str = np.array([0.0736, 0.0125, 0.5959, 0.0059, 0.9078])
        ref_rot_str = np.array([1.7666, 1.7914, -31.0072, 23.7739, -15.7826])
        tol_ene = 1.0e-6
        tol_osc_str = 1.0e-4
        tol_rot_str = 1.0e-2
        self.run_scf_and_rpa(mol, bas, min_bas, xcfun_label, nstates,
                             ref_scf_energy, ref_eigvals, ref_osc_str,
                             ref_rot_str, tol_ene, tol_osc_str, tol_rot_str)

    def test_lda(self):

        here = Path(__file__).parent
        xyz_file = str(here / 'data' / 'noradrenaline.xyz')
        mol = Molecule.read_xyz_file(xyz_file)

        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)
        min_bas = MolecularBasis.read(mol, 'ao-start-guess', ostream=None)

        xcfun_label = 'slda'
        nstates = 5
        ref_scf_energy = -588.3389316250
        ref_eigvals = np.array(
            [0.14289744, 0.15587609, 0.17027188, 0.19101771, 0.19523502])
        ref_osc_str = np.array([0.0018, 0.0120, 0.0618, 0.0015, 0.0665])
        ref_rot_str = np.array([-0.2233, 2.3799, -1.0926, 0.1593, -16.8884])
        tol_ene = 1.0e-6
        tol_osc_str = 1.0e-4
        tol_rot_str = 1.0e-4
        self.run_scf_and_rpa(mol, bas, min_bas, xcfun_label, nstates,
                             ref_scf_energy, ref_eigvals, ref_osc_str,
                             ref_rot_str, tol_ene, tol_osc_str, tol_rot_str)

    def test_gga(self):

        here = Path(__file__).parent
        xyz_file = str(here / 'data' / 'noradrenaline.xyz')
        mol = Molecule.read_xyz_file(xyz_file)

        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)
        min_bas = MolecularBasis.read(mol, 'ao-start-guess', ostream=None)

        xcfun_label = 'b3lyp'
        nstates = 5
        ref_scf_energy = -591.4505576871
        ref_eigvals = np.array(
            [0.18307882, 0.19854738, 0.20962970, 0.21116341, 0.21163050])
        ref_osc_str = np.array([0.0561, 0.0001, 0.0307, 0.0587, 0.0065])
        ref_rot_str = np.array([-0.0400, -1.1835, -9.5903, -5.1921, 21.3354])
        tol_ene = 1.0e-6
        tol_osc_str = 1.0e-4
        tol_rot_str = 1.0e-3
        self.run_scf_and_rpa(mol, bas, min_bas, xcfun_label, nstates,
                             ref_scf_energy, ref_eigvals, ref_osc_str,
                             ref_rot_str, tol_ene, tol_osc_str, tol_rot_str)
