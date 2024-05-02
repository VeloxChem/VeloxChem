import numpy as np

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver


class TestChargeTransfer:

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
            if ref_rot_str is not None:
                assert np.max(np.abs(rot_str - ref_rot_str)) < tol_rot_str

    def get_xyz_string(self):

        return """17
            xyz
            C -4.208800502216 1.293225309841 0.221551598989
            C -2.952270842031 1.128119464149 0.648143185934
            H -4.744457906203 0.490650280784 -0.257353343764
            H -4.719590150518 2.233545340929 0.345269677541
            C -2.250268294891 -0.099106399978 0.500261548074
            C -2.217795218109 2.166079272508 1.283535045385
            N -1.625288817427 3.000502225424 1.795148291861
            N -1.683368268605 -1.085895493656 0.38230399598
            C 12.728353166525482 -3.9292358967565337 -0.3519791087758094
            C 11.456224648876482 -3.520119250924534 -0.6030565924898094
            C 11.391018630170482 -2.172533778849534 -0.17023018196280948
            H 10.666161298975481 -4.097307814747534 -1.037524511553809
            C 12.629152591891481 -1.8760026307035336 0.30709830789219067
            H 10.541888241705482 -1.521886408923534 -0.21175701062380947
            O 13.444872603674481 -2.9369967678605335 0.20002799531919097
            H 13.226820596466483 -4.864143995235533 -0.5145952282898093
            H 13.038745208803482 -0.9816097585085335 0.7328653857981908
        """

    def test_gga(self):

        mol = Molecule.read_xyz_string(self.get_xyz_string())
        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)
        min_bas = MolecularBasis.read(mol, 'ao-start-guess', ostream=None)

        xcfun_label = 'bhandhlyp'
        nstates = 5
        ref_scf_energy = -492.4194736215
        ref_eigvals = np.array(
            [0.19246507, 0.22511728, 0.24181816, 0.25374163, 0.28693381])
        ref_osc_str = np.array([0.0, 0.3206, 0.1589, 0.0, 0.0])
        ref_rot_str = np.array([0.0, -1.2202, 0.7409, 0.0, 0.0])
        tol_ene = 1.0e-6
        tol_osc_str = 1.0e-4
        tol_rot_str = 1.0e-4
        self.run_scf_and_rpa(mol, bas, min_bas, xcfun_label, nstates,
                             ref_scf_energy, ref_eigvals, ref_osc_str,
                             ref_rot_str, tol_ene, tol_osc_str, tol_rot_str)

    def test_rs_gga(self):

        mol = Molecule.read_xyz_string(self.get_xyz_string())
        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)
        min_bas = MolecularBasis.read(mol, 'ao-start-guess', ostream=None)

        xcfun_label = 'lrc-wpbeh'
        nstates = 5
        ref_scf_energy = -492.212402609668
        ref_eigvals = np.array(
            [0.224518580, 0.235169752, 0.243835641, 0.289719169, 0.326796249])
        ref_osc_str = np.array([0.3054478465, 0.0, 0.1610691327, 0.0, 0.0])
        ref_rot_str = None
        tol_ene = 1.0e-4
        tol_osc_str = 1.0e-4
        tol_rot_str = None
        self.run_scf_and_rpa(mol, bas, min_bas, xcfun_label, nstates,
                             ref_scf_energy, ref_eigvals, ref_osc_str,
                             ref_rot_str, tol_ene, tol_osc_str, tol_rot_str)
