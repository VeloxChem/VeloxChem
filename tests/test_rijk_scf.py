import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.rijkfockdriver import RIJKFockDriver


@pytest.mark.solvers
class TestScfDriverWithRIJK:

    def run_scf(self, scf_flag, mol, bas, xcfun_label, ref_scf_energy, tol):

        if scf_flag == 'restricted':
            scf_drv = ScfRestrictedDriver()
        elif scf_flag == 'unrestricted':
            scf_drv = ScfUnrestrictedDriver()

        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_drv.ri_jk = True
        scf_results = scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            assert abs(ref_scf_energy - scf_results['scf_energy']) < tol

    def run_scf_with_modifier(self, scf_flag, mol, bas, ref_scf_energy, tol,
                              configure):
        """
        Runs an RI-JK SCF calculation with a convergence modifier enabled.
        """

        if scf_flag == 'restricted':
            scf_drv = ScfRestrictedDriver()
        elif scf_flag == 'unrestricted':
            scf_drv = ScfUnrestrictedDriver()
        else:
            raise ValueError(f'unknown scf_flag: {scf_flag}')

        scf_drv.ostream.mute()
        scf_drv.ri_jk = True
        configure(scf_drv)
        scf_results = scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            assert scf_drv.is_converged
            assert scf_results is not None
            assert abs(ref_scf_energy -
                       scf_results['scf_energy']) < tol

        return scf_drv, scf_results

    def _record_density_factor_calls(self, monkeypatch):
        """
        Patches the RI-JK screened exchange method to record whether the
        density-factor path is requested.
        """

        calls = []
        original = RIJKFockDriver.compute_screened_k_fock

        def spy(instance, *args, **kwargs):
            calls.append(kwargs.get('use_density_factor', False))
            return original(instance, *args, **kwargs)

        monkeypatch.setattr(RIJKFockDriver, 'compute_screened_k_fock', spy)

        return calls

    def _get_methanol_data(self):

        xyz_string = """6
        xyz
        H      1.2001      0.0363      0.8431
        C      0.7031      0.0083     -0.1305
        H      0.9877      0.8943     -0.7114
        H      1.0155     -0.8918     -0.6742
        O     -0.6582     -0.0067      0.1730
        H     -1.1326     -0.0311     -0.6482
        """

        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)

        return mol, bas

    def _get_methanol_radical_data(self):

        xyz_string = """5
        xyz
        H      1.2001      0.0363      0.8431
        C      0.7031      0.0083     -0.1305
        H      0.9877      0.8943     -0.7114
        H      1.0155     -0.8918     -0.6742
        O     -0.6582     -0.0067      0.1730
        """

        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_multiplicity(2)
        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)

        return mol, bas

    def test_rijk_hf_with_pfon(self, monkeypatch):

        mol, bas = self._get_methanol_data()

        calls = self._record_density_factor_calls(monkeypatch)

        def configure(drv):
            drv.pfon = True
            drv.pfon_temperature = 1000

        scf_drv, _ = self.run_scf_with_modifier(
            'restricted', mol, bas, -114.954012228, 1.0e-8, configure)

        assert scf_drv.pfon_temperature == 0
        assert True in calls

    def test_rijk_hf_with_density_damping(self, monkeypatch):

        mol, bas = self._get_methanol_data()

        calls = self._record_density_factor_calls(monkeypatch)

        def configure(drv):
            drv.density_damping = True
            drv.acc_type = 'diis'

        self.run_scf_with_modifier(
            'restricted', mol, bas, -114.954012228, 1.0e-8, configure)

        assert True in calls

    def test_rijk_hf_openshell_with_pfon(self, monkeypatch):

        mol, bas = self._get_methanol_radical_data()

        calls = self._record_density_factor_calls(monkeypatch)

        def configure(drv):
            drv.pfon = True
            drv.pfon_temperature = 1000

        scf_drv, _ = self.run_scf_with_modifier(
            'unrestricted', mol, bas, -114.331645181, 1.0e-8, configure)

        assert scf_drv.pfon_temperature == 0
        assert True in calls

    def test_rijk_hf_openshell_with_density_damping(self, monkeypatch):

        mol, bas = self._get_methanol_radical_data()

        calls = self._record_density_factor_calls(monkeypatch)

        def configure(drv):
            drv.density_damping = True
            drv.acc_type = 'diis'

        self.run_scf_with_modifier(
            'unrestricted', mol, bas, -114.331645181, 1.0e-8, configure)

        assert True in calls

    def test_rijk_hf(self):

        xyz_string = """6
        xyz
        H      1.2001      0.0363      0.8431
        C      0.7031      0.0083     -0.1305
        H      0.9877      0.8943     -0.7114
        H      1.0155     -0.8918     -0.6742
        O     -0.6582     -0.0067      0.1730
        H     -1.1326     -0.0311     -0.6482
        """

        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)

        self.run_scf('restricted', mol, bas, None, -114.954012228, 1.0e-8)

    def test_rijk_hf_openshell(self):

        xyz_string = """5
        xyz
        H      1.2001      0.0363      0.8431
        C      0.7031      0.0083     -0.1305
        H      0.9877      0.8943     -0.7114
        H      1.0155     -0.8918     -0.6742
        O     -0.6582     -0.0067      0.1730
        """

        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_multiplicity(2)

        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)

        self.run_scf('unrestricted', mol, bas, None, -114.331645181, 1.0e-8)
