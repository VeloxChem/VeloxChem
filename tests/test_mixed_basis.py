import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


@pytest.mark.solvers
class TestMixedBasis:

    def run_scf(self, molstr, molstr_units, basis_label):

        molecule = Molecule.read_molecule_string(molstr, units=molstr_units)
        basis = MolecularBasis.read(molecule, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        if scf_drv.rank == mpi_master():
            scf_energy = scf_results['scf_energy']
        else:
            scf_energy = None
        scf_energy = scf_drv.comm.bcast(scf_energy, root=mpi_master())

        return scf_energy

    def test_mixed_basis(self):

        nh3_str = """
        N   -3.710    3.019   -0.037   def2-svpd
        H   -3.702    4.942    0.059
        H   -4.704    2.415    1.497
        H   -4.780    2.569   -1.573
        """
        nh3_str = nh3_str.strip()

        ch4_str = """
        C   -1.621   -5.080    0.444   def2-svp
        H   -0.819   -6.698   -0.465
        H   -3.412   -4.654   -0.393
        H   -0.381   -3.498    0.222
        H   -1.872   -5.468    2.413
        """
        ch4_str = ch4_str.strip()

        molstr_units = 'bohr'
        basis_label = 'sto-3g'
        tol = 1.0e-7

        e_mixed_basis = self.run_scf(nh3_str + '\n' + ch4_str, molstr_units,
                                     basis_label)

        assert abs(e_mixed_basis - (-96.27701050)) < tol
