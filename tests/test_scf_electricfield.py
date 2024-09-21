import pytest

from veloxchem import mpi_master
from veloxchem import Molecule, MolecularBasis
from veloxchem import ScfRestrictedDriver


@pytest.mark.solvers
class TestScfInElectricField:

    def run_scf_electricfield(self, xcfun_label, ref_scf_energy, electric_field,
                              tol):

        mol_string = """
        O   0.0   0.0   0.0
        H   0.0   1.4   1.1
        H   0.0  -1.4   1.1
        """
        basis_label = 'aug-cc-pvdz'

        mol = Molecule.read_molecule_string(mol_string, units='bohr')
        mol.check_multiplicity()

        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_drv.electric_field = electric_field
        scf_results = scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            assert abs(ref_scf_energy - scf_results['scf_energy']) < tol

    def test_hf(self):

        self.run_scf_electricfield('hf', -76.0688447658, [0, 0, 0.03], 1.0e-6)
