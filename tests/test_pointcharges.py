from pathlib import Path

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


class TestPointCharges:

    @staticmethod
    def get_molecule_and_basis():

        xyz_string = """
        3
        xyz
        O    1.2361419   1.0137761  -0.0612424
        H    0.5104418   0.8944555   0.5514190
        H    1.9926927   1.1973129   0.4956931
        """
        mol = Molecule.read_xyz_string(xyz_string)
        basis_label = 'def2-svp'
        bas = MolecularBasis.read(mol, basis_label, ostream=None)
        return mol, bas

    def test_scf_with_point_charges(self):

        mol, bas = self.get_molecule_and_basis()

        here = Path(__file__).parent
        potfile = str(here / 'data' / 'pe_water.pot')

        scf_drv = ScfRestrictedDriver()
        scf_drv._point_charges = potfile
        scf_drv.ostream.mute()

        scf_results = scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            assert abs(scf_results['scf_energy'] - (-75.9798409316)) < 1e-9
