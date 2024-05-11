from pathlib import Path

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


class TestNoradrenaline:

    def run_scf(self, mol, bas, min_bas, ref_scf_energy, tol_ene):

        scf_drv = ScfRestrictedDriver()
        scf_drv.filename = None
        scf_drv.checkpoint_file = None
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas, min_bas)

        if scf_drv.rank == 0:
            scf_energy = scf_results['scf_energy']
            assert abs(scf_energy - ref_scf_energy) < tol_ene

    def test_hf(self):

        here = Path(__file__).parent
        xyz_file = str(here / 'data' / 'noradrenaline.xyz')
        mol = Molecule.read_xyz_file(xyz_file)
        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)
        min_bas = MolecularBasis.read(mol, 'ao-start-guess', ostream=None)

        ref_scf_energy = -587.9305634780
        tol_ene = 1.0e-6

        self.run_scf(mol, bas, min_bas, ref_scf_energy, tol_ene)
