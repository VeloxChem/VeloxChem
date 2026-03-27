from mpi4py import MPI
from pathlib import Path
import pytest
import sys

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.outputstream import OutputStream
from veloxchem.scfrestdriver import ScfRestrictedDriver


@pytest.mark.solvers
class TestScfEnergyComponents:

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

    def run_scf_energy_components(self,
                                  tmp_path,
                                  mol,
                                  bas,
                                  solvation_model=None,
                                  potfile=None,
                                  vdwfile=None):

        fname = str(tmp_path / 'vlx_energy_components.out')

        comm = MPI.COMM_WORLD
        ostream = OutputStream.create_mpi_ostream(comm, fname)

        scf_drv = ScfRestrictedDriver(comm, ostream)
        scf_drv.solvation_model = solvation_model
        scf_drv.point_charges = potfile
        scf_drv.qm_vdw_params = vdwfile
        scf_results_not_used = scf_drv.compute(mol, bas)

        energy_keys = [
            'Total Energy',
            'Electronic Energy',
            'Nuclear Repulsion Energy',
            'Nuclei-Point Charges Energy',
            'Electrostatic Solvation Energy',
            'SMD Solvation Energy',
            '... ENP contribution',
            '... CDS contribution',
        ]

        energy_components = {}

        if scf_drv.rank == mpi_master():
            with open(fname, 'r') as fh:
                for line in fh:
                    for key in energy_keys:
                        if key in line and ':' in line:
                            content = line.split(':')
                            assert content[0].strip() == key
                            energy_components[key] = float(
                                content[1].split()[0])

            assert 'Total Energy' in energy_components

            total_energy = energy_components['Total Energy']
            sum_of_components = sum([
                val for key, val in energy_components.items()
                if key != 'Total Energy' and not key.startswith('...')
            ])
            assert abs(total_energy - sum_of_components) < 1.0e-9

            if 'SMD Solvation Energy' in energy_components:
                sum_of_smd_components = sum([
                    val for key, val in energy_components.items()
                    if key.startswith('...')
                ])
                assert abs(energy_components['SMD Solvation Energy'] -
                           sum_of_smd_components) < 1.0e-9

    def test_scf_with_point_charges_and_vdw(self, tmp_path):

        mol, bas = self.get_molecule_and_basis()

        here = Path(__file__).parent
        potfile = str(here / 'data' / 'pe_water.pot')
        vdwfile = str(here / 'data' / 'pe_water.qm_vdw_params.txt')

        self.run_scf_energy_components(tmp_path,
                                       mol,
                                       bas,
                                       potfile=potfile,
                                       vdwfile=vdwfile)

    def test_scf_with_point_charges(self, tmp_path):

        mol, bas = self.get_molecule_and_basis()

        here = Path(__file__).parent
        potfile = str(here / 'data' / 'pe_water.pot')

        self.run_scf_energy_components(tmp_path, mol, bas, potfile=potfile)

    def test_scf_with_cpcm(self, tmp_path):

        mol, bas = self.get_molecule_and_basis()

        self.run_scf_energy_components(tmp_path,
                                       mol,
                                       bas,
                                       solvation_model='cpcm')

    @pytest.mark.skipif("rdkit" not in sys.modules,
                        reason="rdkit not available")
    def test_scf_with_smd(self, tmp_path):

        mol, bas = self.get_molecule_and_basis()

        self.run_scf_energy_components(tmp_path,
                                       mol,
                                       bas,
                                       solvation_model='smd')
