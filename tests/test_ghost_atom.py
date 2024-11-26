import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


@pytest.mark.solvers
class TestGhostAtom:

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

    def test_ghost_atom(self):

        nh3_str = """
        N   -3.710    3.019   -0.037
        H   -3.702    4.942    0.059
        H   -4.704    2.415    1.497
        H   -4.780    2.569   -1.573
        """
        nh3_str = nh3_str.strip()

        nh3_bq_str = nh3_str.replace('N', 'N_Bq').replace('H', 'H_Bq')

        ch4_str = """
        C   -1.621   -5.080    0.444
        H   -0.819   -6.698   -0.465
        H   -3.412   -4.654   -0.393
        H   -0.381   -3.498    0.222
        H   -1.872   -5.468    2.413
        """
        ch4_str = ch4_str.strip()

        ch4_bq_str = ch4_str.replace('C', 'C_Bq').replace('H', 'H_Bq')

        ref_energies = {
            'monomer_1': -56.14844212,
            'monomer_2': -40.16809185,
            'm1_in_dimer': -56.14873979,
            'm2_in_dimer': -40.16810402,
            'dimer': -96.31683373,
        }

        molstr_units = 'bohr'
        basis_label = 'def2-svp'
        tol = 1.0e-7

        e_monomer_1 = self.run_scf(nh3_str, molstr_units, basis_label)
        e_monomer_2 = self.run_scf(ch4_str, molstr_units, basis_label)

        e_m1_in_dimer = self.run_scf(nh3_str + '\n' + ch4_bq_str, molstr_units,
                                     basis_label)
        e_m2_in_dimer = self.run_scf(nh3_bq_str + '\n' + ch4_str, molstr_units,
                                     basis_label)

        e_dimer = self.run_scf(nh3_str + '\n' + ch4_str, molstr_units,
                               basis_label)

        calc_energies = {
            'monomer_1': e_monomer_1,
            'monomer_2': e_monomer_2,
            'm1_in_dimer': e_m1_in_dimer,
            'm2_in_dimer': e_m2_in_dimer,
            'dimer': e_dimer,
        }

        for key in calc_energies:
            assert abs(calc_energies[key] - ref_energies[key]) < tol
