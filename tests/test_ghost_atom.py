import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


@pytest.mark.solvers
class TestGhostAtom:

    def run_scf(self, molstr, molstr_units, basis_label, ri):

        molecule = Molecule.read_molecule_string(molstr, units=molstr_units)
        basis = MolecularBasis.read(molecule, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        if scf_drv.rank == mpi_master():
            scf_energy = scf_results['scf_energy']
        else:
            scf_energy = None
        scf_energy = scf_drv.comm.bcast(scf_energy, root=mpi_master())

        return scf_energy

    # NOTE: with RI, Hartree-Fock runs RI-JK and is compared with the
    # four-center reference within the error of the fitting. The basis set
    # labels of the atoms in the mixed basis test once became the fitting basis
    # of those atoms as well, which was off by a quarter of a hartree.
    @pytest.mark.parametrize('ri, tol', [(False, 1.0e-8), (True, 1.0e-4)])
    def test_ghost_atom(self, ri, tol):

        nh3_str = """
        N   -3.710    3.019   -0.037
        H   -3.702    4.942    0.059
        H   -4.704    2.415    1.497
        H   -4.780    2.569   -1.573
        """
        nh3_str = nh3_str.strip()

        nh3_bq_str = nh3_str.replace('N', 'Bq_N').replace('H', 'Bq_H')

        ch4_str = """
        C   -1.621   -5.080    0.444
        H   -0.819   -6.698   -0.465
        H   -3.412   -4.654   -0.393
        H   -0.381   -3.498    0.222
        H   -1.872   -5.468    2.413
        """
        ch4_str = ch4_str.strip()

        ch4_bq_str = ch4_str.replace('C', 'Bq_C').replace('H', 'Bq_H')

        ref_energies = {
            'nh3_in_dimer': -56.1487397917,
            'ch4_in_dimer': -40.1681040195,
        }

        molstr_units = 'bohr'
        basis_label = 'def2-svp'

        e_nh3_in_dimer = self.run_scf(nh3_str + '\n' + ch4_bq_str, molstr_units,
                                      basis_label, ri)
        e_ch4_in_dimer = self.run_scf(nh3_bq_str + '\n' + ch4_str, molstr_units,
                                      basis_label, ri)

        calc_energies = {
            'nh3_in_dimer': e_nh3_in_dimer,
            'ch4_in_dimer': e_ch4_in_dimer,
        }

        for key in calc_energies:
            assert abs(calc_energies[key] - ref_energies[key]) < tol

    @pytest.mark.parametrize('ri, tol', [(False, 1.0e-8), (True, 1.0e-4)])
    def test_ghost_atom_mixed_basis(self, ri, tol):

        nh3_str = """
        N   -3.710    3.019   -0.037  def2-svpd
        H   -3.702    4.942    0.059
        H   -4.704    2.415    1.497
        H   -4.780    2.569   -1.573
        """
        nh3_str = nh3_str.strip()

        nh3_bq_str = nh3_str.replace('N', 'N_Bq').replace('H', 'H_Bq')

        ch4_str = """
        C   -1.621   -5.080    0.444  def2-svp
        H   -0.819   -6.698   -0.465
        H   -3.412   -4.654   -0.393
        H   -0.381   -3.498    0.222
        H   -1.872   -5.468    2.413
        """
        ch4_str = ch4_str.strip()

        ch4_bq_str = ch4_str.replace('C', 'C_Bq').replace('H', 'H_Bq')

        ref_energies = {
            'nh3_in_dimer': -56.1293087216,
            'ch4_in_dimer': -40.1477356376,
        }

        molstr_units = 'bohr'
        basis_label = 'sto-3g'

        e_nh3_in_dimer = self.run_scf(nh3_str + '\n' + ch4_bq_str, molstr_units,
                                      basis_label, ri)
        e_ch4_in_dimer = self.run_scf(nh3_bq_str + '\n' + ch4_str, molstr_units,
                                      basis_label, ri)

        calc_energies = {
            'nh3_in_dimer': e_nh3_in_dimer,
            'ch4_in_dimer': e_ch4_in_dimer,
        }

        for key in calc_energies:
            assert abs(calc_energies[key] - ref_energies[key]) < tol
