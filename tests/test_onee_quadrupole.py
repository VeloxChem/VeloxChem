import numpy as np

from veloxchem.veloxchemlib import compute_quadrupole_integrals
from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


class TestOneElecIntsQuadrupole:

    def get_molecule_and_basis(self):

        mol_string = """
        N         -1.96309        1.59755       -0.01963
        H         -1.95876        2.61528        0.03109
        H         -2.48929        1.27814        0.79244
        H         -2.52930        1.35928       -0.83265
        """
        basis_label = 'def2-tzvpp'

        mol = Molecule.read_molecule_string(mol_string, units='bohr')
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        return mol, bas

    def test_quadrupole_integrals(self):

        mol, bas = self.get_molecule_and_basis()

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            total_density = scf_results['D_alpha'] + scf_results['D_beta']
        else:
            total_density = None
        total_density = scf_drv.comm.bcast(total_density, root=mpi_master())

        coords = mol.get_coordinates_in_bohr()
        nuclear_charges = mol.get_element_ids()

        # center of nuclear charge as origin

        origin = np.sum(coords.T * nuclear_charges,
                        axis=1) / np.sum(nuclear_charges)

        # electronic contribution

        quadrupole_ints = compute_quadrupole_integrals(mol, bas, origin)

        electronic_quadrupole = np.array(
            [np.sum(quadrupole_ints[d] * total_density) for d in range(6)])

        # nuclear contribution

        nuclear_quadrupole = np.zeros(6)

        for iatom in range(mol.number_of_atoms()):
            q = nuclear_charges[iatom]
            x, y, z = coords[iatom] - origin

            nuclear_quadrupole[0] += q * x * x
            nuclear_quadrupole[1] += q * x * y
            nuclear_quadrupole[2] += q * x * z
            nuclear_quadrupole[3] += q * y * y
            nuclear_quadrupole[4] += q * y * z
            nuclear_quadrupole[5] += q * z * z

        # quadrupole moment

        total_quadrupole = nuclear_quadrupole + electronic_quadrupole

        # reference values

        ref_Q_nuc = np.array(
            [0.479107, 0.357448, 0.038664, 1.173405, -0.016335, 1.322784])
        ref_Q_elec = np.array(
            [-5.613973, 0.055904, 0.006044, -5.505387, -0.002561, -5.482035])
        ref_Q_tot = ref_Q_nuc + ref_Q_elec

        assert np.max(np.abs(nuclear_quadrupole - ref_Q_nuc)) < 1.0e-5
        assert np.max(np.abs(electronic_quadrupole - ref_Q_elec)) < 1.0e-5
        assert np.max(np.abs(total_quadrupole - ref_Q_tot)) < 1.0e-5
