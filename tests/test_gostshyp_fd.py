import numpy as np
import pytest

from veloxchem import ElectricDipoleIntegralsDriver
from veloxchem.lrsolver import LinearResponseSolver
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.molecule import Molecule
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.veloxchemlib import mpi_master


@pytest.mark.solvers
class TestGostshypFiniteField:
    """GOSTSHYP finite-field verification for water at 10 GPa.

    The polarizability is obtained in two independent ways: numerically
    from SCF energies and dipole moments at finite electric fields, and
    analytically from a linear response calculation. The two routes must
    agree with each other and with the reference values.
    """

    @pytest.mark.timeconsuming
    def test_water_def2svpd_at_10_gpa(self):

        xyz_string = """3

        O    0.000000000000        0.000000000000        0.000000000000
        H    0.000000000000        0.740848095288        0.582094932012
        H    0.000000000000       -0.740848095288        0.582094932012
        """
        molecule = Molecule.read_xyz_string(xyz_string)
        basis = MolecularBasis.read(molecule, 'def2-svpd', ostream=None)

        pressure = 10  # GPa

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.conv_thresh = 1.0e-8
        scf_drv.pressure_units = 'GPa'
        scf_drv.gostshyp_discretization = 'swig'
        scf_drv.gostshyp_num_lebedev_points = 110
        scf_drv.gostshyp_tssf = 1.2
        scf_drv.gostshyp_switching_thresh = 1.0e-8
        scf_drv.gostshyp_r_ext = 0.0

        # electric-dipole integrals
        dipole_drv = ElectricDipoleIntegralsDriver()
        dipole_mats = dipole_drv.compute(molecule, basis)

        mu_x = -1.0 * dipole_mats.x_to_numpy()
        mu_y = -1.0 * dipole_mats.y_to_numpy()
        mu_z = -1.0 * dipole_mats.z_to_numpy()

        def dipmom(D):
            # electronic part
            mu_e = np.array([
                np.einsum('ab, ab', D, mu_x),
                np.einsum('ab, ab', D, mu_y),
                np.einsum('ab, ab', D, mu_z),
            ])
            # nuclear part
            element_charges = {'H': 1.0, 'O': 8.0}
            mu_n = np.zeros(3)
            for A, label in enumerate(molecule.get_labels()):
                R_A = np.array(molecule.get_atom_coordinates(A))
                mu_n += element_charges[label] * R_A
            # total dipole moment along the field direction (z)
            return (mu_e + mu_n)[2]

        # five-point stencil for the numerical derivatives
        h = 0.001
        field_strengths = np.linspace(-2, 2, 5) * h

        E = np.zeros(len(field_strengths))
        mu = np.zeros(len(field_strengths))

        comm = scf_drv.comm
        rank = comm.Get_rank()

        for F_idx, F in enumerate(field_strengths):
            scf_drv.electric_field = [0.0, 0.0, F]
            scf_drv.pressure = pressure
            scf_results = scf_drv.compute(molecule, basis)

            if rank == mpi_master():
                E[F_idx] = scf_results['scf_energy']
                mu[F_idx] = dipmom(scf_results['D_alpha'] +
                                   scf_results['D_beta'])

            E[F_idx] = comm.bcast(E[F_idx], root=mpi_master())
            mu[F_idx] = comm.bcast(mu[F_idx], root=mpi_master())

        # numerical differentiation (five-point stencil)
        dE = np.zeros(2)
        dmu = np.zeros(2)
        dE[0] = (-E[4] + 8 * E[3] - 8 * E[1] + E[0]) / (12 * h)
        dE[1] = (-E[4] + 16 * E[3] - 30 * E[2] + 16 * E[1] - E[0]) / (12 * h**2)
        dmu[0] = (-mu[4] + 8 * mu[3] - 8 * mu[1] + mu[0]) / (12 * h)
        dmu[1] = (-mu[4] + 16 * mu[3] - 30 * mu[2] + 16 * mu[1] - mu[0]) / (
            12 * h**2)

        # linear response polarizability at zero field
        scf_drv.electric_field = [0.0, 0.0, 0.0]
        scf_results = scf_drv.compute(molecule, basis)

        lrf_drv = LinearResponseSolver()
        lrf_drv.ostream.mute()
        lrf_drv.conv_thresh = 1.0e-5
        lrf_drv.a_operator = 'electric dipole'
        lrf_drv.b_operator = 'electric dipole'
        lrf_drv.a_components = 'z'
        lrf_drv.b_components = 'z'
        lrf_drv.frequencies = [0.0]
        lrf_drv.pressure = pressure
        lrf_drv.pressure_units = 'GPa'
        lrf_drv.gostshyp_discretization = 'swig'
        lrf_drv.gostshyp_num_lebedev_points = 110
        lrf_drv.gostshyp_tssf = 1.2
        lrf_drv.gostshyp_switching_thresh = 1.0e-8
        lrf_drv.gostshyp_r_ext = 0.0

        lrf_results = lrf_drv.compute(molecule, basis, scf_results)

        if rank == mpi_master():
            alpha = -lrf_results['response_functions'][('z', 'z', 0.0)]
        else:
            alpha = None
        alpha = comm.bcast(alpha, root=mpi_master())

        if rank == mpi_master():
            # reference values from the GOSTSHYP finite-field verification
            # (printed to 6 and 5 decimals, respectively)
            ref_mu = 0.801786
            ref_alpha = 7.85045

            # analytical dipole moment and polarizability at 10 GPa
            assert abs(mu[2] - ref_mu) < 1.0e-6
            assert abs(alpha - ref_alpha) < 1.0e-5

            # energy finite-field vs linear response polarizabilities
            assert abs(alpha - (-dE[1])) < 1.0e-6

            # dipole moment finite-field vs linear response polarizabilities
            assert abs(dmu[0] - alpha) < 1.0e-6

            # the two finite-field polarizabilities
            assert abs(dmu[0] - (-dE[1])) < 1.0e-6

            # analytical vs finite-field dipole moments
            assert abs(mu[2] - (-dE[0])) < 1.0e-7
