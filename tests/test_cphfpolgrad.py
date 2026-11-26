from pathlib import Path
import numpy as np
import pytest
import h5py

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lrsolver import LinearResponseSolver
from veloxchem.cppsolver import ComplexResponseSolver
from veloxchem.polorbitalresponse import PolOrbitalResponse
from veloxchem.linearsolver import LinearSolver


class FullDimensionPolOrbitalResponse(PolOrbitalResponse):

    def compute_rhs(self, molecule, basis, scf_tensors, eri_dict, dft_dict,
                    pe_dict, lr_results):
        if self.is_complex:
            return self.compute_rhs_complex(molecule, basis, scf_tensors,
                                            eri_dict, dft_dict, pe_dict,
                                            lr_results)
        return self.compute_rhs_real(molecule, basis, scf_tensors, eri_dict,
                                     dft_dict, pe_dict, lr_results)

    def compute_omega(self, molecule, basis, scf_tensors, lr_results):
        if self.is_complex:
            return self.compute_omega_complex(molecule, basis, scf_tensors,
                                              lr_results)
        return self.compute_omega_real(molecule, basis, scf_tensors, lr_results)


@pytest.mark.solvers
class TestCphfPolgrad:

    def run_cphfpolgrad_real(self,
                             molecule,
                             basis,
                             xcfun=None,
                             label1=None,
                             label2=None,
                             solver_cls=PolOrbitalResponse,
                             *,
                             ri):

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.xcfun = xcfun
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        # linear response
        rsp_settings = {'conv_thresh': 1.0e-5, 'frequencies': (0.0, 0.4)}
        lr_drv = LinearResponseSolver()
        lr_drv.a_operator = "electric dipole"
        lr_drv.b_operator = "electric dipole"
        lr_drv.update_settings(rsp_settings)
        lr_drv.ostream.mute()
        lr_results = lr_drv.compute(molecule, basis, scf_results)

        # test real analytical gradient
        cphfpolgrad_solver = solver_cls()
        cphfpolgrad_settings = {'conv_thresh': 2e-7, 'frequencies': (0.0, 0.4)}
        cphfpolgrad_solver.update_settings(cphfpolgrad_settings)
        cphfpolgrad_solver.ostream.mute()
        cphfpolgrad_solver.compute(molecule, basis, scf_results, lr_results)
        cphfpolgrad_solver.compute_omega(molecule, basis, scf_results, lr_results)

        # get CPHF lambda coefficients from dist array
        dist_cphf_coefficients = cphfpolgrad_solver.cphf_results['dist_cphf_ov']
        dof = len(dist_cphf_coefficients)

        cphfpolgrad_coefficients = []
        for x in range(dof):
            solution_vec = dist_cphf_coefficients[x].get_full_vector(0)
            if scf_drv.rank == mpi_master():
                cphfpolgrad_coefficients.append(solution_vec)

        # get CPHF omega coefficients from dist. array
        dist_omega_coefficients = cphfpolgrad_solver.cphf_results['dist_omega_ao']
        dof = len(dist_omega_coefficients)

        cphfpolgrad_omega_coefficients = []
        for x in range(dof):
            solution_vec = dist_omega_coefficients[x].get_full_vector(0)
            if scf_drv.rank == mpi_master():
                cphfpolgrad_omega_coefficients.append(solution_vec)

        if scf_drv.rank == mpi_master():
            # cphfpolgrad_results = cphfpolgrad_solver.cphf_results
            # cphfpolgrad_omega_coefficients = np.array([cphfpolgrad_results[0.0]['omega_ao'],
            #                                            cphfpolgrad_results[0.4]['omega_ao']])

            here = Path(__file__).parent
            hf_file_name = str(
                here / 'data' / {
                    False: 'cphfpolgrad_coefficients.h5',
                    True: 'cphfpolgrad_coefficients_ri.h5'
                }[ri])
            hf = h5py.File(hf_file_name, 'r')
            cphfpolgrad_reference = np.array(hf.get(label1))
            cphfpolgrad_omega_reference = np.array(hf.get(label2))
            hf.close()

            cphfpolgrad_coefficients = np.array(cphfpolgrad_coefficients)
            cphfpolgrad_coefficients = cphfpolgrad_coefficients.reshape(cphfpolgrad_reference.shape)

            # Here we are comparing the CPHF coefficients in MO basis, so
            # there might be sign differences; we compare absolute values instead.
            assert np.max(np.abs(np.abs(cphfpolgrad_coefficients) - np.abs(cphfpolgrad_reference))) < 1.0e-6
            # NOTE: omega is in the AO basis and does not depend on the phases
            # of the molecular orbitals, so it is compared with its signs.
            assert np.max(np.abs(np.array(cphfpolgrad_omega_coefficients) -
                                 cphfpolgrad_omega_reference)) < 1.0e-6

    def run_cphfpolgrad_complex(self,
                                molecule,
                                basis,
                                xcfun=None,
                                label1=None,
                                label2=None,
                                solver_cls=PolOrbitalResponse,
                                *,
                                ri):

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.xcfun = xcfun
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        # linear response
        rsp_settings = {'conv_thresh': 1.0e-5, 'frequencies': (0.0, 0.4),
                        'damping': 0.5}
        lr_drv = ComplexResponseSolver()
        lr_drv.a_operator = "electric dipole"
        lr_drv.b_operator = "electric dipole"
        lr_drv.update_settings(rsp_settings)
        lr_drv.ostream.mute()
        lr_results = lr_drv.compute(molecule, basis, scf_results)

        # test complex analytical gradient
        cphfpolgrad_solver = solver_cls()
        cphfpolgrad_settings = {'conv_thresh': 2e-7, 'frequencies': (0.0, 0.4),
                                'is_complex': 'yes', 'damping': 0.5}
        cphfpolgrad_solver.update_settings(cphfpolgrad_settings)
        cphfpolgrad_solver.ostream.mute()
        cphfpolgrad_solver.compute(molecule, basis, scf_results, lr_results)
        cphfpolgrad_solver.compute_omega(molecule, basis, scf_results, lr_results)

        # get CPHF lambda coefficients from dist. array
        dist_cphf_coefficients = cphfpolgrad_solver.cphf_results['dist_cphf_ov']
        dof = len(dist_cphf_coefficients)

        cphfpolgrad_coefficients = []
        for x in range(dof):
            solution_vec = dist_cphf_coefficients[x].get_full_vector(0)
            if scf_drv.rank == mpi_master():
                cphfpolgrad_coefficients.append(solution_vec)

        # get CPHF omega coefficients from dist. array
        dist_omega_coefficients = cphfpolgrad_solver.cphf_results['dist_omega_ao']
        dof = len(dist_omega_coefficients)

        cphfpolgrad_omega_coefficients = []
        for x in range(dof):
            solution_vec = dist_omega_coefficients[x].get_full_vector(0)
            if scf_drv.rank == mpi_master():
                cphfpolgrad_omega_coefficients.append(solution_vec)

        if scf_drv.rank == mpi_master():
            # cphfpolgrad_results = cphfpolgrad_solver.cphf_results
            # cphfpolgrad_omega_coefficients = np.array([cphfpolgrad_results[0.0]['omega_ao'],
            #                                            cphfpolgrad_results[0.4]['omega_ao']])

            here = Path(__file__).parent
            hf_file_name = str(
                here / 'data' / {
                    False: 'cphfpolgrad_coefficients.h5',
                    True: 'cphfpolgrad_coefficients_ri.h5'
                }[ri])
            hf = h5py.File(hf_file_name, 'r')
            # no zero frequency: polorbrsp_sanity_check_1 drops it for complex
            cphfpolgrad_reference = np.array(hf.get(label1))
            cphfpolgrad_omega_reference = np.array(hf.get(label2))
            hf.close()

            cphfpolgrad_coefficients = np.array(cphfpolgrad_coefficients)
            cphfpolgrad_coefficients = cphfpolgrad_coefficients.reshape(cphfpolgrad_reference.shape)

            # Here we ar comparing the CPHF coefficients in MO basis, so
            # there might be sign differences; we compare absolute values instead.
            assert np.max(np.abs(np.abs(cphfpolgrad_coefficients) - np.abs(cphfpolgrad_reference))) < 1.0e-6
            # NOTE: omega is in the AO basis and does not depend on the phases
            # of the molecular orbitals, so it is compared with its signs.
            assert np.max(np.abs(np.array(cphfpolgrad_omega_coefficients) -
                                 cphfpolgrad_omega_reference)) < 1.0e-6

    @pytest.mark.parametrize(
        'solver_cls',
        [PolOrbitalResponse, FullDimensionPolOrbitalResponse],
        ids=['reduced', 'full'])
    @pytest.mark.parametrize('ri', [False, True])
    def test_cphfpolgrad_coefficients_real(self, solver_cls, ri):
        h2o_xyz = """3

        O     0.000000    0.000000    0.000000
        H     0.000000    0.504284    0.758602
        H     0.000000   -0.504284    0.758602
        """
        basis_set_label = "sto-3g"

        molecule = Molecule.from_xyz_string(h2o_xyz)
        basis = MolecularBasis.read(molecule, basis_set_label)

        self.run_cphfpolgrad_real(molecule,
                                  basis,
                                  None,
                                  "cphfpolgrad_coefficients_red_real",
                                  "cphfpolgrad_omega_coefficients_red_real",
                                  solver_cls,
                                  ri=ri)

    @pytest.mark.parametrize(
        'solver_cls',
        [PolOrbitalResponse, FullDimensionPolOrbitalResponse],
        ids=['reduced', 'full'])
    @pytest.mark.parametrize('ri', [False, True])
    def test_cphfpolgrad_coefficients_complex(self, solver_cls, ri):
        h2o_xyz = """3

        O     0.000000    0.000000    0.000000
        H     0.000000    0.504284    0.758602
        H     0.000000   -0.504284    0.758602
        """
        basis_set_label = "sto-3g"

        molecule = Molecule.from_xyz_string(h2o_xyz)
        basis = MolecularBasis.read(molecule, basis_set_label)

        self.run_cphfpolgrad_complex(
            molecule,
            basis,
            None,
            "cphfpolgrad_coefficients_red_complex",
            "cphfpolgrad_omega_coefficients_red_complex",
            solver_cls,
            ri=ri)

    @pytest.mark.parametrize(
        'solver_cls',
        [PolOrbitalResponse, FullDimensionPolOrbitalResponse],
        ids=['reduced', 'full'])
    @pytest.mark.parametrize('ri', [False, True])
    def test_cpkspolgrad_coefficients_real(self, solver_cls, ri):
        h2o_xyz = """3

        O     0.000000    0.000000    0.000000
        H     0.000000    0.504284    0.758602
        H     0.000000   -0.504284    0.758602
        """
        basis_set_label = "sto-3g"

        molecule = Molecule.from_xyz_string(h2o_xyz)
        basis = MolecularBasis.read(molecule, basis_set_label)

        self.run_cphfpolgrad_real(
            molecule,
            basis,
            "b3lyp",
            "cpkspolgrad_coefficients_red_b3lyp_real",
            "cpkspolgrad_omega_coefficients_red_b3lyp_real",
            solver_cls,
            ri=ri)

    @pytest.mark.parametrize(
        'solver_cls',
        [PolOrbitalResponse, FullDimensionPolOrbitalResponse],
        ids=['reduced', 'full'])
    @pytest.mark.parametrize('ri', [False, True])
    def test_cpkspolgrad_coefficients_complex(self, solver_cls, ri):
        h2o_xyz = """3

        O     0.000000    0.000000    0.000000
        H     0.000000    0.504284    0.758602
        H     0.000000   -0.504284    0.758602
        """
        basis_set_label = "sto-3g"

        molecule = Molecule.from_xyz_string(h2o_xyz)
        basis = MolecularBasis.read(molecule, basis_set_label)

        self.run_cphfpolgrad_complex(
            molecule,
            basis,
            "b3lyp",
            "cpkspolgrad_coefficients_red_b3lyp_complex",
            "cpkspolgrad_omega_coefficients_red_b3lyp_complex",
            solver_cls,
            ri=ri)

    @pytest.mark.parametrize(
        'solver_cls',
        [PolOrbitalResponse, FullDimensionPolOrbitalResponse],
        ids=['reduced', 'full'])
    @pytest.mark.parametrize('is_complex', [False, True])
    def test_ri_jk_builds_take_factors(self, monkeypatch, solver_cls,
                                       is_complex):
        """With RI-JK every Fock build of the orbital response, those of the
        right-hand side and of omega included, takes the factors of its
        densities, the factors give the densities, and the build agrees with
        a four-center build of them to the accuracy of the fitting."""

        molecule = Molecule.from_xyz_string("""3

        O     0.000000    0.000000    0.000000
        H     0.000000    0.504284    0.758602
        H     0.000000   -0.504284    0.758602
        """)
        basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = 'b3lyp'
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        settings = {'conv_thresh': 1.0e-5, 'frequencies': (0.0, 0.4)}
        if is_complex:
            settings['damping'] = 0.5
            lr_drv = ComplexResponseSolver()
        else:
            lr_drv = LinearResponseSolver()
        lr_drv.update_settings(settings)
        lr_drv.ostream.mute()
        lr_results = lr_drv.compute(molecule, basis, scf_results)

        calls = []
        real_build = LinearSolver._comp_lr_fock

        def recorded(self, dens, *args, dens_factors=None, **kwargs):
            # NOTE: a copy of the list, which some callers empty after the
            # build.
            copied = list(dens) if self.rank == mpi_master() else None

            fock = real_build(self, dens, *args, dens_factors=dens_factors,
                              **kwargs)

            # NOTE: the same densities built again from the four-center
            # integrals, which a build without factors is. Every rank takes
            # part, so this is decided by _uses_ri_jk, which every rank
            # answers alike, and not by the factors, which only the master
            # rank holds.
            if self._uses_ri_jk():
                four_center = real_build(self, copied, *args,
                                         dens_factors=None, **kwargs)
            else:
                four_center = None

            if self.rank == mpi_master():
                calls.append((self._uses_ri_jk(), copied, dens_factors, fock,
                              four_center))

            return fock

        monkeypatch.setattr(LinearSolver, '_comp_lr_fock', recorded)

        solver = solver_cls()
        settings = {'conv_thresh': 2e-7, 'frequencies': (0.0, 0.4)}
        if is_complex:
            settings.update({'is_complex': 'yes', 'damping': 0.5})
        solver.update_settings(settings)
        solver.ostream.mute()
        solver.compute(molecule, basis, scf_results, lr_results)
        solver.compute_omega(molecule, basis, scf_results, lr_results)

        if scf_drv.rank == mpi_master():
            assert calls

            for uses_ri_jk, dens, factors, fock, four_center in calls:
                assert uses_ri_jk
                assert factors is not None

                # NOTE: the RI-JK build against the four-center build of the
                # same densities. Here they differ by 7e-5 to 1.3e-4 of the
                # largest Fock element; factors in a form the build reads
                # differently give a difference of the size of the Fock
                # matrices.
                largest = max(np.max(np.abs(f)) for f in four_center)
                difference = max(
                    np.max(np.abs(f - g)) for f, g in zip(fock, four_center))
                assert difference < 1.0e-2 * largest

                # NOTE: a left factor and the right factors of the densities,
                # or two of them for densities block diagonal in the
                # occupied and the virtual orbitals.
                if len(factors) == 2:
                    left, rights = factors
                    rebuilt = [np.matmul(left, r.T) for r in rights]
                else:
                    left, rights, left_vir, rights_vir = factors
                    rebuilt = [
                        np.matmul(left, r.T) + np.matmul(left_vir, rv.T)
                        for r, rv in zip(rights, rights_vir)
                    ]

                assert len(rebuilt) == len(dens)
                for density, from_factors in zip(dens, rebuilt):
                    assert np.max(np.abs(density - from_factors)) < 1.0e-12
