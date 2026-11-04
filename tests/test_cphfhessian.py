from collections import namedtuple
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import h5py
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.hessianorbitalresponse import HessianOrbitalResponse
from veloxchem.hessianorbitalresponseunrest import (
    UnrestrictedHessianOrbitalResponse)
import veloxchem.batchsize as batchsize
from veloxchem.linearsolver import LinearSolver


@pytest.mark.solvers
class TestCphfSolver:

    def run_cphfsolver(self,
                       molecule,
                       basis,
                       xcfun=None,
                       label=None,
                       use_subcomms=False,
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

        hess_orbrsp_drv = HessianOrbitalResponse()
        orbrsp_settings = {'conv_thresh': 2e-7}
        hess_orbrsp_drv.update_settings(orbrsp_settings)
        hess_orbrsp_drv.ostream.mute()

        hess_orbrsp_drv.use_subcomms = use_subcomms

        # TODO: hess_orbrsp_drv should return cphf_results
        hess_orbrsp_drv.compute(molecule, basis, scf_results)

        dist_cphf_coefficients = hess_orbrsp_drv.cphf_results['dist_cphf_ov']
        dof = len(dist_cphf_coefficients)

        cphf_coefficients = []
        for x in range(dof):
            solution_vec = dist_cphf_coefficients[x].get_full_vector(0)
            if scf_drv.rank == mpi_master():
                cphf_coefficients.append(solution_vec)

        if scf_drv.rank == mpi_master():
            here = Path(__file__).parent
            hf_file_name = str(here / 'data' / {
                False: 'cphf_coefficients.h5',
                True: 'cphf_coefficients_ri.h5'
            }[ri])
            hf = h5py.File(hf_file_name, 'r')
            cphf_reference = np.array(hf.get(label))
            hf.close()

            cphf_coefficients = np.array(cphf_coefficients)
            cphf_coefficients = cphf_coefficients.reshape(cphf_reference.shape)

            # Here we are comparing the CPHF coefficients in MO basis, so
            # there might be sign differences; we compare absolute values instead.
            assert np.max(
                np.abs(np.abs(cphf_coefficients) -
                       np.abs(cphf_reference))) < 1.0e-6

    @pytest.mark.parametrize('ri', [False, True])
    def test_cphf_coefficients(self, ri):
        nh3_xyz = """4

        N     0.000000000     0.000000000     0.000000000
        H    -0.653401663     0.309213352     0.817609879
        H     0.695693936     0.071283622    -0.702632331
        H     0.330018952    -0.826347607     0.052270023
        """
        basis_set_label = "sto-3g"

        molecule = Molecule.from_xyz_string(nh3_xyz)
        basis = MolecularBasis.read(molecule, basis_set_label)

        self.run_cphfsolver(molecule, basis, "hf", "cphf_coefficients", ri=ri)

        self.run_cphfsolver(molecule,
                            basis,
                            "hf",
                            "cphf_coefficients",
                            use_subcomms=True,
                            ri=ri)

    @pytest.mark.parametrize('ri', [False, True])
    def test_cpks_coefficients(self, ri):
        nh3_xyz = """4

        N     0.000000000     0.000000000     0.000000000
        H    -0.653401663     0.309213352     0.817609879
        H     0.695693936     0.071283622    -0.702632331
        H     0.330018952    -0.826347607     0.052270023
        """
        basis_set_label = "sto-3g"

        molecule = Molecule.from_xyz_string(nh3_xyz)
        basis = MolecularBasis.read(molecule, basis_set_label)

        self.run_cphfsolver(molecule,
                            basis,
                            "b3lyp",
                            "cpks_coefficients_b3lyp",
                            ri=ri)

        self.run_cphfsolver(molecule,
                            basis,
                            "b3lyp",
                            "cpks_coefficients_b3lyp",
                            use_subcomms=True,
                            ri=ri)

    @pytest.mark.parametrize('ri', [False, True])
    def test_cphf_trial_orthogonality(self, ri):
        nh3_xyz = """4

        N     0.000000000     0.000000000     0.000000000
        H    -0.653401663     0.309213352     0.817609879
        H     0.695693936     0.071283622    -0.702632331
        H     0.330018952    -0.826347607     0.052270023
        """
        basis_set_label = "sto-3g"

        molecule = Molecule.from_xyz_string(nh3_xyz)
        basis = MolecularBasis.read(molecule, basis_set_label)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        cphf = HessianOrbitalResponse()
        cphf.update_settings({'conv_thresh': 1.0e-6})
        cphf.ostream.mute()
        cphf.compute(molecule, basis, scf_results)

        assert cphf.is_converged

        # matmul_AtB is collective; the Gram matrix is only returned on master.
        gram = cphf.dist_trials.matmul_AtB(cphf.dist_trials)

        if cphf.rank == mpi_master():
            assert gram.shape[0] > 0
            assert np.max(np.abs(gram - np.eye(gram.shape[0]))) < 1.0e-10

    @pytest.mark.parametrize('unrestricted', [False, True])
    def test_little_free_memory(self, monkeypatch, unrestricted):
        """With so little free memory that get_batch_size gives one Fock
        matrix, the exchange-correlation part still takes one atom a batch
        and gives what it gives with the machine's memory."""

        molecule = Molecule.read_molecule_string(
            """
            O  0.0           0.0  0.0
            H   .7586020000  0.0  -.5042840000
            H   .7586020000  0.0   .5042840000
            """,
            units='angstrom')
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)

        scf_drv = (ScfUnrestrictedDriver()
                   if unrestricted else ScfRestrictedDriver())
        scf_drv.xcfun = 'b3lyp'
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        def run():
            cphf = (UnrestrictedHessianOrbitalResponse()
                    if unrestricted else HessianOrbitalResponse())
            cphf.update_settings({'conv_thresh': 1.0e-7})
            cphf.ostream.mute()
            cphf.compute(molecule, basis, scf_results)
            coefficients = cphf.cphf_results['dist_cphf_ov']
            vectors = [c.get_full_vector(0) for c in coefficients]
            return (np.array(vectors)
                    if cphf.rank == mpi_master() else None)

        whole = run()

        # NOTE: only the name psutil in batchsize is replaced, so the RI
        # memory budgets, which read the free memory elsewhere, are not
        # touched.
        memory = namedtuple('memory', ['total', 'available'])
        monkeypatch.setattr(
            batchsize, 'psutil',
            SimpleNamespace(virtual_memory=lambda: memory(64 * 1024**3, 1)))

        little = run()

        if scf_drv.rank == mpi_master():
            assert np.max(np.abs(whole - little)) < 1.0e-10

    @staticmethod
    def rebuilt(factors):
        """The densities the factors give: a left factor and the right
        factors of the densities, or two of them for densities block diagonal
        in the occupied and the virtual orbitals."""

        if len(factors) == 2:
            left, rights = factors
            return [np.matmul(left, r.T) for r in rights]

        left, rights, left_vir, rights_vir = factors
        return [
            np.matmul(left, r.T) + np.matmul(left_vir, rv.T)
            for r, rv in zip(rights, rights_vir)
        ]

    @pytest.mark.parametrize('unrestricted', [False, True])
    def test_ri_jk_builds_take_factors(self, monkeypatch, unrestricted):
        """With RI-JK every Fock build of the Hessian orbital response, those
        of uij on the right-hand side included, takes the factors of its
        densities, the factors give the densities, and the build agrees with
        a four-center build of them to the accuracy of the fitting."""

        molecule = Molecule.read_molecule_string(
            """
            O  0.0           0.0  0.0
            H   .7586020000  0.0  -.5042840000
            H   .7586020000  0.0   .5042840000
            """,
            units='angstrom')
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)

        scf_drv = (ScfUnrestrictedDriver()
                   if unrestricted else ScfRestrictedDriver())
        scf_drv.xcfun = 'b3lyp'
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        name = ('_comp_lr_fock_unrestricted'
                if unrestricted else '_comp_lr_fock')
        real_build = getattr(LinearSolver, name)

        calls = []

        def recorded(self, dens, *args, dens_factors=None, **kwargs):
            # NOTE: a copy of the lists, which some callers empty after the
            # build.
            # NOTE: the unrestricted build takes a pair of lists on every
            # rank, of None on all but the master rank.
            if self.rank == mpi_master():
                copied = ([list(d) for d in dens]
                          if unrestricted else list(dens))
            else:
                copied = (None, None) if unrestricted else None

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

        monkeypatch.setattr(LinearSolver, name, recorded)

        cphf = (UnrestrictedHessianOrbitalResponse()
                if unrestricted else HessianOrbitalResponse())
        cphf.update_settings({'conv_thresh': 1.0e-7})
        cphf.ostream.mute()
        cphf.compute(molecule, basis, scf_results)

        if scf_drv.rank == mpi_master():
            assert calls

            for uses_ri_jk, dens, factors, fock, four_center in calls:
                assert uses_ri_jk
                assert factors is not None

                if unrestricted:
                    spins = list(zip(dens, factors))
                else:
                    spins = [(dens, factors)]

                for spin_dens, spin_factors in spins:
                    from_factors = self.rebuilt(spin_factors)
                    assert len(from_factors) == len(spin_dens)
                    for density, rebuilt in zip(spin_dens, from_factors):
                        assert np.max(np.abs(density - rebuilt)) < 1.0e-12

                # NOTE: the RI-JK build against the four-center build of the
                # same densities, relative to the largest Fock element. Here
                # the uij builds differ by 7e-5 to 1.1e-4 and those of the
                # iterations by up to 3.8e-3; a build which read the factors
                # of an occupied-by-virtual pair as L R^T + R L^T differed by
                # 1.0.
                fock = np.array(fock)
                four_center = np.array(four_center)
                largest = np.max(np.abs(four_center))
                assert np.max(np.abs(fock - four_center)) < 5.0e-2 * largest
