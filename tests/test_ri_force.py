from mpi4py import MPI
import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfdriver import ScfDriver
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfgradientdriver import ScfGradientDriver
from veloxchem.tdaeigensolver import TdaEigenSolver
from veloxchem.tddftgradientdriver import TddftGradientDriver


@pytest.mark.solvers
class TestForceRi:
    """force_ri on the SCF keeps RI in the SCF and in what is formed on it,
    whatever memory they are predicted to need, and the numerical gradient
    keeps the choice the SCF made at the reference geometry."""

    @staticmethod
    def get_water():

        molecule = Molecule.read_str(
            """
            O 0.0000000000    0.0000000000   -0.0254395383
            H 0.0000000000    0.7695699584    0.5948147012
            H 0.0000000000   -0.7695699584    0.5948147012
            """,
            units='angstrom')
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)

        return molecule, basis

    @staticmethod
    def run_scf(molecule, basis, xcfun_label, force_ri=False,
                ri_memory_budget=None):

        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = xcfun_label
        scf_drv.conv_thresh = 1.0e-8
        scf_drv.force_ri = force_ri
        scf_drv.ri_memory_budget = ri_memory_budget
        scf_drv.ostream.mute()
        scf_drv.compute(molecule, basis)

        return scf_drv

    def test_keyword(self):

        scf_drv = ScfRestrictedDriver()
        assert not scf_drv.force_ri

        scf_drv.update_settings({}, {'force_ri': 'yes'})
        assert scf_drv.force_ri

    @pytest.mark.parametrize('xcfun_label, ri_label', [('b3lyp', 'RI-JK'),
                                                       ('blyp', 'RI-J')])
    def test_force_ri_with_budgets_of_1_kb(self, monkeypatch, xcfun_label,
                                           ri_label):

        molecule, basis = self.get_water()

        # NOTE: budgets of 1 kB, which no RI driver fits in: without force_ri
        # every one of them falls back to the four-center integrals.
        plain = self.run_scf(molecule, basis, xcfun_label)
        forced = self.run_scf(molecule, basis, xcfun_label, True, 1.0e-6)

        assert forced._uses_ri_jk() or forced._uses_ri_coulomb()
        assert not forced._ri_fallback

        energies = []
        for scf_drv in [plain, forced]:
            rsp_drv = TdaEigenSolver()
            rsp_drv.nstates = 3
            rsp_drv.conv_thresh = 1.0e-6
            if scf_drv is forced:
                rsp_drv.ri_memory_budget = 1.0e-6
            rsp_drv.ostream.mute()
            rsp_results = rsp_drv.compute(molecule, basis,
                                          scf_drv.scf_results)
            assert rsp_drv._uses_ri_jk() or rsp_drv._uses_ri_coulomb()
            energies.append(rsp_results['eigenvalues'] if rsp_drv.rank ==
                            mpi_master() else None)

        gradients = []
        for scf_drv in [plain, forced]:
            grad_drv = ScfGradientDriver(scf_drv)
            if scf_drv is forced:
                grad_drv.ri_gradient_memory_budget = 1.0e-6
            grad_drv.ostream.mute()
            grad_drv.compute(molecule, basis, scf_drv.scf_results)
            assert not grad_drv._ri_fallback
            gradients.append(grad_drv.get_gradient())

        chosen = []
        real_choose = TddftGradientDriver._choose_ri_gradient

        def recorded(self, *args, **kwargs):
            label = real_choose(self, *args, **kwargs)
            chosen.append(label)
            return label

        monkeypatch.setattr(TddftGradientDriver, '_choose_ri_gradient',
                            recorded)

        rsp_drv = TdaEigenSolver()
        rsp_drv.nstates = 3
        rsp_drv.conv_thresh = 1.0e-5
        rsp_drv.ri_memory_budget = 1.0e-6
        rsp_drv.ostream.mute()
        rsp_results = rsp_drv.compute(molecule, basis, forced.scf_results)

        tddft_grad = TddftGradientDriver(forced)
        tddft_grad.update_settings({}, {'tamm_dancoff': 'yes'},
                                   {'conv_thresh': 1.0e-7},
                                   {'xcfun': xcfun_label})
        tddft_grad.ri_gradient_memory_budget = 1.0e-6
        tddft_grad.ostream.mute()
        tddft_grad.compute(molecule, basis, forced, rsp_drv, rsp_results)

        assert chosen == [ri_label]

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert abs(forced.get_scf_energy() -
                       plain.get_scf_energy()) < 1.0e-8
            assert np.max(np.abs(energies[1] - energies[0])) < 1.0e-6
            assert np.max(np.abs(gradients[1] - gradients[0])) < 1.0e-6

    def numerical_gradient(self, scf_drv, molecule, basis):

        grad_drv = ScfGradientDriver(scf_drv)
        grad_drv.numerical = True
        grad_drv.ostream.mute()
        grad_drv.compute(molecule, basis, scf_drv.scf_results)

        return grad_drv.get_gradient()

    @pytest.mark.timeconsuming
    def test_numerical_gradient_keeps_the_choice_of_the_reference(
            self, monkeypatch):

        molecule, basis = self.get_water()

        # the reference choices with the memory as it is: RI, and with a
        # budget of 1 kB at the reference geometry, the four-center integrals
        with_ri = self.run_scf(molecule, basis, 'hf')
        ref_ri = self.numerical_gradient(with_ri, molecule, basis)

        four_center = self.run_scf(molecule, basis, 'hf')
        four_center.ri = False
        four_center.ri_jk = False
        four_center.ri_coulomb = False
        four_center.compute(molecule, basis)
        ref_four_center = self.numerical_gradient(four_center, molecule,
                                                  basis)

        # NOTE: a reference which used RI, and displaced geometries which
        # would each fall back: they keep RI.
        kept_ri = self.run_scf(molecule, basis, 'hf')
        assert kept_ri._uses_ri_jk()

        with monkeypatch.context() as patch:
            patch.setattr(ScfDriver, '_get_ri_memory_budget',
                          lambda self, ao_basis: 1000)
            grad_kept_ri = self.numerical_gradient(kept_ri, molecule, basis)

        assert not kept_ri.force_ri

        # NOTE: a reference which fell back, and displaced geometries which
        # would each use RI: they keep the four-center integrals.
        kept_fallback = self.run_scf(molecule, basis, 'hf', False, 1.0e-6)
        assert kept_fallback._ri_fallback
        kept_fallback.ri_memory_budget = None

        grad_kept_fallback = self.numerical_gradient(kept_fallback, molecule,
                                                     basis)

        assert kept_fallback.ri_memory_budget is None

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            # NOTE: the two references differ by the fitting error divided by
            # the step, which is what mixing the two would bring in.
            assert np.max(np.abs(ref_ri - ref_four_center)) > 1.0e-5
            assert np.max(np.abs(grad_kept_ri - ref_ri)) < 1.0e-7
            assert np.max(np.abs(grad_kept_fallback -
                                 ref_four_center)) < 1.0e-7
