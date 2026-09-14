from pathlib import Path
import numpy as np
import pytest
from mpi4py import MPI

from veloxchem.errorhandler import VeloxChemError
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lrsolver import LinearResponseSolver
from veloxchem.polarizabilitygradient import PolarizabilityGradient
from veloxchem.veloxchemlib import mpi_master


@pytest.mark.solvers
class TestPolgradEnvironmentGap:
    """The analytical polarizability gradient is vacuum-only: all five
    environment settings are rejected at compute entry, while the numerical
    path stays implicitly consistent (docs/environment-gap.md)."""

    def setup_method(self):

        xyz_string = """3
        xyz
        O   0.000000    0.000000    0.000000
        H   0.758602    0.000000   -0.504284
        H   0.758602    0.000000    0.504284
        """

        self.molecule = Molecule.read_xyz_string(xyz_string)
        self.basis = MolecularBasis.read(self.molecule, 'sto-3g')

    def run_scf_and_lr(self, **env_settings):

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        for key, value in env_settings.items():
            setattr(scf_drv, key, value)
        scf_results = scf_drv.compute(self.molecule, self.basis)

        rsp_settings = {'conv_thresh': 1.0e-5, 'frequencies': (0.0,)}
        lr_drv = LinearResponseSolver()
        lr_drv.a_operator = 'electric dipole'
        lr_drv.b_operator = 'electric dipole'
        lr_drv.update_settings(rsp_settings)
        lr_drv.ostream.mute()
        lr_results = lr_drv.compute(self.molecule, self.basis, scf_results)

        return scf_drv, lr_results

    def run_analytical_polgrad(self, scf_drv, lr_results):

        an_drv = PolarizabilityGradient(scf_drv)
        an_drv.update_settings({'frequencies': (0.0,)},
                               {'conv_thresh': 2.0e-7})
        an_drv.ostream.mute()
        return an_drv.compute(self.molecule, self.basis,
                              scf_drv.scf_results, lr_results)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='pytest.raises only valid in serial')
    def test_compute_analytical_rejects_electric_field(self):

        scf_drv, lr_results = self.run_scf_and_lr(
            electric_field=[0.0, 0.0, 0.001])

        errmsg = 'Analytical polarizability gradient is not supported '
        errmsg += 'with a static electric field'
        with pytest.raises(VeloxChemError, match=errmsg):
            self.run_analytical_polgrad(scf_drv, lr_results)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='pytest.raises only valid in serial')
    def test_compute_analytical_rejects_point_charges(self):

        here = Path(__file__).parent
        potfile = str(here / 'data' / 'pe_water.pot')

        scf_drv, lr_results = self.run_scf_and_lr(point_charges=potfile)

        errmsg = 'Analytical polarizability gradient is not supported '
        errmsg += 'with point charges'
        with pytest.raises(VeloxChemError, match=errmsg):
            self.run_analytical_polgrad(scf_drv, lr_results)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='pytest.raises only valid in serial')
    def test_compute_analytical_rejects_potfile(self):

        pytest.importorskip('pyframe')
        here = Path(__file__).parent
        potfile = str(here / 'data' / 'pe_water.json')

        scf_drv, lr_results = self.run_scf_and_lr(potfile=potfile)

        errmsg = 'Analytical polarizability gradient is not supported '
        errmsg += 'with polarizable embedding'
        with pytest.raises(VeloxChemError, match=errmsg):
            self.run_analytical_polgrad(scf_drv, lr_results)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='pytest.raises only valid in serial')
    def test_compute_analytical_rejects_solvation_model(self):

        scf_drv, lr_results = self.run_scf_and_lr(solvation_model='cpcm',
                                                  cpcm_epsilon=10.0)

        errmsg = 'Analytical polarizability gradient is not supported '
        errmsg += 'with CPCM/SMD solvation'
        with pytest.raises(VeloxChemError, match=errmsg):
            self.run_analytical_polgrad(scf_drv, lr_results)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='pytest.raises only valid in serial')
    def test_compute_analytical_rejects_pressure(self):

        scf_drv, lr_results = self.run_scf_and_lr(pressure=20000.0)

        errmsg = 'Analytical polarizability gradient is not supported '
        errmsg += 'with GOSTSHYP'
        with pytest.raises(VeloxChemError, match=errmsg):
            self.run_analytical_polgrad(scf_drv, lr_results)

    def test_compute_analytical_accepts_vacuum(self):

        scf_drv, lr_results = self.run_scf_and_lr()

        polgrad_results = self.run_analytical_polgrad(scf_drv, lr_results)

        if scf_drv.rank == mpi_master():
            assert 0.0 in polgrad_results
            polgrad = polgrad_results[0.0].reshape(3, 3, 3, 3)
            assert np.max(np.abs(polgrad)) > 0.0

    def test_compute_numerical_accepts_environment(self):
        """The numerical path re-runs SCF + LR per displacement and stays
        implicitly consistent, so it must NOT reject environments."""

        here = Path(__file__).parent
        potfile = str(here / 'data' / 'pe_water.pot')

        scf_drv, _ = self.run_scf_and_lr(point_charges=potfile)

        num_drv = PolarizabilityGradient(scf_drv)
        num_drv.update_settings({'numerical': 'yes', 'frequencies': (0.0,)},
                                {})
        num_drv.ostream.mute()
        polgrad_results = num_drv.compute(self.molecule, self.basis,
                                          scf_drv.scf_results, None)

        if scf_drv.rank == mpi_master():
            assert 0.0 in polgrad_results
            polgrad = polgrad_results[0.0].reshape(3, 3, 3, 3)
            assert np.max(np.abs(polgrad)) > 0.0
