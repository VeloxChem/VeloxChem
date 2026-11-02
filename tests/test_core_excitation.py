from mpi4py import MPI
import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.tdaeigensolver import TdaEigenSolver


@pytest.mark.solvers
class TestCoreExcitation:

    def run_scf(self, ri):

        water_xyz = """
        O   0.0   0.0   0.0
        H   0.0   1.4   1.1
        H   0.0  -1.4   1.1
        """
        molecule = Molecule.read_molecule_string(water_xyz, units='au')
        basis = MolecularBasis.read(molecule, 'def2-svpd', ostream=None)
        xcfun_label = 'b3lyp'

        # ground state
        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.xcfun = xcfun_label
        scf_drv.ostream.mute()
        scf_result = scf_drv.compute(molecule, basis)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_result['ri_jk'] or scf_result['ri_coulomb']
                assert not scf_result['ri_fallback']

        return molecule, basis, xcfun_label, scf_result

    @pytest.mark.parametrize('ri', [False, True])
    def test_core_excitation_rpa(self, ri):

        molecule, basis, xcfun_label, scf_result = self.run_scf(ri)

        # core excitation
        rpa_drv = LinearResponseEigenSolver()
        rpa_drv.xcfun = xcfun_label
        rpa_drv.core_excitation = True
        rpa_drv.num_core_orbitals = 1
        rpa_drv.nstates = 10
        rpa_drv.ostream.mute()
        rpa_result = rpa_drv.compute(molecule, basis, scf_result)

        ref_exc_ene = {
            False: np.array([
                19.09409612, 19.1579597, 19.20965781, 19.21608707, 19.25755683,
                19.28712811, 19.33009438, 19.33571862, 19.36632322, 19.44542822
            ]),
            True: np.array([
                19.09410513, 19.15800754, 19.20970148, 19.21606118, 19.25756739,
                19.28711518, 19.33015103, 19.33577924, 19.36633514, 19.44546557
            ])
        }[ri]
        ref_osc_str = {
            False: np.array([
                0.01, 0.0221, 0.0118, 0.0084, 0.0047, 0.0028, 0., 0.0086,
                0.0044, 0.0185
            ]),
            True: np.array([
                0.010018, 0.022154, 0.011829, 0.008422, 0.004724, 0.002782,
                0.000000, 0.008625, 0.004396, 0.018485
            ])
        }[ri]

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            exc_ene = rpa_result['eigenvalues']
            osc_str = rpa_result['oscillator_strengths']
            assert np.max(np.abs(exc_ene - ref_exc_ene)) < 1.0e-4
            assert np.max(np.abs(osc_str - ref_osc_str)) < 1.0e-4

    @pytest.mark.parametrize('ri', [False, True])
    def test_core_excitation_tda(self, ri):

        molecule, basis, xcfun_label, scf_result = self.run_scf(ri)

        # core excitation
        tda_drv = TdaEigenSolver()
        tda_drv.xcfun = xcfun_label
        tda_drv.core_excitation = True
        tda_drv.num_core_orbitals = 1
        tda_drv.nstates = 10
        tda_drv.ostream.mute()
        tda_result = tda_drv.compute(molecule, basis, scf_result)

        ref_exc_ene = {
            False: np.array([
                19.0941462, 19.15797992, 19.20966688, 19.21610156, 19.25756094,
                19.28713597, 19.33009442, 19.33573228, 19.36632619, 19.44544945
            ]),
            True: np.array([
                19.09415281, 19.15802601, 19.20970930, 19.21607423, 19.25757065,
                19.28712239, 19.33015081, 19.33579136, 19.36633746, 19.44548426
            ])
        }[ri]
        ref_osc_str = {
            False: np.array([
                0.0101, 0.0222, 0.0119, 0.0084, 0.0048, 0.0028, 0., 0.0087,
                0.0044, 0.0186
            ]),
            True: np.array([
                0.010061, 0.022234, 0.011863, 0.008441, 0.004742, 0.002792,
                0.000000, 0.008662, 0.004409, 0.018564
            ])
        }[ri]

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            exc_ene = tda_result['eigenvalues']
            osc_str = tda_result['oscillator_strengths']
            assert np.max(np.abs(exc_ene - ref_exc_ene)) < 1.0e-4
            assert np.max(np.abs(osc_str - ref_osc_str)) < 1.0e-4
