import numpy as np

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.tdaeigensolver import TdaEigenSolver


class TestCoreExcitation:

    def run_scf(self):

        water_xyz = """
        O   0.0   0.0   0.0
        H   0.0   1.4   1.1
        H   0.0  -1.4   1.1
        """
        molecule = Molecule.read_str(water_xyz, units='au')
        basis = MolecularBasis.read(molecule, 'def2-svpd', ostream=None)
        xcfun_label = 'b3lyp'

        # ground state
        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = xcfun_label
        scf_drv.ostream.mute()
        scf_result = scf_drv.compute(molecule, basis)

        return molecule, basis, xcfun_label, scf_result

    def test_core_excitation_rpa(self):

        molecule, basis, xcfun_label, scf_result = self.run_scf()

        # core excitation
        rpa_drv = LinearResponseEigenSolver()
        rpa_drv.xcfun = xcfun_label
        rpa_drv.core_excitation = True
        rpa_drv.num_core_orbitals = 1
        rpa_drv.nstates = 10
        rpa_drv.ostream.mute()
        rpa_result = rpa_drv.compute(molecule, basis, scf_result)

        ref_exc_ene = np.array([
            19.09409611, 19.15795969, 19.20965780, 19.21608706, 19.25755683,
            19.28712810, 19.33009437, 19.33571861, 19.36632321, 19.44738600
        ])
        ref_osc_str = np.array([
            0.0100, 0.0221, 0.0118, 0.0084, 0.0047, 0.0028, 0.0000, 0.0086,
            0.0044, 0.0017
        ])

        if is_mpi_master():
            exc_ene = rpa_result['eigenvalues']
            osc_str = rpa_result['oscillator_strengths']
            assert np.max(np.abs(exc_ene - ref_exc_ene)) < 1.0e-4
            assert np.max(np.abs(osc_str - ref_osc_str)) < 1.0e-4

    def test_core_excitation_tda(self):

        molecule, basis, xcfun_label, scf_result = self.run_scf()

        # core excitation
        tda_drv = TdaEigenSolver()
        tda_drv.xcfun = xcfun_label
        tda_drv.core_excitation = True
        tda_drv.num_core_orbitals = 1
        tda_drv.nstates = 10
        tda_drv.ostream.mute()
        tda_result = tda_drv.compute(molecule, basis, scf_result)

        ref_exc_ene = np.array([
            19.09414620, 19.15797993, 19.20966688, 19.21610156, 19.25756094,
            19.28713597, 19.33009443, 19.33573228, 19.36632620, 19.44744147
        ])
        ref_osc_str = np.array([
            0.0101, 0.0222, 0.0119, 0.0084, 0.0048, 0.0028, 0.0000, 0.0087,
            0.0044, 0.0017
        ])

        if is_mpi_master():
            exc_ene = tda_result['eigenvalues']
            osc_str = tda_result['oscillator_strengths']
            assert np.max(np.abs(exc_ene - ref_exc_ene)) < 1.0e-4
            assert np.max(np.abs(osc_str - ref_osc_str)) < 1.0e-4
