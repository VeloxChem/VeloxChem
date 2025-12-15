import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.cppsolver import ComplexResponse


@pytest.mark.solvers
class TestCppPropertyDensities:

    def run_cpp_prop_densities(self, mol, bas, xcfun_label, cpp_flag,
                               cpp_frequencies, freq_prop_tuple, tol):

        # run SCF

        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = xcfun_label
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        # run CPP

        cpp_drv = ComplexResponse()
        cpp_drv.frequencies = cpp_frequencies
        cpp_drv.cpp_flag = cpp_flag
        cpp_drv.ostream.mute()
        cpp_results = cpp_drv.compute(mol, bas, scf_results)

        # get property densities

        freq, ref_value = freq_prop_tuple

        raw_densitiy_dict = cpp_drv.get_cpp_property_densities(
            mol, bas, scf_results, cpp_results, freq, normalize_densities=False)

        densitiy_dict = cpp_drv.get_cpp_property_densities(
            mol, bas, scf_results, cpp_results, freq, normalize_densities=True)

        # check property densities

        if scf_drv.rank == mpi_master():
            S = scf_results['S']

            prop_dens_D = raw_densitiy_dict['property_density_detachment']
            prop_dens_A = raw_densitiy_dict['property_density_attachment']

            assert abs(-np.sum(prop_dens_D * S) - ref_value) < tol
            assert abs(np.sum(prop_dens_A * S) - ref_value) < tol

            prop_dens_D = densitiy_dict['property_density_detachment']
            prop_dens_A = densitiy_dict['property_density_attachment']

            if ref_value > 0:
                assert abs(-np.sum(prop_dens_D * S) - 1.0) < tol
                assert abs(np.sum(prop_dens_A * S) - 1.0) < tol
            else:
                assert abs(-np.sum(prop_dens_D * S) + 1.0) < tol
                assert abs(np.sum(prop_dens_A * S) + 1.0) < tol

    def test_cpp_prop_dens_absorption(self):

        molstr = """
        C        0.0000        0.0000       -1.2436
        C       -0.0000        0.0000        1.2436
        H        0.0000        1.7272       -2.3230
        H        0.0000       -1.7272       -2.3230
        H       -0.0000        1.7272        2.3230
        H       -0.0000       -1.7272        2.3230
        """
        mol = Molecule.read_molecule_string(molstr, units='bohr')

        bas = MolecularBasis.read(mol, 'sto-3g', verbose=False)

        xcfun_label = 'bhandhlyp'

        cpp_flag = 'absorption'
        cpp_frequencies = np.round(np.arange(0.7, 0.8, 0.005), 4)

        freq_prop_tuple = (0.735, 7.46130779)
        tol = 1e-7

        self.run_cpp_prop_densities(mol, bas, xcfun_label, cpp_flag,
                                    cpp_frequencies, freq_prop_tuple, tol)

    def test_cpp_prop_dens_ecd(self):

        xyzstr = """4
        xyz
        O           -0.65564532 -0.06106286  -0.03621403
        O            0.65564532  0.06106286  -0.03621403
        H           -0.97628735  0.65082652   0.57474201
        H            0.97628735 -0.65082652   0.57474201
        """
        mol = Molecule.read_xyz_string(xyzstr)

        bas = MolecularBasis.read(mol, 'sto-3g', verbose=False)

        xcfun_label = 'bhandhlyp'

        cpp_flag = 'ecd'
        cpp_frequencies = np.round(np.arange(0.3, 0.4, 0.005), 4)

        freq_prop_tuple = (0.380, -33.00280072)
        tol = 1e-7

        self.run_cpp_prop_densities(mol, bas, xcfun_label, cpp_flag,
                                    cpp_frequencies, freq_prop_tuple, tol)
