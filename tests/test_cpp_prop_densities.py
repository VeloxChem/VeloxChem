import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.cppsolver import ComplexResponse


@pytest.mark.solvers
class TestCppPropertyDensities:

    def run_cpp_prop_densities(self, mol, bas, xcfun_label, freq, ref_value,
                               tol):

        # run SCF

        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = xcfun_label
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        # run CPP

        cpp_drv = ComplexResponse()
        cpp_drv.frequencies = np.round(np.arange(0.7, 0.8, 0.005), 3)
        cpp_drv.cpp_flag = 'absorption'
        cpp_drv.ostream.mute()
        cpp_results = cpp_drv.compute(mol, bas, scf_results)

        # get property densities

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

            assert abs(-np.sum(prop_dens_D * S) - 1.0) < tol
            assert abs(np.sum(prop_dens_A * S) - 1.0) < tol

    def test_cpp_prop_densities(self):

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

        freq = 0.735
        ref_value = 7.46130779
        tol = 1e-7

        self.run_cpp_prop_densities(mol, bas, xcfun_label, freq, ref_value, tol)
