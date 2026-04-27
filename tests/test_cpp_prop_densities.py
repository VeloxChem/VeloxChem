from mpi4py import MPI
import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.cppsolver import ComplexResponseSolver
from veloxchem.resultsio import read_results


@pytest.mark.solvers
class TestCppPropertyDensities:

    def run_cpp_prop_densities(self, tmp_path, mol, bas, xcfun_label, cpp_property,
                               cpp_frequencies, freq_prop_tuple, tol):

        comm = MPI.COMM_WORLD
        if comm.Get_rank() == mpi_master():
            filename = str(tmp_path / 'cpp_prop_densities')
        else:
            filename = None
        filename = comm.bcast(filename, root=mpi_master())

        # run SCF

        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = xcfun_label
        scf_drv.acc_type = 'l2_c2diis'
        scf_drv.filename = filename
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        # run CPP

        cpp_drv = ComplexResponseSolver()
        cpp_drv.frequencies = cpp_frequencies
        cpp_drv.property = cpp_property
        cpp_drv.ostream.mute()
        cpp_results_orig = cpp_drv.compute(mol, bas, scf_results)

        if scf_drv.rank == mpi_master():
            cpp_results_read = read_results(filename + '.h5', 'rsp')
        else:
            cpp_results_read = {}

        for cpp_results in [cpp_results_orig, cpp_results_read]:

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

    def test_cpp_prop_dens_absorption(self, tmp_path):

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

        cpp_property = 'absorption'
        cpp_frequencies = np.round(np.arange(0.7, 0.8, 0.005), 4)

        freq_prop_tuple = (0.735, 7.46130779)
        tol = 1e-7

        self.run_cpp_prop_densities(tmp_path, mol, bas, xcfun_label, cpp_property,
                                    cpp_frequencies, freq_prop_tuple, tol)

    def test_cpp_prop_dens_ecd(self, tmp_path):

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

        cpp_property = 'ecd'
        cpp_frequencies = np.round(np.arange(0.3, 0.4, 0.005), 4)

        freq_prop_tuple = (0.380, -33.00280072)
        tol = 1e-7

        self.run_cpp_prop_densities(tmp_path, mol, bas, xcfun_label, cpp_property,
                                    cpp_frequencies, freq_prop_tuple, tol)
