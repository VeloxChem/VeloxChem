import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.cppsolver import ComplexResponse


@pytest.mark.solvers
class TestCPP:

    def run_cpp(self,
                xcfun_label,
                cpp_flag,
                ref_x_data,
                ref_y_data,
                tol,
                use_subcomms=False):

        xyz_string = """6
        xyz
        O -3.42904  1.55532  0.01546
        C -1.99249  1.74379  0.02665
        H -1.74709  2.74160  0.44749
        H -1.59636  1.67836 -1.00861
        H -1.51398  0.95881  0.64937
        H -3.84726  2.33620 -0.34927
        """
        mol = Molecule.read_xyz_string(xyz_string)
        mol.check_multiplicity()

        basis_label = 'def2-svp'
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = ComplexResponse()
        lr_drv.ostream.mute()
        lr_drv.set_cpp_flag(cpp_flag)
        lr_drv.frequencies = list(ref_x_data)
        lr_drv.use_subcomms = use_subcomms
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            lr_spec = lr_drv.get_spectrum(lr_results, 'au')
            assert np.max(
                np.abs(np.array(lr_spec['y_data']) -
                       np.array(ref_y_data))) < tol

    def test_hf_absorption(self):

        # vlxtag: RHF, Absorption, CPP

        xcfun_label = 'hf'
        cpp_flag = 'absorption'
        ref_x_data = [0.39, 0.40, 0.41]
        ref_y_data = [0.11359517, 0.34554765, 1.41940099]

        self.run_cpp(xcfun_label, cpp_flag, ref_x_data, ref_y_data, 1.0e-6)

    def test_hf_ecd(self):

        # vlxtag: RHF, ECD, CPP

        xcfun_label = 'hf'
        cpp_flag = 'ecd'
        ref_x_data = [0.39, 0.40, 0.41]
        ref_y_data = [5.74701958, 35.38200618, -15.89867350]

        self.run_cpp(xcfun_label, cpp_flag, ref_x_data, ref_y_data, 1.0e-6)

        self.run_cpp(xcfun_label,
                     cpp_flag,
                     ref_x_data,
                     ref_y_data,
                     1.0e-6,
                     use_subcomms=True)

    def test_b3lyp_absorption(self):

        # vlxtag: RKS, Absorption, CPP

        xcfun_label = 'b3lyp'
        cpp_flag = 'absorption'
        ref_x_data = [0.39, 0.40, 0.41]
        ref_y_data = [0.20856963, 0.21772352, 0.71434376]

        self.run_cpp(xcfun_label, cpp_flag, ref_x_data, ref_y_data, 1.0e-4)

    def test_b3lyp_ecd(self):

        # vlxtag: RKS, ECD, CPP

        xcfun_label = 'b3lyp'
        cpp_flag = 'ecd'
        ref_x_data = [0.39, 0.40, 0.41]
        ref_y_data = [0.48472627, -1.18394788, -9.81377126]

        self.run_cpp(xcfun_label, cpp_flag, ref_x_data, ref_y_data, 1.0e-4)

        self.run_cpp(xcfun_label,
                     cpp_flag,
                     ref_x_data,
                     ref_y_data,
                     1.0e-4,
                     use_subcomms=True)
