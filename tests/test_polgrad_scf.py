from pathlib import Path
import numpy as np
import pytest
import h5py
import sys

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lrsolver import LinearResponseSolver
from veloxchem.cppsolver import ComplexResponse
from veloxchem.polarizabilitygradient import PolarizabilityGradient


class TestPolgrad:

    def run_polgrad_real(self, molecule, basis, xcfun=None, label=None):
        scf_drv = ScfRestrictedDriver()
        scf_dict = {}
        method_settings = {}
        #if xcfun is not None:
        #    scf_drv._dft = True
        #    scf_drv.xcfun = xcfun
        #    method_settings = {'xcfun': xcfun}

        scf_drv.update_settings(scf_dict, method_settings)
        scf_drv.ostream.mute()
        scf_tensors = scf_drv.compute(molecule, basis)

        # linear response
        rsp_settings = {'conv_thresh': 1.0e-5, 'frequencies': (0.0, 0.4)}
        lr_drv = LinearResponseSolver()
        lr_drv.a_operator = "electric dipole"
        lr_drv.b_operator = "electric dipole"
        lr_drv.update_settings(rsp_settings, method_settings)
        lr_drv.ostream.mute()
        lr_results = lr_drv.compute(molecule, basis, scf_tensors)

        # test real analytical gradient
        an_polgrad_drv = PolarizabilityGradient(scf_drv)
        cphf_settings = {'conv_thresh':2e-7, 'use_subspace_solver': 'no'}
        polgrad_settings = {'frequencies': (0.0, 0.4)}
        an_polgrad_drv.update_settings(polgrad_settings, cphf_settings, method_settings)
        an_polgrad_drv.ostream.mute()
        an_polgrad_drv.compute(molecule, basis, scf_tensors, lr_results)

        if scf_drv.rank == mpi_master():
            polgrad_results = an_polgrad_drv.polgradient
            polgrad_static = polgrad_results[0.0].reshape(3,3,3,3)
            polgrad_dynamic = polgrad_results[0.4].reshape(3,3,3,3)
            np.set_printoptions(suppress=True, precision=10)
            here = Path(__file__).parent
            hf_file_name = str(here /'data'/'polarizabilitygradients.h5')
            hf = h5py.File(hf_file_name, 'r')
            an_label = 'an_' + label
            polgrad_reference = np.array(hf.get(an_label))
            polgrad_static_reference = polgrad_reference[0]
            polgrad_dynamic_reference = polgrad_reference[1]
            hf.close()

            assert np.max(np.abs(polgrad_static) - np.abs(polgrad_static_reference)) < 1.0e-6
            assert np.max(np.abs(polgrad_dynamic) - np.abs(polgrad_dynamic_reference)) < 1.0e-6

        # test real numerical gradient
        num_polgrad_drv = PolarizabilityGradient(scf_drv)
        polgrad_settings = {'numerical': 'yes', 'do_four_point': 'yes', 'frequencies': (0.0, 0.4)}
        cphf_settings = {}
        num_polgrad_drv.update_settings(polgrad_settings, cphf_settings, method_settings)
        num_polgrad_drv.ostream.mute()
        num_polgrad_drv.compute(molecule, basis, scf_tensors, lr_results=None)

        if scf_drv.rank == mpi_master():
            polgrad_results = num_polgrad_drv.polgradient
            polgrad_static = polgrad_results[0.0].reshape(3,3,3,3)
            polgrad_dynamic = polgrad_results[0.4].reshape(3,3,3,3)
            np.set_printoptions(suppress=True, precision=10)
            here = Path(__file__).parent
            hf_file_name = str(here /'data'/'polarizabilitygradients.h5')
            hf = h5py.File(hf_file_name, 'r')
            num_label = 'num_' + label
            polgrad_reference = np.array(hf.get(num_label))
            polgrad_static_reference = polgrad_reference[0]
            polgrad_dynamic_reference = polgrad_reference[1]
            hf.close()

            assert np.max(np.abs(polgrad_dynamic) - np.abs(polgrad_dynamic_reference)) < 1.0e-6

    def run_polgrad_complex(self, molecule, basis, xcfun=None, label=None):
        scf_drv = ScfRestrictedDriver()
        scf_dict = {}
        method_settings = {}
        #if xcfun is not None:
        #    scf_drv._dft = True
        #    scf_drv.xcfun = xcfun
        #    method_settings = {'xcfun': xcfun}

        scf_drv.update_settings(scf_dict, method_settings)
        scf_drv.ostream.mute()
        scf_tensors = scf_drv.compute(molecule, basis)

        # linear response
        rsp_settings = {'conv_thresh': 1.0e-5, 'frequencies': (0.0, 0.4),
                        'damping': 0.5}
        lr_drv = ComplexResponse()
        lr_drv.a_operator = "electric dipole"
        lr_drv.b_operator = "electric dipole"
        lr_drv.update_settings(rsp_settings, method_settings)
        lr_drv.ostream.mute()
        lr_results = lr_drv.compute(molecule, basis, scf_tensors)

        # test complex analytical gradient
        an_polgrad_drv = PolarizabilityGradient(scf_drv)
        cphf_settings = {'conv_thresh':2e-7, 'use_subspace_solver': 'no'}
        polgrad_settings = {'frequencies': (0.0, 0.4), 'is_complex': 'yes',
                            'damping': 0.5}
        an_polgrad_drv.update_settings(polgrad_settings, cphf_settings, method_settings)
        an_polgrad_drv.ostream.mute()
        an_polgrad_drv.compute(molecule, basis, scf_tensors, lr_results)

        if scf_drv.rank == mpi_master():
            polgrad_results = an_polgrad_drv.polgradient
            polgrad_static = polgrad_results[0.0].reshape(3,3,3,3)
            polgrad_dynamic = polgrad_results[0.4].reshape(3,3,3,3)
            np.set_printoptions(suppress=True, precision=10)
            here = Path(__file__).parent
            hf_file_name = str(here /'data'/'polarizabilitygradients.h5')
            hf = h5py.File(hf_file_name, 'r')
            an_label = 'an_' + label
            polgrad_reference = np.array(hf.get(an_label))
            polgrad_static_reference = polgrad_reference[0]
            polgrad_dynamic_reference = polgrad_reference[1]
            hf.close()

            assert np.max(np.abs(polgrad_static) - np.abs(polgrad_static_reference)) < 1.0e-6
            assert np.max(np.abs(polgrad_dynamic) - np.abs(polgrad_dynamic_reference)) < 1.0e-6

        # test complex numerical gradient
        num_polgrad_drv = PolarizabilityGradient(scf_drv)
        polgrad_settings = {'frequencies': (0.0, 0.4), 'is_complex': 'yes',
                            'damping': 0.5, 'numerical': 'yes', 'do_four_point': 'yes'}
        num_polgrad_drv.update_settings(polgrad_settings, cphf_settings, method_settings)
        num_polgrad_drv.ostream.mute()
        num_polgrad_drv.compute(molecule, basis, scf_tensors)

        if scf_drv.rank == mpi_master():
            polgrad_results = num_polgrad_drv.polgradient
            polgrad_static = polgrad_results[0.0].reshape(3,3,3,3)
            polgrad_dynamic = polgrad_results[0.4].reshape(3,3,3,3)
            np.set_printoptions(suppress=True, precision=10)
            here = Path(__file__).parent
            hf_file_name = str(here /'data'/'polarizabilitygradients.h5')
            hf = h5py.File(hf_file_name, 'r')
            num_label = 'num_' + label
            polgrad_reference = np.array(hf.get(num_label))
            polgrad_static_reference = polgrad_reference[0]
            polgrad_dynamic_reference = polgrad_reference[1]
            hf.close()

            assert np.max(np.abs(polgrad_static) - np.abs(polgrad_static_reference)) < 1.0e-6
            assert np.max(np.abs(polgrad_dynamic) - np.abs(polgrad_dynamic_reference)) < 1.0e-6

    def test_hf_polarizabilitygradient_real(self):
        h2o_xyz = """3

        O     0.000000    0.000000    0.000000
        H     0.000000    0.504284    0.758602
        H     0.000000   -0.504284    0.758602 
        """
        basis_set_label = "sto-3g"

        molecule = Molecule.from_xyz_string(h2o_xyz)
        basis = MolecularBasis.read(molecule, basis_set_label)

        self.run_polgrad_real(molecule, basis, None, "polarizabilitygradient_hf_real")

    def test_hf_polarizabilitygradient_complex(self):
        h2o_xyz = """3

        O     0.000000    0.000000    0.000000
        H     0.000000    0.504284    0.758602
        H     0.000000   -0.504284    0.758602 
        """
        basis_set_label = "sto-3g"

        molecule = Molecule.from_xyz_string(h2o_xyz)
        basis = MolecularBasis.read(molecule, basis_set_label)

        self.run_polgrad_complex(molecule, basis, None, "polarizabilitygradient_hf_complex")

