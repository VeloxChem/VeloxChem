from mpi4py import MPI
from types import SimpleNamespace
import numpy as np
import pytest

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.cppsolver import ComplexResponseSolver
from veloxchem.errorhandler import VeloxChemError
from veloxchem.main import select_rsp_property
from veloxchem.outputstream import OutputStream
from veloxchem.rspordspec import OpticalRotatoryDispersionSpectrum
from veloxchem.veloxchemlib import mpi_master


class _MoleculeStub:

    @staticmethod
    def number_of_alpha_occupied_orbitals(basis):
        return 1

    @staticmethod
    def number_of_beta_occupied_orbitals(basis):
        return 1


class _MolecularOrbitalsStub:

    @staticmethod
    def number_of_mos():
        return 2


@pytest.mark.solvers
class TestORD:

    @staticmethod
    def _task_stub():

        return SimpleNamespace(mpi_rank=mpi_master(),
                               mpi_size=1,
                               molecule=_MoleculeStub(),
                               ao_basis=object())

    @pytest.mark.parametrize(
        'property_name',
        [
            'optical rotatory dispersion (cpp)',
            'optical rotatory dispersion(cpp)',
            'optical rotation (cpp)',
            'optical rotation(cpp)',
            'ord (cpp)',
            'ord(cpp)',
        ],
    )
    def test_select_rsp_property_dispatches_ord(self, property_name):

        rsp_dict = {'property': property_name, 'frequencies': '0.1'}
        rsp_prop = select_rsp_property(self._task_stub(),
                                       _MolecularOrbitalsStub(), rsp_dict, {})

        assert isinstance(rsp_prop, OpticalRotatoryDispersionSpectrum)
        assert rsp_prop._rsp_dict['property'] == 'ord (cpp)'
        assert rsp_prop._rsp_dict['a_operator'] == 'electric dipole'
        assert rsp_prop._rsp_dict['b_operator'] == 'magnetic dipole'

    def test_ord_initializes_restricted_cpp_solver(self):

        rsp_prop = OpticalRotatoryDispersionSpectrum({'frequencies': '0.1'})
        rsp_prop.init_driver(ostream=OutputStream(None),
                             method_type='restricted')

        assert isinstance(rsp_prop.rsp_driver, ComplexResponseSolver)
        assert rsp_prop.rsp_driver.property == 'ord'
        assert rsp_prop.rsp_driver.a_operator == 'electric dipole'
        assert rsp_prop.rsp_driver.b_operator == 'magnetic dipole'

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    @pytest.mark.parametrize(
        'rsp_dict, method_type, error_message',
        [
            (
                {'frequencies': '0.1'},
                'unrestricted',
                'only implemented for restricted reference'
            ),
            (
                {'frequencies': '0.1', 'tamm_dancoff': 'yes'},
                'restricted',
                'not implemented with the Tamm-Dancoff approximation'
            ),
        ],
    )
    def test_ord_rejects_unsupported_solver_modes(self, rsp_dict, method_type,
                                                  error_message):

        rsp_prop = OpticalRotatoryDispersionSpectrum(rsp_dict)

        with pytest.raises(VeloxChemError, match=error_message):
            rsp_prop.init_driver(ostream=OutputStream(None),
                                 method_type=method_type)

    def test_ord(self):

        xyz_string = """10
        xyz
        O             -1.009866880244         1.299407071912         1.951947754409
        C              0.031038818173         2.001294498224         1.244371321259
        C              0.339506903623         0.795807201271         2.029917688849
        C              0.600316751832        -0.544462572436         1.394186918859
        H             -0.026285836651         1.930258644799         0.154459855148
        H              0.271804379347         2.990968281608         1.641505800108
        H              0.794935917667         0.948100885444         3.014899839267
        H              0.113110324011        -0.610670447580         0.412641743927
        H              1.681576973773        -0.701096452623         1.264343536564
        H              0.213803265264        -1.354101659046         2.029564064183
        """
        mol = Molecule.read_xyz_string(xyz_string)

        basis_label = '6-31g'
        xcfun_label = 'b3lyp'
        cpp_property = 'ord'

        ref_x_data = [0.3]
        ref_y_data = [12.18]

        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = ComplexResponseSolver()
        lr_drv.ostream.mute()
        lr_drv.property = cpp_property
        lr_drv.frequencies = list(ref_x_data)
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            lr_spec = lr_drv.get_spectrum(lr_results, 'au')
            assert np.max(
                np.abs(np.array(lr_spec['y_data']) -
                       np.array(ref_y_data))) < 0.01
