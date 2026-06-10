from types import SimpleNamespace

import pytest

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


def _task_stub():

    return SimpleNamespace(mpi_rank=mpi_master(),
                           mpi_size=1,
                           molecule=_MoleculeStub(),
                           ao_basis=object())


@pytest.mark.parametrize(
    'property_name',
    [
        'optical rotatory dispersion (cpp)',
        'optical rotation(cpp)',
        'ord(cpp)',
    ],
)
def test_select_rsp_property_dispatches_ord(property_name):

    rsp_dict = {'property': property_name, 'frequencies': '0.1'}
    rsp_prop = select_rsp_property(_task_stub(), _MolecularOrbitalsStub(),
                                   rsp_dict, {})

    assert isinstance(rsp_prop, OpticalRotatoryDispersionSpectrum)
    assert rsp_prop._rsp_dict['property'] == 'ord (cpp)'
    assert rsp_prop._rsp_dict['a_operator'] == 'electric dipole'
    assert rsp_prop._rsp_dict['b_operator'] == 'magnetic dipole'


def test_ord_initializes_restricted_cpp_solver():

    rsp_prop = OpticalRotatoryDispersionSpectrum({'frequencies': '0.1'})
    rsp_prop.init_driver(ostream=OutputStream(None), method_type='restricted')

    assert isinstance(rsp_prop.rsp_driver, ComplexResponseSolver)
    assert rsp_prop.rsp_driver.property == 'ord'
    assert rsp_prop.rsp_driver.a_operator == 'electric dipole'
    assert rsp_prop.rsp_driver.b_operator == 'magnetic dipole'


@pytest.mark.parametrize(
    'rsp_dict, method_type, error_message',
    [
        ({
            'frequencies': '0.1'
        }, 'unrestricted', 'only implemented for restricted reference'),
        ({
            'frequencies': '0.1',
            'tamm_dancoff': 'yes'
        }, 'restricted', 'not implemented with the Tamm-Dancoff approximation'),
    ],
)
def test_ord_rejects_unsupported_solver_modes(rsp_dict, method_type,
                                              error_message):

    rsp_prop = OpticalRotatoryDispersionSpectrum(rsp_dict)

    with pytest.raises(VeloxChemError, match=error_message):
        rsp_prop.init_driver(ostream=OutputStream(None),
                             method_type=method_type)
