from pathlib import Path
import h5py
import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.lrsolver import LinearResponseSolver
from veloxchem.cppsolver import ComplexResponseSolver
from veloxchem.tdacppsolver import ComplexResponseTDA
from veloxchem.c6driver import C6Driver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.tdaeigensolver import TdaEigenSolver
from veloxchem.lrsolverunrest import LinearResponseUnrestrictedSolver
from veloxchem.inputparser import (unparse_input, read_unparsed_input_from_hdf5)


def _get_restricted_water():

    xyz_string = """3
    xyz
    O   -0.1858140  -1.1749469   0.7662596
    H   -0.1285513  -0.8984365   1.6808606
    H   -0.0582782  -0.3702550   0.2638279
    """
    molecule = Molecule.read_xyz_string(xyz_string)
    basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)

    return molecule, basis


def _get_unrestricted_water_cation():

    xyz_string = """3
    xyz
    O   -0.1858140  -1.1749469   0.7662596
    H   -0.1285513  -0.8984365   1.6808606
    H   -0.0582782  -0.3702550   0.2638279
    """
    molecule = Molecule.read_xyz_string(xyz_string)
    molecule.set_charge(1)
    molecule.set_multiplicity(2)
    basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)

    return molecule, basis


def _run_restricted_scf(molecule, basis):

    scf_drv = ScfRestrictedDriver()
    scf_drv.ostream.mute()

    return scf_drv.compute(molecule, basis)


def _run_unrestricted_scf(molecule, basis):

    scf_drv = ScfUnrestrictedDriver()
    scf_drv.ostream.mute()

    return scf_drv.compute(molecule, basis)


def _run_restarted_restricted_scf(molecule, basis, filename):

    scf_drv = ScfRestrictedDriver()
    scf_drv.ostream.mute()
    scf_drv.filename = filename
    scf_drv.compute(molecule, basis)
    scf_drv.restart = True

    return scf_drv.compute(molecule, basis)


def _configure_lr(driver, filename):

    driver.filename = filename
    driver.conv_thresh = 1.0e-5
    driver.max_iter = 120
    driver.a_components = 'xz'
    driver.b_components = 'yz'
    driver.frequencies = (0.0, 0.05)
    driver.non_equilibrium_solv = False
    driver.ri_auxiliary_basis = 'def2-universal-jfit'


def _configure_cpp(driver, filename):

    driver.filename = filename
    driver.conv_thresh = 1.0e-5
    driver.max_iter = 120
    driver.a_components = 'xz'
    driver.b_components = 'yz'
    driver.frequencies = (0.10,)
    driver.damping = 0.02
    driver.non_equilibrium_solv = False
    driver.ri_auxiliary_basis = 'def2-universal-jfit'


def _configure_rpa(driver, filename):

    driver.filename = filename
    driver.conv_thresh = 1.0e-5
    driver.max_iter = 60
    driver.nstates = 2
    driver.initial_guess_multiplier = 2
    driver.max_subspace_dim = 8
    driver.non_equilibrium_solv = False
    driver.ri_auxiliary_basis = 'def2-universal-jfit'


def _configure_tdacpp(driver, filename):

    driver.filename = filename
    driver.conv_thresh = 1.0e-5
    driver.max_iter = 120
    driver.a_components = 'xz'
    driver.b_components = 'yz'
    driver.frequencies = (0.10,)
    driver.damping = 0.02
    driver.non_equilibrium_solv = False
    driver.ri_auxiliary_basis = 'def2-universal-jfit'


def _configure_c6(driver, filename):

    driver.filename = filename
    driver.conv_thresh = 1.0e-5
    driver.max_iter = 120
    driver.n_points = 5
    driver.w0 = 0.4
    driver.non_equilibrium_solv = False
    driver.ri_auxiliary_basis = 'def2-universal-jfit'


def _configure_tda(driver, filename):

    driver.filename = filename
    driver.conv_thresh = 1.0e-5
    driver.max_iter = 60
    driver.nstates = 2
    driver.initial_guess_multiplier = 2
    driver.max_subspace_dim = 8
    driver.non_equilibrium_solv = False
    driver.ri_auxiliary_basis = 'def2-universal-jfit'


def _compare_rsp_results(reference_results, copied_results):

    assert copied_results['response_functions'].keys(
    ) == reference_results['response_functions'].keys()

    for key, value in reference_results['response_functions'].items():
        assert copied_results['response_functions'][key] == pytest.approx(
            value, abs=1.0e-8)


def _compare_eigen_results(reference_results, copied_results):

    assert np.max(
        np.abs(copied_results['eigenvalues'] -
               reference_results['eigenvalues'])) < 1.0e-8


def _compare_c6_results(reference_results, copied_results):

    assert copied_results['c6'] == pytest.approx(reference_results['c6'],
                                                 abs=1.0e-8)
    _compare_rsp_results(reference_results, copied_results)


def _normalize_setting(value):

    if isinstance(value, (list, tuple)):
        return tuple(value)

    return value


def _bcast_path_string(comm, path_string):

    if comm.Get_rank() != mpi_master():
        path_string = None

    return comm.bcast(path_string, root=mpi_master())


@pytest.mark.solvers
@pytest.mark.parametrize(
    'solver_cls, configure_solver, compare_results',
    [
        (LinearResponseSolver, _configure_lr, _compare_rsp_results),
        (ComplexResponseSolver, _configure_cpp, _compare_rsp_results),
        (ComplexResponseTDA, _configure_tdacpp, _compare_rsp_results),
        (C6Driver, _configure_c6, _compare_c6_results),
        (LinearResponseEigenSolver, _configure_rpa, _compare_eigen_results),
        (TdaEigenSolver, _configure_tda, _compare_eigen_results),
    ],
)
def test_restricted_response_settings_roundtrip(tmp_path, solver_cls,
                                                configure_solver,
                                                compare_results):

    molecule, basis = _get_restricted_water()
    scf_results = _run_restricted_scf(molecule, basis)

    driver = solver_cls()
    driver.ostream.mute()

    filename = str(tmp_path / solver_cls.__name__.lower())
    filename = driver.comm.bcast(filename, root=mpi_master())

    configure_solver(driver, filename)

    reference_scf_results = dict(scf_results)
    reference_scf_results['filename'] = filename

    reference_results = driver.compute(molecule, basis, reference_scf_results)
    checkpoint_file = Path(
        _bcast_path_string(driver.comm, driver.checkpoint_file))

    rsp_keywords = {
        key: val[0] for key, val in driver._input_keywords['response'].items()
    }
    method_keywords = {
        key: val[0]
        for key, val in driver._input_keywords['method_settings'].items()
    }

    expected_rsp = unparse_input(driver, rsp_keywords)
    expected_method = unparse_input(driver, method_keywords)

    if driver.rank == mpi_master():
        assert checkpoint_file.is_file()
        checkpoint_rsp_input = read_unparsed_input_from_hdf5(
            str(checkpoint_file), group_name='response_settings')
        checkpoint_method_input = read_unparsed_input_from_hdf5(
            str(checkpoint_file), group_name='method_settings')
    else:
        checkpoint_rsp_input = None
        checkpoint_method_input = None

    checkpoint_rsp_input = driver.comm.bcast(checkpoint_rsp_input,
                                             root=mpi_master())
    checkpoint_method_input = driver.comm.bcast(checkpoint_method_input,
                                                root=mpi_master())

    assert checkpoint_rsp_input == expected_rsp
    assert checkpoint_method_input == expected_method

    second_drv = solver_cls()
    second_drv.ostream.mute()
    second_drv.update_settings(checkpoint_rsp_input, checkpoint_method_input)

    for key in rsp_keywords:
        assert _normalize_setting(getattr(second_drv,
                                          key)) == _normalize_setting(
                                              getattr(driver, key))
    for key in method_keywords:
        assert _normalize_setting(getattr(second_drv,
                                          key)) == _normalize_setting(
                                              getattr(driver, key))

    third_drv = solver_cls()
    third_drv.ostream.mute()
    third_drv.restart = False
    copied_filename = _bcast_path_string(
        driver.comm, str(tmp_path / f'{solver_cls.__name__.lower()}_copy'))
    copied_checkpoint_file = _bcast_path_string(
        driver.comm, str(tmp_path / f'{solver_cls.__name__.lower()}_copy.h5'))
    third_drv.filename = copied_filename
    third_drv.checkpoint_file = copied_checkpoint_file

    third_drv.read_settings(str(checkpoint_file))

    assert third_drv.restart is False
    assert third_drv.filename == copied_filename
    assert third_drv.checkpoint_file == copied_checkpoint_file

    copied_scf_results = dict(scf_results)
    copied_scf_results['filename'] = third_drv.filename

    copied_results = third_drv.compute(molecule, basis, copied_scf_results)

    for key in rsp_keywords:
        if key in ('restart', 'filename', 'checkpoint_file'):
            continue
        assert _normalize_setting(getattr(third_drv,
                                          key)) == _normalize_setting(
                                              getattr(driver, key))
    for key in method_keywords:
        assert _normalize_setting(getattr(third_drv,
                                          key)) == _normalize_setting(
                                              getattr(driver, key))

    if third_drv.rank == mpi_master():
        compare_results(reference_results, copied_results)


@pytest.mark.solvers
def test_unrestricted_response_settings_roundtrip(tmp_path):

    molecule, basis = _get_unrestricted_water_cation()
    scf_results = _run_unrestricted_scf(molecule, basis)

    driver = LinearResponseUnrestrictedSolver()
    driver.ostream.mute()

    filename = str(tmp_path / 'linearresponseunrestrictedsolver')
    filename = driver.comm.bcast(filename, root=mpi_master())

    _configure_lr(driver, filename)

    reference_scf_results = dict(scf_results)
    reference_scf_results['filename'] = filename

    reference_results = driver.compute(molecule, basis, reference_scf_results)
    checkpoint_file = Path(
        _bcast_path_string(driver.comm, driver.checkpoint_file))

    rsp_keywords = {
        key: val[0] for key, val in driver._input_keywords['response'].items()
    }
    method_keywords = {
        key: val[0]
        for key, val in driver._input_keywords['method_settings'].items()
    }

    expected_rsp = unparse_input(driver, rsp_keywords)
    expected_method = unparse_input(driver, method_keywords)

    if driver.rank == mpi_master():
        assert checkpoint_file.is_file()
        checkpoint_rsp_input = read_unparsed_input_from_hdf5(
            str(checkpoint_file), group_name='response_settings')
        checkpoint_method_input = read_unparsed_input_from_hdf5(
            str(checkpoint_file), group_name='method_settings')
    else:
        checkpoint_rsp_input = None
        checkpoint_method_input = None

    checkpoint_rsp_input = driver.comm.bcast(checkpoint_rsp_input,
                                             root=mpi_master())
    checkpoint_method_input = driver.comm.bcast(checkpoint_method_input,
                                                root=mpi_master())

    assert checkpoint_rsp_input == expected_rsp
    assert checkpoint_method_input == expected_method

    second_drv = LinearResponseUnrestrictedSolver()
    second_drv.ostream.mute()
    second_drv.update_settings(checkpoint_rsp_input, checkpoint_method_input)

    for key in rsp_keywords:
        assert _normalize_setting(getattr(second_drv,
                                          key)) == _normalize_setting(
                                              getattr(driver, key))
    for key in method_keywords:
        assert _normalize_setting(getattr(second_drv,
                                          key)) == _normalize_setting(
                                              getattr(driver, key))

    third_drv = LinearResponseUnrestrictedSolver()
    third_drv.ostream.mute()
    third_drv.restart = False
    copied_filename = _bcast_path_string(
        driver.comm, str(tmp_path / 'linearresponseunrestrictedsolver_copy'))
    copied_checkpoint_file = _bcast_path_string(
        driver.comm, str(tmp_path / 'linearresponseunrestrictedsolver_copy.h5'))
    third_drv.filename = copied_filename
    third_drv.checkpoint_file = copied_checkpoint_file

    third_drv.read_settings(str(checkpoint_file))

    assert third_drv.restart is False
    assert third_drv.filename == copied_filename
    assert third_drv.checkpoint_file == copied_checkpoint_file

    copied_scf_results = dict(scf_results)
    copied_scf_results['filename'] = third_drv.filename

    copied_results = third_drv.compute(molecule, basis, copied_scf_results)

    for key in rsp_keywords:
        if key in ('restart', 'filename', 'checkpoint_file'):
            continue
        assert _normalize_setting(getattr(third_drv,
                                          key)) == _normalize_setting(
                                              getattr(driver, key))
    for key in method_keywords:
        assert _normalize_setting(getattr(third_drv,
                                          key)) == _normalize_setting(
                                              getattr(driver, key))

    if third_drv.rank == mpi_master():
        _compare_rsp_results(reference_results, copied_results)


@pytest.mark.solvers
def test_match_settings_rejects_restart_with_mismatched_settings(tmp_path):

    molecule, basis = _get_restricted_water()

    reference_drv = LinearResponseSolver()
    reference_drv.ostream.mute()

    filename = _bcast_path_string(reference_drv.comm,
                                  str(tmp_path / 'lr_match_settings'))
    scf_results = _run_restarted_restricted_scf(molecule, basis, filename)
    reference_scf_results = dict(scf_results)
    reference_scf_results['filename'] = filename

    _configure_lr(reference_drv, filename)
    reference_drv.compute(molecule, basis, reference_scf_results)

    checkpoint_file = _bcast_path_string(reference_drv.comm,
                                         reference_drv.checkpoint_file)

    mismatch_drv = LinearResponseSolver()
    mismatch_drv.ostream.mute()
    _configure_lr(mismatch_drv, filename)
    mismatch_drv.checkpoint_file = checkpoint_file
    mismatch_drv.restart = True
    mismatch_drv.a_components = 'xy'

    assert mismatch_drv.match_settings(checkpoint_file) is False

    mismatch_scf_results = dict(scf_results)
    mismatch_scf_results['filename'] = filename
    mismatch_drv.compute(molecule, basis, mismatch_scf_results)

    assert mismatch_drv.restart is False


@pytest.mark.solvers
def test_match_settings_allows_restart_without_settings_groups(tmp_path):

    molecule, basis = _get_restricted_water()

    reference_drv = LinearResponseSolver()
    reference_drv.ostream.mute()

    filename = _bcast_path_string(reference_drv.comm,
                                  str(tmp_path / 'lr_backward_compatible'))
    scf_results = _run_restarted_restricted_scf(molecule, basis, filename)
    reference_scf_results = dict(scf_results)
    reference_scf_results['filename'] = filename

    _configure_lr(reference_drv, filename)
    reference_drv.compute(molecule, basis, reference_scf_results)

    checkpoint_file = _bcast_path_string(reference_drv.comm,
                                         reference_drv.checkpoint_file)

    if reference_drv.rank == mpi_master():
        with h5py.File(checkpoint_file, 'a') as h5f:
            del h5f['response_settings']
            del h5f['method_settings']

    reference_drv.comm.barrier()

    restarted_drv = LinearResponseSolver()
    restarted_drv.ostream.mute()
    _configure_lr(restarted_drv, filename)
    restarted_drv.checkpoint_file = checkpoint_file
    restarted_drv.restart = True

    assert restarted_drv.match_settings(checkpoint_file) is True

    restarted_scf_results = dict(scf_results)
    restarted_scf_results['filename'] = filename
    restarted_drv.compute(molecule, basis, restarted_scf_results)

    assert restarted_drv.restart is True
