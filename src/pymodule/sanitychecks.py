#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from pathlib import Path
import numpy as np

from .veloxchemlib import parse_xc_func
from .veloxchemlib import mpi_master
from .dftutils import get_default_grid_level
from .errorhandler import assert_msg_critical


def molecule_sanity_check(mol, method_type=None, caller_name=None):
    """
    Checks molecule for charge/multiplicity combination and geometry.

    :param mol:
        The molecule.
    :param method_type:
        The method type (restricted, unrestricted or restricted_openshell).
    :param caller_name:
        The caller name used for method-specific error messages.
    """

    assert_msg_critical(
        mol.check_multiplicity(),
        'Molecule: Incompatible multiplicity and number of electrons')

    if method_type == 'restricted':
        if caller_name is None:
            err_msg = f'Molecule: Invalid multiplicity for {method_type}'
        else:
            err_msg = f'{caller_name}: not implemented for unrestricted case'
        assert_msg_critical(
            mol.get_multiplicity() == 1,
            err_msg)

    assert_msg_critical(mol.check_proximity(0.1), 'Molecule: Atoms too close')


def _warn_if_environment_overwritten(obj, key, scf_value):
    """
    Warns if an environment setting explicitly specified in response
    settings differs from the value used in SCF and will be overwritten
    by it.

    :param obj:
        The object (response driver) that is being updated.
    :param key:
        The environment setting name.
    :param scf_value:
        The value from the SCF results.
    """

    rsp_value = getattr(obj, key, None)

    if rsp_value is None:
        return

    if isinstance(rsp_value, np.ndarray) or isinstance(scf_value, np.ndarray):
        differs = not np.array_equal(rsp_value, scf_value)
    elif (isinstance(rsp_value, (list, tuple)) and
            isinstance(scf_value, (list, tuple))):
        differs = (list(rsp_value) != list(scf_value))
    else:
        differs = (rsp_value != scf_value)

    if differs:
        warn_msg = f"The '{key}' in response settings differs from the "
        warn_msg += 'value used in SCF and will be overwritten: '
        warn_msg += str(scf_value)
        obj.ostream.print_warning(warn_msg)
        obj.ostream.flush()


def _warn_if_environment_absent_from_scf(obj, scf_results):
    """
    Warns about environments active in a calculation but absent from its SCF
    reference.

    :param obj:
        The driver being updated.
    :param scf_results:
        The dictionary containing SCF results.
    """

    environments_absent_from_scf = []

    for key in ('potfile', 'solvation_model'):
        if (getattr(obj, key, None) is not None and
                scf_results.get(key, None) is None):
            environments_absent_from_scf.append(key)

    if (getattr(obj, 'pressure', 0.0) not in (None, 0.0) and
            scf_results.get('pressure', 0.0) in (None, 0.0)):
        environments_absent_from_scf.append('pressure')

    if environments_absent_from_scf:
        warn_msg = ('Environment settings active in the current calculation '
                    'but absent from the SCF reference: ')
        warn_msg += ', '.join(environments_absent_from_scf)
        warn_msg += ('. The SCF reference orbitals and density are not '
                     'relaxed for these settings.')
        obj.ostream.print_warning(warn_msg)
        obj.ostream.flush()


def scf_results_sanity_check(obj, scf_results):
    """
    Checks SCF results and inherits method and environment settings.

    :param obj:
        The object (response driver) that is being updated.
    :param scf_results:
        A dictionary containing SCF results.
    """

    updated_scf_info = {}

    if obj.rank == mpi_master():
        if 'filename' in scf_results:
            # force filename to be the same since we want a single final h5
            # file to store all results
            updated_scf_info['filename'] = scf_results['filename']

        if scf_results.get('eri_thresh', None) is not None:
            updated_scf_info['eri_thresh'] = scf_results['eri_thresh']

        updated_scf_info['ri_coulomb'] = scf_results.get('ri_coulomb', False)
        updated_scf_info['ri_jk'] = scf_results.get('ri_jk', False)
        if 'ri_auxiliary_basis' in scf_results:
            updated_scf_info['ri_auxiliary_basis'] = scf_results[
                'ri_auxiliary_basis']

        if scf_results.get('restart', None) is not None:
            # do not restart if scf is not restarted from checkpoint
            if not scf_results['restart']:
                updated_scf_info['restart'] = scf_results['restart']

        if scf_results.get('xcfun', None) is not None:
            # do not overwrite xcfun if it is already specified: a different
            # functional in response will not be rejected
            if obj.xcfun is None:
                updated_scf_info['xcfun'] = scf_results['xcfun']
                if 'grid_level' in scf_results:
                    updated_scf_info['grid_level'] = scf_results['grid_level']

        if scf_results.get('xc_screening_threshold', None) is not None:
            # the threshold used for the SCF quadrature is also used for the
            # response quadrature, so that the two screen consistently
            updated_scf_info['xc_screening_threshold'] = scf_results[
                'xc_screening_threshold']

        if scf_results.get('potfile', None) is not None:
            _warn_if_environment_overwritten(obj, 'potfile',
                                             scf_results['potfile'])
            updated_scf_info['potfile'] = scf_results['potfile']

        # the embedding method of the SCF reference (None when gas phase) is
        # always updated, so that a reused driver does not keep a stale value
        updated_scf_info['_scf_embedding_method'] = scf_results.get(
            'embedding_method', None)

        if scf_results.get('solvation_model', None) is not None:
            # the environment is inherited from SCF; a different solvation
            # model in response settings is overwritten with a warning
            if hasattr(obj, 'solvation_model'):
                _warn_if_environment_overwritten(
                    obj, 'solvation_model', scf_results['solvation_model'])
                for key in [
                        'solvation_model',
                        'cpcm_epsilon',
                        'cpcm_radii_scaling',
                        'cpcm_grid_per_sphere',
                        'cpcm_cg_thresh',
                        'cpcm_x',
                        'cpcm_custom_vdw_radii',
                ]:
                    if key in scf_results:
                        updated_scf_info[key] = scf_results[key]

        if scf_results.get('pressure', None) is not None:
            # the environment (including GOSTSHYP pressure) is inherited
            # from SCF; response cannot override it
            if hasattr(obj, 'pressure'):
                # warn only for explicitly set pressures: 0.0 is the default
                if obj.pressure != 0.0:
                    _warn_if_environment_overwritten(obj, 'pressure',
                                                     scf_results['pressure'])
                    _warn_if_environment_overwritten(obj, 'pressure_units',
                                                     scf_results['pressure_units'])
                for key in [
                        'pressure',
                        'pressure_units',
                        'gostshyp_num_lebedev_points',
                        'gostshyp_tssf',
                        'gostshyp_discretization',
                        'gostshyp_switching_thresh',
                        'gostshyp_r_ext',
                        'gostshyp_tco_tol',
                ]:
                    if key in scf_results:
                        updated_scf_info[key] = scf_results[key]

        if scf_results.get('electric_field', None) is not None:
            # the environment is inherited from SCF; a different electric
            # field in response settings is overwritten with a warning
            if hasattr(obj, 'electric_field'):
                _warn_if_environment_overwritten(obj, 'electric_field',
                                                 scf_results['electric_field'])
                updated_scf_info['electric_field'] = scf_results['electric_field']

    updated_scf_info = obj.comm.bcast(updated_scf_info, root=mpi_master())

    for key, val in updated_scf_info.items():
        setattr(obj, key, val)

    # double check xcfun in SCF and response

    if obj.rank == mpi_master():
        _warn_if_environment_absent_from_scf(obj, scf_results)

        scf_xcfun_label = scf_results.get('xcfun', 'HF').upper()
        if obj.xcfun is None:
            rsp_xcfun_label = 'HF'
        elif isinstance(obj.xcfun, str):
            rsp_xcfun_label = obj.xcfun.upper()
        else:
            rsp_xcfun_label = obj.xcfun.get_func_label().upper()
        if rsp_xcfun_label != scf_xcfun_label:
            warn_msg = f'{rsp_xcfun_label} will be used in response'
            warn_msg += f' but {scf_xcfun_label} was used in SCF.'
            warn_msg += ' Please double check.'
            obj.ostream.print_warning(warn_msg)
            obj.ostream.flush()


def nonlinear_response_environment_sanity_check(obj):
    """
    Checks environment settings in nonlinear response drivers.

    Solvation, pressure, polarizable embedding and the fika embedding are not
    supported in nonlinear response calculations. The settings are rejected
    here whether they were set explicitly or inherited from the SCF results,
    so that an unsupported environment is never silently ignored.

    :param obj:
        The nonlinear response driver.
    """

    if obj.potfile is not None:
        errmsg = "NonlinearSolver: The 'potfile' keyword is not supported "
        errmsg += 'in nonlinear response calculation.'
        if obj.rank == mpi_master():
            assert_msg_critical(False, errmsg)

    from .fikaembedding import is_fika_embedding

    if (getattr(obj, '_scf_embedding_method', None) == 'fika' or
            is_fika_embedding(getattr(obj, 'embedding', None))):
        errmsg = 'NonlinearSolver: The fika embedding is not supported in '
        errmsg += 'nonlinear response calculation.'
        if obj.rank == mpi_master():
            assert_msg_critical(False, errmsg)

    if obj.solvation_model is not None:
        errmsg = "NonlinearSolver: The 'solvation_model' keyword is not "
        errmsg += 'supported in nonlinear response calculation.'
        if obj.rank == mpi_master():
            assert_msg_critical(False, errmsg)

    if obj.pressure != 0.0:
        errmsg = "NonlinearSolver: The 'pressure' keyword is not supported "
        errmsg += 'in nonlinear response calculation.'
        if obj.rank == mpi_master():
            assert_msg_critical(False, errmsg)

    if obj.electric_field is not None:
        # the electric field is supported: it is a static one-electron term
        # that enters the calculation through the SCF; only validate it
        assert_msg_critical(
            len(obj.electric_field) == 3,
            "NonlinearSolver: Expecting 3 values in 'electric field' "
            'input')

    # point charges are also supported: they enter the calculation through
    # the SCF as a static one-electron term, with no driver attribute


def rsp_results_solvation_sanity_check(obj, rsp_results):
    """
    Checks response results for C-PCM information.

    :param obj:
        The object (TDDFT orbital response driver) that is being updated.
    :param rsp_results:
        A dictionary containing the linear response results.
    """

    updated_rsp_info = {}

    if obj.rank == mpi_master():
        if rsp_results.get('solvation_model', None) is not None:
            if hasattr(obj, 'solvation_model'):
                for key in [
                        'solvation_model',
                        'cpcm_optical_epsilon',
                        'non_equilibrium_solv',
                ]:
                    updated_rsp_info[key] = rsp_results[key]

    updated_rsp_info = obj.comm.bcast(updated_rsp_info, root=mpi_master())

    for key, val in updated_rsp_info.items():
        setattr(obj, key, val)


def dft_sanity_check(obj, method_flag='compute', response_flag='none'):
    """
    Checks DFT settings and updates relevant attributes.

    :param obj:
        The object (SCF or response driver) that is being updated.
    :param method_flag:
        The flag indicating the method in which the sanity check is
        called.
    :param response_flag:
        The flag indicating the type of response calculation in which the
        sanity check is called.
    """

    xcfun_is_none = (obj.xcfun is None)
    xcfun_is_hf = (isinstance(obj.xcfun, str) and obj.xcfun.lower() == 'hf')

    # Hartree-Fock: xcfun is None or 'hf'
    if xcfun_is_none or xcfun_is_hf:
        obj._dft = False

    # DFT: xcfun is functional object or string (other than 'hf')
    else:
        if isinstance(obj.xcfun, str):
            obj.xcfun = parse_xc_func(obj.xcfun.upper())
        assert_msg_critical(not obj.xcfun.is_undefined(),
                            f'{type(obj).__name__}: Undefined XC functional')
        obj._dft = True
        # update leading dimension (max batch size for DFT grid)
        obj.xcfun._set_leading_dimension(obj._xcfun_ldstaging)

    # check grid level
    if obj._dft and obj.grid_level is not None:
        if (obj.grid_level < 1 or obj.grid_level > 8):
            warn_msg = f'Invalid DFT grid level {obj.grid_level}. '
            warn_msg += 'Using default value.'
            obj.ostream.print_warning(warn_msg)
            obj.grid_level = None
        elif (method_flag.lower() == 'compute' and
              obj.grid_level < get_default_grid_level(obj.xcfun)):
            warn_msg = 'DFT grid level is below the recommended value. '
            warn_msg += 'Please double check.'
            obj.ostream.print_warning(warn_msg)
        obj.ostream.flush()

    # check if SCAN family of functional is used in nonliear response
    if obj._dft and response_flag.lower() == 'nonlinear':
        err_msg_scan = f'{type(obj).__name__}: Nonlinear response with '
        err_msg_scan += 'SCAN family of functional is not supported'
        assert_msg_critical('scan' not in obj.xcfun.get_func_label().lower(),
                            err_msg_scan)


def ri_sanity_check(obj):
    """
    Checks RI settings and updates relevant attributes.

    :param obj:
        The object (SCF or response driver) that is being updated.
    """

    ri_coulomb = getattr(obj, 'ri_coulomb', False)
    ri_jk = getattr(obj, 'ri_jk', False)

    if not (ri_coulomb or ri_jk):
        return

    if (ri_coulomb or ri_jk) and hasattr(obj, 'acc_type'):
        # RI uses plain DIIS
        if obj.acc_type.upper() == 'L2_DIIS':
            obj.acc_type = 'DIIS'

    if ri_coulomb and obj.ri_auxiliary_basis == 'def2-universal-jkfit':
        obj.ri_auxiliary_basis = 'def2-universal-jfit'

    if ri_jk and obj.ri_auxiliary_basis == 'def2-universal-jfit':
        obj.ri_auxiliary_basis = 'def2-universal-jkfit'


def polorbrsp_sanity_check_1(obj):
    """
    Checks frequencies in settings for polarizability orbital response.

    :param obj:
        The object (orbital response driver).
    :param method_flag:
        The flag indicating the method in which the sanity check is
        called.
    """
    # check that there is no zero frequency for the complex case
    # this can cause divergence with the subspace solver
    if obj.is_complex:
        if 0.0 in obj.frequencies:
            idx0 = obj.frequencies.index(0.0)
            warn_msg = 'Zero (0.0) frequency in input frequencies for complex'
            warn_msg += ' polarizability gradient/orbital response!\n\n'
            if len(obj.frequencies) == 1:
                warn_msg += 'No other frequencies requested;'
                warn_msg += ' Will continue with zero frequency,'
                warn_msg += ' CPHF solver might diverge.'
            else:
                # converting to a list because "pop()" does not exist for tuples
                freq_list = list(obj.frequencies)
                freq_list.pop(idx0)
                obj.frequencies = freq_list
                warn_msg += 'Zero (0.0) has been removed from the list of frequencies'
                warn_msg += ' due to risk of divergent CPHF solver.\n'
                warn_msg += 'Computations will be carried out for frequencies: '
                warn_msg += str(obj.frequencies)
            obj.ostream.print_warning(warn_msg)


def polorbrsp_sanity_check_2(obj, method_flag, lr_results):
    """
    Checks settings for polarizability orbital response against
    linear response results.

    :param obj:
        The object (orbital response driver).
    :param method_flag:
        The flag indicating the method in which the sanity check is
        called.
    :param lr_results:
        A dictionary containing linear response results.
    """

    if obj.rank == mpi_master():
        # check that frequencies agree with LR
        response_results = lr_results.get('solutions', None)
        for frequency in obj.frequencies:
            if (obj.vector_components[0],
                    frequency) not in response_results.keys():
                error_msg = f'Frequency {frequency:2.3f} in '
                error_msg += method_flag + ' not found in linear response results '
                error_msg += 'for vector compontent ' + obj.vector_components[0]
                raise ValueError(error_msg)


def polgrad_sanity_check_1(obj):
    """
    Checks frequencies in settings for polarizability gradient.

    :param obj:
        The object (polarizability gradient driver).
    :param method_flag:
        The flag indicating the method in which the sanity check is
        called.
    """

    polorbrsp_sanity_check_1(obj)


def polgrad_sanity_check_2(obj, method_flag, lr_results):
    """
    Checks settings for polarizability gradient against
    linear response results.

    :param obj:
        The object (polarizability gradient driver).
    :param method_flag:
        The flag indicating the method in which the sanity check is
        called.
    :param lr_results:
        A dictionary containing linear response results.
    """

    polorbrsp_sanity_check_2(obj, method_flag, lr_results)


# def polgrad_sanity_check(obj, method_flag, lr_results):
#     """
#     Checks settings for polarizability gradient and polarizability
#     orbital response against linear response results.
#
#     :param obj:
#         The object (polarizability gradient or orbital response driver).
#     :param method_flag:
#         The flag indicating the method in which the sanity check is
#         called.
#     :param lr_results:
#         A dictionary containing linear response results.
#     """
#
#     if obj.rank == mpi_master():
#         # check that frequencies agree with LR
#         response_results = lr_results.get('solutions', None)
#         for frequency in obj.frequencies:
#             if (obj.vector_components[0], frequency) not in response_results.keys():
#                 error_msg = f'Frequency {frequency:2.3f} in '
#                 error_msg += method_flag + ' not found in linear response results '
#                 error_msg += 'for vector compontent ' + obj.vector_components[0]
#                 raise ValueError(error_msg)
#
#     # check that there is no zero frequency for the complex case
#     # this can cause divergence with the subspace solver
#     if obj.is_complex:
#         try:
#             idx0 = obj.frequencies.index(0.0)
#             warn_msg = 'Zero (0.0) frequency in input frequencies for complex'
#             warn_msg += ' polarizability gradient/orbital response!\n\n'
#             if len(obj.frequencies) == 1:
#                 warn_msg += 'No other frequencies requested;'
#                 warn_msg += ' Will continue with zero frequency,'
#                 warn_msg += ' CPHF solver might diverge.'
#             else:
#                 # converting to a list because "pop()" does not exist for tuples
#                 freq_list = list(obj.frequencies)
#                 freq_list.pop(idx0)
#                 obj.frequencies = freq_list
#                 warn_msg += 'Zero (0.0) has been removed from the list of frequencies'
#                 warn_msg += ' due to risk of divergent CPHF solver.\n'
#                 warn_msg += 'Computations will be carried out for frequencies: '
#                 warn_msg += str(obj.frequencies)
#             obj.ostream.print_warning(warn_msg)
#         except ValueError:
#             pass


def raman_sanity_check(obj):
    """
    Checks whether both normal and resonance Raman has been requested.
    Print warning message and set normal Raman flag to False.
    The driver will continue with resonance Raman.

    Checks if zero frequency has been input to resonance Raman. If so,
    the value is removed from frequency list.

    :param obj:
        The object (vibrational analysis driver)
    """

    if (obj.do_raman and obj.do_resonance_raman):
        warn_msg = 'Both normal and resonance Raman have been requested, '
        warn_msg += 'but only one can run at a time.\n'
        obj.ostream.print_warning(warn_msg)
        info_msg = 'Will continue with resonance Raman.'
        obj.ostream.print_warning(info_msg)
        obj.ostream.flush()
        obj.do_raman = False

    # This check is due to convergence/singularity issues in the cphf
    # subspace solver for some molecules.
    if obj.do_resonance_raman:
        if 0.0 in obj.frequencies:
            idx0 = obj.frequencies.index(0.0)
            warn_msg = 'Zero frequency in input frequencies for resonance Raman!\n'
            if len(obj.frequencies) == 1:
                warn_msg += 'No other frequencies has been requested, computations'
                warn_msg += 'will continue with normal Raman.'
                obj.do_raman = True
                obj.do_resonance_raman = False
            else:
                # converting to a list because "pop()" does not exist for tuples
                freq_list = list(obj.frequencies)
                freq_list.pop(idx0)
                obj.frequencies = freq_list
                warn_msg += 'Zero (0.0) has been removed from the list of frequencies'
                warn_msg += ' due to risk of divergent CPHF solver.\n'
                warn_msg += 'Resonance Raman will be calculated for frequencies: '
                warn_msg += str(obj.frequencies)
            obj.ostream.print_warning(warn_msg)


def pe_sanity_check(obj, method_dict=None, molecule=None):
    """
    Checks PE settings and updates relevant attributes.

    A ``potfile`` is a fika potential file and takes priority over an
    assigned embedding.

    :param method_dict:
        The dictionary of method settings.
    :param molecule:
        The molecule (kept for backward compatibility).
    """

    if method_dict:
        if 'pe_options' in method_dict:
            obj.pe_options = dict(method_dict['pe_options'])
        else:
            obj.pe_options = {}

    from .fikaembedding import is_fika_embedding, fika_embedding_sanity_check

    if obj.potfile:
        potfile = None
        if obj.rank == mpi_master():
            potfile = obj.potfile
            if not Path(potfile).is_file():
                potfile = str(Path(obj.filename).parent / Path(potfile).name)
            assert_msg_critical(
                Path(potfile).is_file(),
                'PE sanity check: potfile does not exist')

        potfile = obj.comm.bcast(potfile, root=mpi_master())
        obj.embedding = {
            'settings': {
                'embedding_method': 'fika'
            },
            'inputs': {
                'potfile': potfile
            },
        }
        obj.potfile = None

    if getattr(obj, '_scf_embedding_method', None) == 'fika':
        # a response on a fika SCF reference needs the fika embedding as well
        assert_msg_critical(
            is_fika_embedding(obj.embedding),
            f'{type(obj).__name__}: The SCF reference used the fika ' +
            'embedding; set the same fika embedding on the response driver')

    if obj.embedding is not None and not is_fika_embedding(obj.embedding):
        assert_msg_critical(
            False, f'{type(obj).__name__}: only the fika embedding is ' +
            'supported')

    if is_fika_embedding(obj.embedding):
        fika_embedding_sanity_check(obj.embedding)
        obj.pe_options = {}
        obj._pe = True
        return

    obj._pe = False


def gostshyp_sanity_check(obj, basis=None):
    """
    Checks the GOSTSHYP settings and updates relevant attributes.

    :param basis:
        The basis set.
    """

    if obj.pressure is None:
        obj.pressure = 0.0

    obj._gostshyp = (obj.pressure != 0.0)

    if obj._gostshyp:
        assert_msg_critical(obj.pressure > 0.0,
                            'GOSTSHYP: Unphysical negative pressures invalid')

        assert_msg_critical(
            np.isfinite(obj.gostshyp_tssf) and obj.gostshyp_tssf > 0.0,
            'GOSTSHYP: Tessellation sphere scaling factor must be finite and '
            'positive')

        assert_msg_critical(
            obj.solvation_model is None,
            type(obj).__name__ +
            ': GOSTSHYP is incompatible with the solvation model')

        valid_discretizations = ('fixed', 'swig', 'iswig')
        assert_msg_critical(
            isinstance(obj.gostshyp_discretization, str) and
            obj.gostshyp_discretization.lower() in valid_discretizations,
            "GOSTSHYP: Invalid discretization. Valid options are 'fixed', "
            "'swig', and 'iswig'")

        assert_msg_critical(
            obj.gostshyp_tco_tol > 0.0,
            'GOSTSHYP: Three-center overlap integral screening threshold '
            'must be positive')

        from .tessellation import validate_num_lebedev_points
        obj.gostshyp_num_lebedev_points = validate_num_lebedev_points(
            obj.gostshyp_num_lebedev_points,
            obj.gostshyp_discretization,
            obj.ostream)

        if basis is not None:

            # check max angular momentum of basis
            max_am = basis.max_angular_momentum()
            gost_max_am = 3

            # TODO: consider moving this check to gostshyp 1e integral drivers
            if max_am > gost_max_am:
                labels = {0: 's', 1: 'p', 2: 'd', 3: 'f', 4: 'g', 5: 'h'}
                assert_msg_critical(
                    False,
                    f"GOSTSHYP supports basis functions up to f-type; this basis "
                    f"contains {labels.get(max_am, f'l={max_am}')}-type functions.")


def solvation_model_sanity_check(obj):
    """
    Checks solvation model and updates relevant attributes.
    """

    if obj.solvation_model is not None:
        assert_msg_critical(
            not obj._pe,
            type(obj).__name__ +
            ': The \'solvation_model\' option is incompatible with ' +
            'polarizable embedding')

        assert_msg_critical(
            getattr(obj, 'point_charges', None) is None,
            type(obj).__name__ +
            ': The \'solvation_model\' option is incompatible with ' +
            'point charges')

        assert_msg_critical(
            obj.pressure == 0.0,
            type(obj).__name__ +
            ': The \'solvation_model\' option is incompatible with ' +
            'GOSTSHYP (hydrostatic pressure)')

        assert_msg_critical(
            obj.solvation_model.lower() in ['cpcm', 'c-pcm', 'c_pcm', 'smd'],
            type(obj).__name__ +
            ': Only the C-PCM and SMD solvation models are implemented.')

        obj._cpcm = True
        obj._smd = (obj.solvation_model.lower() == 'smd')

    else:
        obj._cpcm = False
        obj._smd = False


def environment_compatibility_sanity_check(obj):
    """
    Checks the pairwise compatibility of the environment settings.

    The environment settings are: ``potfile`` (polarizable embedding),
    ``solvation_model`` (C-PCM / SMD), ``pressure`` (GOSTSHYP),
    ``electric_field`` and ``point_charges``. Each setting can be used
    individually. The supported combinations are ``electric_field``
    together with ``point_charges`` and ``electric_field`` together
    with ``pressure`` (GOSTSHYP): the electric field is a static
    one-electron term with no mutual response loop, and the GOSTSHYP
    pressure potential is converged within the SCF iterations rather
    than in a separate self-consistent loop. All other combinations
    are rejected, because the density-coupled response loops of the
    other environments (induced dipoles, surface charges) are not
    supported in combination.

    The helper is intended to be reusable: it only relies on the driver
    attributes that declare the five settings (all queried with
    ``getattr``), so it can be called from any driver that carries them,
    e.g. after the environment sanity checks in ``ScfDriver``.

    :param obj:
        The driver.
    """

    has_pe = (bool(getattr(obj, '_pe', False)) or
              bool(getattr(obj, 'potfile', None)))
    has_cpcm = getattr(obj, 'solvation_model', None) is not None
    pressure = getattr(obj, 'pressure', 0.0)
    has_gostshyp = (bool(getattr(obj, '_gostshyp', False)) or
                    (pressure is not None and pressure != 0.0))
    has_field = getattr(obj, 'electric_field', None) is not None
    has_pc = getattr(obj, 'point_charges', None) is not None

    if has_field:
        assert_msg_critical(
            len(obj.electric_field) == 3,
            f"{type(obj).__name__}: Expecting 3 values in "
            "'electric field' input")

    names = []
    if has_pe:
        names.append("'potfile' (polarizable embedding)")
    if has_cpcm:
        names.append("'solvation_model'")
    if has_gostshyp:
        names.append("'pressure' (GOSTSHYP)")
    if has_field:
        names.append("'electric_field'")
    if has_pc:
        names.append("'point_charges'")

    if len(names) <= 1:
        return

    if has_field and len(names) == 2 and (has_pc or has_gostshyp):
        # supported combinations: 'electric_field' with 'point_charges'
        # (two static one-electron terms) and 'electric_field' with
        # 'pressure' (GOSTSHYP): the field is a static one-electron term
        # added to the density-dependent GOSTSHYP potential within the
        # SCF iterations, with no mutual response loop
        return

    err_msg = f"{type(obj).__name__}: Incompatible environment settings: "
    err_msg += ', '.join(names) + '. '
    err_msg += "Environment settings can only be used individually; the "
    err_msg += "supported combinations are 'electric_field' with "
    err_msg += "'point_charges' and 'electric_field' with 'pressure' "
    err_msg += "(GOSTSHYP)."
    assert_msg_critical(False, err_msg)


