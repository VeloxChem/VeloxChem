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

from .veloxchemlib import mpi_master
from .veloxchemlib import parse_xc_func
from .dftutils import get_default_grid_level
from .errorhandler import assert_msg_critical


def molecule_sanity_check(mol):
    """
    Checks molecule for charge/multiplicity combination and geometry.

    :param mol:
        The molecule.
    """

    mol.check_multiplicity()
    mol.check_proximity(0.1)


def scf_results_sanity_check(obj, scf_results):
    """
    Checks SCF results for ERI, DFT and PE information.

    :param obj:
        The object (response driver) that is being updated.
    :param scf_results:
        A dictionary containing SCF results.
    """

    updated_scf_info = {}

    if obj.rank == mpi_master():
        if scf_results.get('eri_thresh', None) is not None:
            updated_scf_info['eri_thresh'] = scf_results['eri_thresh']

        if scf_results.get('restart', None) is not None:
            # do not restart if scf is not restarted from checkpoint
            if not scf_results['restart']:
                updated_scf_info['restart'] = scf_results['restart']

        if scf_results.get('xcfun', None) is not None:
            # do not overwrite xcfun if it is already specified
            if obj.xcfun is None:
                updated_scf_info['xcfun'] = scf_results['xcfun']
                if 'grid_level' in scf_results:
                    updated_scf_info['grid_level'] = scf_results['grid_level']

        if scf_results.get('potfile', None) is not None:
            # do not overwrite potfile if it is already specified
            if obj.potfile is None:
                updated_scf_info['potfile'] = scf_results['potfile']

    updated_scf_info = obj.comm.bcast(updated_scf_info, root=mpi_master())

    for key, val in updated_scf_info.items():
        setattr(obj, key, val)

    # double check xcfun in SCF and response

    if obj.rank == mpi_master():
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
