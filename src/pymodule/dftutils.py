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

from .veloxchemlib import xcfun as xcfun_enum
from .veloxchemlib import parse_xc_func
from .veloxchemlib import XCFunctional
from .veloxchemlib import XCIntegrator
from .veloxchemlib import get_number_of_threads
from .errorhandler import assert_msg_critical


def get_default_grid_level(xc_func):
    """
    Gets default grid level for an exchange-correlation functional.

    :param xc_func:
        The exchange-correlation functional.

    :return:
        The default grid level.
    """

    if isinstance(xc_func, str):
        xc_func_obj = parse_xc_func(xc_func.upper())
        func_name = xc_func_obj.get_func_label()
        func_type = xc_func_obj.get_func_type()
    else:
        func_name = xc_func.get_func_label()
        func_type = xc_func.get_func_type()

    # LDA

    if func_type == xcfun_enum.lda:

        return 4

    # GGA

    elif func_type == xcfun_enum.gga:

        if func_name.upper() in [
                'B97',
                'B97-1',
                'B97-2',
                'B97-3',
        ]:
            return 5

        else:
            return 4

    # meta-GGA

    elif func_type == xcfun_enum.mgga:

        if func_name.upper() in [
                'SCAN',
        ]:
            return 7

        elif func_name.upper() in [
                'RSCAN',
                'R2SCAN',
        ]:
            return 6

        elif func_name.upper() in [
                'M05',
                'M05-2X',
                'M06',
                'M06-2X',
                'M06-HF',
                'M06-L',
                'M11-L',
        ]:
            return 6

        else:
            return 5

    else:
        assert_msg_critical(
            False, 'get_default_grid_level: Invalid XC functional type')


def get_optimal_grid_box_size(nbf):
    """
    Gets the optimal DFT grid box size for a given number of basis functions,
    according to an empirically determined nbf-only rule.

    :param nbf:
        The number of basis functions.

    :return:
        The optimal grid box size.
    """

    if nbf < 430:
        return 128
    elif nbf < 1500:
        return 256
    elif nbf < 2900:
        return 512
    else:
        return 1024


def get_response_fock_xc_density_memory(molecule,
                                      basis,
                                      dft_dict,
                                      xcfun,
                                      xc_screening_threshold=None):
    """
    Gets the memory the exchange-correlation part of a response Fock build
    (integrate_fxc_fock, integrate_kxc_fock, integrate_kxclxc_fock) holds on
    this rank for each density it slices, in bytes. Ground-state Fock builds
    and the gradients and Hessians are not what this counts. The first call
    for a dft_dict keeps the answer in it as 'response_fock_xc_density_memory',
    and later calls with the same dft_dict return it.

    :param molecule:
        The molecule.
    :param basis:
        The AO basis set.
    :param dft_dict:
        The dictionary of DFT information a solver's _init_dft returns, whose
        'molgrid' the integration takes.
    :param xcfun:
        The exchange-correlation functional, or None without one.
    :param xc_screening_threshold:
        The screening threshold the solver sets on its XCIntegrator, or None
        where it sets none.

    :return:
        The memory in bytes, 0 without a functional.

    .. note::
        integrate_fxc_fock, integrate_kxc_fock and integrate_kxclxc_fock give
        each grid box to one thread, and a thread copies every density of the
        call into a matrix over the basis functions the box keeps. With
        get_number_of_threads() threads, as many boxes are held at once, so the
        most this takes for one density is the sum of aocount^2 over that many
        of the largest boxes. A caller multiplies this by the densities a build
        slices: two for each density of an open shell (alpha and beta), and
        for the nonlinear builds the first-order densities as well. What a
        call holds besides its sliced densities is
        get_response_fock_xc_call_memory.
    """

    if xcfun is None:
        return 0

    if 'response_fock_xc_density_memory' not in dft_dict:
        xc_drv = XCIntegrator()

        if xc_screening_threshold is not None:
            xc_drv.set_screening_threshold(xc_screening_threshold)

        counts = xc_drv.get_largest_box_ao_counts_for_response_fock_xc(
            molecule, basis, dft_dict['molgrid'], xcfun,
            get_number_of_threads())

        dft_dict['response_fock_xc_density_memory'] = 8 * sum(
            c * c for c in counts)

    return dft_dict['response_fock_xc_density_memory']


# NOTE: the matrices a thread holds for one grid box in integrate_fxc_fock,
# integrate_kxc_fock and integrate_kxclxc_fock besides the densities it
# slices, as (matrices over the box's functions squared, matrices over its
# functions times its points), by functional type and shell:
# - the sliced ground-state density, one for each spin;
# - the GTO values: 1 for LDA; 4 (values and gradient) for GGA and meta-GGA;
# - the matrices one density's contribution is formed in, one density at a
#   time: G, and the partial Fock matrices that G times the GTO values gives.
#   LDA: 1 G and 1 partial for each spin. GGA: 2 G (G, G_gga) and 2 partials
#   for each spin. Meta-GGA: 5 G (G, G_gga, G_gga_x, G_gga_y, G_gga_z) and 5
#   partials.
# The kxc and kxclxc integrations hold the same as the fxc integration of the
# closed shell (integratePartialKxcFock* and integratePartialLxcFock* form the
# same G and partials), and have no open-shell form; nor does the fxc
# integration with a meta-GGA.
_RESPONSE_FOCK_XC_BOX_MATRICES = {
    (xcfun_enum.lda, False): (2, 2),
    (xcfun_enum.gga, False): (3, 6),
    (xcfun_enum.mgga, False): (6, 9),
    (xcfun_enum.lda, True): (4, 3),
    (xcfun_enum.gga, True): (6, 8),
}


def get_response_fock_xc_call_memory(molecule,
                                     basis,
                                     dft_dict,
                                     xcfun,
                                     open_shell=False,
                                     xc_screening_threshold=None):
    """
    Gets the memory the exchange-correlation part of a response Fock build
    (integrate_fxc_fock, integrate_kxc_fock, integrate_kxclxc_fock) holds on
    this rank for each call, whatever the number of densities, in bytes. The
    densities the call slices come on top, from
    get_response_fock_xc_density_memory. The first call for a dft_dict keeps
    the largest boxes in it as 'response_fock_xc_largest_boxes', and later
    calls with the same dft_dict take them from there.

    :param molecule:
        The molecule.
    :param basis:
        The AO basis set.
    :param dft_dict:
        The dictionary of DFT information a solver's _init_dft returns, whose
        'molgrid' the integration takes.
    :param xcfun:
        The exchange-correlation functional, or None without one.
    :param open_shell:
        True where the call integrates open-shell densities.
    :param xc_screening_threshold:
        The screening threshold the solver sets on its XCIntegrator, or None
        where it sets none.

    :return:
        The memory in bytes, 0 without a functional.

    .. note::
        With get_number_of_threads() threads, as many boxes are held at once.
        For each of the boxes that keep the most functions, a thread holds the
        matrices _RESPONSE_FOCK_XC_BOX_MATRICES counts over the box's
        functions squared and over its functions times its points. Arrays of
        one value per point (the density and the functional derivatives on
        the grid) are not counted.
    """

    if xcfun is None:
        return 0

    xcfun_type = xcfun.get_func_type()

    assert_msg_critical(
        (xcfun_type, bool(open_shell)) in _RESPONSE_FOCK_XC_BOX_MATRICES,
        'get_response_fock_xc_call_memory: the response Fock XC integration '
        'is implemented for LDA, GGA and closed-shell meta-GGA')

    if 'response_fock_xc_largest_boxes' not in dft_dict:
        xc_drv = XCIntegrator()

        if xc_screening_threshold is not None:
            xc_drv.set_screening_threshold(xc_screening_threshold)

        dft_dict['response_fock_xc_largest_boxes'] = (
            xc_drv.get_largest_boxes_for_response_fock_xc(
                molecule, basis, dft_dict['molgrid'], xcfun,
                get_number_of_threads()))

    n_square, n_points = _RESPONSE_FOCK_XC_BOX_MATRICES[(xcfun_type,
                                                         bool(open_shell))]

    return 8 * sum(n_square * aocount * aocount + n_points * aocount * npoints
                   for aocount, npoints in
                   dft_dict['response_fock_xc_largest_boxes'])


def print_xc_reference(xcfun, ostream):
    """
    Prints libxc reference.

    :param xcfun:
        The XC functional.
    :param ostream:
        The output stream.
    """

    if isinstance(xcfun, XCFunctional):
        valstr = f'Using the {xcfun.get_func_label()} functional.'
        ostream.print_info(valstr)
        ostream.print_blank()
        printed_refs = []
        for ref in xcfun.get_functional_reference():
            if ref not in printed_refs:
                ostream.print_reference(ref)
                printed_refs.append(ref)
        ostream.print_blank()

        valstr = 'Using the Libxc library '
        valstr += f'(v{xcfun.get_libxc_version()}).'
        ostream.print_info(valstr)
        ostream.print_blank()
        ostream.print_reference(xcfun.get_libxc_reference())
        ostream.print_blank()

        valstr = 'Using the following algorithm for XC numerical integration.'
        ostream.print_info(valstr)
        ostream.print_blank()
        valstr = 'J. Kussmann, H. Laqua and C. Ochsenfeld, '
        valstr += 'J. Chem. Theory Comput. 2021, 17, 1512-1521'
        ostream.print_reference(valstr)
        ostream.print_blank()

        ostream.flush()
