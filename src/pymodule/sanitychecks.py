#
#                              VELOXCHEM
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
#
#  SPDX-License-Identifier: LGPL-3.0-or-later
#
#  This file is part of VeloxChem.
#
#  VeloxChem is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Lesser General Public License as published by the Free
#  Software Foundation, either version 3 of the License, or (at your option)
#  any later version.
#
#  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

from .veloxchemlib import mpi_master


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

        if scf_results.get('qq_type', None) is not None:
            updated_scf_info['qq_type'] = scf_results['qq_type']

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
