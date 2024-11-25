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

from pathlib import Path

from .veloxchemlib import parse_xc_func
from .veloxchemlib import mpi_master
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


def pe_sanity_check(obj, method_dict=None):
    """
    Checks PE settings and updates relevant attributes.

    :param method_dict:
        The dicitonary of method settings.
    """

    if method_dict:
        if 'pe_options' in method_dict:
            obj.pe_options = dict(method_dict['pe_options'])
        else:
            obj.pe_options = {}

    if obj.potfile:
        obj.pe_options['potfile'] = obj.potfile

    obj._pe = (('potfile' in obj.pe_options) or
               (obj.embedding_options is not None))

    if obj._pe:
        if obj.embedding_options is None:
            potfile = None
            if obj.rank == mpi_master():
                potfile = obj.pe_options['potfile']
                if not Path(potfile).is_file():
                    potfile = str(
                        Path(obj.filename).parent / Path(potfile).name)
            potfile = obj.comm.bcast(potfile, root=mpi_master())
            obj.pe_options['potfile'] = potfile
            # TODO: include more options from pe_options
            obj.embedding_options = {
                'settings': {
                    'embedding_method': 'PE',
                },
                'inputs': {
                    'json_file': potfile,
                },
            }
        else:
            potfile = obj.embedding_options['inputs']['json_file']
            obj.pe_options['potfile'] = potfile

        embedding_options_sanity_check(obj.embedding_options)


def embedding_options_sanity_check(options):
    """
    Checks the validity of the given options dictionary.

    :param options: Dictionary with settings and inputs.
        settings:
            - embedding_method: string to set embedding method.
            - vdw: dictionary with keys: method and combination_rule.
                - method: string
                - combination_rule: string
            - induced_dipoles: dictionary with keys: threshold, max_iterations,
                               and solver.
                - solver: string
                - threshold: float
                - max_iterations: integer
            - environment_energy: bool which decides if the environment energy
                                  will be calculated.
        inputs:
            - json_file: string that is the path to the json file that contains
                         the embedding potentials.
            - objects:
                - quantum_subsystem: PyFraME class subsystem.QuantumSubsystem.
                - classical_subsystem: PyFraME class subsystem.ClassicalSubsystem.
    """

    from pyframe.embedding import subsystem

    if not isinstance(options, dict):
        raise TypeError("The options parameter must be a dictionary.")
    if 'settings' not in options:
        raise KeyError("Missing 'settings' key in options dictionary.")

    settings = options['settings']

    if 'embedding_method' not in settings or not isinstance(
            settings['embedding_method'], str):
        raise KeyError("Missing or invalid 'embedding_method' in settings.")

    if 'vdw' in settings:
        if not isinstance(settings['vdw'], dict):
            raise TypeError("'vdw' must be a dictionary.")

        vdw = settings['vdw']

        if 'method' not in vdw or not isinstance(vdw['method'], str):
            raise KeyError("Missing or invalid 'method' in 'vdw'.")
        if 'combination_rule' not in vdw or not isinstance(
                vdw['combination_rule'], str):
            raise KeyError("Missing or invalid 'combination_rule' in 'vdw'.")

    if 'induced_dipoles' in settings:
        if not isinstance(settings['induced_dipoles'], dict):
            raise TypeError("'induced_dipoles' must be a dictionary.")

        dipoles = settings['induced_dipoles']

        if 'solver' in dipoles and not isinstance(dipoles['solver'], str):
            raise KeyError("'solver' must be a string.")
        if 'threshold' in dipoles and not isinstance(dipoles['threshold'],
                                                     (int, float)):
            raise KeyError("'threshold' must be an integer or float.")
        if 'max_iterations' in dipoles and not isinstance(
                dipoles['max_iterations'], int):
            raise KeyError("'max_iterations' must be an integer.")
        if 'max_iterations' in dipoles and dipoles['max_iterations'] <= 0:
            raise ValueError("'max_iterations' must be a positive integer.")

    if 'environment_energy' in settings:
        if not isinstance(settings['environment_energy'], bool):
            raise TypeError("'environment_energy' must be a boolean.")

    if 'inputs' not in options:
        raise KeyError("Missing 'inputs' key in options dictionary.")

    inputs = options['inputs']

    if 'json_file' in inputs and not isinstance(inputs['json_file'], str):
        raise TypeError("'json_file' must be a string.")

    if 'objects' in inputs:
        if not isinstance(inputs['objects'], dict):
            raise TypeError("'objects' must be a dictionary.")

        objects = inputs['objects']

        if 'quantum_subsystem' not in objects:
            raise KeyError("Missing 'quantum_subsystem' in 'objects'.")
        if 'classical_subsystem' not in objects:
            raise KeyError("Missing 'classical_subsystem' in 'objects'.")
        if not isinstance(objects['quantum_subsystem'],
                          subsystem.QuantumSubsystem):
            raise TypeError(
                "'quantum_subsystem' must be an instance of PyFraME class " +
                "subsystem.QuantumSubsystem.")
        if not isinstance(objects['classical_subsystem'],
                          subsystem.ClassicalSubsystem):
            raise TypeError(
                "'classical_subsystem' must be an instance of PyFraME class " +
                "subsystem.ClassicalSubsystem.")

    if 'json_file' not in inputs and 'objects' not in inputs:
        raise KeyError(
            "At least one of 'json_file' or 'objects' must be provided in 'inputs'."
        )
