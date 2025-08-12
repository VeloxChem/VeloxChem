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
import json

from .veloxchemlib import bohr_in_angstroms
from .veloxchemlib import parse_xc_func
from .veloxchemlib import mpi_master
from .inputparser import InputParser
from .dftutils import get_default_grid_level
from .errorhandler import assert_msg_critical


def molecule_sanity_check(mol):
    """
    Checks molecule for charge/multiplicity combination and geometry.

    :param mol:
        The molecule.
    """

    assert_msg_critical(
        mol.check_multiplicity(),
        'Molecule: Incompatible multiplicity and number of electrons')
    assert_msg_critical(mol.check_proximity(0.1), 'Molecule: Atoms too close')


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
        if 'filename' in scf_results:
            # force filename to be the same since we want a single final h5
            # file to store all results
            updated_scf_info['filename'] = scf_results['filename']

        if scf_results.get('eri_thresh', None) is not None:
            updated_scf_info['eri_thresh'] = scf_results['eri_thresh']

        if scf_results.get('ri_coulomb', None) is not None:
            # TODO: consider storing ri_aux_basis_obj in scf_results
            updated_scf_info['ri_coulomb'] = scf_results['ri_coulomb']
            updated_scf_info['ri_auxiliary_basis'] = scf_results[
                'ri_auxiliary_basis']

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

        if scf_results.get('solvation_model', None) is not None:
            # do not overwrite solvation_model if it is already specified
            if hasattr(obj, 'solvation_model') and obj.solvation_model is None:
                for key in [
                        'solvation_model',
                        'cpcm_epsilon',
                        'cpcm_grid_per_sphere',
                        'cpcm_cg_thresh',
                        'cpcm_x',
                        'cpcm_custom_vdw_radii',
                ]:
                    updated_scf_info[key] = scf_results[key]

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
        try:
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
        except ValueError:
            pass


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
            if (obj.vector_components[0], frequency) not in response_results.keys():
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


#def polgrad_sanity_check(obj, method_flag, lr_results):
#    """
#    Checks settings for polarizability gradient and polarizability
#    orbital response against linear response results.
#
#    :param obj:
#        The object (polarizability gradient or orbital response driver).
#    :param method_flag:
#        The flag indicating the method in which the sanity check is
#        called.
#    :param lr_results:
#        A dictionary containing linear response results.
#    """
#
#    if obj.rank == mpi_master():
#        # check that frequencies agree with LR
#        response_results = lr_results.get('solutions', None)
#        for frequency in obj.frequencies:
#            if (obj.vector_components[0], frequency) not in response_results.keys():
#                error_msg = f'Frequency {frequency:2.3f} in '
#                error_msg += method_flag + ' not found in linear response results '
#                error_msg += 'for vector compontent ' + obj.vector_components[0]
#                raise ValueError(error_msg)
#
#    # check that there is no zero frequency for the complex case
#    # this can cause divergence with the subspace solver
#    if obj.is_complex:
#        try:
#            idx0 = obj.frequencies.index(0.0)
#            warn_msg = 'Zero (0.0) frequency in input frequencies for complex'
#            warn_msg += ' polarizability gradient/orbital response!\n\n'
#            if len(obj.frequencies) == 1:
#                warn_msg += 'No other frequencies requested;'
#                warn_msg += ' Will continue with zero frequency,'
#                warn_msg += ' CPHF solver might diverge.'
#            else:
#                # converting to a list because "pop()" does not exist for tuples
#                freq_list = list(obj.frequencies)
#                freq_list.pop(idx0)
#                obj.frequencies = freq_list
#                warn_msg += 'Zero (0.0) has been removed from the list of frequencies'
#                warn_msg += ' due to risk of divergent CPHF solver.\n'
#                warn_msg += 'Computations will be carried out for frequencies: '
#                warn_msg += str(obj.frequencies)
#            obj.ostream.print_warning(warn_msg)
#        except ValueError:
#            pass


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
        try:
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
        except ValueError:
            pass


def pe_sanity_check(obj, method_dict=None, molecule=None):
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

    obj._pe = (('potfile' in obj.pe_options) or (obj.embedding is not None))

    if obj._pe:
        if obj.embedding is None:
            potfile = None
            if obj.rank == mpi_master():
                potfile = obj.pe_options['potfile']
                if not Path(potfile).is_file():
                    potfile = str(
                        Path(obj.filename).parent / Path(potfile).name)
                assert_msg_critical(
                    Path(potfile).is_file(),
                    'PE sanity check: potfile does not exist')

            potfile = obj.comm.bcast(potfile, root=mpi_master())
            obj.pe_options['potfile'] = potfile

            # TODO: include more options from pe_options
            obj.embedding = {
                'settings': {
                    'embedding_method': 'PE',
                    'induced_dipoles': {
                        'solver': 'jacobi',
                        'mic': False,
                        'threshold': 1e-8,
                        'max_iterations': 100,
                    },
                },
                'inputs': {
                    'json_file': potfile,
                },
            }
        else:
            potfile = obj.embedding['inputs']['json_file']
            obj.pe_options['potfile'] = potfile

        # update potfile in case it is not in json format

        if Path(potfile).suffix != '.json' and molecule is not None:
            if obj.rank == mpi_master():
                new_potfile = write_pe_jsonfile(molecule, potfile)
                potfile = new_potfile

            potfile = obj.comm.bcast(potfile, root=mpi_master())
            obj.pe_options['potfile'] = potfile
            obj.embedding['inputs']['json_file'] = potfile

        embedding_sanity_check(obj.embedding)


def embedding_sanity_check(options):
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
            obj.point_charges is None,
            type(obj).__name__ +
            ': The \'solvation_model\' option is incompatible with ' +
            'point charges')

        assert_msg_critical(
            obj.solvation_model.lower() in ['cpcm', 'c-pcm', 'c_pcm', 'smd'],
            type(obj).__name__ +
            ': Only the C-PCM and SMD solvation models are implemented.')

        obj._cpcm = True
        
        if obj.solvation_model.lower() == 'smd':
            obj._smd = True

    else:
        obj._cpcm = False


def write_pe_jsonfile(molecule, potfile):
    """
    Writes potential json file for PyFraME.

    :param molecule:
        The molecule.
    :param potfile:
        The name of the potential file.

    :return:
        The name of the temporary potential json file.
    """

    # read potential file (new format)

    pe_inp = InputParser(potfile).input_dict

    # process units

    if 'units' in pe_inp['environment']:
        units = pe_inp['environment']['units'].lower()
        if len(units) >= 3 and units == 'angstrom'[:len(units)]:
            prefac = 1.0
        elif units in ['au', 'bohr']:
            prefac = bohr_in_angstroms()
        else:
            assert_msg_critical(False, 'potential file: invalid units')
    else:
        prefac = 1.0

    # process residues (a dictionary with residue ID as key)

    residues = {}
    for line in pe_inp['environment']['xyz']:
        val = line.split()
        resname = val[4]
        resid = int(val[5])
        if resid not in residues:
            residues[resid] = {'resname': resname, 'atoms': []}
        residues[resid]['atoms'].append([
            val[0],
            float(val[1]) * prefac,
            float(val[2]) * prefac,
            float(val[3]) * prefac,
        ])

    # process charges (a dictionary with residue name as key)

    charges = {}
    if 'charges' in pe_inp:
        for line in pe_inp['charges']:
            val = line.split()
            resname = val[2]
            if resname not in charges:
                charges[resname] = []
            charges[resname].append(val[1])

    # process polarizabilities (a dictionary with residue name as key)

    polarizabilities = {}
    if 'polarizabilities' in pe_inp:
        for line in pe_inp['polarizabilities']:
            val = line.split()
            resname = val[7]
            if resname not in polarizabilities:
                polarizabilities[resname] = []
            polarizabilities[resname].append(val[1:7])

    # create dictionary for json

    # QM atoms

    qm_nuclei = []

    qm_coords = molecule.get_coordinates_in_bohr()
    qm_labels = molecule.get_labels()
    qm_elem_ids = molecule.get_element_ids()

    for atom_idx in range(molecule.number_of_atoms()):
        # Note: make sure all elements in qm_nuclei are
        # serializable by json
        qm_nuclei.append({
            'index': atom_idx + 1,
            'element': qm_labels[atom_idx].capitalize(),
            'charge': float(qm_elem_ids[atom_idx]),
            'coordinate': list(qm_coords[atom_idx]),
        })

    embedding_json = {
        "quantum_subsystems": [{
            "nuclei": qm_nuclei,
        }],
    }

    # MM atoms

    classical_fragments = []

    mm_atom_count = 0

    for res_count, resid in enumerate(sorted(list(residues.keys()))):

        classical_fragments.append({
            "index": res_count + 1,
            "name": residues[resid]['resname'],
            "atoms": [],
        })

        # charges

        if resname in charges:
            assert len(charges[resname]) == len(residues[resid]['atoms'])
            atom_chgs = [float(x) for x in charges[resname]]
        else:
            atom_chgs = [0.0 for x in range(len(residues[resid]['atoms']))]

        # polarizabilities

        if resname in polarizabilities:
            assert len(polarizabilities[resname]) == len(
                residues[resid]['atoms'])
            atom_pols = [
                [float(p) for p in pol] for pol in polarizabilities[resname]
            ]
        else:
            atom_pols = [[0.0
                          for p in range(6)]
                         for x in range(len(residues[resid]['atoms']))]

        # coordinates

        res_atom_start = mm_atom_count

        for atom_idx, atom in enumerate(residues[resid]['atoms']):
            # Note: make sure all elements in classical_fragments
            # are serializable by json
            classical_fragments[-1]["atoms"].append({
                "index": mm_atom_count + 1,
                "element": atom[0].capitalize(),
                "coordinate": [
                    float(atom[1]) / bohr_in_angstroms(),
                    float(atom[2]) / bohr_in_angstroms(),
                    float(atom[3]) / bohr_in_angstroms(),
                ],
                "multipoles": {
                    "elements": [atom_chgs[atom_idx]],
                },
                "exclusions": list(
                    range(res_atom_start + 1,
                          res_atom_start + 1 + len(residues[resid]['atoms']))),
                "polarizabilities": {
                    "elements": ([0.0, 0.0, 0.0, 0.0] + atom_pols[atom_idx]),
                    "order": [1, 1],
                },
            })

            mm_atom_count += 1

    embedding_json.update({
        "classical_subsystems": [{
            "classical_fragments": classical_fragments,
        }],
    })

    # write json file

    pe_jsonfile = Path(potfile).with_suffix('.json')

    with open(pe_jsonfile, 'w') as fh:
        json.dump(embedding_json, fh, indent=4)

    return str(pe_jsonfile)
