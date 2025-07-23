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
import time as tm
import h5py

from .veloxchemlib import mpi_master
from .distributedarray import DistributedArray
from .errorhandler import assert_msg_critical
from .molecule import Molecule

def create_hdf5(fname, molecule, basis, dft_func_label, potfile_text, orbital_details=None, eigenvals=None, excitation_details=None):
    """
    Creates HDF5 file for a calculation.

    :param fname:
        Name of the HDF5 file.
    :param molecule:
        The molecule.
    :param basis:
        The AO basis set.
    :param dft_func_label:
        The name of DFT functional.
    :param potfile_text:
        The content of potential file for polarizable embedding.
    """

    valid_checkpoint = (fname and isinstance(fname, str))

    if valid_checkpoint:
        hf = h5py.File(fname, 'w')

        hf.create_dataset('nuclear_repulsion',
                          data=np.array([molecule.nuclear_repulsion_energy()]))

        hf.create_dataset('nuclear_charges', data=molecule.get_element_ids())

        hf.create_dataset('atom_coordinates',
                          data=molecule.get_coordinates_in_bohr())

        hf.create_dataset('number_of_atoms',
                          data=np.array([molecule.number_of_atoms()]))

        hf.create_dataset('number_of_alpha_electrons',
                          data=np.array([molecule.number_of_alpha_electrons()]))

        hf.create_dataset('number_of_beta_electrons',
                          data=np.array([molecule.number_of_beta_electrons()]))

        hf.create_dataset('molecular_charge',
                          data=np.array([molecule.get_charge()]))

        hf.create_dataset('spin_multiplicity',
                          data=np.array([molecule.get_multiplicity()]))

        hf.create_dataset('basis_set', data=np.bytes_([basis.get_label()]))

        hf.create_dataset('dft_func_label', data=np.bytes_([dft_func_label]))

        hf.create_dataset('potfile_text', data=np.bytes_([potfile_text]))

        if eigenvals is not None:
            hf.create_dataset('eigenvalues', data=eigenvals)
        
        if orbital_details is not None:
            for key, val in orbital_details.items():
                hf.create_dataset(key, data=np.array([val]))

        if excitation_details is not None:
            grp = hf.create_group('excitation_details')
            for i, item in enumerate(excitation_details):
                grp.create_dataset(str(i), data=np.bytes_(item))

        hf.close()


def write_scf_results_to_hdf5(fname, scf_results, scf_history):
    """
    Writes SCF results to HDF5 file.

    :param fname:
        Name of the HDF5 file.
    :param scf_results:
        The dictionary containing SCF results.
    :param scf_history:
        The list containing SCF history.
    """

    valid_checkpoint = (fname and isinstance(fname, str) and
                        Path(fname).is_file())

    if valid_checkpoint:

        hf = h5py.File(fname, 'a')

        scf_group = hf.create_group('scf')

        # write SCF tensors
        keys = ['S'] + [
            f'{x}_{y}' for x in ['C', 'E', 'occ', 'D', 'F']
            for y in ['alpha', 'beta']
        ]
        for key in keys:
            # TODO: remove this if statement since all keys should be available
            if key in scf_results:
                scf_group.create_dataset(key, data=scf_results[key])

        # write dipole moment
        scf_group.create_dataset('dipole_moment',
                                 data=scf_results['dipole_moment'])

        # write SCF energy
        scf_group.create_dataset('scf_type',
                                 data=np.bytes_([scf_results['scf_type']]))
        scf_group.create_dataset('scf_energy',
                                 data=np.array([scf_results['scf_energy']]))

        # write SCF history
        keys = list(scf_history[0].keys())
        for key in keys:
            data = np.array([step[key] for step in scf_history])
            scf_group.create_dataset(f'scf_history_{key}', data=data)

        hf.close()


def write_rsp_solution(fname, key, vec, group_label='rsp'):
    """
    Writes a response solution vector to HDF5 file.

    :param fname:
        The name of the checkpoint file.
    :param key:
        The key for the solution vector.
    :param vec:
        The solution vector.
	:param group_label:
		The checkpoint file group label.
    """

    if fname and isinstance(fname, str):
        hf = h5py.File(fname, 'a')
        label = group_label + '/' + key
        if label in hf:
            del hf[label]
        hf.create_dataset(label, data=vec)
        hf.close()


def write_rsp_solution_with_multiple_keys(fname, keys, vec, group_label='rsp'):
    """
    Writes a response solution vector with multiple keys to HDF5 file.

    :param fname:
        The name of the checkpoint file.
    :param keys:
        The list of keys for the solution vector.
    :param vec:
        The solution vector.
	:param group_label:
		The checkpoint file group label.
    """

    if fname and isinstance(fname, str):
        hf = h5py.File(fname, 'a')

        label = group_label + '/' + keys[0]
        if label in hf:
            del hf[label]
        dset = hf.create_dataset(label, data=vec)

        for key in keys[1:]:
            label = group_label + '/' + key
            if label in hf:
                del hf[label]
            hf[label] = dset

        hf.close()


def write_lr_rsp_results_to_hdf5(fname, rsp_results, group_label='rsp'):
    """
    Writes the results of a linear response calculation to HDF5 file.

    :param fname:
        Name of the HDF5 file.
    :param rsp_results:
        The dictionary containing the linear response results.
	:param group_label:
		The checkpoint file group label.
    """

    if fname and isinstance(fname, str):

        hf = h5py.File(fname, 'a')

        for key in rsp_results:
            # Do not write the eigenvectors, file names and excitation details
            if "vector" in key or "cube" in key or "file" in key or "details" in key:
                continue

            label = group_label + '/' + key
            if label in hf:
                del hf[label]
            hf.create_dataset(label, data=rsp_results[key])

        hf.close()


def write_rsp_hdf5(fname, arrays, labels, molecule, basis, dft_dict, pe_dict,
                   ostream):
    """
    Writes response vectors to checkpoint file. Nuclear charges and basis
    set can also be written to the checkpoint file.

    :param fname:
        Name of the checkpoint file.
    :param arrays:
        The response vectors.
    :param labels:
        The list of labels for trial vectors and transformed vectors.
    :param molecule:
        The molecule.
    :param basis:
        The AO basis set.
    :param dft_dict:
        The dictionary containing DFT information.
    :param pe_dict:
        The dictionary containing PE information.
    :param ostream:
        The output stream.

    :return:
        True if checkpoint file is written. False if checkpoint file is not
        valid.
    """

    valid_checkpoint = (fname and isinstance(fname, str))

    if not valid_checkpoint:
        return False

    create_hdf5(fname, molecule, basis, dft_dict['dft_func_label'],
                pe_dict['potfile_text'])

    hf = h5py.File(fname, 'a')
    for label, array in zip(labels, arrays):
        hf.create_dataset(label, data=array)
    hf.close()

    checkpoint_text = 'Checkpoint written to file: '
    checkpoint_text += fname
    ostream.print_info(checkpoint_text)
    ostream.print_blank()

    return True


def write_detach_attach_to_hdf5(fname, state_label, dens_detach, dens_attach, group_label='rsp'):
    """
    Writes the detachment and attachment density matrices for a specific
    excited state to the checkpoint file.

    :param fname:
        The checkpoint file name.
    :param state_label:
        The excited state label.
    :param dens_detach:
        The detachment density matrix.
    :param dens_attach:
        The attachment density matrix.
	:param group_label:
		The checkpoint file group label.
    """

    if fname and isinstance(fname, str):

        hf = h5py.File(fname, 'a')

        # add detachment/attachment densities to the rsp group

        detach_label = group_label + "/detach_" + state_label
        if detach_label in hf:
            del hf[detach_label]
        hf.create_dataset(detach_label, data=dens_detach)

        attach_label = group_label + "/attach_" + state_label
        if attach_label in hf:
            del hf[attach_label]
        hf.create_dataset(attach_label, data=dens_attach)

        hf.close()


def read_rsp_hdf5(fname, labels, molecule, basis, dft_dict, pe_dict, ostream):
    """
    Reads response vectors from checkpoint file. Nuclear charges and basis
    set will be used to validate the checkpoint file.

    :param fname:
        Name of the checkpoint file.
    :param labels:
        The list of labels for trial vectors and transformed vectors.
    :param molecule:
        The molecule.
    :param basis:
        The AO basis set.
    :param dft_dict:
        The dictionary containing DFT information.
    :param pe_dict:
        The dictionary containing PE information.
    :param ostream:
        The output stream.

    :return:
        The trial vectors and transformed vectors.
    """

    valid_checkpoint = check_rsp_hdf5(fname, labels, molecule, basis, dft_dict,
                                      pe_dict)

    if not valid_checkpoint:
        return tuple([None for x in range(len(labels))])

    hf = h5py.File(fname, 'r')

    arrays = [None for x in range(len(labels))]

    for i in range(len(labels)):
        if labels[i] in hf.keys():
            arrays[i] = np.array(hf.get(labels[i]))

    hf.close()

    is_empty = [a is None for a in arrays]
    if True not in is_empty:
        checkpoint_text = 'Restarting from checkpoint file: '
        checkpoint_text += fname
        ostream.print_info(checkpoint_text)
        ostream.print_blank()

    return arrays


def check_rsp_hdf5(fname, labels, molecule, basis, dft_dict, pe_dict):
    """
    Checks validity of the checkpoint file. Nuclear charges and basis
    set will be used to validate the checkpoint file.

    :param fname:
        Name of the checkpoint file.
    :param labels:
        The list of labels for trial vectors and transformed vectors.
    :param molecule:
        The molecule.
    :param basis:
        The AO basis set.
    :param dft_dict:
        The dictionary containing DFT information.
    :param pe_dict:
        The dictionary containing PE information.

    :return:
        True if the checkpoint file is valid, False otherwise.
    """

    valid_checkpoint = (fname and isinstance(fname, str) and
                        Path(fname).is_file())

    if not valid_checkpoint:
        return False

    e_nuc = molecule.nuclear_repulsion_energy()
    nuclear_charges = molecule.get_element_ids()
    basis_set = basis.get_label()

    dft_func_label = dft_dict['dft_func_label']
    potfile_text = pe_dict['potfile_text']

    hf = h5py.File(fname, 'r')

    match_labels = True
    for label in labels:
        if label not in hf.keys():
            match_labels = False
            break

    match_nuclear_repulsion = False
    if 'nuclear_repulsion' in hf:
        h5_e_nuc = np.array(hf.get('nuclear_repulsion'))[0]
        if h5_e_nuc > 0.0 and e_nuc > 0.0:
            match_nuclear_repulsion = (abs(1.0 - h5_e_nuc / e_nuc) < 1.0e-13)
        else:
            match_nuclear_repulsion = (h5_e_nuc == e_nuc)

    match_nuclear_charges = False
    if 'nuclear_charges' in hf:
        h5_nuclear_charges = np.array(hf.get('nuclear_charges'))
        if h5_nuclear_charges.shape == nuclear_charges.shape:
            match_nuclear_charges = (
                h5_nuclear_charges == nuclear_charges).all()

    match_basis_set = False
    if 'basis_set' in hf:
        h5_basis_set = hf.get('basis_set')[0].decode('utf-8')
        match_basis_set = (h5_basis_set.upper() == basis_set.upper())

    match_dft_func = False
    if 'dft_func_label' in hf:
        h5_func_label = hf.get('dft_func_label')[0].decode('utf-8')
        match_dft_func = (h5_func_label.upper() == dft_func_label.upper())

    match_potfile = False
    if 'potfile_text' in hf:
        h5_potfile_text = hf.get('potfile_text')[0].decode('utf-8')
        match_potfile = (h5_potfile_text == potfile_text)

    hf.close()

    return (match_labels and match_nuclear_repulsion and
            match_nuclear_charges and match_basis_set and match_dft_func and
            match_potfile)


def write_distributed_focks(fname, dist_focks, key_freq_pairs, comm, ostream):
    """
    Writes distributed Fock matrices to checkpoint file.

    :param fname:
        Name of the checkpoint file.
    :param dist_focks:
        The distributed Fock matrices as DistributedArray.
    :param key_freq_pairs:
        The key-frequency pairs.
    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    :return:
        True if checkpoint file is written. False if checkpoint file is not
        valid.
    """

    valid_checkpoint = (fname and isinstance(fname, str))

    if not valid_checkpoint:
        return False

    t0 = tm.time()

    if comm.Get_rank() == mpi_master():
        hf = h5py.File(fname, 'w')
        str_key_freq_pairs = [str((key, w)) for key, w in key_freq_pairs]
        hf.create_dataset('key_freq_pairs', data=np.bytes_(str_key_freq_pairs))
        hf.close()
    comm.barrier()

    dist_focks.append_to_hdf5_file(fname, 'distributed_focks')

    checkpoint_text = 'Checkpoint written to file: '
    checkpoint_text += fname
    ostream.print_info(checkpoint_text)
    ostream.print_blank()

    checkpoint_text = 'Time spent in writing checkpoint file: '
    checkpoint_text += f'{(tm.time() - t0):.2f} sec'
    ostream.print_info(checkpoint_text)
    ostream.print_blank()

    return True


def read_distributed_focks(fname, comm, ostream):
    """
    Reads distributed Fock matrices from checkpoint file.

    :param fname:
        Name of the checkpoint file.
    :param key_freq_pairs:
        The key-frequency pairs.
    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    :return:
        A dictionary containing the distributed Fock matrices.
    """

    dist_focks = DistributedArray.read_from_hdf5_file(fname,
                                                      'distributed_focks', comm)

    checkpoint_text = 'Restarting from checkpoint file: '
    checkpoint_text += fname
    ostream.print_info(checkpoint_text)
    ostream.print_blank()

    return dist_focks


def check_distributed_focks(fname, key_freq_pairs):
    """
    Checks validity of the checkpoint file for distributed Fock matrices.

    :param fname:
        Name of the checkpoint file.
    :param key_freq_pairs:
        The key-frequency pairs.

    :return:
        True if the checkpoint file is valid, False otherwise.
    """

    valid_checkpoint = (fname and isinstance(fname, str) and
                        Path(fname).is_file())

    if not valid_checkpoint:
        return False

    str_key_freq_pairs = [str((key, w)) for key, w in key_freq_pairs]

    hf = h5py.File(fname, 'r')

    hf_key_freq_pairs = [
        x.decode('utf-8') for x in np.array(hf.get('key_freq_pairs'))
    ]

    valid_checkpoint = (str_key_freq_pairs == hf_key_freq_pairs and
                        'distributed_focks' in hf.keys())

    hf.close()

    return valid_checkpoint


def write_cpcm_charges(fname, cpcm_q):
    """
    Write C-PCM surface charges.

    :param fname:
        Name of the checkpoint file.
    :param cpcm_q:
        The C-PCM surface charges.
    """

    hf = h5py.File(fname, 'a')
    key = 'cpcm_q'
    if key in hf:
        del hf[key]
    hf.create_dataset('cpcm_q', data=cpcm_q)
    hf.close()


def read_cpcm_charges(fname):
    """
    Read C-PCM surface charges.

    :param fname:
        Name of the checkpoint file.

    :return:
        The C-PCM surface charges.
    """

    hf = h5py.File(fname, 'r')
    key = 'cpcm_q'
    if key in hf:
        cpcm_q = np.array(hf.get(key))
    else:
        cpcm_q = None
    hf.close()

    return cpcm_q

def read_results(fname, label):
    """ Read the results dictionary from a checkpoint file.

        :param fname:
            Name of the checkpoint file.
        :param label:
            The response dictionary label (scf, rsp, vib, opt, etc.).

        :return:
            the dictionary of results.
    """

    valid_checkpoint = (fname and isinstance(fname, str) and
                        Path(fname).is_file())

    assert_msg_critical(valid_checkpoint, fname + " is not a valid checkpoint file.") 

    res_dict = {}
    h5f = h5py.File(fname, "r")

    label_found = (label in h5f)

    if not label_found:
        h5f.close()

    assert_msg_critical(label_found, label + " section not found in the checkpoint file.")

    # Always read general information about the molecule, basis, and settings
    for key in h5f:
        if key not in ["vib", "rsp", "scf", "opt"]:
            data = np.array(h5f.get(key))
            # Check if data is a number or an array
            if len(data.shape) == 1 and data.shape[0] == 1:
                res_dict[key] = data[0]
            else:
                res_dict[key] = data

    h5f_dict = h5f[label]

    for key in h5f_dict:
        # For known lists of arrays (normal modes, raman and resonance raman),
        # read the sub-group accordingly
        if "normal_modes" in key or "raman_activities" in key:
            sub_dict = dict(h5f_dict[key])
            sub_dict_results = {}
            for sub_key in sub_dict:
                sub_dict_results[sub_key] = np.array(sub_dict[sub_key])
            res_dict[key] = sub_dict_results
        else:
            data = np.array(h5f_dict[key])
            # Check if data is a number or an array
            if len(data.shape) == 1 and data.shape[0] == 1:
                res_dict[key] = data[0]
            else:
                res_dict[key] = data

    if "opt" in label:
        # Create the list of xyz geometries
        nuclear_charges = np.array(res_dict["nuclear_charges"], dtype=int)
        xyz_geometries = []
        for coords in res_dict["opt_coordinates_au"]:
            molecule = Molecule(nuclear_charges, coords, units="au")
            xyz_geometries.append(molecule.get_xyz_string())
        res_dict["opt_geometries"] = xyz_geometries
    
    # Create molecule xyz
    nuclear_charges = np.array(res_dict["nuclear_charges"], dtype=int)
    coords = res_dict["atom_coordinates"]
    molecule = Molecule(nuclear_charges, coords, units="au")

    xyz_lines = molecule.get_xyz_string().splitlines()

    xyz = []
    for line in xyz_lines[2:]:
        xyz.append(line)

    res_dict["xyz"] = xyz

    if "vib" in label:
        # needed to animate the molecular vibrations
        res_dict["molecule_xyz_string"] = molecule.get_xyz_string()

    h5f.close()
    
    return res_dict
