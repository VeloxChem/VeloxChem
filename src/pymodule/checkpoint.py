#
#                           VELOXCHEM 1.0-RC3
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2022 by VeloxChem developers. All rights reserved.
#  Contact: https://veloxchem.org/contact
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
import numpy as np
import time as tm
import h5py

from .veloxchemlib import mpi_master
from .distributedarray import DistributedArray


def create_hdf5(fname, molecule, basis, dft_func_label, potfile_text):
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

        hf.create_dataset('nuclear_charges', data=molecule.elem_ids_to_numpy())

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

        hf.create_dataset('basis_set', data=np.string_([basis.get_label()]))

        hf.create_dataset('dft_func_label', data=np.string_([dft_func_label]))

        hf.create_dataset('potfile_text', data=np.string_([potfile_text]))

        hf.close()


def write_scf_tensors(fname, scf_tensors):
    """
    Writes SCF tensors to HDF5 file.

    :param fname:
        Name of the HDF5 file.
    :param scf_tensors:
        The dictionary of tensors from converged SCF wavefunction.
    """

    valid_checkpoint = (fname and isinstance(fname, str) and
                        Path(fname).is_file())

    if valid_checkpoint:
        hf = h5py.File(fname, 'a')
        keys = ['S'] + [f'{x}_{y}' for x in 'CEDF' for y in ['alpha', 'beta']]
        for key in keys:
            hf.create_dataset(key, data=scf_tensors[key])
        hf.close()


def write_rsp_solution(fname, key, vec):
    """
    Writes a response solution vector to HDF5 file.

    :param fname:
        The name of the checkpoint file.
    :param key:
        The key for the solution vector.
    :param vec:
        The solution vector.
    """

    valid_checkpoint = (fname and isinstance(fname, str) and
                        Path(fname).is_file())

    if valid_checkpoint:
        hf = h5py.File(fname, 'a')
        hf.create_dataset(key, data=vec)
        hf.close()


def write_rsp_solution_with_multiple_keys(fname, keys, vec):
    """
    Writes a response solution vector with multiple keys to HDF5 file.

    :param fname:
        The name of the checkpoint file.
    :param keys:
        The list of keys for the solution vector.
    :param vec:
        The solution vector.
    """

    valid_checkpoint = (fname and isinstance(fname, str) and
                        Path(fname).is_file())

    if valid_checkpoint:
        hf = h5py.File(fname, 'a')
        dset = hf.create_dataset(keys[0], data=vec)
        for key in keys[1:]:
            hf[key] = dset
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
        The list of labels for trial vecotrs and transformed vectors.
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


def read_rsp_hdf5(fname, labels, molecule, basis, dft_dict, pe_dict, ostream):
    """
    Reads response vectors from checkpoint file. Nuclear charges and basis
    set will be used to validate the checkpoint file.

    :param fname:
        Name of the checkpoint file.
    :param labels:
        The list of labels for trial vecotrs and transformed vectors.
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
        return tuple([None] * len(labels))

    hf = h5py.File(fname, 'r')

    arrays = [None] * len(labels)

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
        The list of labels for trial vecotrs and transformed vectors.
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
    nuclear_charges = molecule.elem_ids_to_numpy()
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
        hf.create_dataset('key_freq_pairs', data=np.string_(str_key_freq_pairs))
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
