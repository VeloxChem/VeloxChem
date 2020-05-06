from os.path import isfile
import numpy as np
import h5py


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

    e_nuc = molecule.nuclear_repulsion_energy()
    nuclear_charges = molecule.elem_ids_to_numpy()
    basis_set = basis.get_label()

    dft_func_label = dft_dict['dft_func_label']
    potfile_text = pe_dict['potfile_text']

    hf = h5py.File(fname, 'w')

    for label, array in zip(labels, arrays):
        hf.create_dataset(label, data=array, compression='gzip')

    hf.create_dataset('nuclear_repulsion',
                      data=np.array([e_nuc]),
                      compression='gzip')

    hf.create_dataset('nuclear_charges',
                      data=nuclear_charges,
                      compression='gzip')

    hf.create_dataset('basis_set',
                      data=np.string_([basis_set]),
                      compression='gzip')

    hf.create_dataset('dft_func_label',
                      data=np.string_([dft_func_label]),
                      compression='gzip')

    hf.create_dataset('potfile_text',
                      data=np.string_([potfile_text]),
                      compression='gzip')

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

    valid_checkpoint = (fname and isinstance(fname, str) and isfile(fname))

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
