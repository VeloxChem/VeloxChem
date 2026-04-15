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
import numbers
import numpy as np
import h5py

from .errorhandler import assert_msg_critical
from .molecule import Molecule
from .molecularbasis import MolecularBasis


def _write_value_to_hdf5(parent, key, value, value_label='HDF5 value'):
    """
    Writes a Python value to HDF5 with metadata for future round-tripping.

    Supported values are ``None``, strings, booleans, numbers, NumPy arrays,
    dictionaries, and list/tuple containers composed of supported values.

    :param parent:
        The parent HDF5 group.
    :param key:
        The dataset/group name.
    :param value:
        The value to write.
    :param value_label:
        A label used in error messages for the serialized value domain.
    """

    string_type = h5py.string_dtype(encoding='utf-8')

    if value is None:
        dset = parent.create_dataset(key, shape=(0,), dtype=string_type)
        dset.attrs['value_type'] = 'none'
        return

    if isinstance(value, str):
        dset = parent.create_dataset(key, data=value, dtype=string_type)
        dset.attrs['value_type'] = 'str'
        return

    if isinstance(value, (bytes, np.bytes_)):
        dset = parent.create_dataset(key,
                                     data=np.bytes_(value.decode('utf-8')
                                                    if isinstance(value, bytes)
                                                    else value.decode('utf-8')))
        dset.attrs['value_type'] = 'str'
        return

    if isinstance(value, dict):
        group = parent.create_group(key)
        group.attrs['value_type'] = 'dict'
        if all(isinstance(child_key, str) for child_key in value):
            group.attrs['dict_storage'] = 'named'
            for child_key, child_value in value.items():
                _write_value_to_hdf5(group, child_key, child_value, value_label)
        else:
            # Fall back to indexed dict entries when keys are not valid HDF5
            # object names, e.g. floats used for response frequencies.
            group.attrs['dict_storage'] = 'entries'
            group.attrs['length'] = len(value)
            for idx, (child_key, child_value) in enumerate(value.items()):
                entry_group = group.create_group(str(idx))
                entry_group.attrs['value_type'] = 'dict_entry'
                _write_value_to_hdf5(entry_group, '__key__', child_key,
                                     value_label)
                _write_value_to_hdf5(entry_group, '__value__', child_value,
                                     value_label)
        return

    if isinstance(value, (list, tuple)):
        group = parent.create_group(key)
        group.attrs['value_type'] = ('list'
                                     if isinstance(value, list) else 'tuple')
        group.attrs['length'] = len(value)
        for idx, item in enumerate(value):
            _write_value_to_hdf5(group, str(idx), item, value_label)
        return

    if isinstance(value, np.ndarray):
        if value.dtype.kind in ('U', 'S', 'O'):
            assert_msg_critical(
                value.dtype.kind != 'O' or all(isinstance(x, str)
                                               for x in value.flat),
                f'Unsupported object array for {value_label} key {key!r}.')
            dset = parent.create_dataset(key, data=value.astype(str),
                                         dtype=string_type)
            dset.attrs['value_type'] = 'ndarray'
            dset.attrs['array_data_type'] = 'str'
        else:
            dset = parent.create_dataset(key, data=value)
            dset.attrs['value_type'] = 'ndarray'
            dset.attrs['array_data_type'] = str(value.dtype)
        return

    if isinstance(value, (np.bool_, bool)):
        dset = parent.create_dataset(key, data=bool(value))
        dset.attrs['value_type'] = 'bool'
        return

    if isinstance(value, numbers.Integral):
        dset = parent.create_dataset(key, data=int(value))
        dset.attrs['value_type'] = 'int'
        return

    if isinstance(value, numbers.Real):
        dset = parent.create_dataset(key, data=float(value))
        dset.attrs['value_type'] = 'float'
        return

    if isinstance(value, numbers.Complex):
        dset = parent.create_dataset(key, data=complex(value))
        dset.attrs['value_type'] = 'complex'
        return

    errmsg = (f'Unsupported {value_label} type for key {key!r}: '
              f'{type(value)!r}')
    assert_msg_critical(False, errmsg)


def _read_value_from_hdf5(h5obj, value_label='HDF5 value'):
    """
    Reads a Python value from an HDF5 dataset/group using stored metadata.

    :param h5obj:
        The HDF5 dataset or group.
    :param value_label:
        A label used in error messages for the serialized value domain.

    :return:
        The reconstructed Python value.
    """

    value_type = h5obj.attrs.get('value_type', None)

    if isinstance(h5obj, h5py.Group):
        if value_type == 'dict':
            dict_storage = h5obj.attrs.get('dict_storage', 'named')
            if dict_storage == 'named':
                return {
                    key: _read_value_from_hdf5(h5obj[key], value_label)
                    for key in h5obj
                    # TODO: Drop this temporary skip once legacy raw nested
                    # groups such as rsp/nto are serialized with value_type
                    # metadata or stored outside the serialized results group.
                    if (not isinstance(h5obj[key], h5py.Group) or
                        h5obj[key].attrs.get('value_type', None) is not None)
                }

            if dict_storage == 'entries':
                length = int(h5obj.attrs['length'])
                values = {}
                for i in range(length):
                    entry_group = h5obj[str(i)]
                    key = _read_value_from_hdf5(entry_group['__key__'],
                                                value_label)
                    value = _read_value_from_hdf5(entry_group['__value__'],
                                                  value_label)
                    values[key] = value
                return values

            errmsg = (f'Unsupported {value_label} dict storage '
                      f'{dict_storage!r}.')
            assert_msg_critical(False, errmsg)

        if value_type in ('list', 'tuple'):
            length = int(h5obj.attrs['length'])
            values = [
                _read_value_from_hdf5(h5obj[str(i)], value_label)
                for i in range(length)
            ]
            return values if value_type == 'list' else tuple(values)

        errmsg = (f'Unsupported {value_label} group value_type '
                  f'{value_type!r}.')
        assert_msg_critical(False, errmsg)

    if value_type == 'none':
        return None

    if value_type == 'str':
        data = h5obj[()]
        if isinstance(data, bytes):
            return data.decode('utf-8')
        return data.astype(str) if isinstance(data, np.ndarray) else str(data)

    if value_type == 'bool':
        return bool(h5obj[()])

    if value_type == 'int':
        return int(h5obj[()])

    if value_type == 'float':
        return float(h5obj[()])

    if value_type == 'complex':
        return complex(h5obj[()])

    if value_type == 'ndarray':
        data = np.array(h5obj)
        if h5obj.attrs.get('array_data_type', None) == 'str':
            return data.astype(str)
        return data

    if value_type is None:
        data = np.array(h5obj)
        if len(data.shape) == 1 and data.shape[0] == 1:
            return data[0]
        return data

    errmsg = (f'Unsupported {value_label} dataset value_type '
              f'{value_type!r}.')
    assert_msg_critical(False, errmsg)


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

        e_nuc = molecule.effective_nuclear_repulsion_energy(basis)

        hf.create_dataset('nuclear_repulsion', data=np.array([e_nuc]))

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

        atom_basis_labels_flattened = []
        for atom_bs_name, elem in molecule.get_atom_basis_labels():
            atom_basis_labels_flattened += [atom_bs_name, elem]

        hf.create_dataset('atom_basis_labels_flattened',
                          data=np.bytes_(atom_basis_labels_flattened))

        hf.create_dataset('dft_func_label', data=np.bytes_([dft_func_label]))

        hf.create_dataset('potfile_text', data=np.bytes_([potfile_text]))

        hf.close()


def write_results_to_hdf5(fname,
                          label,
                          results,
                          value_label='result',
                          replace_group=True):
    """
    Writes a group-local results dictionary to an HDF5 file.

    :param fname:
        Name of the HDF5 file.
    :param label:
        The HDF5 group label.
    :param results:
        The dictionary containing result data.
    :param value_label:
        A label used in serialization error messages.
    :param replace_group:
        If ``True``, replace the whole target group before writing. If
        ``False``, keep existing group members that are not overwritten.
    """

    valid_checkpoint = (fname and isinstance(fname, str) and
                        Path(fname).is_file())

    if valid_checkpoint:

        with h5py.File(fname, 'a') as hf:
            if replace_group and label in hf:
                del hf[label]

            result_group = hf.require_group(label)
            result_group.attrs['value_type'] = 'dict'
            result_group.attrs['dict_storage'] = 'named'

            for key, value in results.items():
                if key in result_group:
                    del result_group[key]
                _write_value_to_hdf5(result_group,
                                     key,
                                     value,
                                     value_label=value_label)


def write_scf_results_to_hdf5(fname, scf_results):
    """
    Writes SCF results to HDF5 file.

    :param fname:
        Name of the HDF5 file.
    :param scf_results:
        The dictionary containing SCF results.
    """

    write_results_to_hdf5(fname,
                          'scf',
                          scf_results,
                          value_label='SCF result')


def write_scf_property_to_hdf5(fname, key, data):
    """
    Writes a single SCF property dataset to the scf group in HDF5 file.

    :param fname:
        Name of the HDF5 file.
    :param key:
        Dataset name inside the scf group.
    :param data:
        Dataset values.
    """

    valid_checkpoint = (fname and isinstance(fname, str) and
                        Path(fname).is_file())

    if valid_checkpoint:
        hf = h5py.File(fname, 'a')

        scf_group = hf.require_group('scf')
        label = f'scf/{key}'
        if label in hf:
            del hf[label]

        scf_group.create_dataset(key, data=data)

        hf.close()


def write_opt_results_to_hdf5(fname, opt_results):
    """
    Writes optimization results to HDF5 file.

    :param fname:
        Name of the HDF5 file.
    :param opt_results:
        The dictionary containing optimization results.
    """

    write_results_to_hdf5(fname,
                          'opt',
                          opt_results,
                          value_label='optimization result')


def write_rsp_solution(fname, key, vec, group_label='rsp'):
    """
    Writes a response solution vector to HDF5 file.

    :param fname:
        The name of the HDF5 file.
    :param key:
        The key for the solution vector.
    :param vec:
        The solution vector.
    :param group_label:
        The HDF5 group label.
    """

    if fname and isinstance(fname, str):
        with h5py.File(fname, 'a') as hf:
            label = group_label + '/' + key
            if label in hf:
                del hf[label]
            hf.create_dataset(label, data=vec)


def write_rsp_solution_with_multiple_keys(fname,
                                          keys,
                                          vec,
                                          group_label='rsp'):
    """
    Writes one response solution vector under multiple HDF5 keys.

    :param fname:
        The name of the HDF5 file.
    :param keys:
        The list of keys for the solution vector.
    :param vec:
        The solution vector.
    :param group_label:
        The HDF5 group label.
    """

    if fname and isinstance(fname, str):
        with h5py.File(fname, 'a') as hf:
            label = group_label + '/' + keys[0]
            if label in hf:
                del hf[label]
            dset = hf.create_dataset(label, data=vec)

            for key in keys[1:]:
                label = group_label + '/' + key
                if label in hf:
                    del hf[label]
                hf[label] = dset


def write_lr_rsp_results_to_hdf5(fname, rsp_results):
    """
    Writes the results of a linear response calculation to HDF5 file.

    :param fname:
        Name of the HDF5 file.
    :param rsp_results:
        The dictionary containing the linear response results.
    """

    write_results_to_hdf5(fname,
                          'rsp',
                          rsp_results,
                          value_label='response result',
                          replace_group=False)


def write_detach_attach_to_hdf5(fname,
                                state_label,
                                dens_detach,
                                dens_attach,
                                chg_detach=None,
                                chg_attach=None,
                                group_label='rsp'):
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
    :param chg_detach:
        The detachment charges.
    :param chg_attach:
        The attachment charges.
    :param group_label:
        The checkpoint file group label.
    """

    if fname and isinstance(fname, str):

        hf = h5py.File(fname, 'a')

        detach_label = group_label + "/detach_attach/detach_" + state_label
        if detach_label in hf:
            del hf[detach_label]
        hf.create_dataset(detach_label, data=dens_detach)

        attach_label = group_label + "/detach_attach/attach_" + state_label
        if attach_label in hf:
            del hf[attach_label]
        hf.create_dataset(attach_label, data=dens_attach)

        if chg_detach is not None:
            detach_label = (group_label + "/detach_attach/detach_charges_" +
                            state_label)
            if detach_label in hf:
                del hf[detach_label]
            hf.create_dataset(detach_label, data=chg_detach)

        if chg_attach is not None:
            attach_label = (group_label + "/detach_attach/attach_charges_" +
                            state_label)
            if attach_label in hf:
                del hf[attach_label]
            hf.create_dataset(attach_label, data=chg_attach)

        hf.close()


def read_results_old(fname, label):
    """ Read the results dictionary from an HDF5 results file.

        :param fname:
            Name of the HDF5 file.
        :param label:
            The response dictionary label (scf, rsp, vib, opt, etc.).

        :return:
            the dictionary of results.
    """

    valid_filename = (fname and isinstance(fname, str))
    assert_msg_critical(valid_filename, f"{fname!r} is not a valid filename.")

    file_exists = Path(fname).is_file()
    assert_msg_critical(file_exists, f"{fname!r} does not exist.")

    res_dict = {}
    h5f = h5py.File(fname, "r")

    label_found = (label in h5f)

    if not label_found:
        h5f.close()

    assert_msg_critical(label_found,
                        label + " section not found in the checkpoint file.")

    for key in h5f:
        if key not in ["vib", "rsp", "scf", "opt"]:
            data = np.array(h5f.get(key))
            if len(data.shape) == 1 and data.shape[0] == 1:
                res_dict[key] = data[0]
            else:
                res_dict[key] = data

    h5f_dict = h5f[label]

    known_keys_for_arrays = [
        'normal_modes',
        'vib_frequencies',
        'force_constants',
        'reduced_masses',
        'ir_intensities',
        'external_frequencies',
        'raman_activities',
    ]

    for key in h5f_dict:
        data = np.array(h5f_dict[key])
        if (len(data.shape) == 1 and data.shape[0] == 1 and
                key not in known_keys_for_arrays):
            res_dict[key] = data[0]
        else:
            res_dict[key] = data

    if "opt" in label:
        nuclear_charges = np.array(res_dict["nuclear_charges"]).astype(int)
        xyz_geometries = []
        for coords in res_dict["opt_coordinates_au"]:
            molecule = Molecule(nuclear_charges, coords, units="au")
            xyz_geometries.append(molecule.get_xyz_string())
        res_dict["opt_geometries"] = xyz_geometries

    molecule, basis = read_molecule_and_basis(fname)

    xyz_lines = molecule.get_xyz_string().splitlines()

    xyz = []
    for line in xyz_lines[2:]:
        xyz.append(line)

    res_dict["xyz"] = xyz

    if "vib" in label:
        res_dict["molecule_xyz_string"] = molecule.get_xyz_string()

    h5f.close()

    return res_dict


def read_results(fname, label):
    """ Read a results dictionary from a specific HDF5 group.

        This reader only reconstructs the content stored under ``label`` and
        does not merge unrelated top-level file metadata into the returned
        dictionary.

        :param fname:
            Name of the HDF5 file.
        :param label:
            The results dictionary label (scf, rsp, vib, opt, etc.).

        :return:
            the dictionary of results stored under the requested label.
    """

    valid_filename = (fname and isinstance(fname, str))
    assert_msg_critical(valid_filename, f"{fname!r} is not a valid filename.")

    file_exists = Path(fname).is_file()
    assert_msg_critical(file_exists, f"{fname!r} does not exist.")

    with h5py.File(fname, "r") as h5f:
        label_found = (label in h5f)
        assert_msg_critical(label_found,
                            label + " section not found in the checkpoint file.")

        return _read_value_from_hdf5(h5f[label], value_label='results')


def read_molecule_and_basis(fname):
    """
    Read molecule and AO basis set from an HDF5 results file.

    :param fname:
        Name of the HDF5 file.

    :return:
        The molecule and AO basis set.
    """

    valid_filename = (fname and isinstance(fname, str))
    assert_msg_critical(valid_filename, f"{fname!r} is not a valid filename.")

    file_exists = Path(fname).is_file()
    assert_msg_critical(file_exists, f"{fname!r} does not exist.")

    molecule = None
    basis = None

    with h5py.File(fname, "r") as h5f:
        nuclear_charges = np.array(h5f.get("nuclear_charges")).astype(int)
        coords = np.array(h5f.get("atom_coordinates"))

        atom_basis_labels = None
        if "atom_basis_labels_flattened" in h5f:
            atom_basis_labels_flattened = [
                x.decode("utf-8")
                for x in np.array(h5f.get("atom_basis_labels_flattened"))
            ]
            atom_basis_names = atom_basis_labels_flattened[0::2]
            atom_basis_elems = atom_basis_labels_flattened[1::2]
            atom_basis_labels = list(zip(atom_basis_names, atom_basis_elems))

        if atom_basis_labels is not None:
            molecule = Molecule(nuclear_charges, coords, "au",
                                atom_basis_labels)
        else:
            molecule = Molecule(nuclear_charges, coords, "au")

        if "molecular_charge" in h5f:
            molecule.set_charge(float(np.array(h5f.get("molecular_charge"))[0]))

        if "spin_multiplicity" in h5f:
            molecule.set_multiplicity(
                int(np.array(h5f.get("spin_multiplicity"))[0]))

        basis_set_label = np.array(h5f.get("basis_set"))[0].decode("utf-8")
        basis = MolecularBasis.read(molecule, basis_set_label, verbose=False)

    return molecule, basis
