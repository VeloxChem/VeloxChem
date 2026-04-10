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
import h5py

from .errorhandler import assert_msg_critical
from .molecule import Molecule
from .molecularbasis import MolecularBasis


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

        keys = ['S'] + [
            f'{x}_{y}' for x in ['C', 'E', 'occ', 'D', 'F']
            for y in ['alpha', 'beta']
        ]
        for key in keys:
            if key in scf_results:
                scf_group.create_dataset(key, data=scf_results[key])

        scf_group.create_dataset('dipole_moment',
                                 data=scf_results['dipole_moment'])

        scf_group.create_dataset('scf_type',
                                 data=np.bytes_([scf_results['scf_type']]))
        scf_group.create_dataset('scf_energy',
                                 data=np.array([scf_results['scf_energy']]))

        keys = list(scf_history[0].keys())
        for key in keys:
            data = np.array([step[key] for step in scf_history])
            scf_group.create_dataset(f'scf_history_{key}', data=data)

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
            if "vector" in key or "cube" in key or "file" in key or "details" in key:
                continue

            if key == "esa_results":
                continue

            label = group_label + '/' + key
            if label in hf:
                del hf[label]
            hf.create_dataset(label, data=rsp_results[key])

        hf.close()


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


def read_results(fname, label):
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
