# import veloxchem as vlx
from .molecule import Molecule
from .veloxchemlib import OverlapIntegralsDriver
from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import DenseMatrix
from .veloxchemlib import ao_matrix_to_veloxchem
from .veloxchemlib import ao_matrix_to_dalton
from .veloxchemlib import ElectricDipoleIntegralsDriver
from .veloxchemlib import xcfun
from .veloxchemlib import XCMolecularHessian
from .veloxchemlib import GridDriver
from .dftutils import get_default_grid_level
import numpy as np
import sys
import h5py

try:
    import pyscf
    # TODO: rename grad to pyscf_grad, hessian, etc.
    from pyscf import grad
    from pyscf import hessian
    from pyscf.grad import rks as pyscf_rks_grad
    # TODO: rename tdrks to pyscf_tdrks
    from pyscf.grad import tdrks
    from pyscf.hessian import rks as pyscf_rks_hessian
except ImportError:
    errtxt = "Please install pyscf.\n"
    errtxt += "$pip3 install pyscf"
    # raise ImportError(errtxt)


def translate_to_pyscf(label):
    """
    Translates the basis set label from veloxchem
    format to pyscf format
    """
    if label == "6-31+G_D,P_":
        return "6-31+G(d,p)"
    elif label == "6-311++G_D,P_":
        return "6-311++G(d,p)"
    elif label == "blyp":
        return "b88,lyp"
    else:
        return label

def get_pyscf_integral_type(int_type):
    """
    Translates the integral type to pyscf moleintor format.
    """
    # TODO: routine which returns all the available integral types.
    int_dict = {
            "overlap"                                       : "int1e_ovlp",
            "kinetic_energy"                                : "int1e_kin",
            "nuclear_attraction"                            : "int1e_nuc",
            "electric_dipole"                               : "int1e_r",
            "electric_quadrupole"                           : "int1e_rr",
            "electric_octupole"                             : "int1e_rrr",
            "electron_repulsion"                            : "int2e",
            "electron_repulsion_erf"                        : "int2e", #_coulerf
            "overlap_derivative"                            : "int1e_ipovlp",
            "kinetic_energy_derivative"                     : "int1e_ipkin",
            "nuclear_attraction_derivative_operator"        : "int1e_iprinv",
            "nuclear_attraction_derivative_orbitals"        : "int1e_ipnuc",
            # electric dipole deriv. calculated numerically
            "numerical_electric_dipole_derivative"          : "int1e",
            "electric_dipole_derivative"                    : "int1e_irp",
            "electron_repulsion_derivative"                 : "int2e_ip1",
            "overlap_second_derivative_2_0"                 : "int1e_ipipovlp",
            "overlap_second_derivative_1_1"                 : "int1e_ipovlpip",
            "kinetic_energy_second_derivative_2_0"          : "int1e_ipipkin",
            "kinetic_energy_second_derivative_1_1"          : "int1e_ipkinip",
            "nuclear_attraction_second_derivative_orbs_2_0" : "int1e_ipipnuc",
            "nuclear_attraction_second_derivative_orbs_1_1" : "int1e_ipnucip",
            "nuclear_attraction_second_derivative_op_2_0"   : "int1e_ipiprinv",
            "nuclear_attraction_second_derivative_op_1_1"   : "int1e_iprinvip",
            "electron_repulsion_second_derivative_2_0_0_0"  : "int2e_ipip1",
            "electron_repulsion_second_derivative_1_1_0_0"  : "int2e_ipvip1",
            "electron_repulsion_second_derivative_1_0_1_0"  : "int2e_ip1ip2",
        }
    if int_type not in int_dict.keys():
        error_text = "Unrecognized integral type: " + int_type +". "
        error_text += "Please use one of the following:"
        raise ValueError(error_text, int_dict.keys())
    else:
        return int_dict[int_type]

def get_sign(pyscf_int_type):
    """
    Returns the sign of the pyscf integrals which corresponds
    to the veloxchem standard. 
    """
    if pyscf_int_type in ["int1e_nuc", "int1e_ipovlp", "int1e_ipkin",
                          "int1e_ipnuc", "int1e_iprinv", "int2e_ip1",
                          "int1e_irp"]:
        return -1
    else:
        return 1

def get_molecule_string(molecule):
    mol_string = ""
    for i in range(molecule.number_of_atoms()):
        mol_string += molecule.get_labels()[i] + " "
        for j in range(3):
            mol_string += str(molecule.get_coordinates()[i][j]) + " "

        mol_string += "\n"

    return mol_string

def write_2d_array_hdf5(fname, arrays, labels=[], atom_index=None):
    """
    Writes the one-electron integral derivatives to the checkpoint file. 

    :param fname:
        Name of the checkpoint file. 
        The file must be created before calling this routine.
    :param arrays:
        The one-electron integral derivatives (with 3 components: x, y, z)
    :param labels:
        The list of labels (x, y, and z)
    :param atom_index:
        The index of the atom with respect to which the derivatives
        are calculated.

    :return:
        True if checkpoint file is written. False if checkpoint file is not
        valid.
    """

    valid_checkpoint = (fname and isinstance(fname, str))

    if not valid_checkpoint:
        return False

    try:
        hf = h5py.File(fname, 'a')
    except:
        return False

    for label, array in zip(labels, arrays):
        if atom_index is not None:
            full_label= str(atom_index) + label
        else:
            full_label = label
        hf.create_dataset(full_label, data=array, compression='gzip')
    hf.close()
    return True

def import_integral(molecule, basis, int_type, atom1, shell1,
                    atom2, shell2, atom3=None, shell3=None,
                    atom4=None, shell4=None, xi1=None, xi2=None,
                    chk_file=None, return_block=True, full_deriv=False,
                    delta_h=1e-5, omega=None):
    """
    Imports integrals and integral derivatives from pyscf and converts 
    them to veloxchem format.
    Specific atoms and shells can be selected.

    :param molecule:
        the vlx molecule object
    :param basis:
        the vlx basis object
    :param int_type:
        the type of integral, or integral derivative 
    :param atom1:
        index of the first atom of interest
    :param shell1:
        list of atomic shells of interest for atom 1
    :param atom2:
        index of the second atom of interest
    :param shell2:
        list of atomic shells of interest for atom 2
    :param atom3:
        index of the third atom of interest (for 2e integrals)
    :param shell3:
        list of atomic shells of interest for atom 3 (for 2e integrals)
    :param atom4:
        index of the fourth atom of interest (for 2e integrals)
    :param shell4:
        list of atomic shells of interest for atom 4 (for 2e integrals)
    :param xi1:
        index of the atom with respect to which derivatives are taken
        (for first and second order derivatives)
    :param xi2:
        index of the second atom with respect to which derivatives are taken
        (for second order derivatives)
    :param chk_file:
        the hdf5 checkpoint file name
    :param return_block:
        return the matrix block, or the full matrix
    :param full_deriv:
        return the full integral derivative (for first and second order
        derivatives)
    :param delta_h:
        step for numerical derivative (required by the electric dipole
        derivatives).


    :return:
        a numpy array corresponding to a specified block of 
        the selected integral or integral derivative matrix.
    """

    pyscf_int_type = get_pyscf_integral_type(int_type)

    if 'int2e' in pyscf_int_type:
        if (atom3 is None) or (shell3 is None) or (atom4 is None) or (shell4 is None):
            errtxt = "Not enough shells defined for a two-electron integral. "
            errtxt += "Please define atom3, shell3, atom4, shell4."
            raise ValueError(errtxt, atom3, shell3, atom4, shell4)

    if 'derivative' in int_type and xi1 is None:
        errtxt = "Please set the index of the atom (xi1) with respect to which"
        errtxt += " the derivative should be taken. "
        raise ValueError(errtxt, xi1)

    if 'second' in int_type:
        if '2_0' not in int_type and xi2 is None:
            errtxt = "Please set the index of the second atom (xi2)"
            errtxt += " with respect to which the derivative should be taken."
            raise ValueError(errtxt, xi2)

    if 'derivative' in int_type:
        if 'second' in int_type:
            if '2_0' in int_type:
                if 'int2e' in pyscf_int_type:
                    return import_2e_integral_derivative(molecule, basis,
                        int_type, atomi=xi1, atom1=atom1, shell1=shell1,
                        atom2=atom2, shell2=shell2, atom3=atom3, shell3=shell3,
                        atom4=atom4, shell4=shell4, chk_file=chk_file,
                        return_block=return_block, full_deriv=full_deriv)
                else:
                    return import_1e_integral_derivative(molecule, basis,
                        int_type, atomi=xi1, atom1=atom1, shell1=shell1,
                        atom2=atom2, shell2=shell2, chk_file=chk_file,
                        return_block=return_block, full_deriv=full_deriv)
            else:
                if 'int2e' in pyscf_int_type:
                    return import_2e_second_order_integral_derivative(molecule,
                        basis, int_type, atomi=xi1, atomj=xi2, atom1=atom1,
                        shell1=shell1, atom2=atom2, shell2=shell2, atom3=atom3,
                        shell3=shell3, atom4=atom4, shell4=shell4,
                        chk_file=chk_file, return_block=return_block,
                        full_deriv=full_deriv)
                else:
                    return import_1e_second_order_integral_derivative(molecule,
                        basis, int_type, atomi=xi1, atomj=xi2, atom1=atom1,
                        shell1=shell1, atom2=atom2, shell2=shell2,
                        chk_file=chk_file, return_block=return_block,
                        full_deriv=full_deriv)

        else:
            if 'int2e' in pyscf_int_type:
                return import_2e_integral_derivative(molecule, basis, int_type,
                        atomi=xi1, atom1=atom1, shell1=shell1, atom2=atom2,
                        shell2=shell2, atom3=atom3, shell3=shell3, atom4=atom4,
                        shell4=shell4, chk_file=chk_file,
                        return_block=return_block, full_deriv=full_deriv)
            else:
                if int_type == "numerical_electric_dipole_derivative":
                    return numerical_electric_dipole_derivatives(molecule,
                        basis, int_type, atomi=xi1, atom1=atom1, shell1=shell1,
                        atom2=atom2, shell2=shell2, delta_h=delta_h,
                        chk_file=chk_file, return_block=return_block) 
                else:
                    return import_1e_integral_derivative(molecule, basis,
                        int_type, 
                        atomi=xi1, atom1=atom1, shell1=shell1,
                        atom2=atom2, shell2=shell2, chk_file=chk_file,
                        return_block=return_block, full_deriv=full_deriv)
    else:
        if 'int2e' in pyscf_int_type:
            return import_2e_integral(molecule, basis, int_type, atom1=atom1,
                       shell1=shell1, atom2=atom2, shell2=shell2, atom3=atom3,
                       shell3=shell3, atom4=atom4, shell4=shell4,
                       chk_file=chk_file, return_block=return_block, omega=omega)
        else:
            return import_1e_integral(molecule, basis, int_type, atom1=atom1,
                       shell1=shell1, atom2=atom2, shell2=shell2,
                       chk_file=chk_file, return_block=return_block)


def import_1e_integral(molecule, basis, int_type, atom1=1, shell1=None,
                       atom2=1, shell2=None, chk_file=None,
                       unit="au", return_block=True):
    """
    Imports one electron integrals from pyscf and converts to veloxchem format.
    Specific atoms and shells can be selected.

    :param molecule:
        the vlx molecule object
    :param basis:
        the vlx basis object
    :param int_type:
        the type of one-electron integral: overlap, kinetic_energy,
                                           nuclear_attraction 
    :param atom1:
        index of the first atom of interest
    :param shell1:
        list of atomic shells of interest for atom 1
    :param atom2:
        index of the second atom of interest
    :param shell2:
        list of atomic shells of interest for atom 2
    :param chk_file:
        the hdf5 checkpoint file name
    :param unit:
        the units to be used for the molecular geometry;
        possible values: "au" (default), "Angstrom"
    :param return_block:
        return the matrix block, or the full matrix.


    :return:
        a numpy array corresponding to a specified block of 
        the selected 1e integral matrix.
    """

    molecule_string = get_molecule_string(molecule)
    basis_set_label = basis.get_label()
    pyscf_basis = translate_to_pyscf(basis_set_label)
    pyscf_molecule = pyscf.gto.M(atom=molecule_string,
                                 basis=pyscf_basis, unit=unit,
                                 charge=molecule.get_charge())

    pyscf_int_type = get_pyscf_integral_type(int_type)
    sign = get_sign(pyscf_int_type)

    pyscf_int = sign * pyscf_molecule.intor(pyscf_int_type, aosym='s1')

    # Transform integral to veloxchem format
    if pyscf_int_type in ["int1e_r", "int1e_rr", "int1e_rrr"]:
        vlx_int = np.zeros_like(pyscf_int)
        dof = vlx_int.shape[0]
        for x in range(dof):
            vlx_int[x] = ao_matrix_to_veloxchem(
                                 DenseMatrix(pyscf_int[x]),
                                 basis, molecule).to_numpy()
    else:
        vlx_int = ao_matrix_to_veloxchem(
                                 DenseMatrix(pyscf_int),
                                 basis, molecule).to_numpy()

    ao_basis_map = basis.get_ao_basis_map(molecule)
    rows = []
    columns = []
    k = 0
    for ao in ao_basis_map:
        parts = ao.split()
        atom = int(parts[0])
        shell = parts[2]
        if atom1 == atom:
            if shell1 is not None:
                for s in shell1:
                    if s in shell:
                        rows.append(k)
            else:
                rows.append(k)
        if atom2 == atom:
            if shell2 is not None:
                for s in shell2:
                    if s in shell:
                        columns.append(k)
            else:
                columns.append(k)
        k += 1
    if rows == []:
        raise ValueError("Atom or shell(s) not found.", atom1, shell1)
    if columns == []:
        raise ValueError("Atom or shell(s) not found.", atom2, shell2)
    

    if pyscf_int_type in ["int1e_r", "int1e_rr", "int1e_rrr"]:
        vlx_int_block = np.zeros((dof, len(rows), len(columns)))
    else:
        vlx_int_block = np.zeros((len(rows), len(columns)))
    
    for i in range(len(rows)):
        for j in range(len(columns)):
            if pyscf_int_type in ["int1e_r", "int1e_rr", "int1e_rrr"]:
                vlx_int_block[:,i,j] = vlx_int[:, rows[i], columns[j]] 
            else:
                vlx_int_block[i,j] = vlx_int[rows[i], columns[j]]
   
    label = int_type+'_atom%d' % (atom1)
    if shell1 is not None:
        for s1 in shell1:
            label += "_%s" % (s1)
    label += "_atom%d" % (atom2)
    if shell2 is not None:
        for s2 in shell2:
            label += "_%s" % (s2)

    print(label) 
    if chk_file is not None:
        write_2d_array_hdf5(chk_file, [vlx_int_block], labels=[label])

    if return_block:
        return vlx_int_block
    else:
        return vlx_int


def import_2e_integral(molecule, basis, int_type, atom1=1, shell1=None,
                       atom2=1, shell2=None, atom3=1, shell3=None,
                       atom4=1, shell4=None, chk_file=None,
                       unit="au", return_block=True, omega=None):
    """
    Imports two electron integrals from pyscf and converts to veloxchem format.
    Specific atoms and shells can be selected.

    :param molecule:
        the vlx molecule object
    :param basis:
        the vlx basis object
    :param int_type:
        the type of one-electron integral: overlap, kinetic_energy,
                                           nuclear_attraction 
    :param atom1:
        index of the first atom of interest
    :param shell1:
        list of atomic shells of interest for atom 1
    :param atom2:
        index of the second atom of interest
    :param shell2:
        list of atomic shells of interest for atom 2
    :param atom3:
        index of the third atom of interest
    :param shell3:
        list of atomic shells of interest for atom 3
    :param atom4:
        index of the fourth atom of interest
    :param shell2:
        list of atomic shells of interest for atom 4
    :param chk_file:
        the hdf5 checkpoint file name
    :param unit:
        the units to be used for the molecular geometry;
        possible values: "au" (default), "Angstrom"
    :param return_block:
        return the matrix block, or the full matrix.
    :param omega:
        Error function parameter for range-separated Coulomb integrals.


    :return:
        a numpy array corresponding to a specified block of
        the selected two electron integral
    """

    molecule_string = get_molecule_string(molecule)
    basis_set_label = basis.get_label()
    pyscf_basis = translate_to_pyscf(basis_set_label)
    pyscf_molecule = pyscf.gto.M(atom=molecule_string,
                                 basis=pyscf_basis, unit=unit,
                                 charge=molecule.get_charge())

    pyscf_int_type = get_pyscf_integral_type(int_type)
    sign = get_sign(pyscf_int_type)
    if 'erf' in int_type:
        with pyscf_molecule.with_range_coulomb(omega):
            pyscf_int = sign * pyscf_molecule.intor(pyscf_int_type, aosym='s1')
    else:
        pyscf_int = sign * pyscf_molecule.intor(pyscf_int_type, aosym='s1')

    nao = pyscf_molecule.nao

    # Transform integral to veloxchem format
    vlx_int = np.zeros_like(pyscf_int)

    basis_set_map = basis.get_index_map(molecule)
    for m in range(nao):
        for n in range(nao):
            for t in range(nao):
                for p in range(nao):
                    vm = basis_set_map[m]
                    vn = basis_set_map[n]
                    vt = basis_set_map[t]
                    vp = basis_set_map[p]
                    vlx_int[vm,vn,vt,vp] = pyscf_int[m,n,t,p]

    ao_basis_map = basis.get_ao_basis_map(molecule)
    index1 = []
    index2 = []
    index3 = []
    index4 = []
    k = 0
    for ao in ao_basis_map:
        parts = ao.split()
        atom = int(parts[0])
        shell = parts[2]
        if atom1 == atom:
            if shell1 is not None:
                for s in shell1:
                    if s in shell:
                        index1.append(k)
            else:
                index1.append(k)
        if atom2 == atom:
            if shell2 is not None:
                for s in shell2:
                    if s in shell:
                        index2.append(k)
            else:
                index2.append(k)
        if atom3 == atom:
            if shell3 is not None:
                for s in shell3:
                    if s in shell:
                        index3.append(k)
            else:
                index3.append(k)
        if atom4 == atom:
            if shell4 is not None:
                for s in shell4:
                    if s in shell:
                        index4.append(k)
            else:
                index4.append(k)
        k += 1
    if index1 == []:
        raise ValueError("Atom or shell(s) not found.", atom1, shell1)
    if index2 == []:
        raise ValueError("Atom or shell(s) not found.", atom2, shell2)
    if index3 == []:
        raise ValueError("Atom or shell(s) not found.", atom3, shell3)
    if index4 == []:
        raise ValueError("Atom or shell(s) not found.", atom4, shell4)
    
    vlx_int_block = np.zeros((len(index1), len(index2),
                              len(index3), len(index4)))
    
    for i in range(len(index1)):
        for j in range(len(index2)):
            for k in range(len(index3)):
                for l in range(len(index4)):
                    vlx_int_block[i,j,k,l] = vlx_int[index1[i], index2[j],
                                                     index3[k], index4[l]]
   
    label = int_type+'_atom%d' % (atom1)
    if shell1 is not None:
        for s1 in shell1:
            label += "_%s" % (s1)

    label += "_atom%d" % (atom2)
    if shell2 is not None:
        for s2 in shell2:
            label += "_%s" % (s2)

    label += "_atom%d" % (atom3)
    if shell3 is not None:
        for s3 in shell3:
            label += "_%s" % (s3)

    label += "_atom%d" % (atom4)
    if shell4 is not None:
        for s4 in shell4:
            label += "_%s" % (s4)

    print(label) 
    if chk_file is not None:
        write_2d_array_hdf5(chk_file, [vlx_int_block], labels=[label])

    if return_block:
        return vlx_int_block
    else:
        return vlx_int

def import_1e_integral_derivative(molecule, basis, int_type, atomi=1, atom1=1,
                                  shell1=None, atom2=1, shell2=None,
                                  chk_file=None, unit="au", return_block=True,
                                  full_deriv=False):
    """
    Imports one electron integral derivatives from pyscf and converts
    to veloxchem format. Specific atoms and shells can be selected.

    :param molecule:
        the vlx molecule object
    :param basis:
        the vlx basis object
    :param int_type:
        the type of one-electron integral: overlap, kinetic_energy,
                                           nuclear_attraction 
    :param atomi:
        the index of the atom with respect to which the derivatives
        are taken (indexing starts at 1)
    :param atom1:
        index of the first atom of interest
    :param shell1:
        list of atomic shells of interest for atom 1
    :param atom2:
        index of the second atom of interest
    :param shell2:
        list of atomic shells of interest for atom 2
    :param chk_file:
        the hdf5 checkpoint file name
    :param unit:
        the units to be used for the molecular geometry;
        possible values: "au" (default), "Angstrom"
    :param return_block:
        return the matrix block, or the full matrix
    :param full_deriv:
        if to calculate the full derivative
        (nabla m | operator | n) + (m | operator | nabla n) (for debugging).


    :return:
        a numpy array corresponding to a specified block of 
        the selected 1e integral derivative matrix.
    """

    molecule_string = get_molecule_string(molecule)
    basis_set_label = basis.get_label()
    pyscf_basis = translate_to_pyscf(basis_set_label)
    pyscf_molecule = pyscf.gto.M(atom=molecule_string,
                                 basis=pyscf_basis, unit=unit,
                                 charge=molecule.get_charge())

    # Get the AO indeces corresponding to atom i
    i = atomi - 1 # to use pyscf indexing

    pyscf_int_type = get_pyscf_integral_type(int_type)
    sign = get_sign(pyscf_int_type)
    if pyscf_int_type in ["int1e_iprinv", "int1e_ipiprinv", "int1e_iprinvip"]:
        # TODO: check if Zilvinas wants it this way
        sign *= molecule.elem_ids_to_numpy()[i] # Z-number
        # (m | nabla_i operator | n)
        with pyscf_molecule.with_rinv_at_nucleus(i): 
            pyscf_int_deriv = sign * pyscf_molecule.intor(pyscf_int_type,
                                                          aosym='s1')
    else:
        pyscf_int_deriv = sign * pyscf_molecule.intor(pyscf_int_type,
                                                      aosym='s1')


    #print("Integral type and sign", sign, pyscf_int_type)
    
    ao_slices = pyscf_molecule.aoslice_by_atom()

    ki, kf = ao_slices[i, 2:]

    pyscf_int_deriv_atom_i = np.zeros(pyscf_int_deriv.shape)

    if pyscf_int_type in ["int1e_iprinv"]:
        # (m | nabla_i operator | n)
        pyscf_int_deriv_atom_i = pyscf_int_deriv
    elif pyscf_int_type in ["int1e_iprinvip"]:
        pyscf_int_deriv_atom_i = -pyscf_int_deriv
        pyscf_int_deriv_atom_i[:,ki:kf] += pyscf_int_deriv[:,ki:kf]
        pyscf_int_deriv_atom_i[:,:,ki:kf] += pyscf_int_deriv[:,:,ki:kf]
    elif pyscf_int_type in ["int1e_ipiprinv"]:
        pyscf_int_deriv_atom_i = -pyscf_int_deriv
        pyscf_int_deriv_atom_i[:,ki:kf] += pyscf_int_deriv[:,ki:kf]
        pyscf_int_deriv_atom_i[:,:,ki:kf] += pyscf_int_deriv[:,ki:kf].transpose(0,2,1)
    elif pyscf_int_type in ["int1e_irp"]:
        pyscf_int_deriv_atom_i[:,:,ki:kf] = pyscf_int_deriv[:,:,ki:kf]
    else:
        pyscf_int_deriv_atom_i[:,ki:kf] = pyscf_int_deriv[:,ki:kf]


    # (nabla m | operator | n) + (m | operator | nabla n)
    # TODO: remove
    if full_deriv:
        pyscf_int_deriv_atom_i += pyscf_int_deriv_atom_i.transpose(0,2,1)

    # Transform integral to veloxchem format
    vlx_int_deriv = np.zeros_like(pyscf_int_deriv_atom_i)

    # number of derivatives (3)
    n_deriv = pyscf_int_deriv_atom_i.shape[0]
    for x in range(n_deriv):
        vlx_int_deriv[x] = ao_matrix_to_veloxchem(
                                 DenseMatrix(pyscf_int_deriv_atom_i[x]),
                                 basis, molecule).to_numpy()

    ao_basis_map = basis.get_ao_basis_map(molecule)
    rows = []
    columns = []
    k = 0
    for ao in ao_basis_map:
        parts = ao.split()
        atom = int(parts[0])
        shell = parts[2]
        if atom1 == atom:
            if shell1 is not None:
                for s in shell1:
                    if s in shell:
                        rows.append(k)
            else:
                rows.append(k)
        if atom2 == atom:
            if shell2 is not None:
                for s in shell2:
                    if s in shell:
                        columns.append(k)
            else:
                columns.append(k)
        k += 1
    if rows == []:
        raise ValueError("Atom or shell(s) not found.", atom1, shell1)
    if columns == []:
        raise ValueError("Atom or shell(s) not found.", atom2, shell2)
    

    vlx_int_block = np.zeros((vlx_int_deriv.shape[0], len(rows), len(columns)))

    for k in range(len(rows)):
        for l in range(len(columns)):
            vlx_int_block[:, k, l] = vlx_int_deriv[:, rows[k], columns[l]]
   
    label = int_type+'_wrt_atom%d_shells_atom%d' % (atomi, atom1)
    if shell1 is not None:
        for s1 in shell1:
            label += "_%s" % (s1)
    label += "_atom%d" % (atom2)
    if shell2 is not None:
        for s2 in shell2:
            label += "_%s" % (s2)

    if chk_file is not None:
        write_2d_array_hdf5(chk_file, [vlx_int_block], labels=[label])

    if return_block:
        return vlx_int_block
    else:
        return vlx_int_deriv

def numerical_electric_dipole_derivatives(molecule, ao_basis, int_type,
                                          atomi=1, atom1=1, shell1=None,
                                          atom2=1, shell2=None, delta_h=1e-5,
                                          chk_file=None, return_block=True,
                                          unit="au"):
    """
    Computes numerical derivatives of the electric dipole integrals.

    :param molecule:
        The molecule.
    :param ao_basis:
        The AO basis set.
    :param int_type:
        The integral type (electric_dipole_derivative)
    :param atomi:
        The atom with respect to which the derivatives are taken
    :param atom1:
        index of the first atom of interest
    :param shell1:
        list of atomic shells of interest for atom 1
    :param atom2:
        index of the second atom of interest
    :param shell2:
        list of atomic shells of interest for atom 2
    :param delta_h:
        The step for the numerical derivative.
    :param chk_file:
        the hdf5 checkpoint file name
    :param return_block:
        return the matrix block, or the full matrix
     :param unit:
        the units to be used for the molecular geometry;
        possible values: "au" (default), "Angstrom"

    :return:
        The selected block of the electric dipole integral derivatives.
    """

    i = atomi - 1
    # atom labels
    labels = molecule.get_labels()

    # number of atoms
    natm = molecule.number_of_atoms()

    # atom coordinates (nx3)
    coords = molecule.get_coordinates()

    # number of atomic orbitals
    molecule_string = get_molecule_string(molecule)
    basis_set_label = ao_basis.get_label()
    pyscf_basis = translate_to_pyscf(basis_set_label)
    pyscf_molecule = pyscf.gto.M(atom=molecule_string,
                                 basis=pyscf_basis, unit=unit,
                                 charge=molecule.get_charge())
    nao = pyscf_molecule.nao

    # Dipole integrals driver
    dipole_drv = ElectricDipoleIntegralsDriver()

    # 3 dipole components x 3 atomic coordinates
    # x No. basis x No. basis
    dipole_integrals_gradient = np.zeros((3, 3, nao, nao))

    for d in range(3):
        coords[i, d] += delta_h
        new_mol = Molecule(labels, coords, units='au')

        dipole_mats_p = dipole_drv.compute(new_mol, ao_basis)
        dipole_ints_p = (dipole_mats_p.x_to_numpy(), 
                         dipole_mats_p.y_to_numpy(),
                         dipole_mats_p.z_to_numpy())

        coords[i, d] -= 2.0 * delta_h
        new_mol = Molecule(labels, coords, units='au')

        dipole_mats_m = dipole_drv.compute(new_mol, ao_basis)
        dipole_ints_m = (dipole_mats_m.x_to_numpy(),
                         dipole_mats_m.y_to_numpy(),
                         dipole_mats_m.z_to_numpy())

        for c in range(3):
            dipole_integrals_gradient[c, d] = (
                ( dipole_ints_p[c] - dipole_ints_m[c] ) / (2.0 * delta_h)
            )

    vlx_int_deriv = dipole_integrals_gradient.reshape(3*3,nao,nao)

    ao_basis_map = ao_basis.get_ao_basis_map(molecule)
    rows = []
    columns = []
    k = 0
    for ao in ao_basis_map:
        parts = ao.split()
        atom = int(parts[0])
        shell = parts[2]
        if atom1 == atom:
            if shell1 is not None:
                for s in shell1:
                    if s in shell:
                        rows.append(k)
            else:
                rows.append(k)
        if atom2 == atom:
            if shell2 is not None:
                for s in shell2:
                    if s in shell:
                        columns.append(k)
            else:
                columns.append(k)
        k += 1
    if rows == []:
        raise ValueError("Atom or shell(s) not found.", atom1, shell1)
    if columns == []:
        raise ValueError("Atom or shell(s) not found.", atom2, shell2)
    

    vlx_int_block = np.zeros((vlx_int_deriv.shape[0], len(rows), len(columns)))

    for k in range(len(rows)):
        for l in range(len(columns)):
            vlx_int_block[:, k, l] = vlx_int_deriv[:, rows[k], columns[l]]
   
    label = int_type+'_wrt_atom%d_shells_atom%d' % (atomi, atom1)
    if shell1 is not None:
        for s1 in shell1:
            label += "_%s" % (s1)
    label += "_atom%d" % (atom2)
    if shell2 is not None:
        for s2 in shell2:
            label += "_%s" % (s2)

    if chk_file is not None:
        write_2d_array_hdf5(chk_file, [vlx_int_block], labels=[label])

    if return_block:
        return vlx_int_block
    else:
        return vlx_int_deriv


def import_1e_second_order_integral_derivative(molecule, basis, int_type,
                                  atomi, atomj, atom1=1, shell1=None,
                                  atom2=1, shell2=None, chk_file=None,
                                  unit="au", return_block=True,
                                  full_deriv=False):
    """
    Imports one electron integral derivatives from pyscf and converts
    to veloxchem format. Specific atoms and shells can be selected.

    :param molecule:
        the vlx molecule object
    :param basis:
        the vlx basis object
    :param int_type:
        the type of one-electron integral: overlap, kinetic_energy,
                                           nuclear_attraction 
    :param atomi:
        the index of the atom with respect to which the derivatives
        are taken (indexing starts at 1)
    :param atomj:
        the index of the second atom with respect to which the derivatives
        are taken (indexing starts at 1)
    :param atom1:
        index of the first atom of interest
    :param shell1:
        list of atomic shells of interest for atom 1
    :param atom2:
        index of the second atom of interest
    :param shell2:
        list of atomic shells of interest for atom 2
    :param chk_file:
        the hdf5 checkpoint file name
    :param unit:
        the units to be used for the molecular geometry;
        possible values: "au" (default), "Angstrom"
    :param return_block:
        return the matrix block, or the full matrix
    :param full_deriv:
        if to calculate the full derivative
        (nabla m | operator | n) + (m | operator | nabla n) (for debugging).


    :return:
        a numpy array corresponding to a specified block of 
        the selected 1e integral derivative matrix.
    """

    molecule_string = get_molecule_string(molecule)
    basis_set_label = basis.get_label()
    pyscf_basis = translate_to_pyscf(basis_set_label)
    pyscf_molecule = pyscf.gto.M(atom=molecule_string,
                                 basis=pyscf_basis, unit=unit,
                                 charge=molecule.get_charge())

    # Get the AO indeces corresponding to atom i
    i = atomi - 1 # to use pyscf indexing
    j = atomj - 1

    pyscf_int_type = get_pyscf_integral_type(int_type)
    sign = get_sign(pyscf_int_type)

    # TODO: modify to include all rinv-related integral derivatives
    if pyscf_int_type in ["int1e_iprinvip", "int1e_ipiprinv"]:
        # TODO: check if Zilvinas wants it this way
        sign *= molecule.elem_ids_to_numpy()[j] # Z-number
        # (nabla_i m | nabla_j operator | n)
        with pyscf_molecule.with_rinv_at_nucleus(j): 
            pyscf_int_deriv = sign * pyscf_molecule.intor(pyscf_int_type,
                                                          aosym='s1')
    else:
        pyscf_int_deriv = sign * pyscf_molecule.intor(pyscf_int_type,
                                                      aosym='s1')


    #print("Integral type and sign", sign, pyscf_int_type)
    
    ao_slices = pyscf_molecule.aoslice_by_atom()

    ki, kf = ao_slices[i, 2:]
    kj, kg = ao_slices[j, 2:]

    pyscf_int_deriv_atoms_ij = np.zeros(pyscf_int_deriv.shape)

    if pyscf_int_type in ["int1e_iprinvip", "int1e_ipiprinv"]:
        # (nabla_i m | nabla_j operator | n)
        pyscf_int_deriv_atoms_ij[:,ki:kf,:] = pyscf_int_deriv[:,ki:kf,:]
    else:
        pyscf_int_deriv_atoms_ij[:,ki:kf,kj:kg] = pyscf_int_deriv[:,ki:kf,kj:kg]


    # (nabla**2 m | operator | n) + (m | operator | nabla**2 n)
    # or (nabla m | operator | nabla n) + (nabla n | operator | nabla m)
    # TODO: remove
    nao = pyscf_molecule.nao
    if full_deriv:
        pyscf_int_deriv_atoms_ij += pyscf_int_deriv_atoms_ij.transpose(0,2,1)
        if pyscf_int_type in ["int1e_iprinvip"]:
            sign = get_sign(pyscf_int_type) * molecule.elem_ids_to_numpy()[i]
            # (nabla_j m | nabla_i operator | n)
            with pyscf_molecule.with_rinv_at_nucleus(i): 
                new_pyscf_int_deriv = sign * pyscf_molecule.intor(pyscf_int_type,
                                                aosym='s1').reshape((3,3,nao,nao))
            transposed_new_pyscf_int_deriv = new_pyscf_int_deriv.transpose(1,0,2,3)
            reshaped_new_pyscf_int_deriv = transposed_new_pyscf_int_deriv.reshape((3*3,nao,nao))
            new_pyscf_int_deriv_atoms_ij = np.zeros(pyscf_int_deriv.shape)
            new_pyscf_int_deriv_atoms_ij[:,kj:kg] += reshaped_new_pyscf_int_deriv[:,kj:kg]
            new_pyscf_int_deriv_atoms_ij += new_pyscf_int_deriv_atoms_ij.transpose(0,2,1)
            pyscf_int_deriv_atoms_ij += new_pyscf_int_deriv_atoms_ij
        elif pyscf_int_type in ["int1e_ipiprinv"]:
            sign = get_sign(pyscf_int_type) * molecule.elem_ids_to_numpy()[i]
            # (nabla_j m | nabla_i operator | n)
            with pyscf_molecule.with_rinv_at_nucleus(i): 
                new_pyscf_int_deriv = sign * pyscf_molecule.intor(pyscf_int_type,
                                                aosym='s1')
            new_pyscf_int_deriv_atoms_ij = np.zeros(pyscf_int_deriv.shape)
            new_pyscf_int_deriv_atoms_ij[:,kj:kg] += new_pyscf_int_deriv[:,kj:kg]
            new_pyscf_int_deriv_atoms_ij += new_pyscf_int_deriv_atoms_ij.transpose(0,2,1)
            pyscf_int_deriv_atoms_ij += new_pyscf_int_deriv_atoms_ij

    # Transform integral to veloxchem format
    vlx_int_deriv = np.zeros_like(pyscf_int_deriv_atoms_ij)

    # number of derivatives (3)
    n_deriv = pyscf_int_deriv_atoms_ij.shape[0]
    for x in range(n_deriv):
        vlx_int_deriv[x] = ao_matrix_to_veloxchem(
                                 DenseMatrix(pyscf_int_deriv_atoms_ij[x]),
                                 basis, molecule).to_numpy()

    ao_basis_map = basis.get_ao_basis_map(molecule)
    rows = []
    columns = []
    k = 0
    for ao in ao_basis_map:
        parts = ao.split()
        atom = int(parts[0])
        shell = parts[2]
        if atom1 == atom:
            if shell1 is not None:
                for s in shell1:
                    if s in shell:
                        rows.append(k)
            else:
                rows.append(k)
        if atom2 == atom:
            if shell2 is not None:
                for s in shell2:
                    if s in shell:
                        columns.append(k)
            else:
                columns.append(k)
        k += 1
    if rows == []:
        raise ValueError("Atom or shell(s) not found.", atom1, shell1)
    if columns == []:
        raise ValueError("Atom or shell(s) not found.", atom2, shell2)
    

    vlx_int_block = np.zeros((n_deriv, len(rows), len(columns)))

    for k in range(len(rows)):
        for l in range(len(columns)):
            vlx_int_block[:, k, l] = vlx_int_deriv[:, rows[k], columns[l]]
   
    label = int_type+'_wrt_atom%d_atom%d_shells_atom%d' % (atomi, atomj, atom1)
    if shell1 is not None:
        for s1 in shell1:
            label += "_%s" % (s1)
    label += "_atom%d" % (atom2)
    if shell2 is not None:
        for s2 in shell2:
            label += "_%s" % (s2)

    if chk_file is not None:
        write_2d_array_hdf5(chk_file, [vlx_int_block], labels=[label])

    if return_block:
        return vlx_int_block
    else:
        return vlx_int_deriv

def import_2e_integral_derivative(molecule, basis, int_type, atomi, atom1=1,
                        shell1=None, atom2=1, shell2=None, atom3=1, shell3=None,
                        atom4=1, shell4=None, chk_file=None,
                        unit="au", return_block=True, full_deriv=False):
    """
    Imports two electron integral derivatives from pyscf and converts 
    them to veloxchem format.
    Specific atoms and shells can be selected.

    :param molecule:
        the vlx molecule object
    :param basis:
        the vlx basis object
    :param int_type:
        the type of one-electron integral: overlap, kinetic_energy,
                                           nuclear_attraction
    :param atomi:
        the index of the atom with respect to which the derivatives
        are taken (indexing starts at 1) 
    :param atom1:
        index of the first atom of interest
    :param shell1:
        list of atomic shells of interest for atom 1
    :param atom2:
        index of the second atom of interest
    :param shell2:
        list of atomic shells of interest for atom 2
    :param atom3:
        index of the third atom of interest
    :param shell3:
        list of atomic shells of interest for atom 3
    :param atom4:
        index of the fourth atom of interest
    :param shell2:
        list of atomic shells of interest for atom 4
    :param chk_file:
        the hdf5 checkpoint file name
    :param unit:
        the units to be used for the molecular geometry;
        possible values: "au" (default), "Angstrom"
    :param return_block:
        return the matrix block, or the full matrix.


    :return:
        a numpy array corresponding to a specified block of 
        the selected two electron integral derivative matrix.
    """

    molecule_string = get_molecule_string(molecule)
    basis_set_label = basis.get_label()
    pyscf_basis = translate_to_pyscf(basis_set_label)
    pyscf_molecule = pyscf.gto.M(atom=molecule_string,
                                 basis=pyscf_basis, unit=unit,
                                 charge=molecule.get_charge())

    # Get the AO indeces corresponding to atom i
    i = atomi - 1 # to use pyscf indexing

    pyscf_int_type = get_pyscf_integral_type(int_type)
    sign = get_sign(pyscf_int_type)
    pyscf_int_deriv = sign * pyscf_molecule.intor(pyscf_int_type, aosym='s1')

    nao = pyscf_molecule.nao

    ao_slices = pyscf_molecule.aoslice_by_atom()

    # Get the AO indeces corresponding to atom i
    ki, kf = ao_slices[i, 2:]

    pyscf_int_deriv_atom_i = np.zeros(pyscf_int_deriv.shape)

    pyscf_int_deriv_atom_i[:, ki:kf] = pyscf_int_deriv[:, ki:kf]

    if full_deriv:
        #   (nabla m n | t p) + ( m nabla n | t p)
        # + (m n | nabla t p) + (m n | t nabla p)
        pyscf_int_deriv_atom_i += ( pyscf_int_deriv_atom_i.transpose(0,2,1,4,3)
                             + pyscf_int_deriv_atom_i.transpose(0,3,4,1,2)
                             + pyscf_int_deriv_atom_i.transpose(0,4,3,2,1) )

    # Transform integral to veloxchem format
    vlx_int_deriv = np.zeros_like(pyscf_int_deriv_atom_i)

    basis_set_map = basis.get_index_map(molecule)
    for m in range(nao):
        for n in range(nao):
            for t in range(nao):
                for p in range(nao):
                    vm = basis_set_map[m]
                    vn = basis_set_map[n]
                    vt = basis_set_map[t]
                    vp = basis_set_map[p]
                    vlx_int_deriv[:,vm,vn,vt,vp] = (
                                    pyscf_int_deriv_atom_i[:,m,n,t,p] )

    ao_basis_map = basis.get_ao_basis_map(molecule)
    index1 = []
    index2 = []
    index3 = []
    index4 = []
    k = 0
    for ao in ao_basis_map:
        parts = ao.split()
        atom = int(parts[0])
        shell = parts[2]
        if atom1 == atom:
            if shell1 is not None:
                for s in shell1:
                    if s in shell:
                        index1.append(k)
            else:
                index1.append(k)
        if atom2 == atom:
            if shell2 is not None:
                for s in shell2:
                    if s in shell:
                        index2.append(k)
            else:
                index2.append(k)
        if atom3 == atom:
            if shell3 is not None:
                for s in shell3:
                    if s in shell:
                        index3.append(k)
            else:
                index3.append(k)
        if atom4 == atom:
            if shell4 is not None:
                for s in shell4:
                    if s in shell:
                        index4.append(k)
            else:
                index4.append(k)
        k += 1
    if index1 == []:
        raise ValueError("Atom or shell(s) not found.", atom1, shell1)
    if index2 == []:
        raise ValueError("Atom or shell(s) not found.", atom2, shell2)
    if index3 == []:
        raise ValueError("Atom or shell(s) not found.", atom3, shell3)
    if index4 == []:
        raise ValueError("Atom or shell(s) not found.", atom4, shell4)
   
    n_deriv = pyscf_int_deriv_atom_i.shape[0]
    vlx_int_block = np.zeros((n_deriv, len(index1), len(index2),
                              len(index3), len(index4)))
    
    for k in range(len(index1)):
        for l in range(len(index2)):
            for m in range(len(index3)):
                for n in range(len(index4)):
                    vlx_int_block[:,k,l,m,n] = vlx_int_deriv[:,index1[k],
                                                       index2[l], index3[m],
                                                       index4[n]]
   
    label = int_type+'_atom%d_shells_atom%d' % (atomi, atom1)
    if shell1 is not None:
        for s1 in shell1:
            label += "_%s" % (s1)

    label += "_atom%d" % (atom2)
    if shell2 is not None:
        for s2 in shell2:
            label += "_%s" % (s2)

    label += "_atom%d" % (atom3)
    if shell3 is not None:
        for s3 in shell3:
            label += "_%s" % (s3)

    label += "_atom%d" % (atom4)
    if shell4 is not None:
        for s4 in shell4:
            label += "_%s" % (s4)

    if chk_file is not None:
        write_2d_array_hdf5(chk_file, [vlx_int_block], labels=[label])

    if return_block:
        return vlx_int_block
    else:
        return vlx_int_deriv

def import_2e_second_order_integral_derivative(molecule, basis, int_type,
                        atomi, atomj, atom1=1, shell1=None, atom2=1,
                        shell2=None, atom3=1, shell3=None, atom4=1,
                        shell4=None, chk_file=None, unit="au",
                        return_block=True, full_deriv=False):
    """
    Imports two electron integral derivatives from pyscf and converts 
    them to veloxchem format.
    Specific atoms and shells can be selected.

    :param molecule:
        the vlx molecule object
    :param basis:
        the vlx basis object
    :param int_type:
        the type of one-electron integral: overlap, kinetic_energy,
                                           nuclear_attraction
    :param atomi:
        the index of the first atom with respect to which the derivatives
        are taken (indexing starts at 1)
    :param atomj:
        the index of the second atom with respect to which the derivatives
        are taken (indexing starts at 1) 
    :param atom1:
        index of the first atom of interest
    :param shell1:
        list of atomic shells of interest for atom 1
    :param atom2:
        index of the second atom of interest
    :param shell2:
        list of atomic shells of interest for atom 2
    :param atom3:
        index of the third atom of interest
    :param shell3:
        list of atomic shells of interest for atom 3
    :param atom4:
        index of the fourth atom of interest
    :param shell2:
        list of atomic shells of interest for atom 4
    :param chk_file:
        the hdf5 checkpoint file name
    :param unit:
        the units to be used for the molecular geometry;
        possible values: "au" (default), "Angstrom"
    :param return_block:
        return the matrix block, or the full matrix.


    :return:
        a numpy array corresponding to a specified block of 
        the selected two electron integral derivative matrix.
    """

    molecule_string = get_molecule_string(molecule)
    basis_set_label = basis.get_label()
    pyscf_basis = translate_to_pyscf(basis_set_label)
    pyscf_molecule = pyscf.gto.M(atom=molecule_string,
                                 basis=pyscf_basis, unit=unit,
                                 charge=molecule.get_charge())

    # Get the indices corresponding to atom i, j
    i = atomi - 1 # to use pyscf indexing
    j = atomj - 1

    pyscf_int_type = get_pyscf_integral_type(int_type)
    sign = get_sign(pyscf_int_type)

    pyscf_int_deriv = sign * pyscf_molecule.intor(pyscf_int_type, aosym='s1')

    nao = pyscf_molecule.nao

    ao_slices = pyscf_molecule.aoslice_by_atom()

    # Get the AO indeces corresponding to atom i
    ki, kf = ao_slices[i, 2:]
    kj, kg = ao_slices[j, 2:]

    pyscf_int_deriv_atoms_ij = np.zeros(pyscf_int_deriv.shape)

    if pyscf_int_type in ["int2e_ipvip1"]:
        # ( nabla m nabla n | p q ) 
        pyscf_int_deriv_atoms_ij[:, ki:kf, kj:kg] = pyscf_int_deriv[:, ki:kf, kj:kg]
    else:
        # ( nabla m n | nabla p q )
        pyscf_int_deriv_atoms_ij[:, ki:kf, :, kj:kg] = pyscf_int_deriv[:, ki:kf, :, kj:kg]

    if full_deriv:
        if pyscf_int_type in ["int2e_ipvip1"]:
            #   (nabla_i m nabla_j n | p t ) + ( nabla_j m nabla_i n | p t )
            # + (m n | nabla_i p  nabla_j t) + ( m n | nabla_j p nabla_i t)
            pyscf_int_deriv_atoms_ij += ( pyscf_int_deriv_atoms_ij.transpose(0,2,1,4,3)
                             + pyscf_int_deriv_atoms_ij.transpose(0,3,4,1,2)
                             + pyscf_int_deriv_atoms_ij.transpose(0,4,3,2,1) )

        elif pyscf_int_type in ["int2e_ip1ip2"]:
            #   (nabla_i m n | nabla_j p t ) + ( m nabla_i n | nabla_j p t )
            # + (nabla_i m n | p  nabla_j t) + ( m nabla_i n | p nabla_j t)
            #   (nabla_j m n | nabla_i p t ) + ( m nabla_j n | nabla_i p t )
            # + (nabla_j m n | p  nabla_i t) + ( m nabla_j n | p nabla_i t)
            pyscf_int_deriv_atoms_ij += ( pyscf_int_deriv_atoms_ij.transpose(0,2,1,3,4)
                         + pyscf_int_deriv_atoms_ij.transpose(0,1,2,4,3)
                         + pyscf_int_deriv_atoms_ij.transpose(0,2,1,4,3)
                         + pyscf_int_deriv_atoms_ij.transpose(0,3,4,1,2)
                         + pyscf_int_deriv_atoms_ij.transpose(0,4,3,1,2)
                         + pyscf_int_deriv_atoms_ij.transpose(0,3,4,2,1)
                         + pyscf_int_deriv_atoms_ij.transpose(0,4,3,2,1) )

    # Transform integral to veloxchem format
    vlx_int_deriv = np.zeros_like(pyscf_int_deriv_atoms_ij)

    basis_set_map = basis.get_index_map(molecule)
    for m in range(nao):
        for n in range(nao):
            for t in range(nao):
                for p in range(nao):
                    vm = basis_set_map[m]
                    vn = basis_set_map[n]
                    vt = basis_set_map[t]
                    vp = basis_set_map[p]
                    vlx_int_deriv[:,vm,vn,vt,vp] = (
                                    pyscf_int_deriv_atoms_ij[:,m,n,t,p] )

    ao_basis_map = basis.get_ao_basis_map(molecule)
    index1 = []
    index2 = []
    index3 = []
    index4 = []
    k = 0
    for ao in ao_basis_map:
        parts = ao.split()
        atom = int(parts[0])
        shell = parts[2]
        if atom1 == atom:
            if shell1 is not None:
                for s in shell1:
                    if s in shell:
                        index1.append(k)
            else:
                index1.append(k)
        if atom2 == atom:
            if shell2 is not None:
                for s in shell2:
                    if s in shell:
                        index2.append(k)
            else:
                index2.append(k)
        if atom3 == atom:
            if shell3 is not None:
                for s in shell3:
                    if s in shell:
                        index3.append(k)
            else:
                index3.append(k)
        if atom4 == atom:
            if shell4 is not None:
                for s in shell4:
                    if s in shell:
                        index4.append(k)
            else:
                index4.append(k)
        k += 1
    if index1 == []:
        raise ValueError("Atom or shell(s) not found.", atom1, shell1)
    if index2 == []:
        raise ValueError("Atom or shell(s) not found.", atom2, shell2)
    if index3 == []:
        raise ValueError("Atom or shell(s) not found.", atom3, shell3)
    if index4 == []:
        raise ValueError("Atom or shell(s) not found.", atom4, shell4)
   
    n_deriv = pyscf_int_deriv_atoms_ij.shape[0]
    vlx_int_block = np.zeros((n_deriv, len(index1), len(index2),
                              len(index3), len(index4)))
    
    for k in range(len(index1)):
        for l in range(len(index2)):
            for m in range(len(index3)):
                for n in range(len(index4)):
                    vlx_int_block[:,k,l,m,n] = vlx_int_deriv[:,index1[k],
                                                       index2[l], index3[m],
                                                       index4[n]]
   
    label = int_type+'_atom%d_atom%d_shells_atom%d' % (atomi, atomj, atom1)
    if shell1 is not None:
        for s1 in shell1:
            label += "_%s" % (s1)

    label += "_atom%d" % (atom2)
    if shell2 is not None:
        for s2 in shell2:
            label += "_%s" % (s2)

    label += "_atom%d" % (atom3)
    if shell3 is not None:
        for s3 in shell3:
            label += "_%s" % (s3)

    label += "_atom%d" % (atom4)
    if shell4 is not None:
        for s4 in shell4:
            label += "_%s" % (s4)

    if chk_file is not None:
        write_2d_array_hdf5(chk_file, [vlx_int_block], labels=[label])

    if return_block:
        return vlx_int_block
    else:
        return vlx_int_deriv

def import_xc_contrib_tddft(molecule, basis, scfdrv, xmy_ao, rel_dm_ao,
                            tda=False, vxc_deriv_only=False, unit="au"):
    """
    Imports the xc contribution to the TDDFT gradient from pyscf.

    :param molecule:
        the vlx molecule object.
    :param basis:
        the vlx basis set object.
    :param scfdrv:
        the vlx SCF driver.
    :param xmy_ao:
        the perturbed density matrix in AO basis.
    :param rel_dm_ao:
        the relaxed one-particle density matrix in AO basis.
    :param tda:
        flag to use the Tamm-Dancoff approximation or not.
    :param vxc_deriv_only:
        Return only the vxc derivative
    """
    # create pyscf objects: molecule, scf_driver, 
    molecule_string = get_molecule_string(molecule)
    basis_set_label = basis.get_label()
    pyscf_basis = translate_to_pyscf(basis_set_label)
    pyscf_molecule = pyscf.gto.M(atom=molecule_string,
                                 basis=pyscf_basis, unit=unit,
                                 charge=molecule.get_charge())

    pyscf_scf = pyscf.scf.RKS(pyscf_molecule)

    # set functional, convergence threshold and grid level 
    pyscf_scf.xc = translate_to_pyscf(scfdrv.xcfun.get_func_label())
    pyscf_scf.conv_tol = scfdrv.conv_thresh
    if scfdrv.grid_level is None:
        scfdrv.grid_level = get_default_grid_level(scfdrv.xcfun)
    if scfdrv.grid_level == 6 or scfdrv.grid_level == 7:
        pyscf_scf.grids.level = 9
    else:
        pyscf_scf.grids.level = scfdrv.grid_level

    pyscf_scf.kernel()

    # create a TDDFT_Gradient object in pyscf
    if tda:
        pyscf_tddft = pyscf_scf.TDA()
    else:
        pyscf_tddft = pyscf_scf.TDDFT()

    pyscf_tddft_grad = pyscf_tddft.Gradients()
    mf = pyscf_tddft_grad.base._scf

    gs_dm = 2*scfdrv.scf_tensors['D_alpha']

    # transform densities to pyscf format:
    pyscf_xmy = ao_matrix_to_dalton(DenseMatrix(xmy_ao),
                                    basis, molecule).to_numpy()
    pyscf_rel_dm = ao_matrix_to_dalton(DenseMatrix(rel_dm_ao),
                                    basis, molecule).to_numpy()
    pyscf_gs_dm = ao_matrix_to_dalton(DenseMatrix(gs_dm),
                                    basis, molecule).to_numpy()

    # calculate xc derivatives:
    fxc_xmy, fxc_rel, vxc, gxc = (
            tdrks._contract_xc_kernel(pyscf_tddft_grad, mf.xc, pyscf_xmy,
                                   pyscf_rel_dm, True, True, True, 2000) )

    # return vxc derivatives in pyscf orbital order
    if vxc_deriv_only:
        return vxc, fxc_rel

    # contract to get the final xc gradient contributions
    natm = molecule.number_of_atoms()
    vxc_contrib = np.zeros((natm, 3))
    vxc_contrib_2 = np.zeros((natm, 3))
    fxc_contrib = np.zeros((natm, 3))
    fxc_contrib_2 = np.zeros((natm, 3))

    atmlst = range(natm)
    offsetdic = pyscf_molecule.offset_nr_by_atom()

    nao = scfdrv.scf_tensors['C_alpha'].shape[0]
    
    for k, ia in enumerate(atmlst):
        shl0, shl1, p0, p1 = offsetdic[ia]

        veff = np.zeros((3, nao, nao))
        veff[:,p0:p1] = vxc[1:,p0:p1]
        veff[:,:,p0:p1] += vxc[1:,p0:p1].transpose(0,2,1)

        vxc_contrib[k] = np.einsum('xpq,pq->x', veff, pyscf_rel_dm)

        vxc_contrib_2[k] = np.einsum('xij,ij->x', fxc_rel[1:,p0:p1],
                                      pyscf_gs_dm[p0:p1])
        fxc_contrib[k] = ( 2.0 * np.einsum('xij,ij->x', fxc_xmy[1:,p0:p1], 
                                            pyscf_xmy[p0:p1,:])
                          + 2.0 * np.einsum('xji,ij->x', fxc_xmy[1:,p0:p1], 
                                            pyscf_xmy[:,p0:p1])
                          )
        fxc_contrib_2[k] = np.einsum('xij,ij->x', gxc[1:,p0:p1],
                                     pyscf_gs_dm[p0:p1])

    return vxc_contrib, vxc_contrib_2, fxc_contrib, fxc_contrib_2 

def overlap_deriv(molecule, basis, i=0, full_deriv=True, unit="au",
                  chk_file=None):
    """
    Imports the derivatives of the overlap matrix
    from pyscf and converts it to veloxchem format

    :param molecule:
        the vlx molecule object
    :param basis:
        the vlx basis object
    :param i:
        the index of the atom for which the derivatives
        are computed.
    :param full_deriv:
        True, to compute ( nabla m | n ) + ( m | nabla n) 
        False, to compute (nabla m | n) only.
    :param unit:
        the units to be used for the molecular geometry;
        possible values: "au" (default), "Angstrom"
    :param chk_file:
        the hdf5 checkpoint file name.

    :return:
        a numpy array of shape 3 x nao x nao
        (nao = number of atomic orbitals)
        corresponding to the derivative of the overlap matrix
        with repsect to the x, y and z coords. of atom i.
    """

    molecule_string = get_molecule_string(molecule)
    basis_set_label = basis.get_label()
    pyscf_basis = translate_to_pyscf(basis_set_label)
    pyscf_molecule = pyscf.gto.M(atom=molecule_string,
                                 basis=pyscf_basis, unit=unit,
                                 charge=molecule.get_charge())

    # The "-" sign is due to the fact that pyscf computes -(nabla m | n)
    pyscf_ovlp_deriv = - pyscf_molecule.intor('int1e_ipovlp', aosym='s1')
    ao_slices = pyscf_molecule.aoslice_by_atom()

    # Get the AO indeces corresponding to atom i
    ki, kf = ao_slices[i, 2:]

    overlap_deriv_atom_i = np.zeros(pyscf_ovlp_deriv.shape)

    overlap_deriv_atom_i[:,ki:kf] = pyscf_ovlp_deriv[:,ki:kf]

    if full_deriv:
        # (nabla m | n) + (m | nabla n)
        overlap_deriv_atom_i += overlap_deriv_atom_i.transpose(0,2,1)

    # Transform the oo block to MO basis
    #pyscf_scf = pyscf.scf.RHF(pyscf_molecule)
    #pyscf_scf.conv_tol = 1e-10
    #pyscf_scf.kernel()
    #mo_occ = pyscf_scf.mo_coeff[:,pyscf_scf.mo_occ>0]
    #print("mo_occ obtained from pyscf:\n", mo_occ)
    #s1oo = np.einsum('mi,xmn,nj->xij', mo_occ, overlap_deriv_atom_i, mo_occ)
    #print("\n\ns1oo transformed with pyscf data:\n", s1oo)

    vlx_ovlp_deriv_atom_i = np.zeros(overlap_deriv_atom_i.shape)


    # Transform each component (x,y,z) to veloxchem format
    vlx_ovlp_deriv_atom_i[0] = ( ao_matrix_to_veloxchem(
                                 DenseMatrix(overlap_deriv_atom_i[0]),
                                 basis, molecule).to_numpy()
                                )

    vlx_ovlp_deriv_atom_i[1] = ( ao_matrix_to_veloxchem(
                                 DenseMatrix(overlap_deriv_atom_i[1]),
                                 basis, molecule).to_numpy()
                                )

    vlx_ovlp_deriv_atom_i[2] = ( ao_matrix_to_veloxchem(
                                 DenseMatrix(overlap_deriv_atom_i[2]),
                                 basis, molecule).to_numpy()
                                )

    if chk_file is not None:
        write_2d_array_hdf5(chk_file, vlx_ovlp_deriv_atom_i,
                            labels=['overlap_x', 'overlap_y', 'overlap_z'],
                            atom_index=i)

    return vlx_ovlp_deriv_atom_i

def compute_dipole_integral_derivatives(molecule, ao_basis, i, unit="au",
                                        delta_h=1e-5, chk_file=None):
    """
    Computes numerical derivatives of dipole integrals.

    :param molecule:
        The molecule.
    :param ao_basis:
        The AO basis set.
    :param i:
        The atom index.
    :param delta_h:
        The step for the numerical derivative.
    :chk_file:
        The name of the checkpoint file. 

    :return:
        The dipole integral derivatives.
    """

    # atom labels
    labels = molecule.get_labels()

    # number of atoms
    natm = molecule.number_of_atoms()

    # atom coordinates (nx3)
    coords = molecule.get_coordinates()

    # number of atomic orbitals
    molecule_string = get_molecule_string(molecule)
    basis_set_label = ao_basis.get_label()
    pyscf_basis = translate_to_pyscf(basis_set_label)
    pyscf_molecule = pyscf.gto.M(atom=molecule_string,
                                 basis=pyscf_basis, unit=unit,
                                 charge=molecule.get_charge())
    nao = pyscf_molecule.nao

    # Dipole integrals driver
    dipole_drv = ElectricDipoleIntegralsDriver()

    # 3 dipole components x 3 atomic coordinates
    # x No. basis x No. basis
    dipole_integrals_gradient = np.zeros((3, 3, nao, nao))

    for d in range(3):
        coords[i, d] += delta_h
        new_mol = Molecule(labels, coords, units='au')

        dipole_mats_p = dipole_drv.compute(new_mol, ao_basis)
        dipole_ints_p = (dipole_mats_p.x_to_numpy(), 
                         dipole_mats_p.y_to_numpy(),
                         dipole_mats_p.z_to_numpy())

        coords[i, d] -= 2.0 * delta_h
        new_mol = Molecule(labels, coords, units='au')

        dipole_mats_m = dipole_drv.compute(new_mol, ao_basis)
        dipole_ints_m = (dipole_mats_m.x_to_numpy(),
                         dipole_mats_m.y_to_numpy(),
                         dipole_mats_m.z_to_numpy())

        for c in range(3):
            dipole_integrals_gradient[c, d] = (
                ( dipole_ints_p[c] - dipole_ints_m[c] ) / (2.0 * delta_h)
            )
    if chk_file is not None:
        write_2d_array_hdf5(chk_file,
                dipole_integrals_gradient.reshape(3*3,nao,nao),
                labels=['xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz'],
                atom_index=i)        

    return dipole_integrals_gradient


def kinen_deriv(molecule, basis, i=0, full_deriv=True, unit="au",
                chk_file=None):
    """
    Imports the derivatives of the kinetic energy matrix
    from pyscf and converts it to veloxchem format
    kinetic energy matrix element: ( m | 0.5 (nabla_e)**2 | n )
    (where nabla_e acts on electronic coordinates)

    :param molecule:
        the vlx molecule object
    :param basis:
        the vlx basis object
    :param i:
        the index of the atom for which the derivatives
        are computed.
    :param full_deriv:
        True, to compute ( nabla m | 0.5 (nabla_e)**2 | n ) 
                       + ( m | 0.5 (nabla_e)**2 | nabla n ) 
        False, to compute (nabla m | nabla_e | n) only.
    :param unit:
        the units to be used for the molecular geometry;
        possible values: "au" (default), "Angstrom"
    :param chk_file:
        the hdf5 checkpoint file name.

    :return:
        a numpy array of shape 3 x nao x nao
        (nao = number of atomic orbitals)
        corresponding to the derivative of the kinetic energy matrix
        with repsect to the x, y and z coords. of atom i.
    """

    molecule_string = get_molecule_string(molecule)
    basis_set_label = basis.get_label()
    pyscf_basis = translate_to_pyscf(basis_set_label)
    pyscf_molecule = pyscf.gto.M(atom=molecule_string,
                                 basis=pyscf_basis, unit=unit,
                                 charge=molecule.get_charge())

    # TODO:check sign; correct mistakes in comments.
    pyscf_kinen_deriv = -pyscf_molecule.intor('int1e_ipkin', aosym='s1')
    ao_slices = pyscf_molecule.aoslice_by_atom()

    # Get the AO indeces corresponding to atom i
    ki, kf = ao_slices[i, 2:]

    kinen_deriv_atom_i = np.zeros(pyscf_kinen_deriv.shape)

    kinen_deriv_atom_i[:,ki:kf] = pyscf_kinen_deriv[:,ki:kf]

    if full_deriv:
        # (nabla m | n) + (m | nabla n)
        kinen_deriv_atom_i += kinen_deriv_atom_i.transpose(0,2,1)

    vlx_kinen_deriv_atom_i = np.zeros(kinen_deriv_atom_i.shape)


    # Transform each component (x,y,z) to veloxchem format
    vlx_kinen_deriv_atom_i[0] = ( ao_matrix_to_veloxchem(
                                 DenseMatrix(kinen_deriv_atom_i[0]),
                                 basis, molecule).to_numpy()
                                )

    vlx_kinen_deriv_atom_i[1] = ( ao_matrix_to_veloxchem(
                                 DenseMatrix(kinen_deriv_atom_i[1]),
                                 basis, molecule).to_numpy()
                                )

    vlx_kinen_deriv_atom_i[2] = ( ao_matrix_to_veloxchem(
                                 DenseMatrix(kinen_deriv_atom_i[2]),
                                 basis, molecule).to_numpy()
                                )

    if chk_file is not None:
        write_2d_array_hdf5(chk_file, vlx_kinen_deriv_atom_i, labels='xyz',
                            atom_index=i)

    return vlx_kinen_deriv_atom_i

def pyscf_hcore_deriv_generator(molecule, basis, i=0, unit="au"):
    """
    Copy of the pyscf hcore_generator which computes the Core Hamiltonian
    derivatives.
    """
    molecule_string = get_molecule_string(molecule)
    basis_set_label = basis.get_label()
    pyscf_basis = translate_to_pyscf(basis_set_label)
    pyscf_molecule = pyscf.gto.M(atom=molecule_string,
                                 basis=pyscf_basis, unit=unit,
                                 charge=molecule.get_charge())

    # The atomic charge of atom i -- required for the derivative of 1/RA
    Z = molecule.elem_ids_to_numpy()[i]

    h = pyscf_molecule.intor('int1e_ipkin', comp=3)
    h+= pyscf_molecule.intor('int1e_ipnuc', comp=3)

    ao_slices = pyscf_molecule.aoslice_by_atom()

    shl0, shl1, p0, p1 = ao_slices[i] 
    with pyscf_molecule.with_rinv_at_nucleus(i):
        vrinv = pyscf_molecule.intor('int1e_iprinv', comp=3) # <\nabla|1/r|>
        vrinv *= -pyscf_molecule.atom_charge(i)
    vrinv[:,p0:p1] -= h[:,p0:p1]

    hcore_deriv_i = vrinv + vrinv.transpose(0,2,1)

    vlx_hcore_deriv_atom_i = np.zeros(hcore_deriv_i.shape)


    # Transform each compnent (x,y,z) to veloxchem format
    vlx_hcore_deriv_atom_i[0] = ( ao_matrix_to_veloxchem(
                                 DenseMatrix(hcore_deriv_i[0]),
                                 basis, molecule).to_numpy()
                                )

    vlx_hcore_deriv_atom_i[1] = ( ao_matrix_to_veloxchem(
                                 DenseMatrix(hcore_deriv_i[1]),
                                 basis, molecule).to_numpy()
                                )

    vlx_hcore_deriv_atom_i[2] = ( ao_matrix_to_veloxchem(
                                 DenseMatrix(hcore_deriv_i[2]),
                                 basis, molecule).to_numpy()
                                )

    return vlx_hcore_deriv_atom_i
    



def coulomb_attract_deriv(molecule, basis, i=0, full_deriv=True, unit="au",
                          chk_file=None):
    """
    Imports the derivatives of the nuclear-electron Coulomb energy matrix
    from pyscf and converts it to veloxchem format
    Coulomb energy matrix element: ( m | 1/RA | n )

    :param molecule:
        the vlx molecule object
    :param basis:
        the vlx basis object
    :param i:
        the index of the atom for which the derivatives
        are computed.
    :param full_deriv:
        True, to compute ( nabla m | 1/RA | n ) 
                       + ( m | 1/RA | nabla n )
                       + ( m | nabla 1/RA | n ) 
        False, to compute (nabla m | 1/RA | n) + (m | nabla 1/RA | n) only.
    :param unit:
        the units to be used for the molecular geometry;
        possible values: "au" (default), "Angstrom"
    :param chk_file:
        the hdf5 checkpoint file name.

    :return:
        a numpy array of shape 3 x nao x nao
        (nao = number of atomic orbitals)
        corresponding to the derivative of the kinetic energy matrix
        with repsect to the x, y and z coords. of atom i.
    """

    molecule_string = get_molecule_string(molecule)
    basis_set_label = basis.get_label()
    pyscf_basis = translate_to_pyscf(basis_set_label)
    pyscf_molecule = pyscf.gto.M(atom=molecule_string,
                                 basis=pyscf_basis, unit=unit,
                                 charge=molecule.get_charge())

    # The atomic charge of atom i -- required for the derivative of 1/RA
    Z = molecule.elem_ids_to_numpy()[i]

    # TODO:check sign; correct mistakes in comments.
    pyscf_ipnuc_deriv = -pyscf_molecule.intor('int1e_ipnuc', aosym='s1')

    with pyscf_molecule.with_rinv_at_nucleus(i):
        pyscf_iprinv_deriv = pyscf_molecule.intor('int1e_iprinv', aosym='s1')
        pyscf_iprinv_deriv *= -Z 

    ao_slices = pyscf_molecule.aoslice_by_atom()

    # Get the AO indeces corresponding to atom i
    ki, kf = ao_slices[i, 2:]

    ipnuc_deriv_atom_i = np.zeros(pyscf_ipnuc_deriv.shape)

    ipnuc_deriv_atom_i[:,ki:kf] = pyscf_ipnuc_deriv[:,ki:kf]

    if full_deriv:
        # (nabla m | 1/RA | n) + (m | 1/RA | nabla n)
        ipnuc_deriv_atom_i += ipnuc_deriv_atom_i.transpose(0,2,1)
        # TODO: write term explicitly and figure out why the transpose has to
        # be summed up here
        pyscf_iprinv_deriv += pyscf_iprinv_deriv.transpose(0,2,1)

    vlx_coulomb_att_deriv_atom_i = np.zeros(ipnuc_deriv_atom_i.shape)


    # Transform each component (x,y,z) to veloxchem format
    vlx_coulomb_att_deriv_atom_i[0] = ( ao_matrix_to_veloxchem(
                DenseMatrix(pyscf_iprinv_deriv[0] + ipnuc_deriv_atom_i[0]),
                                 basis, molecule).to_numpy()
                                )

    vlx_coulomb_att_deriv_atom_i[1] = ( ao_matrix_to_veloxchem(
                 DenseMatrix(pyscf_iprinv_deriv[1] + ipnuc_deriv_atom_i[1]),
                                 basis, molecule).to_numpy()
                                )

    vlx_coulomb_att_deriv_atom_i[2] = ( ao_matrix_to_veloxchem(
               DenseMatrix(pyscf_iprinv_deriv[2] + ipnuc_deriv_atom_i[2]),
                                 basis, molecule).to_numpy()
                                )

    if chk_file is not None:
        write_2d_array_hdf5(chk_file, vlx_coulomb_att_deriv_atom_i, labels='xyz',
                            atom_index=i)

    return vlx_coulomb_att_deriv_atom_i


def hcore_deriv(molecule, basis, i=0, unit="au"):
    """
    Imports the derivatives of the Fock matrix
    from pyscf and converts it to veloxchem format

    :param molecule:
        the vlx molecule object
    :param basis:
        the vlx basis object
    :param i:
        the index of the atom for which the derivatives
        are computed.
    :param unit:
        the units to be used for the molecular geometry;
        possible values: "au" (default), "Angstrom"

    :return:
        a numpy array of shape 3 x nao x nao
        (nao = number of atomic orbitals)
        corresponding to the derivative of the core Hamiltonian matrix
        with respect to the x, y and z coords. of atom i.
    """
    molecule_string = get_molecule_string(molecule)
    basis_set_label = basis.get_label()
    pyscf_basis = translate_to_pyscf(basis_set_label)
    pyscf_molecule = pyscf.gto.M(atom=molecule_string,
                                 basis=pyscf_basis, unit=unit,
                                 charge=molecule.get_charge())
    pyscf_scf = pyscf.scf.RHF(pyscf_molecule)
    nao = pyscf_molecule.nao
    pyscf_grad = grad.RHF(pyscf_scf)

    ao_slices = pyscf_molecule.aoslice_by_atom()
    hcore_generator = pyscf_grad.hcore_generator(pyscf_molecule)

    # Get the AO indices corresponding to atom i
    ki, kf = ao_slices[i, 2:]

    hcore_deriv_atom_i = np.zeros((3,nao,nao))

    hcore_deriv_atom_i = hcore_generator(i)

    vlx_hcore_deriv_atom_i = np.zeros(hcore_deriv_atom_i.shape)


    # Transform each compnent (x,y,z) to veloxchem format
    vlx_hcore_deriv_atom_i[0] = ( ao_matrix_to_veloxchem(
                                 DenseMatrix(hcore_deriv_atom_i[0]),
                                 basis, molecule).to_numpy()
                                )

    vlx_hcore_deriv_atom_i[1] = ( ao_matrix_to_veloxchem(
                                 DenseMatrix(hcore_deriv_atom_i[1]),
                                 basis, molecule).to_numpy()
                                )

    vlx_hcore_deriv_atom_i[2] = ( ao_matrix_to_veloxchem(
                                 DenseMatrix(hcore_deriv_atom_i[2]),
                                 basis, molecule).to_numpy()
                                )

    return vlx_hcore_deriv_atom_i

# TODO: MPI does not work efficiently here. XC part is executed
# using MPI, while the rest is run on each node separately.
# Either fix, or replace with new integrals code.
def fock_deriv(molecule, basis, density, i=0, scfdrv=None, unit="au"):
    """
    Imports the derivatives of the Fock matrix
    from pyscf and converts it to veloxchem format

    :param molecule:
        the vlx molecule object
    :param basis:
        the vlx basis object
    :param density:
        The SCF density matrix (alpha part) in AO basis
    :param i:
        the index of the atom for which the derivatives
        are computed.
    :param unit:
        the units to be used for the molecular geometry;
        possible values: "au" (default), "Angstrom"
    :param scfdrv:
        The SCF driver (needed for xc derivatives in DFT)

    :return:
        a numpy array of shape 3 x nao x nao
        (nao = number of atomic orbitals)
        corresponding to the derivative of the Fock matrix
        with respect to the x, y and z coords. of atom i.
    """
    molecule_string = get_molecule_string(molecule)
    basis_set_label = basis.get_label()
    pyscf_basis = translate_to_pyscf(basis_set_label)
    pyscf_molecule = pyscf.gto.M(atom=molecule_string,
                                 basis=pyscf_basis, unit=unit,
                                 charge=molecule.get_charge())

    gs_dm = ao_matrix_to_dalton(DenseMatrix(density), basis, molecule).to_numpy()
    nao = density.shape[0]
    vlx_vxc_deriv_atom_i = np.zeros((3, nao, nao))

    if scfdrv is not None:
        if scfdrv._dft:
            pyscf_scf = pyscf.scf.RKS(pyscf_molecule)
            # set functional, convergence threshold and grid level 
            pyscf_scf.xc = translate_to_pyscf(scfdrv.xcfun.get_func_label())
            pyscf_scf.conv_tol = scfdrv.conv_thresh
            if scfdrv.grid_level is None:
                scfdrv.grid_level = get_default_grid_level(scfdrv.xcfun)
            if scfdrv.grid_level == 6 or scfdrv.grid_level == 7:
                pyscf_scf.grids.level = 9
            else:
                pyscf_scf.grids.level = scfdrv.grid_level

            #if scfdrv.xcfun.get_func_type() == xcfun.lda:
            grid_drv = GridDriver()
            if scfdrv.grid_level is not None:
                grid_drv.set_level(scfdrv.grid_level)
            else:
                grid_level = get_default_grid_level(scfdrv.xcfun)
            xc_mol_hess = XCMolecularHessian()
            mol_grid = grid_drv.generate(molecule)
            vlx_density = scfdrv.density
            vlx_vxc_deriv_atom_i = xc_mol_hess.integrate_vxc_fock_gradient(
                                    molecule, basis, vlx_density, mol_grid,
                                    scfdrv.xcfun.get_func_label(), i)
            #else:
            #    pyscf_scf.kernel()
            #
            #    pyscf_hessian = pyscf.hessian.rks.Hessian(pyscf_scf)
        else:
            pyscf_scf = pyscf.scf.RHF(pyscf_molecule)
    else:
        pyscf_scf = pyscf.scf.RHF(pyscf_molecule)
    
    if nao != pyscf_molecule.nao:
        error_text = "vlx and pyscf number of atomic orbitals are different!"
        error_text +="\nCheck if the basis sets are defined the same way."
        raise ValueError(error_text)    

    pyscf_grad = grad.RHF(pyscf_scf)

    # The "-" sign is due to the fact that pyscf computes -(nabla m n | t p)
    pyscf_eri_deriv = - pyscf_molecule.intor('int2e_ip1', aosym='s1')
    ao_slices = pyscf_molecule.aoslice_by_atom()
    hcore_generator = pyscf_grad.hcore_generator(pyscf_molecule)

    # Get the AO indices corresponding to atom i
    ki, kf = ao_slices[i, 2:]

    eri_deriv_atom_i = np.zeros(pyscf_eri_deriv.shape)
    fock_deriv_atom_i = np.zeros((3,nao,nao))

    #   (nabla m n | t p) + ( m nabla n | t p)
    # + (m n | nabla t p) + (m n | t nabla p)
    eri_deriv_atom_i[:, ki:kf] = pyscf_eri_deriv[:, ki:kf]
    eri_deriv_atom_i += ( eri_deriv_atom_i.transpose(0,2,1,4,3)
                         + eri_deriv_atom_i.transpose(0,3,4,1,2)
                         + eri_deriv_atom_i.transpose(0,4,3,2,1) )

    fock_deriv_atom_i = hcore_generator(i)

    # fraction of exact exchange
    x_frac = 1.0

    if scfdrv is not None:
        if scfdrv._dft:
            # This one works!

            #if scfdrv.xcfun.get_func_type() != xcfun.lda:
            #    pyscf_mo_coeff = pyscf_scf.mo_coeff
            #    pyscf_mo_occ = pyscf_scf.mo_occ
            #    max_memory = 2000
            #    vxc_deriv1 = pyscf_rks_hessian._get_vxc_deriv1(pyscf_hessian,
            #                                              pyscf_mo_coeff,
            #                                              pyscf_mo_occ,
            #                                              max_memory)

            #    #print(vxc_deriv1)
            #    # select derivatives wrt. atom i
            #    vxc_deriv_atom_i = vxc_deriv1[i]
            #    fock_deriv_atom_i += vxc_deriv_atom_i

            # check fraction of exact exchange
            if scfdrv.xcfun.is_hybrid():
                x_frac = scfdrv.xcfun.get_frac_exact_exchange()
            else:
                x_frac = 0.0
            

    fock_deriv_atom_i += ( 2 * np.einsum('xmntf,tf->xmn',
                                          eri_deriv_atom_i, gs_dm)
                         - x_frac * np.einsum('xmtnf,tf->xmn',
                                      eri_deriv_atom_i, gs_dm)
                        )

    vlx_fock_deriv_atom_i = np.zeros(fock_deriv_atom_i.shape)


    # Transform each compnent (x,y,z) to veloxchem format
    vlx_fock_deriv_atom_i[0] = ( ao_matrix_to_veloxchem(
                                 DenseMatrix(fock_deriv_atom_i[0]),
                                 basis, molecule).to_numpy()
                                )

    vlx_fock_deriv_atom_i[1] = ( ao_matrix_to_veloxchem(
                                 DenseMatrix(fock_deriv_atom_i[1]),
                                 basis, molecule).to_numpy()
                                )

    vlx_fock_deriv_atom_i[2] = ( ao_matrix_to_veloxchem(
                                 DenseMatrix(fock_deriv_atom_i[2]),
                                 basis, molecule).to_numpy()
                                )

    #if scfdrv.xcfun.get_func_type() == xcfun.lda:
    return vlx_fock_deriv_atom_i + vlx_vxc_deriv_atom_i 
    #else:
    #    return vlx_fock_deriv_atom_i


def vxc_deriv(molecule, basis, scfdrv, unit="au"):
    """
    Imports the derivatives of the Vxc matrix
    from pyscf and converts it to veloxchem format

    :param molecule:
        the vlx molecule object
    :param basis:
        the vlx basis object
    :param density:
        the SCF density matrix (alpha part) in AO basis
    :param scfdrv:
        the ScfDriver
    :param unit:
        the units to be used for the molecular geometry;
        possible values: "au" (default), "Angstrom"

    :return:
        a numpy array of shape natm x 3 x nao x nao
        (nao = number of atomic orbitals)
        corresponding to the derivative of the Vxc matrix
        with respect to the x, y and z coords. of all atoms.
    """
    molecule_string = get_molecule_string(molecule)
    natm = molecule.number_of_atoms()
    basis_set_label = basis.get_label()
    pyscf_basis = translate_to_pyscf(basis_set_label)
    pyscf_molecule = pyscf.gto.M(atom=molecule_string,
                                 basis=pyscf_basis, unit=unit,
                                 charge=molecule.get_charge())
    pyscf_scf = pyscf.scf.RKS(pyscf_molecule)
    pyscf_scf.conv_tol = scfdrv.conv_thresh
    # some functional are parametrized differently in PYSCF
    pyscf_scf.xc = translate_to_pyscf(scfdrv.xcfun.get_func_label())
    if scfdrv.grid_level is None:
        scfdrv.grid_level = get_default_grid_level(scfdrv.xcfun)
    if scfdrv.grid_level == 6 or scfdrv.grid_level == 7:
        pyscf_scf.grids.level = 9
    else:
        pyscf_scf.grids.level = scfdrv.grid_level
    pyscf_scf.kernel()
    pyscf_hessian = pyscf.hessian.rks.Hessian(pyscf_scf)

    density = scfdrv.scf_tensors['D_alpha']
    gs_dm = ao_matrix_to_dalton(DenseMatrix(density), basis, molecule).to_numpy()
    nao = density.shape[0]
    
    if nao != pyscf_molecule.nao:
        error_text = "vlx and pyscf number of atomic orbitals are different!"
        error_text +="\nCheck if the basis sets are defined the same way."
        raise ValueError(error_text)    

    pyscf_mo_coeff = pyscf_scf.mo_coeff
    pyscf_mo_occ = pyscf_scf.mo_occ
    max_memory = 2000
    vxc_deriv1 = pyscf_rks_hessian._get_vxc_deriv1(pyscf_hessian,
                                                   pyscf_mo_coeff,
                                                   pyscf_mo_occ,
                                                   max_memory)

    vlx_vxc_deriv = np.zeros(vxc_deriv1.shape)

    # Transform each compnent (x,y,z) to veloxchem format
    for i in range(natm):
        vlx_vxc_deriv[i,0] = ( ao_matrix_to_veloxchem(
                                       DenseMatrix(vxc_deriv1[i,0]),
                                       basis, molecule).to_numpy()
                                    )

        vlx_vxc_deriv[i,1] = ( ao_matrix_to_veloxchem(
                                     DenseMatrix(vxc_deriv1[i,1]),
                                     basis, molecule).to_numpy()
                                    )

        vlx_vxc_deriv[i,2] = ( ao_matrix_to_veloxchem(
                                     DenseMatrix(vxc_deriv1[i,2]),
                                     basis, molecule).to_numpy()
                                    )

    return vlx_vxc_deriv


def eri_deriv(molecule, basis, i=0, full_deriv=True, unit="au",
              chk_file=None):
    """
    Imports the derivatives of the electron repulsion integrals
    from pyscf and converts them to veloxchem format

    :param molecule:
        the vlx molecule object
    :param basis:
        the vlx basis object
    :param i:
        the index of the atom for which the derivatives
        are computed.
    :param full_deriv:
        True, to compute ( nabla m n | p q )
                       + ( m nabla n | p q )
                       + ( m n | nabla p q )
                       + ( m n | p nabla q )
        False, to compute ( nabla m n | p q ) only.
    :param unit:
        the units to be used for the molecular geometry;
        possible values: "au" (default), "Angstrom"
    :param chk_file:
        the hdf5 checkpoint file name.

    :return:
        a numpy array with the derivatives for atom i
    """
    molecule_string = get_molecule_string(molecule)
    basis_set_label = basis.get_label()

    basis_set_map = basis.get_index_map(molecule)

    pyscf_basis = translate_to_pyscf(basis_set_label)
    pyscf_molecule = pyscf.gto.M(atom=molecule_string,
                                 basis=pyscf_basis, unit=unit,
                                 charge=molecule.get_charge())

    nao = pyscf_molecule.nao

    pyscf_eri_deriv = - pyscf_molecule.intor('int2e_ip1', aosym='s1')
    ao_slices = pyscf_molecule.aoslice_by_atom()

    # Get the AO indeces corresponding to atom i
    ki, kf = ao_slices[i, 2:]

    eri_deriv_atom_i = np.zeros(pyscf_eri_deriv.shape)

    eri_deriv_atom_i[:, ki:kf] = pyscf_eri_deriv[:, ki:kf]

    if full_deriv:
        #   (nabla m n | t p) + ( m nabla n | t p)
        # + (m n | nabla t p) + (m n | t nabla p)
        eri_deriv_atom_i += ( eri_deriv_atom_i.transpose(0,2,1,4,3)
                             + eri_deriv_atom_i.transpose(0,3,4,1,2)
                             + eri_deriv_atom_i.transpose(0,4,3,2,1) )

    vlx_eri_deriv_atom_i = np.zeros(eri_deriv_atom_i.shape)

    for m in range(nao):
        for n in range(nao):
            for t in range(nao):
                for p in range(nao):
                    vm = basis_set_map[m]
                    vn = basis_set_map[n]
                    vt = basis_set_map[t]
                    vp = basis_set_map[p]
                    vlx_eri_deriv_atom_i[:,vm,vn,vt,vp] = eri_deriv_atom_i[:,m,n,t,p]

    if chk_file is not None:
        write_2d_array_hdf5(chk_file, vlx_eri_deriv_atom_i, labels='xyz',
                            atom_index=i) 

    return vlx_eri_deriv_atom_i

def dipole_deriv(molecule, basis, i=0, unit="au"):
    """
    Imports the derivatives of the dipole moment integrals
    from pyscf and converts it to veloxchem format

    :param molecule:
        the vlx molecule object
    :param basis:
        the vlx basis object
    :param i:
        the index of the atom for which the derivatives
        are computed.
    :param unit:
        the units to be used for the molecular geometry;
        possible values: "au" (default), "Angstrom"

    :return:
        a numpy array of shape 3 x 3 x nao x nao
        (nao = number of atomic orbitals)
        corresponding to the derivative of the x, y, and z
        coordinates of the dipole moment integrals
        with respect to the x, y and z coords. of atom i.
    """

    molecule_string = get_molecule_string(molecule)
    basis_set_label = basis.get_label()
    pyscf_basis = translate_to_pyscf(basis_set_label)
    pyscf_molecule = pyscf.gto.M(atom=molecule_string,
                                 basis=pyscf_basis, unit=unit,
                                 charge=molecule.get_charge())
    nao = pyscf_molecule.nao


    pyscf_dipole_deriv = -pyscf_molecule.intor('int1e_irp', aosym='s1').reshape(3, 3, nao, nao)
    ao_slices = pyscf_molecule.aoslice_by_atom()

    # Get the AO indeces corresponding to atom i
    ki, kf = ao_slices[i, 2:]

    dipole_deriv_atom_i = np.zeros(pyscf_dipole_deriv.shape)

    dipole_deriv_atom_i[:,:,:,ki:kf] = pyscf_dipole_deriv[:,:,:,ki:kf]

    # (nabla m | n) + (m | nabla n)
    dipole_deriv_atom_i += dipole_deriv_atom_i.transpose(0,1,3,2)

    vlx_dipole_deriv_atom_i = np.zeros(dipole_deriv_atom_i.shape)


    # Transform each component to veloxchem format
    for c in range(3):
        for x in range(3):
            vlx_dipole_deriv_atom_i[c,x] = ( ao_matrix_to_veloxchem(
                                         DenseMatrix(dipole_deriv_atom_i[c,x]),
                                         basis, molecule).to_numpy()
                                         )

    return vlx_dipole_deriv_atom_i


### Second derivatives

def overlap_second_deriv(molecule, basis, i=0, j=0, unit="au"):
    """
    Imports the second derivatives of the overlap matrix
    from pyscf and converts it to veloxchem format

    :param molecule:
        the vlx molecule object
    :param basis:
        the vlx basis object
    :param i, j:
        the indices of the atoms for which the derivatives
        are computed.
    :param unit:
        the units to be used for the molecular geometry;
        possible values: "au" (default), "Angstrom"

    :return:
        two numpy arrays corresponding to (nabla**2 m | n) and
        (nabla m | nabla n) of shape 3 x 3 x nao x nao
        (nao = number of atomic orbitals)
        corresponding to the derivative of the overlap matrix
        with respect to the x, y and z coords. of atoms i and j.
    """

    molecule_string = get_molecule_string(molecule)
    basis_set_label = basis.get_label()
    pyscf_basis = translate_to_pyscf(basis_set_label)
    pyscf_molecule = pyscf.gto.M(atom=molecule_string,
                                 basis=pyscf_basis, unit=unit,
                                 charge=molecule.get_charge())
    # number of atomic orbitals
    nao = pyscf_molecule.nao

    # (nabla**2 m | n)
    pyscf_ovlp_deriv_aa = pyscf_molecule.intor('int1e_ipipovlp', aosym='s1').reshape(3, 3, nao, nao)
    # (nabla m | nabla n)
    pyscf_ovlp_deriv_ab = pyscf_molecule.intor('int1e_ipovlpip', aosym='s1').reshape(3, 3, nao, nao)

    ao_slices = pyscf_molecule.aoslice_by_atom()

    # Get the AO indeces corresponding to atoms i and j
    ki, kf = ao_slices[i, 2:]
    kj, kg = ao_slices[j, 2:]

    overlap_deriv_atoms_ii = np.zeros(pyscf_ovlp_deriv_aa.shape)
    overlap_deriv_atoms_ij = np.zeros(pyscf_ovlp_deriv_ab.shape)

    overlap_deriv_atoms_ii[:,:,ki:kf] = pyscf_ovlp_deriv_aa[:,:,ki:kf]
    overlap_deriv_atoms_ij[:,:,ki:kf,kj:kg] = pyscf_ovlp_deriv_ab[:,:,ki:kf,kj:kg]

    # add the transpose
    overlap_deriv_atoms_ii += overlap_deriv_atoms_ii.transpose(0,1,3,2)
    overlap_deriv_atoms_ij += overlap_deriv_atoms_ij.transpose(0,1,3,2)

    vlx_ovlp_deriv_atoms_ii = np.zeros(overlap_deriv_atoms_ii.shape)
    vlx_ovlp_deriv_atoms_ij = np.zeros(overlap_deriv_atoms_ij.shape)


    # Transform each component (x,y,z),(x,y,z) to veloxchem format
    for x in range(3):
        for y in range(3):
            vlx_ovlp_deriv_atoms_ii[x,y] = ( ao_matrix_to_veloxchem(
                                         DenseMatrix(overlap_deriv_atoms_ii[x,y]),
                                         basis, molecule).to_numpy()
                                        )
            vlx_ovlp_deriv_atoms_ij[x,y] = ( ao_matrix_to_veloxchem(
                                         DenseMatrix(overlap_deriv_atoms_ij[x,y]),
                                         basis, molecule).to_numpy()
                                        )

    return vlx_ovlp_deriv_atoms_ii, vlx_ovlp_deriv_atoms_ij


def hcore_second_deriv(molecule, basis, i=0, j=0, unit="au"):
    """
    Imports the second derivatives of the core Hamiltonian
    from pyscf and converts it to veloxchem format

    :param molecule:
        the vlx molecule object
    :param basis:
        the vlx basis object
    :param i, j:
        the indices of the atoms for which the derivatives
        are computed.
    :param unit:
        the units to be used for the molecular geometry;
        possible values: "au" (default), "Angstrom"

    :return:
        two numpy arrays corresponding to (nabla**2 m | n) and
        (nabla m | nabla n) of shape 3 x 3 x nao x nao
        (nao = number of atomic orbitals)
        corresponding to the derivative of the overlap matrix
        with respect to the x, y and z coords. of atoms i and j.
    """

    molecule_string = get_molecule_string(molecule)
    basis_set_label = basis.get_label()
    pyscf_basis = translate_to_pyscf(basis_set_label)
    pyscf_molecule = pyscf.gto.M(atom=molecule_string,
                                 basis=pyscf_basis, unit=unit,
                                 charge=molecule.get_charge())
    # number of atomic orbitals
    nao = pyscf_molecule.nao

    # PYSCF SCF and Hessian drivers
    pyscf_scf_drv = pyscf.scf.RHF(pyscf_molecule)
    pyscf_hessian_drv = hessian.rhf.Hessian(pyscf_scf_drv)
    pyscf_hcore_generator = pyscf_hessian_drv.hcore_generator(pyscf_molecule)

    # Call hcore_generator for atoms i and j
    pyscf_hcore_2nd_deriv = pyscf_hcore_generator(i, j)


    vlx_hcore_deriv_atoms_ij = np.zeros(pyscf_hcore_2nd_deriv.shape)


    # Transform each component (x,y,z),(x,y,z) to veloxchem format
    for x in range(3):
        for y in range(3):
            vlx_hcore_deriv_atoms_ij[x,y] = ( ao_matrix_to_veloxchem(
                                         DenseMatrix(pyscf_hcore_2nd_deriv[x,y]),
                                         basis, molecule).to_numpy()
                                        )

    return vlx_hcore_deriv_atoms_ij


def dft_xc_second_deriv(molecule, basis, scf_drv, unit="au"):
    """
    Imports the second derivatives of the exchange and correlation
    energy from pyscf.

    :param molecule:
        the vlx molecule object
    :param basis:
        the vlx basis object
    :param scf_drv:
        the vlx Scf Driver.
    :param unit:
        the units to be used for the molecular geometry;
        possible values: "au" (default), "Angstrom"

    :return:
        a numpy array of shape (natm, natm, 3, 3) corresponding to the
        partial derivative of Exc with respect to the x, y and z coords.
        of all atoms.
    """
    molecule_string = get_molecule_string(molecule)
    basis_set_label = basis.get_label()
    pyscf_basis = translate_to_pyscf(basis_set_label)
    pyscf_molecule = pyscf.gto.M(atom=molecule_string,
                                 basis=pyscf_basis, unit=unit,
                                 charge=molecule.get_charge())

    pyscf_scf = pyscf.scf.RKS(pyscf_molecule)
    # set functional, convergence threshold and grid level 
    pyscf_scf.xc = translate_to_pyscf(scf_drv.xcfun.get_func_label())
    pyscf_scf.conv_tol = scf_drv.conv_thresh
    if scf_drv.grid_level == 6 or scf_drv.grid_level == 7:
        pyscf_scf.grids.level = 9
    else:
        pyscf_scf.grids.level = scf_drv.grid_level

    pyscf_scf.kernel()

    pyscf_hessian = pyscf.hessian.rks.Hessian(pyscf_scf)

    pyscf_mo_coeff = pyscf_scf.mo_coeff
    pyscf_mo_occ = pyscf_scf.mo_occ
    max_memory = 2000
    pyscf_vxc_deriv2 = pyscf_rks_hessian._get_vxc_deriv2(pyscf_hessian,
                                                pyscf_mo_coeff,
                                                pyscf_mo_occ, max_memory)
    pyscf_vxc_diag = pyscf_rks_hessian._get_vxc_diag(pyscf_hessian,
                                                pyscf_mo_coeff,
                                                pyscf_mo_occ, max_memory)

    natm = pyscf_molecule.natm
    exc_hessian_contrib = np.zeros((natm, natm, 3, 3))
    ao_slices = pyscf_molecule.aoslice_by_atom()

    density = 2 * scf_drv.scf_tensors['D_alpha']
    gs_dm = ao_matrix_to_dalton(DenseMatrix(density),
                                basis, molecule).to_numpy()

    #print("Shape of vxc_deriv2: ", pyscf_vxc_deriv2.shape)
    for i in range(natm):
        # Get the AO indeces corresponding to atoms i and j
        ki, kf = ao_slices[i, 2:]

        exc_hessian_contrib[i, i] += 2 * np.einsum('xypq,pq->xy',
                                        pyscf_vxc_diag[:,:,ki:kf],
                                        gs_dm[ki:kf])

        for j in range(i+1):
            li, lf = ao_slices[j, 2:]
            exc_hessian_contrib[i, j] += 2 * np.einsum('xypq,pq->xy',
                                        pyscf_vxc_deriv2[i,:,:,li:lf],
                                        gs_dm[li:lf])

        for j in range(i):
            exc_hessian_contrib[j, i] = exc_hessian_contrib[i, j].T

    return exc_hessian_contrib 


def eri_second_deriv(molecule, basis, i=0, j=0, unit="au"):
    """
    Imports the second derivatives of the ERI tensor
    from pyscf and converts it to veloxchem format

    :param molecule:
        the vlx molecule object
    :param basis:
        the vlx basis object
    :param i, j:
        the indices of the atoms for which the derivatives
        are computed.
    :param unit:
        the units to be used for the molecular geometry;
        possible values: "au" (default), "Angstrom"

    :return:
        two numpy arrays corresponding to (nabla**2 m | n) and
        (nabla m | nabla n) of shape 3 x 3 x nao x nao
        (nao = number of atomic orbitals)
        corresponding to the derivative of the overlap matrix
        with respect to the x, y and z coords. of atoms i and j.
    """

    molecule_string = get_molecule_string(molecule)
    basis_set_label = basis.get_label()
    pyscf_basis = translate_to_pyscf(basis_set_label)
    pyscf_molecule = pyscf.gto.M(atom=molecule_string,
                                 basis=pyscf_basis, unit=unit,
                                 charge=molecule.get_charge())
    # number of atomic orbitals
    nao = pyscf_molecule.nao

    basis_set_map = basis.get_index_map(molecule)

    # terms with nabla squared
    pyscf_eri_deriv_ii = pyscf_molecule.intor('int2e_ipip1', aosym='s1').reshape(
                                3, 3, nao, nao, nao, nao)

    # terms with two nablas on same side
    pyscf_eri_deriv_ij = pyscf_molecule.intor('int2e_ipvip1', aosym='s1').reshape(
                                3, 3, nao, nao, nao, nao)

    # terms with two nablas on different sides
    pyscf_eri_deriv_12 = pyscf_molecule.intor('int2e_ip1ip2', aosym='s1').reshape(
                                3, 3, nao, nao, nao, nao)

    ao_slices = pyscf_molecule.aoslice_by_atom()

    # Get the AO indeces corresponding to atoms i and j
    ki, kf = ao_slices[i, 2:]
    kj, kg = ao_slices[j, 2:]

    eri_deriv_atom_ii = np.zeros(pyscf_eri_deriv_ii.shape)
    eri_deriv_atom_ij = np.zeros(pyscf_eri_deriv_ij.shape)
    eri_deriv_atom_12 = np.zeros(pyscf_eri_deriv_12.shape)

    #   (nabla**2 m n | t p) + ( m nabla**2 n | t p)
    # + (m n | nabla**2 t p) + (m n | t nabla**2 p)
    eri_deriv_atom_ii[:, :, ki:kf] = pyscf_eri_deriv_ii[:, :, ki:kf]
    eri_deriv_atom_ii += ( eri_deriv_atom_ii.transpose(0,1,3,2,5,4)
                         + eri_deriv_atom_ii.transpose(0,1,4,5,2,3)
                         + eri_deriv_atom_ii.transpose(0,1,5,4,3,2) )

    #   (nabla_i m nabla_j n | t p) + (nabla_j m nabla_i n | t p)
    # + (m n | nabla_i t nabla_j p) + (m n | nabla_j t nabla_i p)
    eri_deriv_atom_ij[:, :, ki:kf, kj:kg] = pyscf_eri_deriv_ij[:, :, ki:kf, kj:kg]
    eri_deriv_atom_ij += ( eri_deriv_atom_ij.transpose(0,1,3,2,5,4)
                         + eri_deriv_atom_ij.transpose(0,1,4,5,2,3)
                         + eri_deriv_atom_ij.transpose(0,1,5,4,3,2) )

    #   (nabla_i m n | nabla_j t p) + (m nabla_i n | nabla_j t p)
    # + (nabla_i m n | t nabla_j p) + (m nabla_i n | t nabla_j p)
    # + (nabla_j m n | nabla_i t p) + (m nabla_j n | nabla_i t p)
    # + (nabla_j m n | t nabla_i p) + (m nabla_j n | t nabla_i p)
    eri_deriv_atom_12[:, :, ki:kf, :, kj:kg] = pyscf_eri_deriv_12[:, :, ki:kf, :, kj:kg]
    eri_deriv_atom_12 += ( eri_deriv_atom_12.transpose(0,1,3,2,4,5)
                         + eri_deriv_atom_12.transpose(0,1,2,3,5,4)
                         + eri_deriv_atom_12.transpose(0,1,3,2,5,4)
                         + eri_deriv_atom_12.transpose(0,1,4,5,2,3)
                         + eri_deriv_atom_12.transpose(0,1,5,4,2,3)
                         + eri_deriv_atom_12.transpose(0,1,4,5,3,2)
                         + eri_deriv_atom_12.transpose(0,1,5,4,3,2) )

    vlx_eri_deriv_atom_ii = np.zeros(eri_deriv_atom_ii.shape)
    vlx_eri_deriv_atom_ij = np.zeros(eri_deriv_atom_ij.shape)
    vlx_eri_deriv_atom_12 = np.zeros(eri_deriv_atom_12.shape)


    # Transform each component (x,y,z),(x,y,z) to veloxchem format
    for m in range(nao):
        for n in range(nao):
            for t in range(nao):
                for p in range(nao):
                    vm = basis_set_map[m]
                    vn = basis_set_map[n]
                    vt = basis_set_map[t]
                    vp = basis_set_map[p]
                    vlx_eri_deriv_atom_ii[:,:,vm,vn,vt,vp] = eri_deriv_atom_ii[:,:,m,n,t,p]
                    vlx_eri_deriv_atom_ij[:,:,vm,vn,vt,vp] = eri_deriv_atom_ij[:,:,m,n,t,p]
                    vlx_eri_deriv_atom_12[:,:,vm,vn,vt,vp] = eri_deriv_atom_12[:,:,m,n,t,p]

    return vlx_eri_deriv_atom_ii, vlx_eri_deriv_atom_ij + vlx_eri_deriv_atom_12

