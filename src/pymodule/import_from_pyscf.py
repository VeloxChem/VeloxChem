# import veloxchem as vlx
from .veloxchemlib import OverlapIntegralsDriver
from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import DenseMatrix
from .veloxchemlib import ao_matrix_to_veloxchem 
import numpy as np
import sys

try:
    import pyscf
    from pyscf import grad
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
    else:
        return label

def overlap_deriv(molecule, basis, molecule_string, basis_set_label,
                  i=0, unit="au"):
    """
    Imports the derivatives of the overlap matrix
    from pyscf and converts it to veloxchem format
   
    :param molecule:
        the vlx molecule object
    :param basis:
        the vlx basis object 
    :param molecule_string:
        the string which defines the molecular structure
        in xyz format
    :param basis_set_label:
        the label describing the basis set;
    :param i:
        the index of the atom for which the derivatives
        are computed.
    :param unit:
        the units to be used for the molecular geometry;
        possible values: "au" (default), "Angstrom"

    :return:
        a numpy array of shape 3 x nao x nao
        (nao = number of atomic orbitals)
        corresponding to the derivative of the overlap matrix
        with repsect to the x, y and z coords. of atom i.
    """

    pyscf_basis = translate_to_pyscf(basis_set_label)
    pyscf_molecule = pyscf.gto.M(atom=molecule_string,
                                 basis=pyscf_basis, unit=unit)
    
    # The "-" sign is due to the fact that pyscf computes -(nabla m | n)
    pyscf_ovlp_deriv = - pyscf_molecule.intor('int1e_ipovlp', aosym='s1')
    ao_slices = pyscf_molecule.aoslice_by_atom()
    
    # Get the AO indeces corresponding to atom i 
    ki, kf = ao_slices[i, 2:]

    overlap_deriv_atom_i = np.zeros(pyscf_ovlp_deriv.shape)

    overlap_deriv_atom_i[:,ki:kf] = pyscf_ovlp_deriv[:,ki:kf]
    
    # (nabla m | n) + (m | nabla m)
    overlap_deriv_atom_i += overlap_deriv_atom_i.transpose(0,2,1)

    vlx_ovlp_deriv_atom_i = np.zeros(overlap_deriv_atom_i.shape)


    # Transform each compnent (x,y,z) to veloxchem format    
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

    return vlx_ovlp_deriv_atom_i

    

def fock_deriv(molecule, basis, molecule_string, basis_set_label,
               i=0, unit="au"):
    """
    Imports the derivatives of the Fock matrix
    from pyscf and converts it to veloxchem format
   
    :param molecule:
        the vlx molecule object
    :param basis:
        the vlx basis object 
    :param molecule_string:
        the string which defines the molecular structure
        in xyz format
    :param basis_set_label:
        the label describing the basis set;
    :param i:
        the index of the atom for which the derivatives
        are computed.
    :param unit:
        the units to be used for the molecular geometry;
        possible values: "au" (default), "Angstrom"

    :return:
        a numpy array of shape 3 x nao x nao
        (nao = number of atomic orbitals)
        corresponding to the derivative of the Fock matrix
        with repsect to the x, y and z coords. of atom i.
    """
    pyscf_basis = translate_to_pyscf(basis_set_label)
    pyscf_molecule = pyscf.gto.M(atom=molecule_string,
                                 basis=pyscf_basis, unit=unit)
    pyscf_scf = pyscf.scf.RHF(pyscf_molecule)
    pyscf_scf.kernel()
    nocc = pyscf_molecule.nelec[0]
    gs_dm = np.einsum('mi,ni->mn', pyscf_scf.mo_coeff[:,:nocc],
                       pyscf_scf.mo_coeff[:,:nocc])
    nao = pyscf_scf.mo_coeff.shape[0]
    pyscf_grad = grad.RHF(pyscf_scf)
    
    # The "-" sign is due to the fact that pyscf computes -(nabla m n | t p)
    pyscf_eri_deriv = - pyscf_molecule.intor('int2e_ip1', aosym='s1')
    ao_slices = pyscf_molecule.aoslice_by_atom()
    hcore_generator = pyscf_grad.hcore_generator(pyscf_molecule)
    
    # Get the AO indeces corresponding to atom i 
    ki, kf = ao_slices[i, 2:]

    eri_deriv_atom_i = np.zeros(pyscf_eri_deriv.shape)
    fock_deriv_atom_i = np.zeros((3,nao,nao))

    #   (nabla m n | t p) + ( m nabla n | t p) 
    # + (m n | nabla t p) + (m n | t nabla p)
    eri_deriv_atom_i[:, ki:kf] = pyscf_eri_deriv[:, ki:kf]
    eri_deriv_atom_i += ( eri_deriv_atom_i.transpose(0,2,1,4,3) 
                         + eri_deriv_atom_i.transpose(0,3,4,1,2)
                         + eri_deriv_atom_i.transpose(0,4,3,2,1) )
    
    fock_deriv_atom_i = (  hcore_generator(i) 
                         + 2 * np.einsum('xmntf,tf->xmn',
                                          eri_deriv_atom_i, gs_dm) 
                         - np.einsum('xmtnf,tf->xmn',
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

    return vlx_fock_deriv_atom_i


def eri_deriv(molecule, basis, molecule_string, basis_set_label,
               i=0, unit="au"):
    """
    Imports the derivatives of the electron repulsion integrals
    from pyscf and converts them to veloxchem format
   
    :param molecule:
        the vlx molecule object
    :param basis:
        the vlx basis object 
    :param molecule_string:
        the string which defines the molecular structure
        in xyz format
    :param basis_set_label:
        the label describing the basis set;
    :param i:
        the index of the atom for which the derivatives
        are computed.
    :param unit:
        the units to be used for the molecular geometry;
        possible values: "au" (default), "Angstrom"

    :return:
        a numpy array with the derivatives for atom i
    """

    basis_set_map = basis.get_index_map(molecule)

    pyscf_basis = translate_to_pyscf(basis_set_label)
    pyscf_molecule = pyscf.gto.M(atom=molecule_string,
                                 basis=pyscf_basis, unit=unit)

    nao = pyscf_molecule.nao

    pyscf_eri_deriv = - pyscf_molecule.intor('int2e_ip1', aosym='s1')
    ao_slices = pyscf_molecule.aoslice_by_atom()
    
    # Get the AO indeces corresponding to atom i 
    ki, kf = ao_slices[i, 2:]

    eri_deriv_atom_i = np.zeros(pyscf_eri_deriv.shape)

    #   (nabla m n | t p) + ( m nabla n | t p) 
    # + (m n | nabla t p) + (m n | t nabla p)
    eri_deriv_atom_i[:, ki:kf] = pyscf_eri_deriv[:, ki:kf]
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
    
    return vlx_eri_deriv_atom_i



