# import veloxchem as vlx
from .molecule import Molecule
from .veloxchemlib import OverlapIntegralsDriver
from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import DenseMatrix
from .veloxchemlib import ao_matrix_to_veloxchem
from .veloxchemlib import ao_matrix_to_dalton
from .veloxchemlib import ElectricDipoleIntegralsDriver
import numpy as np
import sys
import h5py

try:
    import pyscf
    from pyscf import grad
    from pyscf import hessian
    from pyscf.grad import rks as pyscf_rks_grad
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

def get_molecule_string(molecule):
    mol_string = ""
    for i in range(molecule.number_of_atoms()):
        mol_string += molecule.get_labels()[i] + " "
        for j in range(3):
            mol_string += str(molecule.get_coordinates()[i][j]) + " "

        mol_string += "\n"

    return mol_string

def write_2d_array_hdf5(fname, arrays, labels='xyz', atom_index=0):
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
        full_label= str(atom_index) + label
        hf.create_dataset(full_label, data=array, compression='gzip')
    hf.close()
    return True


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
                                 basis=pyscf_basis, unit=unit)

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
        write_2d_array_hdf5(chk_file, vlx_ovlp_deriv_atom_i, labels='xyz',
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
                                 basis=pyscf_basis, unit=unit)
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
                                 basis=pyscf_basis, unit=unit)

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
                                 basis=pyscf_basis, unit=unit)

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
                                 basis=pyscf_basis, unit=unit)

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
                                 basis=pyscf_basis, unit=unit)
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


def fock_deriv(molecule, basis, density, i=0, unit="au"):
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
                                 basis=pyscf_basis, unit=unit)
    pyscf_scf = pyscf.scf.RHF(pyscf_molecule)
    #pyscf_scf.kernel()
    #nocc = pyscf_molecule.nelec[0]
    gs_dm = ao_matrix_to_dalton(DenseMatrix(density), basis, molecule).to_numpy()
    #gs_dm = np.einsum('mi,ni->mn', pyscf_scf.mo_coeff[:,:nocc],
    #                   pyscf_scf.mo_coeff[:,:nocc])
    nao = density.shape[0]
    
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


def vxc_deriv(molecule, basis, density, xcfun, i=0, unit="au"):
    """
    Imports the derivatives of the Vxc matrix
    from pyscf and converts it to veloxchem format

    :param molecule:
        the vlx molecule object
    :param basis:
        the vlx basis object
    :param density:
        the SCF density matrix (alpha part) in AO basis
    :param xcfun:
        the label of the exchange-correlation functional
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
        with respect to the x, y and z coords. of atom i.
    """
    molecule_string = get_molecule_string(molecule)
    basis_set_label = basis.get_label()
    pyscf_basis = translate_to_pyscf(basis_set_label)
    pyscf_molecule = pyscf.gto.M(atom=molecule_string,
                                 basis=pyscf_basis, unit=unit)
    pyscf_scf = pyscf.scf.RKS(pyscf_molecule)

    # TODO: some functional are parametrized differently in PYSCF
    # TODO: adjust the grid level -- default used right now.
    pyscf_scf.xc = translate_to_pyscf(xcfun)
    #pyscf_scf.kernel()
    #nocc = pyscf_molecule.nelec[0]
    gs_dm = ao_matrix_to_dalton(DenseMatrix(density), basis, molecule).to_numpy()
    #gs_dm = np.einsum('mi,ni->mn', pyscf_scf.mo_coeff[:,:nocc],
    #                   pyscf_scf.mo_coeff[:,:nocc])
    nao = density.shape[0]
    
    if nao != pyscf_molecule.nao:
        error_text = "vlx and pyscf number of atomic orbitals are different!"
        error_text +="\nCheck if the basis sets are defined the same way."
        raise ValueError(error_text)    

    pyscf_grad = pyscf_scf.Gradients()
    pyscf_mf = pyscf_grad.base 
    pyscf_ni = pyscf_mf._numint # Pyscf object for numerical integration
    if pyscf_grad.grids is not None:
        pyscf_grids = pyscf_grad.grids
    else:
        pyscf_grids = pyscf_mf.grids
    if pyscf_grids.coords is None:
        pyscf_grids.build(with_non0tab=True)

    pyscf_exc, pyscf_vxc_deriv = pyscf_rks_grad.get_vxc(pyscf_ni,
                                             pyscf_molecule,
                                             pyscf_grids, pyscf_mf.xc,
                                             2*gs_dm)

    ao_slices = pyscf_molecule.aoslice_by_atom()

    # Get the AO indices corresponding to atom i
    ki, kf = ao_slices[i, 2:]

    vxc_deriv_atom_i = np.zeros((3,nao,nao))

    vxc_deriv_atom_i[:,ki:kf] = pyscf_vxc_deriv[:,ki:kf]

    # (nabla m | vxc | n) + (m | vxc |nabla n)
    vxc_deriv_atom_i += vxc_deriv_atom_i.transpose(0,2,1)

    vlx_vxc_deriv_atom_i = np.zeros(vxc_deriv_atom_i.shape)


    # Transform each compnent (x,y,z) to veloxchem format
    vlx_vxc_deriv_atom_i[0] = ( ao_matrix_to_veloxchem(
                                 DenseMatrix(vxc_deriv_atom_i[0]),
                                 basis, molecule).to_numpy()
                                )

    vlx_vxc_deriv_atom_i[1] = ( ao_matrix_to_veloxchem(
                                 DenseMatrix(vxc_deriv_atom_i[1]),
                                 basis, molecule).to_numpy()
                                )

    vlx_vxc_deriv_atom_i[2] = ( ao_matrix_to_veloxchem(
                                 DenseMatrix(vxc_deriv_atom_i[2]),
                                 basis, molecule).to_numpy()
                                )

    return vlx_vxc_deriv_atom_i


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
                                 basis=pyscf_basis, unit=unit)

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

## Not working correctly
##def dipole_deriv(molecule, basis, i=0, unit="au"):
##    """
##    Imports the derivatives of the dipole moment integrals
##    from pyscf and converts it to veloxchem format
##
##    :param molecule:
##        the vlx molecule object
##    :param basis:
##        the vlx basis object
##    :param i:
##        the index of the atom for which the derivatives
##        are computed.
##    :param unit:
##        the units to be used for the molecular geometry;
##        possible values: "au" (default), "Angstrom"
##
##    :return:
##        a numpy array of shape 3 x 3 x nao x nao
##        (nao = number of atomic orbitals)
##        corresponding to the derivative of the x, y, and z
##        coordinates of the dipole moment integrals
##        with respect to the x, y and z coords. of atom i.
##    """
##
##    molecule_string = get_molecule_string(molecule)
##    basis_set_label = basis.get_label()
##    pyscf_basis = translate_to_pyscf(basis_set_label)
##    pyscf_molecule = pyscf.gto.M(atom=molecule_string,
##                                 basis=pyscf_basis, unit=unit)
##    nao = pyscf_molecule.nao
##
##
##    pyscf_dipole_deriv = pyscf_molecule.intor('int1e_irp', aosym='s1').reshape(3, 3, nao, nao)
##    ao_slices = pyscf_molecule.aoslice_by_atom()
##
##    # Get the AO indeces corresponding to atom i
##    #ki, kf = ao_slices[i, 2:]
##
##    dipole_deriv_atom_i = np.zeros(pyscf_dipole_deriv.shape)
##
##    #dipole_deriv_atom_i[:,:,ki:kf] = pyscf_dipole_deriv[:,:,ki:kf]
##    dipole_deriv_atom_i = pyscf_dipole_deriv
##
##    # (nabla m | n) + (m | nabla n)
##    dipole_deriv_atom_i += dipole_deriv_atom_i.transpose(1,0,3,2)
##
##    #vlx_dipole_deriv_atom_i = np.zeros_like(dipole_deriv_atom_i.reshape(3,3,nao,nao))
##    vlx_dipole_deriv_atom_i = np.zeros(dipole_deriv_atom_i.shape)
##
##
##    # Transform each component to veloxchem format
##    for c in range(3):
##        for x in range(3):
##            vlx_dipole_deriv_atom_i[c,x] = ( ao_matrix_to_veloxchem(
##                                         DenseMatrix(dipole_deriv_atom_i[c,x]),
##                                         basis, molecule).to_numpy()
##                                         )
##
##    return vlx_dipole_deriv_atom_i


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
                                 basis=pyscf_basis, unit=unit)
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
                                 basis=pyscf_basis, unit=unit)
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
                                 basis=pyscf_basis, unit=unit)
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

