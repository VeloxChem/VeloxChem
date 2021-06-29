#
#                           VELOXCHEM 1.0-RC2
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
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

from mpi4py import MPI
import numpy as np
import h5py
import sys
import os
from collections import Counter

from .lrsolver import LinearResponseSolver
from .inputparser import InputParser
from .outputstream import OutputStream
from .veloxchemlib import get_basis_function_indices_for_atom
from .veloxchemlib import ElectricDipoleIntegralsDriver
from .veloxchemlib import mpi_master
from .veloxchemlib import bohr_in_angstroms
from .errorhandler import assert_msg_critical


class LoPropDriver:
    """
    Implements the LoProp driver.
    
    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes the LoProp driver.
        """

        if comm is None:
            comm = MPI.COMM_WORLD
        if ostream is None:
            ostream = OutputStream(sys.stdout)

        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()
        self.ostream = ostream

    def compute(self, molecule, basis, scf_tensors):
        """
        calculate the loprop transformation matrix T         
        patial charge (Qab) and localised polarisabilities
                
        :param molecule:
            the molecule
        :param basis:
            the bais functions 
        :param scf_tensors:
            The tensors from the converged SCF calculation.
            
        :return:
            A dictionary containing localized charges and localized
            polarizabilities.
        """

        natoms = molecule.number_of_atoms()

        S = scf_tensors['S']
        C = scf_tensors['C_alpha']
        D = scf_tensors['D_alpha'] + scf_tensors['D_beta']

        # number of orbitals
        n_ao = C.shape[0]
        n_mo = C.shape[1]

        # re-arrange AOs according to principal quantum number
        re_arranged_indices = []
        for i in range(molecule.number_of_atoms()):
            indices, angmoms = get_basis_function_indices_for_atom(
                molecule, basis, i)
            re_arranged_indices.append(indices)
        re_arranged_indices = np.concatenate(re_arranged_indices)

        # obtain occupied & virtual orbital lists
        ao_per_atom, ao_occ, ao_vir = self.get_ao_indices(molecule, basis)

        # TO transforation: re-arrange S
        S0 = S[re_arranged_indices, :][:, re_arranged_indices]

        # T1 transformation: Gram-Schmidt
        T1 = np.zeros((n_ao, n_ao))
        ri = 0
        for a in range(natoms):
            rf = ri + ao_per_atom[a]
            s = S0[ri:rf, ri:rf]
            L = np.linalg.cholesky(s)
            t = np.linalg.inv(L.T)
            T1[ri:rf, ri:rf] = t
            ri += ao_per_atom[a]
        S1 = np.linalg.multi_dot([T1.T, S0, T1])

        # T2 transformation: Lodwin occupied and virtual
        T2 = np.zeros((n_ao, n_ao))
        n_occ_ao = len(ao_occ)
        n_vir_ao = len(ao_vir)
        # occupied
        s = S1[ao_occ, :][:, ao_occ]
        t = self.lowdin_orthonormalize(s)
        for r in range(n_occ_ao):
            for c in range(n_occ_ao):
                T2[ao_occ[r], ao_occ[c]] = t[r, c]
        # virtural
        s = S1[ao_vir, :][:, ao_vir]
        t = self.lowdin_orthonormalize(s)
        for r in range(n_vir_ao):
            for c in range(n_vir_ao):
                T2[ao_vir[r], ao_vir[c]] = t[r, c]
        S2 = np.linalg.multi_dot([T2.T, S1, T2])

        # T3: projection to virtual
        T3 = np.identity(n_ao)
        # selected the overlap between occupied and virtual
        T3_ov = S2[ao_occ, :][:, ao_vir]
        # projection
        for ao_o_ind, ao_o in enumerate(ao_occ):
            for ao_v_ind, ao_v in enumerate(ao_vir):
                T3[ao_o, ao_v] = -T3_ov[ao_o_ind, ao_v_ind]
        S3 = np.linalg.multi_dot([T3.T, S2, T3])

        # T4: lodwin virtural
        T4_virtual = S3[ao_vir, :][:, ao_vir]
        T4_virtual = self.lowdin_orthonormalize(T4_virtual)
        T4 = np.identity(n_ao)
        for ao_v_ind, ao_v in enumerate(ao_vir):
            for ao_v_1_ind, ao_v_1 in enumerate(ao_vir):
                T4[ao_v, ao_v_1] = T4_virtual[ao_v_ind, ao_v_1_ind]

        # total transformation T becomes:
        T = np.linalg.multi_dot([T1, T2, T3, T4])

        # obtain density matrix D in loprop basis set
        T_inv = np.linalg.pinv(T)
        D = D[re_arranged_indices, :][:, re_arranged_indices]
        D_loprop = np.linalg.multi_dot([T_inv, D, T_inv.T])

        # calculated localised charges
        Qab = np.zeros((natoms, natoms))
        start_point = 0
        for atom in range(natoms):
            select = np.arange(start_point, start_point + ao_per_atom[atom])
            Qab[atom, atom] = -np.trace(D_loprop[select, :][:, select])
            start_point += ao_per_atom[atom]

        # nuclear charge
        nuclear_charge = molecule.elem_ids_to_numpy()
        for a in range(natoms):
            Qab[a, a] = Qab[a, a] + nuclear_charge[a]

        # solve linear response
        lrs_drv = LinearResponseSolver(self.comm, self.ostream)
        lrs_out = lrs_drv.compute(molecule, basis, scf_tensors)

        # obtain response vectors
        Nx = lrs_out['solutions'][('x', 0)]
        Ny = lrs_out['solutions'][('y', 0)]
        Nz = lrs_out['solutions'][('z', 0)]

        # unpact response vectors to matrix form
        nocc = molecule.number_of_alpha_electrons()
        norb = n_mo
        kappa_x = self.lrvec2mat(Nx, nocc, norb)
        kappa_y = self.lrvec2mat(Ny, nocc, norb)
        kappa_z = self.lrvec2mat(Nz, nocc, norb)

        # perturbed densities
        # factor of 2 from spin-adapted excitation vectors
        Dx = 2 * np.linalg.multi_dot([C, kappa_x, C.T])
        Dy = 2 * np.linalg.multi_dot([C, kappa_y, C.T])
        Dz = 2 * np.linalg.multi_dot([C, kappa_z, C.T])

        # convert purterbed density to loprop basis
        Dx = Dx[re_arranged_indices, :][:, re_arranged_indices]
        Dy = Dy[re_arranged_indices, :][:, re_arranged_indices]
        Dz = Dz[re_arranged_indices, :][:, re_arranged_indices]
        Dx_loprop = np.linalg.multi_dot([T_inv, Dx, T_inv.T])
        Dy_loprop = np.linalg.multi_dot([T_inv, Dy, T_inv.T])
        Dz_loprop = np.linalg.multi_dot([T_inv, Dz, T_inv.T])
        Dk_loprop = np.concatenate(([Dx_loprop], [Dy_loprop], [Dz_loprop]))

        # dipole
        dipole_drv = ElectricDipoleIntegralsDriver(self.comm)
        dipole_mats = dipole_drv.compute(molecule, basis)
        x_ao = dipole_mats.x_to_numpy()
        y_ao = dipole_mats.y_to_numpy()
        z_ao = dipole_mats.z_to_numpy()

        # convert to loprop basis
        x_ao = x_ao[re_arranged_indices, :][:, re_arranged_indices]
        y_ao = y_ao[re_arranged_indices, :][:, re_arranged_indices]
        z_ao = z_ao[re_arranged_indices, :][:, re_arranged_indices]
        x_ao_loprop = np.linalg.multi_dot([T.T, x_ao, T])
        y_ao_loprop = np.linalg.multi_dot([T.T, y_ao, T])
        z_ao_loprop = np.linalg.multi_dot([T.T, z_ao, T])
        dipole = np.concatenate(([x_ao_loprop], [y_ao_loprop], [z_ao_loprop]))

        # localised polarisabilities:
        #
        #    aAB = -rAB delta DAB + dQab (Ra-Rb)
        #
        #    where dQab is obtained by solving the following Lagragian with
        #    minimal charge transfer, here a penaty function is also introduced
        #
        #    contribution from localised dipoles: rAB*delta DAB = Aab
        #    bond contribution:  dQab (Ra-Rb)

        # dQa:charge shift per atom
        dQa = np.zeros((natoms, 3))
        # dQa array as [[Ax,Ay,Az],[Bx,By,Bz]]
        it = 0
        for a in range(natoms):
            select = np.arange(it, it + ao_per_atom[a])
            dQa[a][0] = np.trace(Dx_loprop[select, :][:, select])
            dQa[a][1] = np.trace(Dy_loprop[select, :][:, select])
            dQa[a][2] = np.trace(Dz_loprop[select, :][:, select])
            it += ao_per_atom[a]

        # coord_matrix, the rab matrix in equation above
        molecule_coord = molecule.get_coordinates()
        coord_matrix = np.zeros((natoms, natoms, 3))
        for i in range(natoms):
            for j in range(natoms):
                if i == j:
                    # a==b: rab=ra
                    coord_matrix[i][j] = molecule_coord[i]
                else:
                    # a!=b: rab = (ra-rb)/2
                    a = np.abs(molecule_coord[i] - molecule_coord[j])
                    coord_matrix[i][j] = a / 2

        # contribution from localised dipole
        loc_dipole = np.zeros((3, 3, natoms, natoms))
        for i in range(3):
            for j in range(3):
                it_a = 0
                for a in range(natoms):
                    select_a = np.arange(it_a, it_a + ao_per_atom[a])
                    it_b = 0
                    for b in range(natoms):
                        # select the subblock[a][b] region in dks_lp
                        select_b = np.arange(it_b, it_b + ao_per_atom[b])
                        # selected the lp basis for subblock[A][B] in purterbed density matrix
                        D_AB = Dk_loprop[i][select_a, :][:, select_b]
                        # selected the dipole matrice for subblock[A][B] in purterbed density matrix
                        dipole_select = dipole[j][select_a, :][:, select_b]
                        dipole_select = dipole_select.transpose()

                        loc_dipole[i, j, a, b] += np.trace(
                            np.matmul(dipole_select, D_AB))
                        it_b += ao_per_atom[b]

                    loc_dipole[i, j, a, a] -= dQa[a, j] * coord_matrix[a, a, i]
                    it_a += ao_per_atom[a]

        # Lagragian
        Fab = np.zeros((natoms, natoms))
        for a in range(natoms):
            za = nuclear_charge[a]
            Ra = molecule_coord[a]
            for b in range(a):
                zb = nuclear_charge[b]
                Rb = molecule_coord[b]
                Fab[a, b] = self.penalty_fc(za, Ra, zb, Rb)
                Fab[b, a] = Fab[a][b]
            for a in range(natoms):
                Fab[a, a] += -sum(Fab[a, :])

        Lab = Fab + 2 * np.max(np.abs(Fab))

        #dQa = np.swapaxes(dQa, 0, 1)
        lagragian = [np.linalg.solve(Lab, rhs) for rhs in dQa.T]

        #dQab
        dQab = np.zeros((natoms, natoms, 3))
        for i in range(3):
            for a in range(natoms):
                za = nuclear_charge[a]
                Ra = molecule_coord[a]
                for b in range(a):
                    zb = nuclear_charge[b]
                    Rb = molecule_coord[b]
                    dQab[a, b, i] = -(lagragian[i][a] -
                                      lagragian[i][b]) * self.penalty_fc(
                                          za, Ra, zb, Rb)
                    dQab[b, a, i] -= dQab[a, b, i]

        # dRab matrix: mid point of Ra-Rb
        dRab = np.zeros((natoms, natoms, 3))

        for a in range(natoms):
            for b in range(natoms):
                for i in range(3):
                    dRab[a][b][i] = (molecule_coord[a][i] -
                                     molecule_coord[b][i]) / 2
                    dRab[b][a][i] = -dRab[a][b][i]

        # charge transfer from bond polarisability
        bond_polarisabilities = np.zeros((3, 3, natoms, natoms))
        for a in range(natoms):
            for b in range(natoms):
                for i in range(3):
                    for j in range(3):
                        bond_polarisabilities[
                            i, j, a, b] = dRab[a, b, i] * dQab[a, b, j] + dRab[
                                a, b, j] * dQab[a, b, i]

        local_polarisabilities = loc_dipole + 0.5 * bond_polarisabilities

        # molecular polarisabilities
        molecule_polarisabilities = (loc_dipole +
                                     0.5 * bond_polarisabilities).sum(
                                         axis=3).sum(axis=2)

        # obtain atom polarisabilities
        atom_polarisabilities = np.zeros((natoms, 6))
        for i in range(natoms):
            # in the order xx, xy,xz, yy,yz, zz
            atom_polarisabilities[i, 0] = np.sum(
                local_polarisabilities[0, 0, i, :] +
                local_polarisabilities[0, 0, :, i] +
                (natoms - 2) * local_polarisabilities[0, 0, i, i]) / (natoms)

            atom_polarisabilities[i, 1] = np.sum(
                local_polarisabilities[0, 1, i, :] +
                local_polarisabilities[0, 1, :, i] +
                (natoms - 2) * local_polarisabilities[0, 1, i, i]) / (natoms)

            atom_polarisabilities[i, 2] = np.sum(
                local_polarisabilities[0, 2, i, :] +
                local_polarisabilities[0, 2, :, i] +
                (natoms - 2) * local_polarisabilities[0, 2, i, i]) / (natoms)

            atom_polarisabilities[i, 3] = np.sum(
                local_polarisabilities[1, 1, i, :] +
                local_polarisabilities[1, 1, :, i] +
                (natoms - 2) * local_polarisabilities[1, 1, i, i]) / (natoms)

            atom_polarisabilities[i, 4] = np.sum(
                local_polarisabilities[1, 2, i, :] +
                local_polarisabilities[1, 2, :, i] +
                (natoms - 2) * local_polarisabilities[1, 2, i, i]) / (natoms)

            atom_polarisabilities[i, 5] = np.sum(
                local_polarisabilities[2, 2, i, :] +
                local_polarisabilities[2, 2, :, i] +
                (natoms - 2) * local_polarisabilities[2, 2, i, i]) / (natoms)

        self.print_results(molecule, natoms, Qab, local_polarisabilities,
                           molecule_polarisabilities, atom_polarisabilities)

        if self.rank == mpi_master():
            ret_dict = {
                'localised_charges': Qab,
                'localised_polarisabilities': local_polarisabilities
            }
            return ret_dict

    def penalty_fc(self, za, Ra, zb, Rb):
        """
        penalty function        

        :param za:
            atomic number for atom a
        :param Ra:
            position vector for atom a
        :param zb:
            atomic number for atom b
        :param Rb:
            position vector for atom b
             
        :return:
            The penalty function
        """

        #Ra/Rb are the position vectors
        AUG2AU_factor = bohr_in_angstroms()
        AUG2AU = 1.0 / AUG2AU_factor
        RBS = np.array([
            0, 0.25, 0.25, 1.45, 1.05, 0.85, 0.7, 0.65, 0.5, 0.43, 1.8, 1.5,
            1.25, 1.1, 1.0, 1.0, 1.0, 1.0
        ]) * AUG2AU
        assert_msg_critical('za>17 or zb>17',
                            'Response: we currently support up to Cl')
        ra = RBS[za]
        rb = RBS[zb]
        #Ra and Rb taking from r_mat: 0,1,2 represents x y z directions
        rab2 = (Ra[0] - Rb[0])**2 + (Ra[1] - Rb[1])**2 + (Ra[2] - Rb[2])**2
        #alpha set to 2
        f = 0.5 * np.exp(-2 * (rab2 / (ra + rb)**2))
        return f

    #unpact response vector to matrix form
    def lrvec2mat(self, vec, nocc, norb):
        """
         unpact response vector to matrix form
         param: vec
             the response vector
         param: nocc
             number of contracted orbitals (alpha electrons)
         param: norb
             number of orbitals 
             
         return: kappa
              unpacked response vector
         """
        kappa = np.zeros((norb, norb))
        i = 0
        for r in range(nocc):
            for c in range(nocc, norb):
                kappa[r, c] = vec[i]
                i += 1
        kappa = kappa + kappa.T
        return kappa

    def lowdin_orthonormalize(self, A):
        """
        Lodwin orthonormaliztaion.

        :param A:
            The input matrix.

        :return:
            The orthonormalised vector.
        """

        eigs, U = np.linalg.eigh(A)
        X = np.linalg.multi_dot([U, np.diag(1.0 / np.sqrt(eigs)), U.T])
        return X

    def get_ao_indices(self, molecule, basis):
        """
        Gets information about occupied/virtual orbitals by counting the basis
        set functions.
        
        :param molecule:
            The molecule.
        :param basis:
            The molecular basis set.
        
        :return:
            
        """

        elements = molecule.elem_ids_to_numpy()
        basis = basis.get_label()

        #get basis file
        local_name = basis
        basis_file = []
        lib_member = os.path.join(os.environ['VLXBASISPATH'], basis)
        if os.path.exists(local_name):
            basis_file = os.path.abspath(local_name)
        elif os.path.exists(lib_member):
            basis_file = lib_member

        bp = InputParser(basis_file)
        basis = bp.get_dict()
        keys = list(basis.keys())

        ao_occ = []
        ao_vir = []
        iterr = 0
        multiplicity = dict(S=1, P=3, D=5, F=7, G=9, H=11, I=13)
        ao_per_atom = []
        for e in elements:
            k = keys[e - 1]
            atoms_data = basis[k]
            c = Counter(
                line[0] for line in atoms_data if line and line[0] in "SPDFGHI")

            assert_msg_critical('e>10',
                                'Response: we currently support up to Ne')
            # For H and He: 1s
            ao_occ.append(iterr)

            # For Li-Ne: + 2s 2p
            if e >= 3:
                ao_occ.append(iterr + 1)
                offset_p = c['S'] + iterr
                ao_occ.append(offset_p + 0)
                ao_occ.append(offset_p + 1)
                ao_occ.append(offset_p + 2)

            #sum no. orbitals. used to calculate virtual orbitals later
            orb = sum(multiplicity[k] * v for k, v in c.items())
            iterr = iterr + orb
            ao_per_atom.append(orb)

        total_orb_array = np.arange(iterr)

        for i in total_orb_array:
            if i not in ao_occ:
                ao_vir.append(i)

        return ao_per_atom, ao_occ, ao_vir

    def print_results(self, molecule, natoms, Qab, local_polarisabilities, Am,
                      atom_polarisabilities):
        """
        Prints results for loprop calculation.
        
        param: molecule
            the molecule
        param: natoms
            number of atoms
        param: Qab
            localised charges
        param: local_polarisabilities
            localised polarisabilities
        param: Am
            molecular polarisabilities
        """
        element_names = molecule.get_labels()
        width = 92

        self.ostream.print_blank()
        self.ostream.print_header('Local Properties (LoProp) Calculations')
        self.ostream.print_header(19 * '=')
        self.ostream.print_blank()

        output_header = '*** Molecular Polarisabilities *** '
        self.ostream.print_header(output_header.ljust(width))

        direction = ["x", "y", "z"]
        for a in range(3):
            output_iter = '<<{};{}>>: {:15.8f} '.format(direction[a],
                                                        direction[a], Am[a, a])
            self.ostream.print_header(output_iter.ljust(width))
        self.ostream.print_blank()

        #print localised chagres
        output_header = '*** Loprop localised charges *** '
        self.ostream.print_header(output_header.ljust(width))
        for a in range(natoms):
            output_iter = '{:<5s}: {:15.8f} '.format(element_names[a], Qab[a,
                                                                           a])
            self.ostream.print_header(output_iter.ljust(width))
        self.ostream.print_blank()

        output_header = '*** Localised polarisabilities *** '
        self.ostream.print_header(output_header.ljust(width))
        for a in range(natoms):
            output_iter = '{:<5s}: {:10.4f}  {:10.4f}  {:10.4f}  {:10.4f} {:10.4f} {:10.4f}'.format(
                element_names[a], atom_polarisabilities[a][0],
                atom_polarisabilities[a][1], atom_polarisabilities[a][2],
                atom_polarisabilities[a][3], atom_polarisabilities[a][4],
                atom_polarisabilities[a][5])
            self.ostream.print_header(output_iter.ljust(width))
        self.ostream.print_blank()

        self.ostream.flush()
