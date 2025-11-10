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


from mpi4py import MPI
import numpy as np
import sys

from .veloxchem import (mpi_master,
                        compute_kinetic_energy_integrals,
                        compute_nuclear_potential_integrals,
                        compute_overlap_integrals,
                        FockDriver
                        )

from .outputstream import OutputStream


class SOBasisIntegrals:

    def __init__(self, comm=None, ostream=None):

        """
        Initializes the SOIntegrals class.

        :param comm:
            The MPI communicator. 
        :param ostream:
            The output stream.

        Instance variables
            - molecule: The molecular structure.
            - basis_set: The basis set used for calculations.
            - S: Overlap integrals matrix.
            - T: Kinetic energy integrals matrix.
            - V: Nuclear attraction integrals matrix.
            - h_MO: One-electron Hamiltonian in MO basis.
            - g_MO: Two-electron integrals in MO basis.
            - g_MO_physicist: Two-electron integrals in physicist's notation.
            - h_SO: One-electron Hamiltonian in SO basis.
            - g_SO: Two-electron integrals in SO basis.
            - threshold: Threshold for zeroing small integrals.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)
        
        self.comm = comm
        self.ostream = ostream

        # Molecule and basis set
        self.molecule = None
        self.basis_set = None
        
        # AO integrals
        self.S = None  # Overlap integrals
        self.T = None  # Kinetic energy integrals
        self.V = None  # Nuclear attraction integrals
        self.h_AO = None  # One-electron integrals in AO basis
        self.g_AO = None  # Two-electron integrals in AO basis

        # MO integrals
        self.h_MO = None  # One-electron integrals in MO basis
        self.g_MO = None  # Two-electron integrals in MO basis
        self.g_MO_physicist = None  # Two-electron integrals in physicist's notation

        # SO integrals
        self.h_SO = None  # One-electron integrals in SO basis
        self.g_SO = None  # Two-electron integrals in SO basis

        self.threshold = 1e-9  # Threshold for zeroing small integrals


    def make_zeroes(self, matrix):
        """
        Set to zero elements smaller than the threshold.
        
        :param matrix:
            The input matrix.
        :return:
            The matrix with small elements set to zero.
        """
        matrix[np.abs(matrix) < self.threshold] = 0.0

        return matrix

    def compute_integrals(self, molecule, basis_set): 
        """
        Compute one- and two-electron integrals in SO basis.

        :param molecule:
            The molecular structure.
        :param basis_set:
            The basis set used for calculations.
        """

        # Store molecule and basis set for export
        self.molecule = molecule
        self.basis_set = basis_set

        # One electron integrals in AO basis

        self.T = compute_kinetic_energy_integrals(molecule, basis_set)
        self.V = compute_nuclear_potential_integrals(molecule, basis_set)
        self.S = compute_overlap_integrals(molecule, basis_set)

        self.h_AO = self.T + self.V
        self.h_AO = self.make_zeroes(self.h_AO)

        # Transform to MO basis
        eigvals, eigvecs = np.linalg.eigh(self.S)
        S_half_inv = eigvecs @ np.diag(1.0 / np.sqrt(eigvals)) @ eigvecs.T
        T_ortho = S_half_inv.T @ self.T @ S_half_inv
        V_ortho = S_half_inv.T @ self.V @ S_half_inv
        h_ortho = T_ortho + V_ortho
        eigvals_h, C_MO = np.linalg.eigh(h_ortho)
        self.h_MO = C_MO.T @ h_ortho @ C_MO

        # Two electron integrals in AO basis
        fock_driver = FockDriver()
        self.g_AO = fock_driver.compute_eri(molecule, basis_set)

        # Transform to MO basis
        self.g_MO = np.einsum('pqrs,pi,qj,rk,sl->ijkl',
                         self.g_AO, C_MO, C_MO, C_MO, C_MO,
                         optimize=True)
        self.g_MO = self.make_zeroes(self.g_MO)
        self.g_MO_physicist = self.g_MO.transpose(0, 2, 1, 3)

        # Transform to SO basis
        # Get number of MOs and SOs
        n_MO = self.h_MO.shape[0]
        n_SO = 2 * n_MO

        # One-electron integrals in SO basis
        self.h_SO = np.zeros((n_SO, n_SO))
        for p in range(n_MO):
            for q in range(n_MO):
                self.h_SO[2*p, 2*q] = self.h_MO[p, q]
                self.h_SO[2*p+1, 2*q+1] = self.h_MO[p, q]

        # Two-electron integrals in SO basis
        self.g_SO = np.zeros((n_SO, n_SO, n_SO, n_SO))
        for p in range(n_MO):
            for q in range(n_MO):
                for r in range(n_MO):
                    for s in range(n_MO):
                        self.g_SO[2*p, 2*q + 1, 2*r + 1, 2*s] = self.g_MO_physicist[p, q, r, s]
                        self.g_SO[2*p+1, 2*q, 2*r, 2*s+1] = self.g_MO_physicist[p, q, r, s]
                        self.g_SO[2*p, 2*q, 2*r, 2*s] = self.g_MO_physicist[p, q, r, s]
                        self.g_SO[2*p+1, 2*q + 1, 2*r+1, 2*s+1] = self.g_MO_physicist[p, q, r, s]

        
        self.ostream.print_info("One- and two-electron integrals in SO basis computed.")
        self.ostream.print_info("Verifying symmetries of two-electron integrals in SO basis")
        self.ostream.flush()

        if self.verify_two_electron_integral_symmetries(self.g_SO):
            self.ostream.print_info("The two-electron integrals in SO basis satisfy all symmetries.")
            self.ostream.flush()
        else:
            self.ostream.print_warning("The two-electron integrals in SO basis do NOT satisfy all symmetries!")
            self.ostream.flush()
   
    def verify_two_electron_integral_symmetries(self, g_SO, tol=1e-9):
        """
        Verify the symmetries of the two-electron integrals in SO basis.

        :param g_SO:
            The two-electron integrals in SO basis.
        :param tol:
            The tolerance for symmetry checks.
        :return:
            True if all symmetries are satisfied within the tolerance, False otherwise.
        """

        n_SO = g_SO.shape[0]

        for p in range(n_SO):
            for q in range(n_SO):
                for r in range(n_SO):
                    for s in range(n_SO):
                        test1 = np.abs(g_SO[p, q, r, s] - g_SO[q, p, r, s]) < tol
                        test2 = np.abs(g_SO[p, q, r, s] - g_SO[p, q, s, r]) < tol
                        test3 = np.abs(g_SO[p, q, r, s] - g_SO[r, s, p, q]) < tol
                        if test1 is False or test2 is False or test3 is False:
                            return False
        return True
    
    # Helper function to integrate with Qrisp
    
    def to_qrisp_data(self):
        """
        Convert integrals to Qrisp data format to be processed by the function 
        create_electronic_hamiltonian. 
        https://github.com/eclipse-qrisp/Qrisp/blob/main/src/qrisp/algorithms/vqe/problems/electronic_structure.py

        :return:
            A dictionary containing one- and two-electron integrals in SO basis.
        """

        qrisp_data = {
            "mol": self.molecule,
            "one_int": self.h_SO,
            "two_int": self.g_SO,
            "num_orb": self.h_SO.shape[0],
            "num_elec": self.molecule.number_of_alpha_electrons(),
            "energy_nuc": self.molecule.nuclear_repulsion_energy()
        }

        return qrisp_data

