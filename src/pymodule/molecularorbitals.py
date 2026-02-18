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

from enum import Enum
import numpy as np
import h5py
import math
import sys

from .veloxchemlib import tensor_order
from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical


class molorb(Enum):
    rest = 1
    unrest = 2
    restopen = 3


class MolecularOrbitals:
    """
    Implements the molecular orbitals class.

    Instance variables
        - _orbitals: The molecular orbitals.
        - _energies: The molecular orbital energies.
        - _occupations: The molecular orbital occupation numbers.
        - _orbital_type: The type of molecular orbitals.
    """

    def __init__(self, orbs=None, enes=None, occs=None, orbs_type=None):
        """
        Initializes the molecular orbitals object.
        """

        if (orbs is not None and enes is not None and occs is not None and
                orbs_type is not None):
            # Note: all numpy arrays are deep-copied
            self._orbitals = []
            for orb in orbs:
                self._orbitals.append(orb.copy())
            self._energies = []
            for ene in enes:
                self._energies.append(ene.copy())
            self._occupations = []
            for occ in occs:
                self._occupations.append(occ.copy())
            self._orbitals_type = orbs_type

            assert_msg_critical(
                self._orbitals_type
                in [molorb.rest, molorb.unrest, molorb.restopen],
                'MolecularOrbitals: Invalid orbitals type')

            err_num = 'MolecularOrbitals: Inconsistent orbitals, energies or'
            err_num += ' occupation numbers'
            if self._orbitals_type == molorb.rest:
                assert_msg_critical(
                    len(self._orbitals) == 1 and len(self._energies) == 1 and
                    len(self._occupations) == 1, err_num)
            elif self._orbitals_type == molorb.unrest:
                assert_msg_critical(
                    len(self._orbitals) == 2 and len(self._energies) == 2 and
                    len(self._occupations) == 2, err_num)
            elif self._orbitals_type == molorb.restopen:
                assert_msg_critical(
                    len(self._orbitals) == 1 and len(self._energies) == 1 and
                    len(self._occupations) == 2, err_num)

        else:
            self._orbitals = None
            self._energies = None
            self._occupations = None
            self._orbitals_type = None

    def get_orbitals_type(self):

        return self._orbitals_type

    def alpha_to_numpy(self):

        return self._orbitals[0].copy()

    def beta_to_numpy(self):

        if self._orbitals_type in [molorb.rest, molorb.restopen]:
            return self._orbitals[0].copy()
        else:
            return self._orbitals[1].copy()

    def ea_to_numpy(self):

        return self._energies[0].copy()

    def eb_to_numpy(self):

        if self._orbitals_type in [molorb.rest, molorb.restopen]:
            return self._energies[0].copy()
        else:
            return self._energies[1].copy()

    def occa_to_numpy(self):

        return self._occupations[0].copy()

    def occb_to_numpy(self):

        if self._orbitals_type == molorb.rest:
            return self._occupations[0].copy()
        else:
            return self._occupations[1].copy()

    def number_aos(self):

        # for backward compatibility
        return self._orbitals[0].shape[0]

    def number_of_aos(self):

        return self._orbitals[0].shape[0]

    def number_mos(self):

        # for backward compatibility
        return self._orbitals[0].shape[1]

    def number_of_mos(self):

        return self._orbitals[0].shape[1]

    def print_orbitals(self, molecule, basis, orb_inds=None, ostream=None):
        """
        Prints molecular orbitals to output stream.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param orb_inds:
            The starting and ending indices of orbitals.
        :param ostream:
            The output stream.
        """

        if ostream is None:
            ostream = OutputStream(sys.stdout)

        norb = self.number_of_mos()

        ao_map = basis.get_ao_basis_map(molecule)

        if self._orbitals_type in [molorb.rest, molorb.restopen]:

            ostream.print_blank()
            ostream.print_header("Spin Restricted Orbitals")
            ostream.print_header("------------------------")

            nocc = molecule.number_of_electrons() // 2

            if isinstance(orb_inds, (np.ndarray, tuple, list)):
                assert_msg_critical(
                    len(orb_inds) == 2, 'MolecularOrbitals.print_orbitals: ' +
                    'Expecting starting and ending indices of orbitals')
                nstart, nend = max(0, orb_inds[0]), min(norb, orb_inds[1])
            else:
                if not orb_inds:
                    nstart, nend = max(0, nocc - 5), min(norb, nocc + 5)
                else:
                    nstart, nend = 0, norb

            rvecs = self.alpha_to_numpy()
            reigs = self.ea_to_numpy()
            rnocc = self.occa_to_numpy() + self.occb_to_numpy()

            for i in range(nstart, nend):
                self._print_coefficients(reigs[i], rnocc[i], i, rvecs[:, i],
                                         ao_map, 0.15, ostream)

            ostream.print_blank()
            ostream.flush()

        elif self._orbitals_type == molorb.unrest:

            ostream.print_blank()
            ostream.print_header("Spin Unrestricted Alpha Orbitals")
            ostream.print_header("--------------------------------")

            nalpha = molecule.number_of_alpha_electrons()

            if isinstance(orb_inds, (np.ndarray, tuple, list)):
                assert_msg_critical(
                    len(orb_inds) == 2, 'MolecularOrbitals.print_orbitals: ' +
                    'Expecting starting and ending indices of orbitals')
                nstart, nend = max(0, orb_inds[0]), min(norb, orb_inds[1])
            else:
                if not orb_inds:
                    nstart, nend = max(0, nalpha - 5), min(norb, nalpha + 5)
                else:
                    nstart, nend = 0, norb

            uvecs = self.alpha_to_numpy()
            ueigs = self.ea_to_numpy()
            unocc = self.occa_to_numpy()

            for i in range(nstart, nend):
                self._print_coefficients(ueigs[i], unocc[i], i, uvecs[:, i],
                                         ao_map, 0.15, ostream)

            ostream.print_blank()
            ostream.print_header("Spin Unrestricted Beta Orbitals")
            ostream.print_header("-------------------------------")

            nbeta = molecule.number_of_beta_electrons()

            if isinstance(orb_inds, (np.ndarray, tuple, list)):
                assert_msg_critical(
                    len(orb_inds) == 2, 'MolecularOrbitals.print_orbitals: ' +
                    'Expecting starting and ending indices of orbitals')
                nstart, nend = max(0, orb_inds[0]), min(norb, orb_inds[1])
            else:
                if not orb_inds:
                    nstart, nend = max(0, nbeta - 5), min(norb, nbeta + 5)
                else:
                    nstart, nend = 0, norb

            uvecs = self.beta_to_numpy()
            ueigs = self.eb_to_numpy()
            unocc = self.occb_to_numpy()

            for i in range(nstart, nend):
                self._print_coefficients(ueigs[i], unocc[i], i, uvecs[:, i],
                                         ao_map, 0.15, ostream)

            ostream.print_blank()
            ostream.flush()

        else:

            errmsg = "MolecularOrbitals.print_orbitals:"
            errmsg += " Invalid molecular orbitals type"
            assert_msg_critical(False, errmsg)

    @staticmethod
    def _print_coefficients(eigval, focc, iorb, coeffs, ao_map, thresh,
                            ostream):
        """
        Prints molecular orbital coefficients to output stream.

        :param eigval:
            The eigenvalue (orbital energy).
        :param focc:
            The occupation number of the orbital.
        :param iorb:
            The index (0-based) of the orbital.
        :param coeffs:
            The AO coefficients of the orbital.
        :param ao_map:
            The string representation map of basis functions.
        :param thresh:
            The threshold for priting AO coefficients.
        :param ostream:
            The output stream.
        """

        ostream.print_blank()

        valstr = "Molecular Orbital No.{:4d}:".format(iorb + 1)
        ostream.print_header(valstr.ljust(92))
        valstr = 26 * "-"
        ostream.print_header(valstr.ljust(92))

        valstr = "Occupation: {:.3f} Energy: {:10.5f} a.u.".format(focc, eigval)
        ostream.print_header(valstr.ljust(92))

        tuplist = []

        for i in range(coeffs.shape[0]):

            if math.fabs(coeffs[i]) > thresh:
                atomidx = int(ao_map[i].split()[0])
                anglmom = tensor_order(ao_map[i][-3].upper())
                valstr = "(" + ao_map[i] + ": {:8.2f}".format(coeffs[i]) + ") "
                tuplist.append((atomidx, anglmom, valstr))

        valstr = ""
        curidx = 0

        for t in sorted(tuplist):
            valstr += t[-1]
            curidx += 1

            if curidx == 3:
                ostream.print_header(valstr.ljust(92))
                valstr = ""
                curidx = 0

        if curidx > 0:
            ostream.print_header(valstr.ljust(92))

    def get_density(self, molecule, scf_type=None):
        """
        Gets AO density matrix from molecular orbitals.

        :param molecule:
            The molecule.
        :param scf_type:
            The type of SCF calculation (restricted, unrestricted, or
            restricted_openshell).

        :return:
            The AO density matrix.
        """

        if self._orbitals_type == molorb.rest and (scf_type is None or
                                                   scf_type == 'restricted'):

            mo = self.alpha_to_numpy()
            occ = self.occa_to_numpy()

            occ_mo = occ * mo
            dens = np.matmul(occ_mo, occ_mo.T)

            # Note: return a tuple
            return (dens,)

        elif self._orbitals_type == molorb.unrest and (
                scf_type is None or scf_type == 'unrestricted'):

            mo_a = self.alpha_to_numpy()
            mo_b = self.beta_to_numpy()

            occ_a = self.occa_to_numpy()
            occ_b = self.occb_to_numpy()

            occ_mo_a = occ_a * mo_a
            occ_mo_b = occ_b * mo_b

            dens_a = np.matmul(occ_mo_a, occ_mo_a.T)
            dens_b = np.matmul(occ_mo_b, occ_mo_b.T)

            return (dens_a, dens_b)

        elif self._orbitals_type == molorb.restopen and (
                scf_type is None or scf_type == 'restricted_openshell'):

            mo = self.alpha_to_numpy()

            occ_a = self.occa_to_numpy()
            occ_b = self.occb_to_numpy()

            occ_mo_a = occ_a * mo
            occ_mo_b = occ_b * mo

            dens_a = np.matmul(occ_mo_a, occ_mo_a.T)
            dens_b = np.matmul(occ_mo_b, occ_mo_b.T)

            return (dens_a, dens_b)

        else:

            errmsg = "MolecularOrbitals.get_density: "
            errmsg += " Invalid molecular orbitals type"
            assert_msg_critical(False, errmsg)

    def is_empty(self):
        """
        Checks if the molecular orbitals object is empty.
        """

        return (self._orbitals is None)

    def broadcast(self, comm, root=mpi_master()):
        """
        Broadcasts the molecular orbitals object.

        :param comm:
            The MPI communicator.
        :param root:
            The root rank to broadcast from.

        :return:
            The molecular orbitals object.
        """

        mo_occs = []
        mo_coefs = []
        mo_enes = []

        mo_type = comm.bcast(self._orbitals_type, root=root)

        if comm.Get_rank() == root:
            mo_occs.append(self._occupations[0])
            mo_coefs.append(self._orbitals[0])
            mo_enes.append(self._energies[0])

            if mo_type != molorb.rest:
                mo_occs.append(self._occupations[1])
                if mo_type == molorb.unrest:
                    mo_coefs.append(self._orbitals[1])
                    mo_enes.append(self._energies[1])

        else:
            mo_occs.append(None)
            mo_coefs.append(None)
            mo_enes.append(None)

            if mo_type != molorb.rest:
                mo_occs.append(None)
                if mo_type == molorb.unrest:
                    mo_coefs.append(None)
                    mo_enes.append(None)

        mo_coefs[0] = comm.bcast(mo_coefs[0], root=root)
        mo_enes[0] = comm.bcast(mo_enes[0], root=root)
        mo_occs[0] = comm.bcast(mo_occs[0], root=root)

        if mo_type != molorb.rest:
            mo_occs[1] = comm.bcast(mo_occs[1], root=root)
            if mo_type == molorb.unrest:
                mo_coefs[1] = comm.bcast(mo_coefs[1], root=root)
                mo_enes[1] = comm.bcast(mo_enes[1], root=root)

        return MolecularOrbitals(mo_coefs, mo_enes, mo_occs, mo_type)

    def write_hdf5(self, fname, nuclear_charges=None, basis_set=None, label=''):
        """
        Writes molecular orbitals to hdf5 file.

        :param fname:
            The name of the hdf5 file.
        :param nuclear_charges:
            The nuclear charges.
        :param basis_set:
            Name of the basis set.
        """

        if label and isinstance(label, str):
            prefix = label + '_'
        else:
            prefix = ''

        hf = h5py.File(fname, 'a')

        for key in [
                prefix + 'alpha_orbitals',
                prefix + 'alpha_energies',
                prefix + 'alpha_occupations',
                prefix + 'beta_orbitals',
                prefix + 'beta_energies',
                prefix + 'beta_occupations',
                prefix + 'nuclear_charges',
                prefix + 'basis_set',
        ]:
            if key in hf:
                del hf[key]

        hf.create_dataset(prefix + 'alpha_orbitals', data=self.alpha_to_numpy())
        hf.create_dataset(prefix + 'alpha_energies', data=self.ea_to_numpy())
        hf.create_dataset(prefix + 'alpha_occupations',
                          data=self.occa_to_numpy())

        if self._orbitals_type == molorb.unrest:
            hf.create_dataset(prefix + 'beta_orbitals',
                              data=self.beta_to_numpy())
            hf.create_dataset(prefix + 'beta_energies', data=self.eb_to_numpy())
            hf.create_dataset(prefix + 'beta_occupations',
                              data=self.occb_to_numpy())

        elif self._orbitals_type == molorb.restopen:
            hf.create_dataset(prefix + 'beta_occupations',
                              data=self.occb_to_numpy())

        if nuclear_charges is not None:
            hf.create_dataset(prefix + 'nuclear_charges', data=nuclear_charges)

        if basis_set is not None:
            hf.create_dataset(prefix + 'basis_set', data=np.bytes_([basis_set]))

        hf.close()

    @staticmethod
    def check_label_validity(fname, label=''):
        """
        Checks label validity for molecular orbitals hdf5 file.

        :param fname:
            The name of the hdf5 file.

        :return:
            True if the label is valid, False otherwise.
        """

        if label and isinstance(label, str):
            prefix = label + '_'
        else:
            prefix = ''

        hf = h5py.File(fname, 'r')

        is_valid = True

        for key in ['alpha_orbitals', 'alpha_energies', 'alpha_occupations']:
            key_in_hf = ((prefix + key) in hf)
            is_valid = (is_valid and key_in_hf)

        hf.close()

        return is_valid

    @staticmethod
    def read_hdf5(fname, label=''):
        """
        Reads molecular orbitals from hdf5 file.

        :param fname:
            The name of the hdf5 file.

        :return:
            The molecular orbitals.
        """

        if label and isinstance(label, str):
            prefix = label + '_'
        else:
            prefix = ''

        hf = h5py.File(fname, 'r')

        orbs_type = molorb.rest

        for key in ['alpha_orbitals', 'alpha_energies', 'alpha_occupations']:
            assert_msg_critical((prefix + key) in hf,
                                f'MolecularOrbitals.read_hdf5: {key} not found')

        if 'beta_orbitals' in hf or 'beta_energies' in hf:
            orbs_type = molorb.unrest

            for key in ['beta_orbitals', 'beta_energies', 'beta_occupations']:
                assert_msg_critical(
                    (prefix + key) in hf,
                    f'MolecularOrbitals.read_hdf5: {key} not found')

        elif 'beta_occupations' in hf:
            orbs_type = molorb.restopen

        orbs = []
        enes = []
        occs = []

        orbs.append(np.array(hf.get(prefix + 'alpha_orbitals')))
        enes.append(np.array(hf.get(prefix + 'alpha_energies')))
        occs.append(np.array(hf.get(prefix + 'alpha_occupations')))

        if orbs_type == molorb.unrest:
            orbs.append(np.array(hf.get(prefix + 'beta_orbitals')))
            enes.append(np.array(hf.get(prefix + 'beta_energies')))
            occs.append(np.array(hf.get(prefix + 'beta_occupations')))

        elif orbs_type == molorb.restopen:
            occs.append(np.array(hf.get(prefix + 'beta_occupations')))

        hf.close()

        return MolecularOrbitals(orbs, enes, occs, orbs_type)

    @staticmethod
    def match_hdf5(fname, nuclear_charges, basis_set, scf_type):
        """
        Checks if the hdf5 file matches the given nuclear charges and basis set.

        :param fname:
            The name of the hdf5 file.
        :param nuclear_charges:
            The nuclear charges.
        :param basis_set:
            Name of the basis set.
        :param scf_type:
            The type of SCF calculation (restricted, unrestricted, or
            restricted_openshell).

        :return:
            Whether the hdf5 file matches the given nuclear charges and basis set.
        """

        hf = h5py.File(fname, 'r')

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

        if 'beta_orbitals' in hf or 'beta_energies' in hf:
            h5_scf_type = 'unrestricted'
        elif 'beta_occupations' in hf:
            h5_scf_type = 'restricted_openshell'
        else:
            h5_scf_type = 'restricted'
        match_scf_type = (h5_scf_type == scf_type)

        hf.close()

        return (match_nuclear_charges and match_basis_set and match_scf_type)

    @staticmethod
    def create_nto(nto_orbitals_list, nto_lambdas_list, mo_type):
        """
        Creates molecular orbitals object for natural transiton orbitals (NTO).

        :param nto_orbs_list:
            The list of NTO orbital coefficients.
        :param nto_lambdas_list:
            The list of NTO lambda's.
        :param mo_type:
            The molecular orbital type.

        :return:
            The molecular orbitals object for NTO.
        """

        # sanity checks

        for nto_orbitals, nto_lambdas in zip(nto_orbitals_list,
                                             nto_lambdas_list):

            assert_msg_critical(
                nto_orbitals.shape[1] == nto_lambdas.shape[0],
                'MolecularOrbitals.create_nto: Inconsistent size')

            negative_lambdas = [x for x in nto_lambdas if x < 0.0]
            positive_lambdas = [x for x in nto_lambdas if x > 0.0]

            assert_msg_critical(
                len(negative_lambdas) == len(positive_lambdas),
                'MolecularOrbitals.create_nto: Inconsistent number of lambda values'
            )

            for m, p in zip(negative_lambdas[::-1], positive_lambdas):
                assert_msg_critical(
                    abs(m + p) < 1.0e-6,
                    'MolecularOrbitals.create_nto: Inconsistent lambda values')

        nto_energies_list = [
            np.zeros(nto_lambdas_list[idx].shape[0])
            for idx in range(len(nto_lambdas_list))
        ]

        return MolecularOrbitals(nto_orbitals_list, nto_energies_list,
                                 nto_lambdas_list, mo_type)

    def is_nto(self):
        """
        Checks if this molecular orbitals object is natural transiton orbitals
        (NTO).

        :return:
            True if this molecular orbitals object is NTO, False otherwise.
        """

        assert_msg_critical(
            self._orbitals_type == molorb.rest,
            'MolecularOrbitals.is_nto: Only restricted case is implemented')

        nto_lambdas = self.occa_to_numpy()

        negative_lambdas = [x for x in nto_lambdas if x < 0.0]
        positive_lambdas = [x for x in nto_lambdas if x > 0.0]

        if len(negative_lambdas) != len(positive_lambdas):
            return False

        for m, p in zip(negative_lambdas[::-1], positive_lambdas):
            if abs(m + p) >= 1.0e-6:
                return False

        nto_energies = self.ea_to_numpy()

        for e in nto_energies:
            if e != 0.0:
                return False

        return True

    def __deepcopy__(self, memo):
        """
        Implements deepcopy.

        :param memo:
            The memo dictionary for deepcopy.

        :return:
            A deepcopy of self.
        """

        return MolecularOrbitals(self)
