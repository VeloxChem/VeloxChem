#
#                              VELOXCHEM
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2023 by VeloxChem developers. All rights reserved.
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

from enum import Enum
import numpy as np
import h5py
import math
import sys

from .veloxchemlib import to_angular_momentum
from .veloxchemlib import denmat
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import matmul_gpu
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical


class molorb(Enum):
    rest = 1
    unrest = 2
    restopen = 3


class MolecularOrbitals:
    """
    TODO
    """

    def __init__(self, orbs=None, enes=None, occs=None, orbs_type=None):
        """
        TODO
        """

        if (orbs is not None and enes is not None and occs is not None and
                orbs_type is not None):
            self._orbitals = list(orbs)
            self._energies = list(enes)
            self._occupations = list(occs)
            self._orbitals_type = orbs_type
        else:
            self._orbitals = None
            self._energies = None
            self._occupations = None
            self._orbitals_type = None

    def get_orbitals_type(self):

        return self._orbitals_type

    def alpha_to_numpy(self):

        return self._orbitals[0]

    def ea_to_numpy(self):

        return self._energies[0]

    def occa_to_numpy(self):

        return self._occupations[0]

    def number_of_aos(self):

        return self._orbitals[0].shape[0]

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

        if self.get_orbitals_type() in [molorb.rest, molorb.restopen]:

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

        elif self.get_orbitals_type() == molorb.unrest:

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

        else:

            errmsg = "MolecularOrbitals.print_orbitals:"
            errmsg += " Invalid molecular orbitals type"
            assert_msg_critical(False, errmsg)

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
                anglmom = to_angular_momentum(ao_map[i][-3].upper())
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

        if self.get_orbitals_type() == molorb.rest and (
                scf_type is None or scf_type == 'restricted'):

            mo = self.alpha_to_numpy()

            occ = self.occa_to_numpy()

            # TODO: use hipblas
            occ_mo = occ * mo

            dens = matmul_gpu(occ_mo, occ_mo.T)

            return AODensityMatrix([dens], denmat.rest)

        elif self.get_orbitals_type() == molorb.unrest and (
                scf_type is None or scf_type == 'unrestricted'):

            mo_a = self.alpha_to_numpy()
            mo_b = self.beta_to_numpy()

            occ_a = self.occa_to_numpy()
            occ_b = self.occb_to_numpy()

            occ_mo_a = occ_a * mo_a
            occ_mo_b = occ_b * mo_b

            dens_a = np.matmul(occ_mo_a, occ_mo_a.T)
            dens_b = np.matmul(occ_mo_b, occ_mo_b.T)

            return AODensityMatrix([dens_a, dens_b], denmat.unrest)

        elif self.get_orbitals_type() == molorb.restopen and (
                scf_type is None or scf_type == 'restricted_openshell'):

            mo = self.alpha_to_numpy()

            occ_a = self.occa_to_numpy()
            occ_b = self.occb_to_numpy()

            occ_mo_a = occ_a * mo
            occ_mo_b = occ_b * mo

            dens_a = np.matmul(occ_mo_a, occ_mo_a.T)
            dens_b = np.matmul(occ_mo_b, occ_mo_b.T)

            return AODensityMatrix([dens_a, dens_b], denmat.unrest)

        else:

            errmsg = "MolecularOrbitals.get_density: "
            errmsg += " Invalid molecular orbitals type"
            assert_msg_critical(False, errmsg)

    def write_hdf5(self, fname, nuclear_charges=None, basis_set=None):
        """
        Writes molecular orbitals to hdf5 file.

        :param fname:
            The name of the hdf5 file.
        :param nuclear_charges:
            The nuclear charges.
        :param basis_set:
            Name of the basis set.
        """

        hf = h5py.File(fname, 'w')

        hf.create_dataset('alpha_orbitals', data=self.alpha_to_numpy())
        hf.create_dataset('alpha_energies', data=self.ea_to_numpy())
        hf.create_dataset('alpha_occupations', data=self.occa_to_numpy())

        if self.get_orbitals_type() == molorb.unrest:
            hf.create_dataset('beta_orbitals', data=self.beta_to_numpy())
            hf.create_dataset('beta_energies', data=self.eb_to_numpy())
            hf.create_dataset('beta_occupations', data=self.occb_to_numpy())

        elif self.get_orbitals_type() == molorb.restopen:
            hf.create_dataset('beta_occupations', data=self.occb_to_numpy())

        if nuclear_charges is not None:
            hf.create_dataset('nuclear_charges', data=nuclear_charges)

        if basis_set is not None:
            hf.create_dataset('basis_set', data=np.string_([basis_set]))

        hf.close()

    @staticmethod
    def read_hdf5(fname):
        """
        Reads molecular orbitals from hdf5 file.

        :param fname:
            The name of the hdf5 file.

        :return:
            The molecular orbitals.
        """

        hf = h5py.File(fname, 'r')

        orbs_type = molorb.rest

        assert_msg_critical(
            'alpha_orbitals' in hf and 'alpha_energies' in hf and
            'alpha_occupations' in hf, 'MolecularOrbitals.read_hdf5: ' +
            'alpha orbitals/energies/occupations not found')

        if 'beta_orbitals' in hf or 'beta_energies' in hf:
            orbs_type = molorb.unrest

            assert_msg_critical(
                'beta_orbitals' in hf and 'beta_energies' in hf and
                'beta_occupations' in hf, 'MolecularOrbitals.read_hdf5: ' +
                'beta orbitals/energies/occupations not found')

        elif 'beta_occupations' in hf:
            orbs_type = molorb.restopen

        orbs = []
        enes = []
        occs = []

        orbs.append(np.array(hf.get('alpha_orbitals')))
        enes.append(np.array(hf.get('alpha_energies')))
        occs.append(np.array(hf.get('alpha_occupations')))

        if orbs_type == molorb.unrest:
            orbs.append(np.array(hf.get('beta_orbitals')))
            enes.append(np.array(hf.get('beta_energies')))
            occs.append(np.array(hf.get('beta_occupations')))

        elif orbs_type == molorb.restopen:
            occs.append(np.array(hf.get('beta_occupations')))

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

        assert_msg_critical(
            len(nto_orbitals_list) == 1 and len(nto_lambdas_list) == 1 and
            mo_type == molorb.rest,
            'MolecularOrbitals.create_nto: Only restricted case is implemented')

        nto_orbitals = nto_orbitals_list[0]
        nto_lambdas = nto_lambdas_list[0]

        assert_msg_critical(nto_orbitals.shape[1] == nto_lambdas.shape[0],
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

        nto_energies = np.zeros(nto_lambdas.shape[0])

        return MolecularOrbitals([nto_orbitals], [nto_energies], [nto_lambdas],
                                 mo_type)

    def is_nto(self):
        """
        Checks if this molecular orbitals object is natural transiton orbitals
        (NTO).

        :return:
            True if this molecular orbitals object is NTO, False otherwise.
        """

        assert_msg_critical(
            self.get_orbitals_type() == molorb.rest,
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
