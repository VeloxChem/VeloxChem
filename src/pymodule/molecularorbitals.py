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

import numpy as np
import h5py
import math
import sys

from .veloxchemlib import MolecularOrbitals
from .veloxchemlib import molorb
from .veloxchemlib import to_angular_momentum
from .veloxchemlib import denmat
from .aodensitymatrix import AODensityMatrix
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical


def _MolecularOrbitals_print_orbitals(self,
                                      molecule,
                                      ao_basis,
                                      all_orbs=False,
                                      ostream=None):
    """
    Prints molecular orbitals to output stream.

    :param molecule:
        The molecule.
    :param ao_basis:
        The AO basis set.
    :param all_orbs:
        The flag to print all orbitals.
    :param ostream:
        The output stream.
    """

    if ostream is None:
        ostream = OutputStream(sys.stdout)

    norb = self.number_mos()

    ao_map = ao_basis.get_ao_basis_map(molecule)

    if self.get_orbitals_type() == molorb.rest:

        ostream.print_blank()

        ostream.print_header("Spin Restricted Orbitals")
        ostream.print_header("------------------------")

        nocc = molecule.number_of_electrons() // 2

        if all_orbs:
            nstart, nend = 0, norb
        else:
            nstart, nend = max(0, nocc - 5), min(norb, nocc + 5)

        rvecs = self.alpha_to_numpy()
        reigs = self.ea_to_numpy()
        rnocc = [2.0 if x < nocc else 0.0 for x in range(norb)]

        for i in range(nstart, nend):
            _MolecularOrbitals_print_coefficients(reigs[i], rnocc[i], i,
                                                  rvecs[:, i], ao_map, 0.15,
                                                  ostream)

        ostream.print_blank()

    elif self.get_orbitals_type() == molorb.unrest:

        ostream.print_blank()

        ostream.print_header("Spin Unrestricted Alpha Orbitals")
        ostream.print_header("--------------------------------")

        nalpha = molecule.number_of_alpha_electrons()

        if all_orbs:
            nstart, nend = 0, norb
        else:
            nstart, nend = max(0, nalpha - 5), min(norb, nalpha + 5)

        uvecs = self.alpha_to_numpy()
        ueigs = self.ea_to_numpy()
        unocc = [1.0 if x < nalpha else 0.0 for x in range(norb)]

        for i in range(nstart, nend):
            _MolecularOrbitals_print_coefficients(ueigs[i], unocc[i], i,
                                                  uvecs[:, i], ao_map, 0.15,
                                                  ostream)

        ostream.print_blank()

        ostream.print_header("Spin Unrestricted Beta Orbitals")
        ostream.print_header("-------------------------------")

        nbeta = molecule.number_of_beta_electrons()

        if all_orbs:
            nstart, nend = 0, norb
        else:
            nstart, nend = max(0, nbeta - 5), min(norb, nbeta + 5)

        uvecs = self.beta_to_numpy()
        ueigs = self.eb_to_numpy()
        unocc = [1.0 if x < nbeta else 0.0 for x in range(norb)]

        for i in range(nstart, nend):
            _MolecularOrbitals_print_coefficients(ueigs[i], unocc[i], i,
                                                  uvecs[:, i], ao_map, 0.15,
                                                  ostream)

        ostream.print_blank()

    else:

        errmsg = "MolecularOrbitals.get_density:"
        errmsg += " Invalid molecular orbitals type"
        assert_msg_critical(False, errmsg)


def _MolecularOrbitals_print_coefficients(eigval, focc, iorb, coeffs, ao_map,
                                          thresh, ostream):
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

    valstr = "Occupation: {:.1f} Energy: {:10.5f} a.u.".format(focc, eigval)
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


def _MolecularOrbitals_get_density(self, molecule, scf_type):
    """
    Gets AO density matrix from molecular orbitals.

    :param molecule:
        The molecule.
    :param scf_type:
        The type of SCF calculation (restricted, unrestricted, or restricted
        open-shell).

    :return:
        The AO density matrix.
    """

    nalpha = molecule.number_of_alpha_electrons()
    nbeta = molecule.number_of_beta_electrons()

    if (self.get_orbitals_type() == molorb.rest and scf_type == 'restricted'):
        return self.get_ao_density(nalpha + nbeta)

    elif (self.get_orbitals_type() == molorb.unrest and
          scf_type == 'unrestricted'):
        return self.get_ao_density(nalpha, nbeta)

    elif (self.get_orbitals_type() == molorb.rest and
          scf_type == 'restricted_openshell'):
        mo_coef = self.alpha_to_numpy()
        mo_occ_alpha = mo_coef[:, :nalpha]
        mo_occ_beta = mo_coef[:, :nbeta]
        dalpha = np.matmul(mo_occ_alpha, mo_occ_alpha.T)
        dbeta = np.matmul(mo_occ_beta, mo_occ_beta.T)
        return AODensityMatrix([dalpha, dbeta], denmat.unrest)

    else:
        errmsg = "MolecularOrbitals.get_density:"
        errmsg += " Invalid molecular orbitals type"
        assert_msg_critical(False, errmsg)


def _MolecularOrbitals_write_hdf5(self,
                                  fname,
                                  nuclear_charges=None,
                                  basis_set=None):
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

    hf.create_dataset('alpha_orbitals',
                      data=self.alpha_to_numpy(),
                      compression='gzip')

    hf.create_dataset('alpha_energies',
                      data=self.ea_to_numpy(),
                      compression='gzip')

    if self.get_orbitals_type() == molorb.unrest:

        hf.create_dataset('beta_orbitals',
                          data=self.beta_to_numpy(),
                          compression='gzip')

        hf.create_dataset('beta_energies',
                          data=self.eb_to_numpy(),
                          compression='gzip')

    if nuclear_charges is not None:
        hf.create_dataset('nuclear_charges',
                          data=nuclear_charges,
                          compression='gzip')

    if basis_set is not None:
        hf.create_dataset('basis_set',
                          data=np.string_([basis_set]),
                          compression='gzip')

    hf.close()


@staticmethod
def _MolecularOrbitals_read_hdf5(fname):
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
        'alpha_orbitals' in hf and 'alpha_energies' in hf,
        'MolecularOrbitals.read_hdf5: alpha orbitals/energies not found')

    if 'beta_orbitals' in hf or 'beta_energies' in hf:
        orbs_type = molorb.unrest

        assert_msg_critical(
            'beta_orbitals' in hf and 'beta_energies' in hf,
            'MolecularOrbitals.read_hdf5: beta orbitals/energies not found')

    orbs = []
    enes = []

    orbs.append(np.array(hf.get('alpha_orbitals')))
    enes.append(np.array(hf.get('alpha_energies')))

    if orbs_type == molorb.unrest:
        orbs.append(np.array(hf.get('beta_orbitals')))
        enes.append(np.array(hf.get('beta_energies')))

    hf.close()

    return MolecularOrbitals(orbs, enes, orbs_type)


@staticmethod
def _MolecularOrbitals_match_hdf5(fname, nuclear_charges, basis_set, scf_type):
    """
    Checks if the hdf5 file matches the given nuclear charges and basis set.

    :param fname:
        The name of the hdf5 file.
    :param nuclear_charges:
        The nuclear charges.
    :param basis_set:
        Name of the basis set.
    :param scf_type:
        The type of SCF calculation (restricted, unrestricted, or restricted
        open-shell).

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

    h5_restricted = ('beta_orbitals' not in hf and 'beta_energies' not in hf)
    restricted = (scf_type in ['restricted', 'restricted_openshell'])
    match_restricted = (h5_restricted == restricted)

    hf.close()

    return (match_nuclear_charges and match_basis_set and match_restricted)


MolecularOrbitals.print_orbitals = _MolecularOrbitals_print_orbitals
MolecularOrbitals.get_density = _MolecularOrbitals_get_density
MolecularOrbitals.write_hdf5 = _MolecularOrbitals_write_hdf5
MolecularOrbitals.read_hdf5 = _MolecularOrbitals_read_hdf5
MolecularOrbitals.match_hdf5 = _MolecularOrbitals_match_hdf5
