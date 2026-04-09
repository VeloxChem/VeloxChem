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

import numpy as np

from .veloxchemlib import OverlapDriver
from .veloxchemlib import OverlapGeom100Driver
from .veloxchemlib import KineticEnergyDriver
from .veloxchemlib import KineticEnergyGeom100Driver
from .veloxchemlib import NuclearPotentialDriver
from .veloxchemlib import NuclearPotentialGeom100Driver
from .veloxchemlib import NuclearPotentialGeom010Driver
from .veloxchemlib import ECPGradientDriver
from .veloxchemlib import ElectricDipoleMomentDriver
from .veloxchemlib import NuclearPotentialGeom200Driver
from .veloxchemlib import NuclearPotentialGeom101Driver
from .veloxchemlib import compute_quadrupole_integrals
from .veloxchemlib import compute_linear_momentum_integrals
from .veloxchemlib import compute_angular_momentum_integrals
from .veloxchemlib import compute_electric_field_integrals
from .veloxchemlib import compute_electric_field_values
from .veloxchemlib import compute_electric_field_potential_gradient
from .veloxchemlib import compute_electric_field_fock_gradient
from .veloxchemlib import compute_electric_field_potential_gradient_for_mm
from .veloxchemlib import compute_electric_field_potential_hessian
from .matrices import Matrices


def compute_overlap_integrals(molecule, basis):
    """
    Computes overlap integrals.

    :param molecule:
        The molecule.
    :param basis:
        The molecular basis set.

    :return:
        The overlap integral matrix.
    """

    ovl_drv = OverlapDriver()
    ovl_mat = ovl_drv.compute(molecule, basis)

    return ovl_mat.to_numpy()


def compute_kinetic_energy_integrals(molecule, basis):
    """
    Computes kinetic energy integrals.

    :param molecule:
        The molecule.
    :param basis:
        The molecular basis set.

    :return:
        The kinetic energy integral matrix.
    """

    kin_drv = KineticEnergyDriver()
    kin_mat = kin_drv.compute(molecule, basis)

    return kin_mat.to_numpy()


def compute_nuclear_potential_integrals(molecule,
                                        basis,
                                        charges=None,
                                        coordinates=None):
    """
    Computes nuclear potential integrals.

    :param molecule:
        The molecule.
    :param basis:
        The molecular basis set.

    :return:
        The nuclear potential integral matrix.
    """

    npot_drv = NuclearPotentialDriver()

    if charges is None and coordinates is None:
        npot_mat = npot_drv.compute(molecule, basis)
    else:
        npot_mat = npot_drv.compute(molecule, basis, charges, coordinates)

    # Note: factor -1.0 for electron charge
    return -1.0 * npot_mat.to_numpy()


def compute_nuclear_potential_gradient(molecule,
                                       basis,
                                       D_total,
                                       atom_list=None):
    """
    Computes nuclear potential gradient.

    :param molecule:
        The molecule.
    :param basis:
        The molecular basis set.
    :param D_total:
        The sum of spin-alpha and spin-beta density matrix.
    :param atom_list:
        The list of atoms to be computed.

    :return:
        The nuclear potential gradient.
    """

    natoms = molecule.number_of_atoms()

    gradient = np.zeros((natoms, 3))

    if atom_list is None:
        atom_list = list(range(natoms))

    npot_grad_100_drv = NuclearPotentialGeom100Driver()
    npot_grad_010_drv = NuclearPotentialGeom010Driver()

    mol_charges = molecule.get_effective_nuclear_charges(basis)
    mol_coords = molecule.get_coordinates_in_bohr()

    for iatom in atom_list:
        gmats_100 = npot_grad_100_drv.compute(molecule, basis, iatom,
                                              mol_coords, mol_charges)
        gmats_010 = npot_grad_010_drv.compute(molecule, basis, iatom,
                                              mol_charges[iatom])

        for icart, label in enumerate(['X', 'Y', 'Z']):
            gmat_100 = gmats_100.matrix_to_numpy(label)
            gmat_010 = gmats_010.matrix_to_numpy(label)

            gradient[iatom, icart] += np.sum(
                (gmat_100 + gmat_100.T + gmat_010) * D_total)

    # Note: factor -1.0 for electron charge
    return -1.0 * gradient


def compute_point_charge_gradient(molecule,
                                  basis,
                                  D_total,
                                  charges,
                                  coordinates,
                                  atom_list=None):
    """
    Computes point charge contribution to molecular gradient.

    :param molecule:
        The molecule.
    :param basis:
        The molecular basis set.
    :param D_total:
        The sum of spin-alpha and spin-beta density matrix.
    :param charges:
        The point charges.
    :param coordinates:
        The coordinates of point charges.
    :param atom_list:
        The list of atoms to be computed.

    :return:
        The point charge contribution to molecular gradient.
    """

    natoms = molecule.number_of_atoms()

    gradient = np.zeros((natoms, 3))

    if atom_list is None:
        atom_list = list(range(natoms))

    npot_grad_100_drv = NuclearPotentialGeom100Driver()

    for iatom in atom_list:
        gmats_100 = npot_grad_100_drv.compute(molecule, basis, iatom,
                                              coordinates, charges)

        for icart, label in enumerate(['X', 'Y', 'Z']):
            gmat_100 = gmats_100.matrix_to_numpy(label)
            gradient[iatom, icart] += np.sum((gmat_100 + gmat_100.T) *
                                             D_total)

        gmats_100 = Matrices()

    # Note: factor -1.0 for electron charge
    return -1.0 * gradient


def compute_ecp_gradient(molecule, basis, D_total, atom_list=None):
    """
    Computes ECP gradient.

    :param molecule:
        The molecule.
    :param basis:
        The molecular basis set.
    :param D_total:
        The sum of spin-alpha and spin-beta density matrix.
    :param atom_list:
        The list of atoms to be computed.

    :return:
        The nuclear potential gradient.
    """

    natoms = molecule.number_of_atoms()

    gradient = np.zeros((natoms, 3))

    if atom_list is None:
        atom_list = list(range(natoms))

    ecp_grad_drv = ECPGradientDriver()

    core_electrons = basis.get_number_of_ecp_core_electrons()
    ecp_atom_indices = [
        idx for idx, nelec in enumerate(core_electrons) if nelec > 0
    ]

    for iatom in atom_list:
        gmats_100 = ecp_grad_drv.compute_bra_grad(molecule, basis,
                                                  ecp_atom_indices, iatom)
        if iatom in ecp_atom_indices:
            gmats_010 = ecp_grad_drv.compute_pot_grad(molecule, basis, iatom)

        for i, label in enumerate(['X', 'Y', 'Z']):
            gmat_100 = gmats_100.matrix_to_numpy(label)
            if iatom in ecp_atom_indices:
                gmat_010 = gmats_010.matrix_to_numpy(label)

            # Note different signs for 100 and 010 contributions
            gradient[iatom, i] += np.sum((gmat_100 + gmat_100.T) * D_total)
            if iatom in ecp_atom_indices:
                gradient[iatom, i] -= np.sum(gmat_010 * D_total)

    return gradient


def compute_overlap_gradient(molecule, basis, W_total, atom_list=None):
    """
    Computes overlap contribution to molecular gradient.

    :param molecule:
        The molecule.
    :param basis:
        The molecular basis set.
    :param W_total:
        The total energy-weighted density matrix.
    :param atom_list:
        The list of atoms to be computed.

    :return:
        The overlap contribution to molecular gradient.
    """

    natoms = molecule.number_of_atoms()

    gradient = np.zeros((natoms, 3))

    if atom_list is None:
        atom_list = list(range(natoms))

    ovl_grad_drv = OverlapGeom100Driver()

    for iatom in atom_list:
        gmats = ovl_grad_drv.compute(molecule, basis, iatom)

        for icart, label in enumerate(['X', 'Y', 'Z']):
            gmat = gmats.matrix_to_numpy(label)
            gradient[iatom, icart] += np.sum((gmat + gmat.T) * W_total)

        gmats = Matrices()

    # Note: minus sign for energy-weighted density contribution
    return -1.0 * gradient


def compute_kinetic_energy_gradient(molecule, basis, D_total, atom_list=None):
    """
    Computes kinetic energy contribution to molecular gradient.

    :param molecule:
        The molecule.
    :param basis:
        The molecular basis set.
    :param D_total:
        The sum of spin-alpha and spin-beta density matrix.
    :param atom_list:
        The list of atoms to be computed.

    :return:
        The kinetic energy contribution to molecular gradient.
    """

    natoms = molecule.number_of_atoms()

    gradient = np.zeros((natoms, 3))

    if atom_list is None:
        atom_list = list(range(natoms))

    kin_grad_drv = KineticEnergyGeom100Driver()

    for iatom in atom_list:
        gmats = kin_grad_drv.compute(molecule, basis, iatom)

        for icart, label in enumerate(['X', 'Y', 'Z']):
            gmat = gmats.matrix_to_numpy(label)
            gradient[iatom, icart] += np.sum((gmat + gmat.T) * D_total)

        gmats = Matrices()

    return gradient


def compute_electric_dipole_integrals(molecule, basis, origin=(0.0, 0.0, 0.0)):
    """
    Computes electric dipole integrals.

    :param molecule:
        The molecule.
    :param basis:
        The molecular basis set.

    :return:
        A tuple containint the electric dipole integral matrices.
    """

    dip_drv = ElectricDipoleMomentDriver()
    dip_mats = dip_drv.compute(molecule, basis, list(origin))

    # Note: factor -1.0 for electron charge
    return tuple([
        -1.0 * dip_mats.matrix_to_numpy('X'),
        -1.0 * dip_mats.matrix_to_numpy('Y'),
        -1.0 * dip_mats.matrix_to_numpy('Z'),
    ])


def compute_nuclear_potential_gradient_bfs(molecule, basis, charges,
                                           coordinates, density):
    """
    Computes nuclear potential integrals contribution from point charges to
    molecular gradient.

    :param molecule:
        The molecule.
    :param basis:
        The molecular basis set.
    :param charges:
        The point charges.
    :param coordinates:
        The coordinates of point charges.

    :return:
        The nuclear potential integral contribution to molecular gradient.
    """

    natoms = molecule.number_of_atoms()

    grad = np.zeros((natoms, 3))

    npot_grad_100_drv = NuclearPotentialGeom100Driver()

    for iatom in range(natoms):
        gmats_100 = npot_grad_100_drv.compute(molecule, basis, iatom,
                                              coordinates, charges)

        for i, label in enumerate(['X', 'Y', 'Z']):
            gmat_100 = gmats_100.matrix_to_numpy(label)

            grad[iatom, i] += np.sum((gmat_100 + gmat_100.T) * density)

        gmats_100 = Matrices()

    # Note: factor -1.0 for electron charge
    return -1.0 * grad


def compute_electrostatic_potential_hessian(molecule, basis, mm_charges,
                                            mm_coordinates, density,
                                            qm_atom_index_i, qm_atom_index_j):

    hess = np.zeros((3, 3))

    i, j = qm_atom_index_i, qm_atom_index_j

    if i == j:
        npot_hess_200_drv = NuclearPotentialGeom200Driver()

        hmats_200 = npot_hess_200_drv.compute(molecule, basis, i,
                                              mm_coordinates, mm_charges)

        for x, label_x in enumerate('XYZ'):
            for y, label_y in enumerate('XYZ'):
                npot_label = label_x + label_y if x <= y else label_y + label_x
                npot_200_iixy = hmats_200.matrix_to_numpy(npot_label)
                hess[x,
                     y] += np.sum(density * (npot_200_iixy + npot_200_iixy.T))

        hmats_200 = Matrices()

    npot_hess_101_drv = NuclearPotentialGeom101Driver()

    hmats_101 = npot_hess_101_drv.compute(molecule, basis, i, j, mm_coordinates,
                                          mm_charges)

    for x, label_x in enumerate('XYZ'):
        for y, label_y in enumerate('XYZ'):
            npot_xy_label = f'{label_x}_{label_y}'
            npot_101_ijxy = hmats_101.matrix_to_numpy(npot_xy_label)
            hess[x, y] += np.sum(density * (npot_101_ijxy + npot_101_ijxy.T))

    hmats_101 = Matrices()

    # Note: factor -1.0 for electron charge
    return -1.0 * hess


def compute_electrostatic_integrals_gradient(molecule, basis, mm_charges,
                                             mm_coordinates, qm_atom_index_i):

    naos = basis.get_dimensions_of_basis()

    ints_grad = np.zeros((3, naos, naos))

    i = qm_atom_index_i

    npot_grad_100_drv = NuclearPotentialGeom100Driver()

    gmats_100 = npot_grad_100_drv.compute(molecule, basis, i, mm_coordinates,
                                          mm_charges)

    for x, label in enumerate(['X', 'Y', 'Z']):
        gmat_100 = gmats_100.matrix_to_numpy(label)
        # Note: factor -1.0 for electron charge
        ints_grad[x] -= gmat_100 + gmat_100.T

    gmats_100 = Matrices()

    return ints_grad
