#
#                           VELOXCHEM 1.0-RC3
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2022 by VeloxChem developers. All rights reserved.
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

from pathlib import Path

import numpy as np
from mpi4py import MPI

try:
    import cppe
except ImportError:
    raise ImportError('Unable to import cppe. Please install cppe via ' +
                      '\'python3 -m pip install cppe\'')

from .errorhandler import assert_msg_critical
from .inputparser import InputParser
from .subcommunicators import SubCommunicators
from .veloxchemlib import (ElectricFieldIntegralsDriver,
                           NuclearPotentialIntegralsDriver, bohr_in_angstrom,
                           mpi_master)


class PolEmbed:
    """
    Implements interface to the CPPE library.

    :param molecule:
        The molecule.
    :param basis:
        The AO basis set.
    :param pe_dict:
        The dictionary with options for CPPE.
    :param comm:
        The MPI communicator.

    Instance variables
        - molecule: The molecule.
        - basis: The AO basis set.
        - comm: The MPI communicator.
        - rank: The rank of MPI process.
        - nodes: The number of MPI processes.
        - V_es: The multipole contribution to Fock matrix.
        - options: The dictionary with options for CPPE.
        - cppe_state: The CPPE state object.
        - polarizable_coords: The coordinates of the polarizable sites.
    """

    def __init__(self, molecule, basis, pe_dict, comm=None):
        """
        Initializes interface to the CPPE library.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        self.molecule = molecule
        self.basis = basis
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()
        self.output = ''
        self.V_es = None

        self.update_options(pe_dict)

        cppe_mol = cppe.Molecule()
        coords = np.vstack(
            (self.molecule.x_to_numpy(), self.molecule.y_to_numpy(),
             self.molecule.z_to_numpy())).T
        charges = self.molecule.elem_ids_to_numpy()

        for z, coord in zip(charges, coords):
            cppe_mol.append(cppe.Atom(int(z), *coord))
        self.cppe_state = cppe.CppeState(self.options, cppe_mol,
                                         self.print_callback)
        self.cppe_state.calculate_static_energies_and_fields()
        self._enable_induction = False
        if self.cppe_state.get_polarizable_site_number():
            self._enable_induction = True
            coords = np.array([
                site.position
                for site in self.cppe_state.potentials
                if site.is_polarizable
            ])
            self.polarizable_coords = coords

    def get_cppe_version(self):
        """
        Gets the version of CPPE.

        :return:
            The version of CPPE.
        """

        return str(cppe.__version__)

    def update_options(self, pe_dict):
        """
        Updates PE options based on the input dictionary.

        :param pe_dict:
            The dictionary with options for CPPE.
        """

        keytypes = {
            'potfile': 'string',
            'iso_pol': 'bool',
            'induced_thresh': 'float',
            'maxiter': 'int',
            'damp_induced': 'bool',
            'damp_multipole': 'bool',
            'damping_factor_induced': 'float',
            'damping_factor_multipole': 'float',
            'pe_border': 'bool',
            'border_type': 'string',
            'border_rmin': 'float',
            'border_nredist': 'int',
            'border_redist_order': 'int',
            'border_redist_pol': 'bool',
        }

        self.options = {}

        for key in pe_dict:
            if keytypes[key] == 'string':
                self.options[key] = str(pe_dict[key])
            elif keytypes[key] == 'float':
                self.options[key] = float(pe_dict[key])
            elif keytypes[key] == 'int':
                self.options[key] = int(pe_dict[key])
            elif keytypes[key] == 'bool':
                boolstr = pe_dict[key].lower()
                self.options[key] = (boolstr in ['yes', 'y'])

    def print_callback(self, output):
        """
        Handles the output from the CPPE library.

        :param output:
            The output from the CPPE library.
        """

        self.output += output

    def get_pe_contribution(self, dm, elec_only=False):
        """
        Computes contributions from polarizable embedding.

        :param dm:
            The density matrix.
        :param elec_only:
            The flag for computing electronic contribution only.

        :return:
            The polarizable embedding contributions to energy and Fock matrix.
        """

        if self.V_es is None:
            self.V_es = self.compute_multipole_potential_integrals()

        if not elec_only:
            if self.rank == mpi_master():
                e_el = np.sum(self.V_es * dm)
                self.cppe_state.energies["Electrostatic"]["Electronic"] = e_el

        V_ind = np.zeros(self.V_es.shape)

        if self._enable_induction:
            elec_fields = self.compute_electric_field_value(dm)
            # solve induced moments
            if self.rank == mpi_master():
                self.cppe_state.update_induced_moments(elec_fields.flatten(),
                                                       elec_only)
                induced_moments = np.array(
                    self.cppe_state.get_induced_moments()).reshape(
                        self.polarizable_coords.shape)
            else:
                induced_moments = None
            induced_moments = self.comm.bcast(induced_moments,
                                              root=mpi_master())
            V_ind = self.compute_induction_operator(induced_moments)

        if self.rank == mpi_master():
            if not elec_only:
                vmat = self.V_es + V_ind
                e = self.cppe_state.total_energy
            else:
                vmat = V_ind
                e = self.cppe_state.energies["Polarization"]["Electronic"]
            return e, vmat
        else:
            return 0.0, None

    def compute_multipole_potential_integrals(self):
        """
        Computes contribution from multipoles.

        :return:
            The multipole contribution to Fock matrix.
        """

        sites = []
        dipole_sites = []
        charges = []
        dipoles = []

        for p in self.cppe_state.potentials:
            site = p.position
            for m in p.multipoles:
                # zeroth order
                if m.k == 0:
                    charges.append(m.values[0])
                    sites.append(site)
                # first order
                elif m.k == 1:
                    dipoles.append(m.values)
                    dipole_sites.append(site)
                else:
                    raise NotImplementedError(
                        "PE electrostatics only implemented through first order"
                    )

        V_es = np.zeros(0)

        npot_drv = NuclearPotentialIntegralsDriver(self.comm)
        if self.rank == mpi_master():
            V_es = -1.0 * npot_drv.compute(self.molecule, self.basis,
                                           np.array(charges),
                                           np.array(sites)).to_numpy()

        if dipole_sites:
            ef_driver = ElectricFieldIntegralsDriver(self.comm)
            ret = ef_driver.compute(self.molecule, self.basis,
                                    np.array(dipoles), np.array(dipole_sites))
            if self.rank == mpi_master():
                V_es += -1.0 * (ret.x_to_numpy() + ret.y_to_numpy() +
                                ret.z_to_numpy())

        return V_es

    def compute_induction_operator(self, moments):
        """
        Computes contribution from induction operator.

        :param moments:
            The induced moments.

        :return:
            The induction contribution to Fock matrix.
        """

        node_grps = [p for p in range(self.nodes)]
        subcomm = SubCommunicators(self.comm, node_grps)
        local_comm = subcomm.local_comm
        cross_comm = subcomm.cross_comm

        ave, res = divmod(self.polarizable_coords.shape[0], self.nodes)
        counts = [ave + 1 if p < res else ave for p in range(self.nodes)]

        start = sum(counts[:self.rank])
        end = sum(counts[:self.rank + 1])

        ef_driver = ElectricFieldIntegralsDriver(local_comm)
        if start < end:
            ret = ef_driver.compute(self.molecule, self.basis,
                                    moments[start:end, :],
                                    self.polarizable_coords[start:end, :])
            V_size = ret.x_to_numpy().shape[0]
        else:
            V_size = 0

        if self.polarizable_coords.shape[0] < self.nodes:
            V_size = cross_comm.bcast(V_size, root=mpi_master())

        if local_comm.Get_rank() == mpi_master():
            if start < end:
                V_ind = -1.0 * (ret.x_to_numpy() + ret.y_to_numpy() +
                                ret.z_to_numpy())
            else:
                V_ind = np.zeros((V_size, V_size))
            V_ind = cross_comm.reduce(V_ind, root=mpi_master())
        else:
            V_ind = np.zeros(0)

        return V_ind

    def compute_electric_field_value(self, dm):
        """
        Computes electric field at the polarizable sites.

        :param dm:
            The density matrix.

        :return:
            The electric field.
        """

        node_grps = [p for p in range(self.nodes)]
        subcomm = SubCommunicators(self.comm, node_grps)
        local_comm = subcomm.local_comm
        cross_comm = subcomm.cross_comm

        ave, res = divmod(self.polarizable_coords.shape[0], self.nodes)
        counts = [ave + 1 if p < res else ave for p in range(self.nodes)]

        start = sum(counts[:self.rank])
        end = sum(counts[:self.rank + 1])

        ef_driver = ElectricFieldIntegralsDriver(local_comm)
        elec_field = []
        for i in range(start, end):
            coord = self.polarizable_coords[i]
            ret = ef_driver.compute(self.molecule, self.basis, *coord)
            if local_comm.Get_rank() == mpi_master():
                elec_field.append(np.sum(dm * ret.x_to_numpy()))
                elec_field.append(np.sum(dm * ret.y_to_numpy()))
                elec_field.append(np.sum(dm * ret.z_to_numpy()))

        if local_comm.Get_rank() == mpi_master():
            elec_field = cross_comm.gather(elec_field, root=mpi_master())
        else:
            elec_field = []

        if self.rank == mpi_master():
            elec_field = np.array([xyz for ef in elec_field for xyz in ef])
            elec_field = elec_field.reshape(self.polarizable_coords.shape)

        return elec_field

    @staticmethod
    def write_cppe_potfile(potfile):
        """
        Writes potential file for the CPPE module.

        :param potfile:
            The name of the potential file.

        :return:
            The name of the temporary potential file.
        """

        with open(potfile) as f_pot:
            first_line = f_pot.readline().strip()
            while (not first_line) or first_line.startswith('!'):
                first_line = f_pot.readline().strip()

        if first_line.split()[0] == '@COORDINATES':
            return potfile

        # read potential file (new format)

        pe_inp = InputParser(potfile).input_dict
        natoms = len(pe_inp['environment']['xyz'])

        # process units

        if 'units' in pe_inp['environment']:
            units = pe_inp['environment']['units'].lower()
            if len(units) >= 3 and units == 'angstrom'[:len(units)]:
                prefac = 1.0
            elif units in ['au', 'bohr']:
                prefac = bohr_in_angstrom()
            else:
                assert_msg_critical(False, 'potential file: invalid units')
        else:
            prefac = 1.0

        # process residues (a dictionary with residue ID as key)

        residues = {}
        for line in pe_inp['environment']['xyz']:
            val = line.split()
            resname = val[4]
            resid = int(val[5])
            if resid not in residues:
                residues[resid] = {'resname': resname, 'atoms': []}
            residues[resid]['atoms'].append([
                val[0],
                float(val[1]) * prefac,
                float(val[2]) * prefac,
                float(val[3]) * prefac,
            ])

        # process charges (a dictionary with residue name as key)

        charges = {}
        for line in pe_inp['charges']:
            val = line.split()
            resname = val[2]
            if resname not in charges:
                charges[resname] = []
            charges[resname].append(val[1])

        # process polarizabilities (a dictionary with residue name as key)

        polarizabilities = {}
        for line in pe_inp.get('polarizabilities', ''):
            val = line.split()
            resname = val[7]
            if resname not in polarizabilities:
                polarizabilities[resname] = []
            polarizabilities[resname].append(val[1:7])

        # write potential file (old format)

        cppe_potfile = Path(potfile).with_suffix('.cppe.pot')

        with open(str(cppe_potfile), 'w') as f_tmp:

            # coordinates

            f_tmp.write('@COORDINATES\n')
            f_tmp.write('{:d}\n'.format(natoms))
            f_tmp.write('AA\n')

            for resid in sorted(list(residues.keys())):
                for atom in residues[resid]['atoms']:
                    f_tmp.write('{:<2s}{:12.7f}{:12.7f}{:12.7f}\n'.format(
                        atom[0], atom[1], atom[2], atom[3]))

            # charges

            q_lines = ''
            q_num = 0

            ind = 0
            for resid in sorted(list(residues.keys())):
                resname = residues[resid]['resname']
                if resname in charges:
                    for k, q in enumerate(charges[resname]):
                        q_lines += '{:<5d}{:13.8f}\n'.format(
                            ind + k + 1, float(q))
                        q_num += 1
                ind += len(residues[resid]['atoms'])

            f_tmp.write('@MULTIPOLES\n')
            f_tmp.write('ORDER 0\n')
            f_tmp.write('{:d}\n'.format(q_num))
            f_tmp.write(q_lines)

            # polarizabilities

            alpha_lines = ''
            alpha_num = 0

            ind = 0
            for resid in sorted(list(residues.keys())):
                resname = residues[resid]['resname']
                if resname in polarizabilities:
                    for k, alpha in enumerate(polarizabilities[resname]):
                        alpha_lines += '{:<5d}'.format(ind + k + 1)
                        for a in alpha:
                            alpha_lines += '{:13.8f}'.format(float(a))
                        alpha_lines += '\n'
                        alpha_num += 1
                ind += len(residues[resid]['atoms'])

            if alpha_num > 0:
                f_tmp.write('@POLARIZABILITIES\n')
                f_tmp.write('ORDER 1 1\n')
                f_tmp.write('{:d}\n'.format(alpha_num))
                f_tmp.write(alpha_lines)

            # exclusions

            max_natoms_res = 0
            for resid in sorted(list(residues.keys())):
                if max_natoms_res < len(residues[resid]['atoms']):
                    max_natoms_res = len(residues[resid]['atoms'])

            f_tmp.write('EXCLISTS\n')
            f_tmp.write('{:<6d}{:<6d}\n'.format(natoms, max_natoms_res))

            ind = 0
            for resid in sorted(list(residues.keys())):
                natoms_res = len(residues[resid]['atoms'])
                for k in range(natoms_res):
                    f_tmp.write('{:<6d}'.format(ind + k + 1))
                    for j in range(natoms_res):
                        if j != k:
                            f_tmp.write('{:<6d}'.format(ind + j + 1))
                    for j in range(max_natoms_res - natoms_res):
                        f_tmp.write('{:<6d}'.format(0))
                    f_tmp.write('\n')
                ind += natoms_res

        return str(cppe_potfile)
