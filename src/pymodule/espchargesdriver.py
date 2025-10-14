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
from pathlib import Path
import numpy as np
import sys

from .veloxchemlib import compute_nuclear_potential_values
from .veloxchemlib import (mpi_master, bohr_in_angstrom, hartree_in_ev,
                           boltzmann_in_evperkelvin)
from .outputstream import OutputStream
from .molecule import Molecule
from .molecularbasis import MolecularBasis
from .scfrestdriver import ScfRestrictedDriver
from .scfunrestdriver import ScfUnrestrictedDriver
from .inputparser import parse_input, print_keywords
from .errorhandler import assert_msg_critical, safe_solve


class EspChargesDriver:
    """
    Implements ESP charges driver.

    # vlxtag: RHF, ESP_charges
    # vlxtag: RHF, ESP_on_points

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables:
        - comm: The MPI communicator.
        - rank: The MPI rank.
        - nodes: The number of MPI processes.
        - ostream: The output stream.
        - number_layers: The number of layers of scaled van der Waals surfaces.
        - density: The density of grid points in points per square Angstrom.
        - equal_charges: The charges that are constrained to be equal.
        - weights: The weight factors of different conformers.
        - energies: The energies of different conformers for Boltzmann weight factors.
        - temperature: The temperature for Boltzmann weight factors.
        - net_charge: The charge of the molecule.
        - multiplicity: The multiplicity of the molecule.
        - method_dict: The dictionary of method settings.
        - filename: The filename for the calculation.
        - xyz_file: The xyz file containing the conformers.
        - fitting_points: The fitting point on which charges are calculated.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes ESP charges driver.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        # mpi information
        self.comm = comm
        self.rank = comm.Get_rank()
        self.nodes = comm.Get_size()

        # outputstream
        self.ostream = ostream

        # filename
        self.filename = None

        # method
        self.restart = True
        self.conv_thresh = 1.0e-6
        self.xcfun = None
        self.grid_level = None
        self.solvation_model = None

        # conformers
        self.xyz_file = None
        self.net_charge = 0.0
        self.multiplicity = 1
        self.weights = None
        self.energies = None
        self.temperature = 298.15

        # grid information
        self.grid_type = 'mk'
        self.custom_mk_radii = None

        # MK grid settings (density in Angstrom^-2)
        self.number_layers = 4
        self.density = 1.0

        # CHELPG grid settings (in Angstrom)
        self.chelpg_spacing = 0.3
        self.chelpg_margin = 2.8

        # esp fitting to user specified points
        self.fitting_points = None

        # esp fitting
        self.equal_charges = None

        # input keywords
        self._input_keywords = {
            'esp_charges': {
                'restart': ('bool', 'restart flag for SCF driver'),
                'conv_thresh':
                    ('float', 'convergence threshold for SCF driver'),
                'grid_type': ('str_lower', 'type of grid (mk or chelpg)'),
                'number_layers':
                    ('int', 'number of layers of scaled vdW surfaces'),
                'density': ('float', 'density of points in each layer'),
                'chelpg_spacing': ('float', 'step size for CHELPG grid'),
                'chelpg_margin': ('float', 'maximum radius for CHELPG grid'),
                'equal_charges': ('str', 'constraints for equal charges'),
                'custom_mk_radii':
                    ('seq_fixed_str', 'custom MK radii for ESP charges'),
                'xyz_file': ('str', 'xyz file containing the conformers'),
                'net_charge': ('float', 'net charge of the molecule'),
                'multiplicity': ('int', 'spin multiplicity of the molecule'),
                'weights': ('seq_fixed', 'weight factors of conformers'),
                'energies': ('seq_fixed', 'energies of conformers'),
                'temperature':
                    ('float', 'temperature for Boltzmann weight factor'),
                'fitting_points':
                    ('list', 'points on which ESP charges are computed'),
            },
            'method_settings': {
                'xcfun': ('str_upper', 'exchange-correlation functional'),
                'grid_level': ('int', 'accuracy level of DFT grid (1-8)'),
                'solvation_model': ('str', 'solvation model'),
            },
        }

    def print_keywords(self):
        """
        Prints input keywords in SCF driver.
        """

        print_keywords(self._input_keywords, self.ostream)

    def update_settings(self, esp_dict, method_dict=None):
        """
        Updates settings in ESP charges driver.

        :param esp_dict:
            The dictionary of ESP charges input.
        :param method_dict:
            The dicitonary of method settings.
        """

        esp_keywords = {
            key: val[0]
            for key, val in self._input_keywords['esp_charges'].items()
        }

        parse_input(self, esp_keywords, esp_dict)

        method_keywords = {
            key: val[0]
            for key, val in self._input_keywords['method_settings'].items()
        }

        parse_input(self, method_keywords, method_dict)

        if 'filename' in esp_dict:
            self.filename = esp_dict['filename']

    def compute(self, molecule, basis=None, scf_results=None, flag=None):
        """
        Computes ESP charges.

        :param molecule:
            The molecule.
        :param basis:
            The molecular basis set.
        :param scf_results:
            The SCF results.

        :return:
            The charges.
        """

        # backward compatibility
        assert_msg_critical(
            flag is None or flag.lower() == 'esp',
            f'{type(self).__name__}.compute: Use of resp/esp flag is ' +
            'deprecated. Please use either RespChargesDriver or ' +
            'EspChargesDriver')

        if isinstance(molecule, list):
            # conformer-weighted esp charges
            if scf_results is not None:
                cur_str = 'Ignoring scf_results since multiple molecules were '
                cur_str += 'provided.'
                self.ostream.print_warning(cur_str)
            molecule_list = molecule
            return self.compute_multiple_molecules(molecule_list, basis)
        else:
            # single molecule esp charges
            return self.compute_single_molecule(molecule, basis, scf_results)

    def compute_single_molecule(self, molecule, basis=None, scf_results=None):
        """
        Computes ESP charges for a single molecule.

        :param molecule:
            The molecule.
        :param basis:
            The molecular basis set.
        :param scf_results:
            The SCF results.

        :return:
            The ESP charges.
        """

        # sanity check for equal_charges
        if self.equal_charges is not None:
            if isinstance(self.equal_charges, str):
                eq_chgs = []
                for q in self.equal_charges.split(','):
                    if q:
                        eq_chgs.append(list(map(int, q.split('='))))
                self.equal_charges = eq_chgs
            assert_msg_critical(
                isinstance(self.equal_charges, (list, tuple)),
                f'{type(self).__name__}.compute: Invalid equal_charges')

        # check basis and scf_results
        if self.rank == mpi_master():
            need_basis = (basis is None)
            need_scf = (scf_results is None)
        else:
            need_basis = None
            need_scf = None
        need_basis, need_scf = self.comm.bcast((need_basis, need_scf),
                                               root=mpi_master())

        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()

        # default basis is 6-31g*
        if need_basis:
            if self.rank == mpi_master():
                basis = MolecularBasis.read(molecule,
                                            '6-31g*',
                                            ostream=self.ostream)
            basis = self.comm.bcast(basis, root=mpi_master())

        # run SCF if needed
        if need_scf:
            if nalpha == nbeta:
                scf_drv = ScfRestrictedDriver(self.comm, self.ostream)
            else:
                scf_drv = ScfUnrestrictedDriver(self.comm, self.ostream)
            if self.filename is not None:
                scf_drv.filename = self.filename
            scf_drv.restart = self.restart
            scf_drv.conv_thresh = self.conv_thresh
            scf_drv.xcfun = self.xcfun
            scf_drv.grid_level = self.grid_level
            scf_drv.solvation_model = self.solvation_model
            scf_results = scf_drv.compute(molecule, basis)

        # sanity check
        if (self.rank == mpi_master()) and (nalpha != nbeta):
            assert_msg_critical(
                scf_results['scf_type'] != 'restricted',
                f'{type(self).__name__}.compute: Molecule is open-shell ' +
                'while SCF is closed-shell')

        grid_m = self.get_grid_points(molecule)
        esp_m = self.get_electrostatic_potential(grid_m, molecule, basis,
                                                 scf_results)

        q = None

        if self.rank == mpi_master():
            if self.fitting_points is not None:
                q = self.compute_esp_charges_on_fitting_points([molecule],
                                                               [grid_m],
                                                               [esp_m])
            else:
                q = self.compute_esp_charges([molecule], [grid_m], [esp_m],
                                             [1.0])

        return q

    def compute_multiple_molecules(self, molecule_list, basis_set_list=None):
        """
        Computes conformer-weighted ESP charges.

        :param molecule_list:
            The list of molecules.
        :param basis_set_list:
            The list of molecular basis set.

        :return:
            The charges.
        """

        # sanity check for fitting_points
        assert_msg_critical(
            self.fitting_points is None,
            f'{type(self).__name__}.compute: Cannot compute ESP charges on fitting '
            + 'points for multiple conformers')

        # sanity check for equal_charges
        if self.equal_charges is not None:
            if isinstance(self.equal_charges, str):
                eq_chgs = []
                for q in self.equal_charges.split(','):
                    if q:
                        eq_chgs.append(list(map(int, q.split('='))))
                self.equal_charges = eq_chgs
            assert_msg_critical(
                isinstance(self.equal_charges, (list, tuple)),
                f'{type(self).__name__}.compute: Invalid equal_charges')

        molecules = molecule_list
        basis_sets = basis_set_list

        if basis_sets is None:
            if self.rank == mpi_master():
                basis_sets = [
                    MolecularBasis.read(mol, '6-31g*', verbose=False)
                    for mol in molecules
                ]
            basis_sets = self.comm.bcast(basis_sets, root=mpi_master())

        if self.rank == mpi_master():
            if self.weights is not None and self.energies is not None:
                # avoid conflict between weights and energies
                errmsg = f'{type(self).__name__}: Please do not specify weights and '
                errmsg += 'energies at the same time.'
                assert_msg_critical(False, errmsg)

            elif self.weights is not None:
                errmsg = f'{type(self).__name__}: Number of weights does not match '
                errmsg += 'number of conformers.'
                assert_msg_critical(len(self.weights) == len(molecules), errmsg)
                # normalize weights
                sum_of_weights = sum(self.weights)
                self.weights = [w / sum_of_weights for w in self.weights]

            elif self.energies is not None:
                errmsg = f'{type(self).__name__}: Number of energies does not match '
                errmsg += 'number of conformers.'
                assert_msg_critical(
                    len(self.energies) == len(molecules), errmsg)

        grids = []
        esp = []
        scf_energies = []

        for ind, (mol, bas) in enumerate(zip(molecules, basis_sets)):
            info_text = f'Processing conformer {ind + 1}...'
            self.ostream.print_info(info_text)
            self.ostream.print_blank()
            self.ostream.print_block(mol.get_string())
            self.ostream.print_block(mol.more_info())
            self.ostream.print_blank()
            self.ostream.print_block(bas.get_string('Atomic Basis'))
            self.ostream.flush()

            nalpha = mol.number_of_alpha_electrons()
            nbeta = mol.number_of_beta_electrons()

            # run SCF
            if self.filename is not None:
                output_dir = Path(self.filename + '_files')
            else:
                output_dir = None

            if self.rank == mpi_master():
                if output_dir is not None:
                    output_dir.mkdir(parents=True, exist_ok=True)
            self.comm.barrier()

            if output_dir is not None:
                ostream_fname = str(output_dir / f'conformer_{ind + 1}')
                ostream = OutputStream(ostream_fname + '.out')
            else:
                ostream = self.ostream

            # select SCF driver
            if nalpha == nbeta:
                scf_drv = ScfRestrictedDriver(self.comm, ostream)
            else:
                scf_drv = ScfUnrestrictedDriver(self.comm, ostream)
            if self.filename is not None:
                scf_drv.filename = self.filename
            scf_drv.restart = self.restart
            scf_drv.conv_thresh = self.conv_thresh
            scf_drv.xcfun = self.xcfun
            scf_drv.grid_level = self.grid_level
            scf_drv.solvation_model = self.solvation_model
            scf_results = scf_drv.compute(mol, bas)

            if output_dir is not None:
                ostream.close()

            grid_m = self.get_grid_points(mol)
            esp_m = self.get_electrostatic_potential(grid_m, mol, bas,
                                                     scf_results)
            scf_energy_m = scf_drv.get_scf_energy()

            if self.rank == mpi_master():
                grids.append(grid_m)
                esp.append(esp_m)
                scf_energies.append(scf_energy_m)

        if self.rank == mpi_master():
            if self.weights is not None:
                weights = self.weights
            elif self.energies is not None:
                weights = self.get_boltzmann_weights(self.energies,
                                                     self.temperature)
            else:
                weights = self.get_boltzmann_weights(scf_energies,
                                                     self.temperature)
            q = self.compute_esp_charges(molecules, grids, esp, weights)
        else:
            q = None

        return q

    def compute_esp_charges(self, molecules, grids, esp, weights):
        """
        Computes Merz-Kollman or CHELPG ESP charges.

        :param molecules:
            The list of molecules.
        :param grids:
            The list of grid points.
        :param esp:
            The QM ESP on grid points.
        :param weights:
            The weight factors of different conformers.

        :return:
            The ESP charges.
        """

        n_conf = len(molecules)
        n_points = sum([esp[m].size for m in range(n_conf)])
        self.print_header(n_conf, n_points)
        self.print_esp_header()

        if self.equal_charges is not None:
            constr = [0] * molecules[0].number_of_atoms()
            for points in self.equal_charges:
                for i in range(len(points) - 1):
                    constr[points[i + 1] - 1] = points[i]
        else:
            constr = None

        # generate and solve equation system (ESP fit)
        a, b = self.generate_equation_system(molecules, grids, esp, weights,
                                             None, constr)
        q = safe_solve(a, b)

        n_atoms = molecules[0].number_of_atoms()

        # print results
        for i in range(n_atoms):
            cur_str = '{:3d}     {:2s}     {:11.6f}    '.format(
                i + 1, molecules[0].get_labels()[i], q[i])
            self.ostream.print_header(cur_str)

        q_tot = np.sum(q[:n_atoms])
        q_tot = round(q_tot * 1e+6) / 1e+6
        self.ostream.print_header(31 * '-')
        cur_str = 'Total Charge  : {:9.6f}   '.format(q_tot)
        self.ostream.print_header(cur_str)
        self.ostream.print_blank()

        self.print_fit_quality(self.get_rrms(molecules, grids, esp, weights, q))
        self.ostream.print_blank()
        self.print_references(len(molecules))
        self.ostream.flush()

        return q[:n_atoms]

    def compute_esp_charges_on_fitting_points(self, molecules, grids, esp):
        """
        Computes ESP charges on fitting points.

        :param molecules:
            The list of molecules.
        :param grids:
            The list of grid points.
        :param esp:
            The QM ESP on grid points.

        :return:
            The ESP charges on fitting points.
        """

        # get coordinates and constraints from input
        coords = []
        for string in self.fitting_points:
            coords.append([float(j) for j in string.split()])
        coords = np.array(coords)

        n_points = coords.shape[0]

        constr = [0] * n_points
        if self.equal_charges is not None:
            for points in self.equal_charges:
                for i in range(len(points) - 1):
                    constr[points[i + 1] - 1] = points[i]

        # number of constraints
        n_c = 0
        for i in range(n_points):
            if constr[i] > 0:
                n_c += 1

        self.print_header(1, esp[0].size, coords)
        self.print_esp_fitting_points_header(len(molecules))

        # calculate matrix a and vector b
        a = np.zeros((n_points + 1 + n_c, n_points + 1 + n_c))
        b = np.zeros(n_points + 1 + n_c)

        # total charge constraint
        a[:n_points, n_points] = 1
        a[n_points, :n_points] = 1
        b[n_points] = molecules[0].get_charge()

        # eq.10 & 11,  J. Am. Chem. Soc. 1992, 114, 9075-9079
        for i in range(n_points):
            for p in range(esp[0].size):
                r_pi = np.linalg.norm(grids[0][p] - coords[i])
                b[i] += esp[0][p] / r_pi
                for j in range(n_points):
                    a[i, j] += 1.0 / (r_pi *
                                      np.linalg.norm(grids[0][p] - coords[j]))

        # equivalent charge constraints
        j = 0
        for i in range(n_points):
            if constr[i] > 0:
                a[n_points + 1 + j, constr[i] - 1] = 1
                a[n_points + 1 + j, i] = -1
                a[constr[i] - 1, n_points + 1 + j] = 1
                a[i, n_points + 1 + j] = -1
                j += 1

        q = safe_solve(a, b)

        # print results
        for i in range(n_points):
            str_constr = ' '
            if constr[i] > 0:
                str_constr = str(constr[i])
            cur_str = '{:3d}     {:>3s}     {:11.6f}    '.format(
                i + 1, str_constr, q[i])
            self.ostream.print_header(cur_str)

        q_tot = np.sum(q[:n_points])
        q_tot = round(q_tot * 1e+6) / 1e+6
        self.ostream.print_header(31 * '-')
        cur_str = 'Total Charge  : {:9.6f}   '.format(q_tot)
        self.ostream.print_header(cur_str)
        self.ostream.print_blank()

        # calculate RRMS
        chi_square = 0.0
        norm = 0.0
        for i in range(esp[0].size):
            # partial charge ESP
            pot = 0.0
            for j in range(n_points):
                pot += q[j] / np.linalg.norm(grids[0][i] - np.array(coords[j]))
            chi_square += (esp[0][i] - pot)**2
            norm += esp[0][i]**2

        self.print_fit_quality(np.sqrt(chi_square / norm))
        self.ostream.print_blank()
        self.ostream.flush()

        return q[:n_points]

    def generate_equation_system(self,
                                 molecules,
                                 grids,
                                 esp,
                                 weights,
                                 q0=None,
                                 constr=None):
        """
        Generates the equation system for charge fitting to QM ESP.

        :param molecules:
            The list of molecules.
        :param grids:
            The list of grid points.
        :param esp:
            The QM ESP on grid points.
        :param weights:
            The weight factors of different conformers.
        :param q0:
            The initial charges.
        :param constr:
            The constraints of the charges.

        :return:
            The matrix a and vector b.
        """

        n_atoms = molecules[0].number_of_atoms()

        if q0 is None:
            q0 = np.zeros(n_atoms)
        if constr is None:
            constr = [0] * n_atoms

        # number of Lagrangian constraints
        n_c = 0
        for i in range(n_atoms):
            if constr[i] > 0:
                n_c += 1

        # calculate unrestrained matrix a and vector b
        a = np.zeros((n_atoms + 1 + n_c, n_atoms + 1 + n_c))
        b = np.zeros(n_atoms + 1 + n_c)

        # total charge constraint
        a[:n_atoms, n_atoms] = 1
        a[n_atoms, :n_atoms] = 1
        b[n_atoms] = molecules[0].get_charge()

        # eq.10 & 11,  J. Am. Chem. Soc. 1992, 114, 9075-9079
        for m in range(len(molecules)):
            a_m = np.zeros(a.shape)
            b_m = np.zeros(b.shape)

            coords = molecules[m].get_coordinates_in_bohr()

            for i in range(n_atoms):
                if constr[i] != -1:

                    inv_r_pi = 1.0 / np.sqrt(
                        np.sum((grids[m] - coords[i])**2, axis=1))

                    b_m[i] += np.sum(esp[m] * inv_r_pi)

                    for j in range(n_atoms):

                        r_pj = np.sqrt(np.sum((grids[m] - coords[j])**2,
                                              axis=1))

                        a_m[i, j] += np.sum(inv_r_pi / r_pj)

            a += weights[m] * a_m
            b += weights[m] * b_m

        # equivalent charge constraints and frozen/constant charges
        j = 0
        for i in range(n_atoms):
            if constr[i] == -1:
                a[i, i] = 1
                a[i, n_atoms] = 0
                b[i] = q0[i]
            elif constr[i] > 0:
                a[n_atoms + 1 + j, constr[i] - 1] = 1
                a[n_atoms + 1 + j, i] = -1
                a[constr[i] - 1, n_atoms + 1 + j] = 1
                a[i, n_atoms + 1 + j] = -1
                j += 1

        return a, b

    def get_grid_points(self, molecule):
        """
        Gets grid points for (R)ESP fitting.

        :param molecule:
            The molecule.

        :return:
            The coordinates of grid points as an array.
        """

        if self.grid_type == 'mk':
            return self.get_grid_points_mk(molecule)

        elif self.grid_type == 'chelpg':
            return self.get_grid_points_chelpg(molecule)

        else:
            assert_msg_critical(
                False, f'{type(self).__name__}.compute: Invalid grid_type ' +
                self.grid_type)

    def get_grid_points_mk(self, molecule):
        """
        Gets MK grid points.

        :param molecule:
            The molecule.

        :return:
            The coordinates of grid points as an array.
        """

        grid = []
        coords = molecule.get_coordinates_in_bohr()

        mol_mk_radii = molecule.mk_radii_to_numpy()

        if self.custom_mk_radii is not None:
            assert_msg_critical(
                len(self.custom_mk_radii) % 2 == 0,
                f'{type(self).__name__}: expecting even number of entries for '
                + 'user-defined MK radii')

            keys = self.custom_mk_radii[0::2]
            vals = self.custom_mk_radii[1::2]

            self.ostream.print_blank()

            for key, val in zip(keys, vals):
                val_au = float(val) / bohr_in_angstrom()
                try:
                    idx = int(key) - 1
                    assert_msg_critical(
                        0 <= idx and idx < molecule.number_of_atoms(),
                        f'{type(self).__name__}: invalid atom index for ' +
                        'user-defined MK radii')
                    mol_mk_radii[idx] = val_au
                    self.ostream.print_info(
                        f'Applying user-defined MK radius {val} for atom {key}')
                except ValueError:
                    elem_found = False
                    for idx, label in enumerate(molecule.get_labels()):
                        if label.upper() == key.upper():
                            mol_mk_radii[idx] = val_au
                            elem_found = True
                    if elem_found:
                        self.ostream.print_info(
                            f'Applying user-defined MK radius {val} for atom {key}'
                        )

            self.ostream.print_blank()

        for layer in range(self.number_layers):

            # MK radii with layer-dependent scaling factor
            scaling_factor = 1.4 + layer * 0.4 / np.sqrt(self.number_layers)
            r = scaling_factor * mol_mk_radii

            for atom in range(molecule.number_of_atoms()):
                # number of points fitting on the equator
                n_eq = int(2.0 * np.pi * r[atom] * np.sqrt(self.density) *
                           bohr_in_angstrom())

                # number of latitudes
                n_i = n_eq // 2

                for i in range(n_i + 1):
                    phi = np.pi * i / n_i

                    # number of points on the latitude
                    if i == 0 or i == n_i:
                        n_j = 1
                    else:
                        n_j = int(n_eq * np.sin(phi) + 1e-10)

                    # generate new point and check for overlap
                    for j in range(n_j):
                        theta = 2 * np.pi * j / n_j
                        point = coords[atom] + r[atom] * np.array([
                            np.sin(phi) * np.cos(theta),
                            np.sin(phi) * np.sin(theta),
                            np.cos(phi)
                        ])

                        save = True
                        for atom2 in range(molecule.number_of_atoms()):
                            if atom != atom2 and np.linalg.norm(
                                    point - coords[atom2]) < r[atom2]:
                                save = False
                                break
                        if save:
                            grid.append(point)

        return np.array(grid)

    def get_grid_points_chelpg(self, molecule):
        """
        Gets CHELPG grid points.

        :param molecule:
            The molecule.

        :return:
            The coordinates of grid points as an array.
        """

        chelpg_spacing_au = self.chelpg_spacing / bohr_in_angstrom()
        chelpg_margin_au = self.chelpg_margin / bohr_in_angstrom()

        coords = molecule.get_coordinates_in_bohr()

        x_max = np.max(coords[:, 0] + chelpg_margin_au)
        x_min = np.min(coords[:, 0] - chelpg_margin_au)

        y_max = np.max(coords[:, 1] + chelpg_margin_au)
        y_min = np.min(coords[:, 1] - chelpg_margin_au)

        z_max = np.max(coords[:, 2] + chelpg_margin_au)
        z_min = np.min(coords[:, 2] - chelpg_margin_au)

        n_x_half_float = 0.5 * (x_max - x_min) / chelpg_spacing_au
        n_y_half_float = 0.5 * (y_max - y_min) / chelpg_spacing_au
        n_z_half_float = 0.5 * (z_max - z_min) / chelpg_spacing_au

        n_x_half = int(n_x_half_float)
        n_y_half = int(n_y_half_float)
        n_z_half = int(n_z_half_float)

        x_min = 0.5 * (x_min + x_max) - chelpg_spacing_au * n_x_half
        y_min = 0.5 * (y_min + y_max) - chelpg_spacing_au * n_y_half
        z_min = 0.5 * (z_min + z_max) - chelpg_spacing_au * n_z_half

        x_max = 0.5 * (x_min + x_max) + chelpg_spacing_au * n_x_half
        y_max = 0.5 * (y_min + y_max) + chelpg_spacing_au * n_y_half
        z_max = 0.5 * (z_min + z_max) + chelpg_spacing_au * n_z_half

        # try to fill the box with fewer grid points

        if n_x_half_float - n_x_half <= 0.5:
            x_min += 0.5 * chelpg_spacing_au
            x_max -= 0.5 * chelpg_spacing_au

        if n_y_half_float - n_y_half <= 0.5:
            y_min += 0.5 * chelpg_spacing_au
            y_max -= 0.5 * chelpg_spacing_au

        if n_z_half_float - n_z_half <= 0.5:
            z_min += 0.5 * chelpg_spacing_au
            z_max -= 0.5 * chelpg_spacing_au

        natoms = molecule.number_of_atoms()
        atom_radii = molecule.chelpg_radii_to_numpy()

        margin_squared = chelpg_margin_au**2
        atom_radii_squared = atom_radii * atom_radii

        grid_points = []

        x_grid = np.arange(x_min, x_max + 0.01 * chelpg_spacing_au,
                           chelpg_spacing_au)
        y_grid = np.arange(y_min, y_max + 0.01 * chelpg_spacing_au,
                           chelpg_spacing_au)
        z_grid = np.arange(z_min, z_max + 0.01 * chelpg_spacing_au,
                           chelpg_spacing_au)

        for x in x_grid:
            for y in y_grid:
                for z in z_grid:

                    # criteria for including the grid point:
                    # 1. grid point is outside CHELPG radii
                    # 2. grid point is inside maximum radius (margin)

                    within_vdw_radii = False
                    min_r2 = None

                    for a in range(natoms):
                        ax, ay, az = coords[a]
                        r2 = (x - ax)**2 + (y - ay)**2 + (z - az)**2
                        if r2 < atom_radii_squared[a]:
                            within_vdw_radii = True
                            break
                        if min_r2 is None or min_r2 > r2:
                            min_r2 = r2

                    if not (within_vdw_radii or min_r2 > margin_squared):
                        grid_points.append([x, y, z])

        return np.array(grid_points)

    def get_electrostatic_potential(self, grid, molecule, basis, scf_results):
        """
        Gets the QM ESP on the grid points.

        :param grid:
            The grid points.
        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_results:
            The dictionary containing SCF results.

        :return:
            The ESP at each grid point listed in an array.
        """

        if self.rank == mpi_master():
            D = scf_results['D_alpha'] + scf_results['D_beta']
        else:
            D = None
        D = self.comm.bcast(D, root=mpi_master())

        ave, res = divmod(grid.shape[0], self.nodes)
        counts = [ave + 1 if p < res else ave for p in range(self.nodes)]

        start = sum(counts[:self.rank])
        end = sum(counts[:self.rank + 1])

        esp = np.zeros(end - start)

        # nuclear contribution

        coords = molecule.get_coordinates_in_bohr()
        elem_ids = molecule.get_element_ids()

        for a in range(molecule.number_of_atoms()):
            esp += elem_ids[a] / np.sqrt(
                np.sum((grid[start:end, :] - coords[a])**2, axis=1))

        # electrostatic contribution

        esp += compute_nuclear_potential_values(molecule, basis,
                                                grid[start:end, :], D)

        gathered_esp = self.comm.gather(esp, root=mpi_master())

        if self.rank == mpi_master():
            return np.hstack(gathered_esp)
        else:
            return None

    def get_rrms(self, molecules, grids, esp, weights, q):
        """
        Gets the relative root-mean-square (RRMS) error.

        :param molecules:
            The list of molecules.
        :param grids:
            The list of grid points.
        :param esp:
            The QM ESP on grid points.
        :param weights:
            The weight factors of different conformers.
        :param q:
            The atom-centered charges.

        :return:
            The RRMS.
        """

        rrms = 0.0

        for m in range(len(molecules)):
            coords = molecules[m].get_coordinates_in_bohr()
            chi_square = 0.0
            norm = 0.0
            for i in range(esp[m].size):
                # partial charge ESP
                pot = 0.0
                for j in range(molecules[m].number_of_atoms()):
                    pot += q[j] / np.linalg.norm(grids[m][i] - coords[j])
                chi_square += (esp[m][i] - pot)**2
                norm += esp[m][i]**2
            rrms += weights[m] * np.sqrt(chi_square / norm)
        return rrms

    def get_boltzmann_weights(self, energies, temperature):
        """
        Gets Boltzmann weight factors for conformers.

        :param energies:
            The energies.
        :param temperature:
            The temperature.
        """

        k_b = boltzmann_in_evperkelvin() / hartree_in_ev()

        weights = []
        for i in range(len(energies)):
            weights.append(
                np.exp(
                    (-1) * (energies[i] - min(energies)) / (k_b * temperature)))

        # normalize weights
        sum_of_weights = sum(weights)
        weights = [w / sum_of_weights for w in weights]

        return weights

    def print_header(self, n_conf, n_points, fitting_coords=None):
        """
        Prints header for the RESP charges driver.

        :param n_conf:
            The number of conformers.
        :param n_points:
            The number of grid points.
        :param fitting_coords:
            The coordinates of the fitting points.
        """

        self.ostream.print_blank()

        if self.grid_type == 'mk':
            title = 'RESP Charges Driver Setup'
        elif self.grid_type == 'chelpg':
            title = 'ESP Charges Driver Setup'
        self.ostream.print_header(title)
        self.ostream.print_header('=' * (len(title) + 2))
        self.ostream.print_blank()

        str_width = 40
        cur_str = 'Number of Conformers         :  ' + str(n_conf)
        self.ostream.print_header(cur_str.ljust(str_width))
        if n_conf > 1:
            if self.weights is not None:
                cur_str = 'Weights of Conformers        :  {:4f}'.format(
                    self.weights[0])
                self.ostream.print_header(cur_str.ljust(str_width))
                for i in range(1, n_conf):
                    cur_str = '                                {:4f}'.format(
                        self.weights[i])
                    self.ostream.print_header(cur_str.ljust(str_width))

        if self.grid_type == 'mk':
            cur_str = 'Number of Layers             :  ' + str(
                self.number_layers)
            self.ostream.print_header(cur_str.ljust(str_width))
            cur_str = 'Points per Square Angstrom   :  ' + str(self.density)
            self.ostream.print_header(cur_str.ljust(str_width))

        elif self.grid_type == 'chelpg':
            cur_str = 'Grid Spacing in Angstrom     :  ' + str(
                self.chelpg_spacing)
            self.ostream.print_header(cur_str.ljust(str_width))
            cur_str = 'Grid Margin in Angstrom      :  ' + str(
                self.chelpg_margin)
            self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = 'Total Number of Grid Points  :  ' + str(n_points)
        self.ostream.print_header(cur_str.ljust(str_width))
        if fitting_coords is not None:
            cur_str = 'Fitting Points  :'
            self.ostream.print_header(cur_str.ljust(str_width))
            for coords in fitting_coords:
                cur_str = '                 {:7.3f}{:7.3f}{:7.3f}'.format(
                    coords[0], coords[1], coords[2])
                self.ostream.print_header(cur_str.ljust(str_width))

        self.ostream.flush()

    def print_esp_header(self):
        """
        Prints header for ESP charges calculation.

        """

        self.ostream.print_blank()

        if self.grid_type == 'mk':
            self.ostream.print_header('Merz-Kollman ESP Charges')
        elif self.grid_type == 'chelpg':
            self.ostream.print_header('CHELPG ESP Charges')

        self.ostream.print_header(26 * '-')
        self.ostream.print_blank()
        cur_str = '{}   {}      {}'.format('No.', 'Atom', 'Charge (a.u.)')
        self.ostream.print_header(cur_str)
        self.ostream.print_header(31 * '-')
        self.ostream.flush()

    def print_esp_fitting_points_header(self, n_conf):
        """
        Prints header for ESP charges calculation.

        :param n_conf:
            The number of conformers.
        """

        self.ostream.print_blank()
        if n_conf > 1:
            cur_str = 'Multiple conformers are given, but only one '
            cur_str += 'can be processed!'
            self.ostream.print_warning(cur_str)

        self.ostream.print_blank()
        self.ostream.print_header('ESP Charges on Fitting Points')
        self.ostream.print_header(31 * '-')
        self.ostream.print_blank()
        cur_str = '{} {}  {}'.format('No.', 'Constraints', 'Charge (a.u.)')
        self.ostream.print_header(cur_str)
        self.ostream.print_header(31 * '-')
        self.ostream.flush()

    def print_fit_quality(self, rrms):
        """
        Prints fit quality.

        :param rrms:
            The relative root-mean-square error.
        """

        self.ostream.print_header('Fit Quality')
        self.ostream.print_header('-' * 13)
        cur_str = 'Relative Root-Mean-Square Error  : {:9.6f}'.format(rrms)
        self.ostream.print_header(cur_str)
        self.ostream.flush()

    def print_references(self, n_conf):
        """
        Prints references.

        :param n_conf:
            The number of conformers.
        """

        str_width = 44

        title = 'Reference: '
        self.ostream.print_header(title.ljust(str_width))

        if self.grid_type == 'mk':
            title = 'J. Comput. Chem. 1984, 5, 129-145.'
            self.ostream.print_header(title.ljust(str_width))
        elif self.grid_type == 'chelpg':
            title = 'J. Comput. Chem. 1990, 11, 361-373.'
            self.ostream.print_header(title.ljust(str_width))

        if n_conf > 1:
            title = 'J. Am. Chem. Soc. 1992, 114, 9075-9079.'
            self.ostream.print_header(title.ljust(str_width))

        self.ostream.print_blank()
        self.ostream.flush()

    def read_multiple_molecules(self):
        """
        Reads multiple molecules from xyz file.

        :return:
            A list of molecules.
        """

        assert_msg_critical(
            self.xyz_file is not None,
            f'{type(self).__name__}.compute: xyz_file not specified')

        if self.rank == mpi_master():
            with open(self.xyz_file, 'r') as f_xyz:
                xyz_lines = f_xyz.readlines()
            n_atoms = int(xyz_lines[0].split()[0])
            n_geoms = len(xyz_lines) // (n_atoms + 2)
        else:
            n_geoms = None
        n_geoms = self.comm.bcast(n_geoms, root=mpi_master())

        molecules = []

        for i_geom in range(n_geoms):
            if self.rank == mpi_master():
                i_start = i_geom * (n_atoms + 2)
                i_end = i_start + (n_atoms + 2)

                assert_msg_critical(
                    int(xyz_lines[i_start].split()[0]) == n_atoms,
                    f'{type(self).__name__}.compute: inconsistent number of atoms'
                )

                mol_str = ''.join(xyz_lines[i_start + 2:i_end])

                mol = Molecule.read_molecule_string(mol_str, units='angstrom')
                mol.set_charge(self.net_charge)
                mol.set_multiplicity(self.multiplicity)
            else:
                mol = Molecule()

            mol = self.comm.bcast(mol, root=mpi_master())
            molecules.append(mol)

        info_text = f'Found {len(molecules):d} conformers.'
        self.ostream.print_info(info_text)
        self.ostream.print_blank()
        self.ostream.flush()

        return molecules
