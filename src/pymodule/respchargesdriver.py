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
from .inputparser import parse_input, print_keywords, get_random_string_parallel
from .errorhandler import assert_msg_critical, safe_solve


class RespChargesDriver:
    """
    Implements ESP and RESP charges.

    # vlxtag: RHF, RESP_charges
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
        - restrain_hydrogen: If hydrogens should be restrained or not.
        - weak_restraint: The strength of the restraint in first stage of RESP fit.
        - strong_restraint: The strength of the restraint in second stage of RESP fit.
        - weights: The weight factors of different conformers.
        - energies: The energies of different conformers for Boltzmann weight factors.
        - temperature: The temperature for Boltzmann weight factors.
        - net_charge: The charge of the molecule.
        - multiplicity: The multiplicity of the molecule.
        - method_dict: The dictionary of method settings.
        - max_iter: The maximum number of iterations of the RESP fit.
        - threshold: The convergence threshold of the RESP fit.
        - filename: The filename for the calculation.
        - xyz_file: The xyz file containing the conformers.
        - fitting_points: The fitting point on which charges are calculated.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes RESP charges driver.
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
        self.xcfun = None
        self.restart = True

        # conformers
        self.xyz_file = None
        self.net_charge = 0.0
        self.multiplicity = 1
        self.method_dict = None
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

        # resp fitting
        self.equal_charges = None
        self.restrain_hydrogen = False
        self.weak_restraint = 0.0005
        self.strong_restraint = 0.001
        self.max_iter = 50
        self.threshold = 1.0e-6

        # input keywords
        self.input_keywords = {
            'resp_charges': {
                'xcfun': ('str_upper', 'exchange-correlation functional'),
                'restart': ('bool', 'restart from checkpoint file'),
                'grid_type': ('str_lower', 'type of grid (mk or chelpg)'),
                'number_layers':
                    ('int', 'number of layers of scaled vdW surfaces'),
                'density': ('float', 'density of points in each layer'),
                'chelpg_spacing': ('float', 'step size for CHELPG grid'),
                'chelpg_margin': ('float', 'maximum radius for CHELPG grid'),
                'restrain_hydrogen': ('bool', 'restrain hydrogen atoms'),
                'weak_restraint':
                    ('float', 'strength of restraint in 1st RESP stage'),
                'strong_restraint':
                    ('float', 'strength of restraint in 2nd RESP stage'),
                'max_iter': ('int', 'maximum iterations in RESP fit'),
                'threshold': ('float', 'convergence threshold of RESP fit'),
                'equal_charges': ('str', 'constraints for equal charges'),
                'custom_mk_radii':
                    ('seq_fixed_str', 'custom MK radii for RESP charges'),
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
        }

    def print_keywords(self):
        """
        Prints input keywords in SCF driver.
        """

        print_keywords(self.input_keywords, self.ostream)

    def update_settings(self, resp_dict, method_dict=None):
        """
        Updates settings in RESP charges driver.

        :param resp_dict:
            The dictionary of RESP charges input.
        :param method_dict:
            The dicitonary of method settings.
        """

        resp_keywords = {
            key: val[0]
            for key, val in self.input_keywords['resp_charges'].items()
        }

        parse_input(self, resp_keywords, resp_dict)

        if 'filename' in resp_dict:
            self.filename = resp_dict['filename']

        if method_dict is not None:
            self.method_dict = dict(method_dict)

    def compute(self, molecule, basis, *args):
        """
        Computes RESP or ESP charges.

        :param molecule:
            The molecule.
        :param basis:
            The molecular basis set.
        :param scf_results:
            The SCF results.
        :param chg_type:
            The type of charges ('resp' or 'esp').

        :return:
            The charges.
        """

        # for backward compatibility, where three arguments were given: (mol,
        # bas, chg_type)
        if len(args) == 1 and isinstance(args[0], str):
            chg_type = args[0]
            return self.compute_weighted_charges(molecule, basis, chg_type)

        # we expect scf_results and chg_type in addition to molecule and basis
        assert_msg_critical(
            len(args) == 2,
            'RespChargesDriver.compute: Expecting scf_results and chg_type')
        scf_results, chg_type = args

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
                'RespChargesDriver.compute: Invalid equal_charges')

        use_631gs_basis = (basis.get_label() in ['6-31G*', '6-31G_D_'])

        # sanity check for RESP
        if chg_type.lower() == 'resp':
            assert_msg_critical(
                self.grid_type == 'mk',
                'RespChargesDriver.compute: For RESP charges, ' +
                'grid_type must be \'mk\'')
            if not use_631gs_basis:
                cur_str = 'Recommended basis set 6-31G* is not used!'
                self.ostream.print_warning(cur_str)
                self.ostream.print_blank()

        if self.rank == mpi_master():
            nalpha = molecule.number_of_alpha_electrons()
            nbeta = molecule.number_of_beta_electrons()
            if nalpha != nbeta:
                assert_msg_critical(
                    scf_results['scf_type'] != 'restricted',
                    'RespChargesDriver.compute: Molecule is open-shell while ' +
                    'SCF is closed-shell')

        grid_m = self.get_grid_points(molecule)
        esp_m = self.get_electrostatic_potential(grid_m, molecule, basis,
                                                 scf_results)

        if self.rank == mpi_master():
            scf_energy_m = scf_results['scf_energy']
            if self.energies is None:
                self.energies = [scf_energy_m]
            if self.weights is None:
                self.weights = self.get_boltzmann_weights(
                    self.energies, self.temperature)

        q = None

        if self.rank == mpi_master():
            if self.fitting_points is not None:
                q = self.compute_esp_charges_on_fitting_points([molecule],
                                                               [grid_m],
                                                               [esp_m])
            elif chg_type.lower() == 'resp':
                q = self.compute_resp_charges([molecule], [grid_m], [esp_m],
                                              self.weights)
            elif chg_type.lower() == 'esp':
                q = self.compute_esp_charges([molecule], [grid_m], [esp_m],
                                             self.weights)

        return q

    def compute_weighted_charges(self, molecule, basis, flag):
        """
        Computes RESP or ESP charges.

        :param molecule:
            The molecule.
        :param basis:
            The molecular basis set.
        :param flag:
            The flag for the task ("resp" or "esp").

        :return:
            The charges.
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
                'RespChargesDriver.compute: Invalid equal_charges')

        molecules = []
        basis_sets = []

        use_xyz_file = (molecule.number_of_atoms() == 0)

        if not use_xyz_file:
            molecules.append(molecule)
            basis_sets.append(basis)

        else:
            errmsg = 'RespChargesDriver: The \'xyz_file\' keyword is not '
            errmsg += 'specified.'
            assert_msg_critical(self.xyz_file is not None, errmsg)

            xyz_lines = None
            if self.rank == mpi_master():
                with open(self.xyz_file, 'r') as f_xyz:
                    xyz_lines = f_xyz.readlines()
            xyz_lines = self.comm.bcast(xyz_lines, root=mpi_master())

            n_atoms = int(xyz_lines[0].split()[0])
            n_geoms = len(xyz_lines) // (n_atoms + 2)

            for i_geom in range(n_geoms):
                i_start = i_geom * (n_atoms + 2)
                i_end = i_start + (n_atoms + 2)

                assert_msg_critical(
                    int(xyz_lines[i_start].split()[0]) == n_atoms,
                    'RespChargesDriver.compute: inconsistent number of atoms')

                mol_str = ''.join(xyz_lines[i_start + 2:i_end])

                if self.rank == mpi_master():
                    # create molecule
                    mol = Molecule.read_molecule_string(mol_str,
                                                        units='angstrom')
                    mol.set_charge(self.net_charge)
                    mol.set_multiplicity(self.multiplicity)
                    # create basis set
                    basis_path = '.'
                    if 'basis_path' in self.method_dict:
                        basis_path = self.method_dict['basis_path']
                    basis_name = self.method_dict['basis'].upper()
                    bas = MolecularBasis.read(mol,
                                              basis_name,
                                              basis_path,
                                              ostream=None)
                else:
                    mol = Molecule()
                    bas = MolecularBasis()

                # broadcast molecule and basis set
                mol = self.comm.bcast(mol, root=mpi_master())
                bas = self.comm.bcast(bas, root=mpi_master())

                molecules.append(mol)
                basis_sets.append(bas)

            info_text = f'Found {len(molecules):d} conformers.'
            self.ostream.print_info(info_text)
            self.ostream.print_blank()
            self.ostream.flush()

        if self.rank == mpi_master():
            if self.weights is not None and self.energies is not None:
                # avoid conflict between weights and energies
                errmsg = 'RespChargesDriver: Please do not specify weights and '
                errmsg += 'energies at the same time.'
                assert_msg_critical(False, errmsg)

            elif self.weights is not None:
                errmsg = 'RespChargesDriver: Number of weights does not match '
                errmsg += 'number of conformers.'
                assert_msg_critical(len(self.weights) == len(molecules), errmsg)
                # normalize weights
                sum_of_weights = sum(self.weights)
                self.weights = [w / sum_of_weights for w in self.weights]

            elif self.energies is not None:
                errmsg = 'RespChargesDriver: Number of energies does not match '
                errmsg += 'number of conformers.'
                assert_msg_critical(
                    len(self.energies) == len(molecules), errmsg)

        use_631gs_basis = all(
            [bas.get_label() in ['6-31G*', '6-31G_D_'] for bas in basis_sets])

        # sanity check for RESP
        if flag.lower() == 'resp':
            assert_msg_critical(
                self.grid_type == 'mk',
                'RespChargesDriver.compute: For RESP charges, ' +
                'grid_type must be \'mk\'')
            if not use_631gs_basis:
                cur_str = 'Recommended basis set 6-31G* is not used!'
                self.ostream.print_warning(cur_str)
                self.ostream.print_blank()

        grids = []
        esp = []
        scf_energies = []

        for ind, (mol, bas) in enumerate(zip(molecules, basis_sets)):
            if use_xyz_file:
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
            if not use_xyz_file:
                # select SCF driver
                if nalpha == nbeta:
                    scf_drv = ScfRestrictedDriver(self.comm, self.ostream)
                else:
                    scf_drv = ScfUnrestrictedDriver(self.comm, self.ostream)
                if self.filename is not None:
                    scf_drv.filename = self.filename
                scf_drv.restart = self.restart
                if (self.method_dict is None or
                        'xcfun' not in self.method_dict):
                    scf_drv.xcfun = self.xcfun
                scf_results = scf_drv.compute(mol, bas)

            else:
                if self.filename is not None:
                    output_dir = Path(self.filename + '_files')
                else:
                    name_string = get_random_string_parallel(self.comm)
                    output_dir = Path('vlx_' + name_string + '_files')

                if self.rank == mpi_master():
                    output_dir.mkdir(parents=True, exist_ok=True)
                self.comm.barrier()

                ostream_fname = str(output_dir / f'conformer_{ind + 1}')
                ostream = OutputStream(ostream_fname + '.out')

                # select SCF driver
                if nalpha == nbeta:
                    scf_drv = ScfRestrictedDriver(self.comm, ostream)
                else:
                    scf_drv = ScfUnrestrictedDriver(self.comm, ostream)
                if self.filename is not None:
                    scf_drv.filename = self.filename
                scf_drv.restart = self.restart
                if (self.method_dict is None or
                        'xcfun' not in self.method_dict):
                    scf_drv.xcfun = self.xcfun
                scf_results = scf_drv.compute(mol, bas)
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
            if self.energies is None:
                self.energies = scf_energies
            if self.weights is None:
                self.weights = self.get_boltzmann_weights(
                    self.energies, self.temperature)

        q = None

        if self.rank == mpi_master():
            if self.fitting_points is not None:
                q = self.compute_esp_charges_on_fitting_points(
                    molecules, grids, esp)
            elif flag.lower() == 'resp':
                q = self.compute_resp_charges(molecules, grids, esp,
                                              self.weights)
            elif flag.lower() == 'esp':
                q = self.compute_esp_charges(molecules, grids, esp,
                                             self.weights)

        return q

    def compute_resp_charges(self, molecules, grids, esp, weights):
        """
        Computes RESP charges.

        :param molecules:
            The list of molecules.
        :param grids:
            The list of grid points.
        :param esp:
            The QM ESP on grid points.
        :param weights:
            The weight factors of different conformers.

        :return:
            The RESP charges.
        """

        n_conf = len(molecules)
        n_points = sum([esp[m].size for m in range(n_conf)])
        self.print_header(n_conf, n_points)

        constr_1, constr_2 = self.generate_constraints(molecules[0])

        # first stage of resp fit
        self.print_resp_stage_header('first')
        q0 = np.zeros(molecules[0].number_of_atoms())
        q = self.optimize_charges(molecules, grids, esp, weights, q0, constr_1,
                                  self.weak_restraint)

        self.print_resp_stage_results('first', molecules[0], q, constr_1)
        self.print_fit_quality(self.get_rrms(molecules, grids, esp, weights, q))
        self.ostream.print_blank()

        # second stage of RESP fit
        if constr_2 == [-1] * molecules[0].number_of_atoms():
            cur_str = '*** No refitting in second stage needed.'
            self.ostream.print_header(cur_str.ljust(40))
        else:
            self.print_resp_stage_header('second')
            q = self.optimize_charges(molecules, grids, esp, weights, q,
                                      constr_2, self.strong_restraint)

            self.print_resp_stage_results('second', molecules[0], q, constr_2)
            self.print_fit_quality(
                self.get_rrms(molecules, grids, esp, weights, q))

        self.ostream.print_blank()
        self.print_references('resp', len(molecules))
        self.ostream.flush()

        return q

    def compute_esp_charges(self, molecules, grids, esp, weights):
        """
        Computes Merz-Kollman ESP charges.

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
        self.print_references('esp', len(molecules))
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

    def optimize_charges(self, molecules, grids, esp, weights, q0, constr,
                         restraint_strength):
        """
        Iterative procedure to optimize charges for a single stage of the RESP fit.

        :param molecules:
            The list of molecules.
        :param grids:
            The list of grid points.
        :param esp:
            The QM ESP on grid points.
        :param weights:
            The weight factors of different conformers.
        :param q0:
            The intial charges.
        :param constr:
            The constraints of the charges.
        :param restraint_strength:
            The strength of the hyperbolic penalty function.

        :return:
            The optimized charges.
        """

        n_atoms = molecules[0].number_of_atoms()

        # initializes equation system
        a, b = self.generate_equation_system(molecules, grids, esp, weights, q0,
                                             constr)
        q_new = np.concatenate((q0, np.zeros(b.size - n_atoms)), axis=None)

        fitting_converged = False
        current_iteration = 0

        # iterative RESP fitting procedure
        for iteration in range(self.max_iter):
            q_old = q_new

            # add restraint to matrix a
            rstr = np.zeros(b.size)
            for i in range(n_atoms):
                if (not self.restrain_hydrogen and
                        molecules[0].get_labels()[i] == 'H'):
                    continue
                elif constr[i] == -1:
                    continue
                else:
                    # eq.10, J. Phys. Chem. 1993, 97, 10269-10280
                    rstr[i] += restraint_strength / np.sqrt(q_old[i]**2 +
                                                            0.1**2)
            a_tot = a + np.diag(rstr)

            q_new = safe_solve(a_tot, b)
            dq_norm = np.linalg.norm(q_new - q_old)

            current_iteration = iteration
            if dq_norm < self.threshold:
                fitting_converged = True
                break

        errmsg = 'RespChargesDriver: Charge fitting is not converged!'
        assert_msg_critical(fitting_converged, errmsg)

        cur_str = '*** Charge fitting converged in '
        cur_str += f'{current_iteration + 1} iterations.'
        self.ostream.print_header(cur_str.ljust(40))
        self.ostream.print_blank()

        return q_new[:n_atoms]

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

                    inv_r_pi = 1.0 / np.sqrt(np.sum(
                        (grids[m] - coords[i])**2, axis=1))

                    b_m[i] += np.sum(esp[m] * inv_r_pi)

                    for j in range(n_atoms):

                        r_pj = np.sqrt(np.sum(
                            (grids[m] - coords[j])**2, axis=1))

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

    def generate_constraints(self, molecule):
        """
        Determines Lagrangian constraints and frozen charges for both RESP fit stages.

        :param molecule:
            The molecule.

        return:
            Two lists (constr_1 and constr_2) for variations in the two stages of
            the RESP fitting procedure for which -1 tells a charge to be frozen,
            and x > 0 that the atom has the same charge like atom x.
        """

        n_atoms = molecule.number_of_atoms()
        coords = molecule.get_coordinates_in_bohr()
        labels = molecule.get_labels()
        covalent_radii = molecule.covalent_radii_to_numpy()

        connectivity = np.full((n_atoms, n_atoms), False, dtype='bool')

        # connect atoms with distances close to sum of covalent radii with
        # tolerance of 0.4 Angstrom
        tolerance = 0.4 / bohr_in_angstrom()
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                r_ij = np.linalg.norm(coords[i] - coords[j])
                if abs(r_ij - covalent_radii[i] -
                       covalent_radii[j]) < tolerance:
                    connectivity[i, j] = True
                    connectivity[j, i] = True

        # identify C(sp3)-H bonds from connectivity
        # c_nh: number of H bonded to atom if atom is C(sp3), else 0
        # h_cbond: atom number of C(sp3) to which H is bonded
        c_nh = [0] * n_atoms
        h_cbond = [0] * n_atoms

        for i in range(n_atoms):
            if labels[i] == 'C' and np.sum(connectivity[i]) == 4:
                for j in range(n_atoms):
                    if labels[j] == 'H' and connectivity[i, j]:
                        c_nh[i] += 1
                        h_cbond[j] = i + 1

        constr_1 = [0] * n_atoms
        constr_2 = [-1] * n_atoms

        # automatic constraints (refitting methyl and methylene in the 2nd stage)
        if self.equal_charges is None:
            for i in range(n_atoms):
                if c_nh[i] > 1:
                    constr_2[i] = 0
                    j = 0
                    for k in range(n_atoms):
                        if h_cbond[k] == (i + 1) and j == 0:
                            constr_2[k] = 0
                            j = k + 1
                        elif h_cbond[k] == (i + 1) and j > 0:
                            constr_2[k] = j
                            j = k + 1

        # use constraints from input
        else:
            # put all constraints to constr_1
            for atoms in self.equal_charges:
                for i in range(len(atoms) - 1):
                    constr_1[atoms[i + 1] - 1] = atoms[i]

            # move methyl and methylene constraints to constr_2
            for i in range(n_atoms):
                if c_nh[i] > 1:
                    constr_2[i] = constr_1[i]
                    constr_1[i] = 0
                    for j in range(n_atoms):
                        if h_cbond[j] == (i + 1):
                            constr_2[j] = constr_1[j]
                            constr_1[j] = 0

        return constr_1, constr_2

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
                False, 'RespChargesDriver.compute: Invalid grid_type ' +
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
                'RespChargesDriver: expecting even number of entries for ' +
                'user-defined MK radii')

            keys = self.custom_mk_radii[0::2]
            vals = self.custom_mk_radii[1::2]

            self.ostream.print_blank()

            for key, val in zip(keys, vals):
                val_au = float(val) / bohr_in_angstrom()
                try:
                    idx = int(key) - 1
                    assert_msg_critical(
                        0 <= idx and idx < molecule.number_of_atoms(),
                        'RespChargesDriver: invalid atom index for ' +
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

    def print_resp_stage_header(self, stage):
        """
        Prints header for a stage of the RESP fit.

        :param stage:
            The stage of the RESP fit.
        """

        str_width = 40
        self.ostream.print_blank()

        if stage.lower() == 'first':
            if not self.recommended_resp_parameters():
                cur_str = 'Parameters for RESP fitting differ from '
                cur_str += 'recommended choices!'
                self.ostream.print_warning(cur_str)
            self.ostream.print_blank()
            self.ostream.print_header('First Stage Fit')
            self.ostream.print_header(17 * '-')
            self.ostream.print_blank()
            cur_str = 'Restraint Strength           :  ' + str(
                self.weak_restraint)
        elif stage.lower() == 'second':
            self.ostream.print_header('Second Stage Fit')
            self.ostream.print_header(18 * '-')
            self.ostream.print_blank()
            cur_str = 'Restraint Strength           :  ' + str(
                self.strong_restraint)

        self.ostream.print_header(cur_str.ljust(str_width))
        str_restrain_hydrogen = 'Yes' if self.restrain_hydrogen else 'No'
        cur_str = 'Restrained Hydrogens         :  ' + str_restrain_hydrogen
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Max. Number of Iterations    :  ' + str(self.max_iter)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Convergence Threshold (a.u.) :  ' + str(self.threshold)
        self.ostream.print_header(cur_str.ljust(str_width))
        self.ostream.print_blank()
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

    def print_resp_stage_results(self, stage, molecule, q, constr):
        """
        Prints results for a stage of the RESP fit.

        :param stage:
            The stage of the RESP fit.
        :param molecule:
            The molecule.
        :param q:
            The final charges of the stage.
        :param constr:
            The variations and constraints for the stage.
        """

        if stage.lower() == 'first':
            str_variation = ''
            width = 44
        elif stage.lower() == 'second':
            str_variation = 'Frozen |'
            width = 52

        cur_str = '{} | {} | {} {} | {}'.format('No.', 'Atom', str_variation,
                                                'Constraints', 'Charges (a.u.)')
        self.ostream.print_header(cur_str)
        self.ostream.print_header(width * '-')

        for i in range(q.size):
            str_constraint = ''

            if constr[i] < 0:
                str_variation = 'Yes     '
            elif constr[i] == 0:
                str_variation = ' No     '
            else:
                str_variation = ' No     '
                str_constraint = str(constr[i])

            if stage.lower() == 'first':
                str_variation = ''

            cur_str = '{:3d}     {:2s}     {}    {:3.6s}     {:12.6f}   '.format(
                i + 1,
                molecule.get_labels()[i], str_variation, str_constraint, q[i])
            self.ostream.print_header(cur_str)

        self.ostream.print_header(width * '-')

        q_tot = np.sum(q)
        q_tot = round(q_tot * 1e+6) / 1e+6
        cur_str = 'Total Charge  : {:9.6f}   '.format(q_tot)
        self.ostream.print_header(cur_str)
        self.ostream.print_blank()
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

    def print_references(self, flag, n_conf):
        """
        Prints references.

        :param flag:
            The flag for the task ("resp" or "esp").
        :param n_conf:
            The number of conformers.
        """

        str_width = 44

        title = 'Reference: '
        self.ostream.print_header(title.ljust(str_width))

        if flag.lower() == 'resp':
            title = 'J. Phys. Chem. 1993, 97, 10269-10280.'
            self.ostream.print_header(title.ljust(str_width))

        elif flag.lower() == 'esp':
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

    def recommended_resp_parameters(self):
        """
        Checks if the RESP parameters are recommended.

        :return:
            True if the RESP parameters are recommended, False otherwise.
        """

        if self.number_layers != 4:
            return False
        if self.density != 1.0:
            return False
        if self.weak_restraint != 0.0005:
            return False
        if self.strong_restraint != 0.001:
            return False
        return True
