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

from .veloxchemlib import mpi_master, bohr_in_angstrom
from .outputstream import OutputStream
from .molecularbasis import MolecularBasis
from .scfrestdriver import ScfRestrictedDriver
from .scfunrestdriver import ScfUnrestrictedDriver
from .espchargesdriver import EspChargesDriver
from .inputparser import parse_input, print_keywords
from .errorhandler import assert_msg_critical, safe_solve


class RespChargesDriver(EspChargesDriver):
    """
    Implements RESP charges driver.

    # vlxtag: RHF, RESP_charges

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
        - fitting_max_iter: The maximum number of iterations of the RESP fit.
        - fitting_threshold: The convergence threshold of the RESP fit.
        - filename: The filename for the calculation.
        - xyz_file: The xyz file containing the conformers.
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

        # resp fitting
        self.equal_charges = None
        self.restrain_hydrogen = False
        self.weak_restraint = 0.0005
        self.strong_restraint = 0.001
        self.fitting_max_iter = 50
        self.fitting_threshold = 1.0e-6

        # input keywords
        self._input_keywords = {
            'resp_charges': {
                'restart': ('bool', 'restart flag for SCF driver'),
                'conv_thresh':
                    ('float', 'convergence threshold for SCF driver'),
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
                'fitting_max_iter': ('int', 'maximum iterations in RESP fit'),
                'fitting_threshold':
                    ('float', 'convergence threshold of RESP fit'),
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
            for key, val in self._input_keywords['resp_charges'].items()
        }

        parse_input(self, resp_keywords, resp_dict)

        method_keywords = {
            key: val[0]
            for key, val in self._input_keywords['method_settings'].items()
        }

        parse_input(self, method_keywords, method_dict)

        if 'filename' in resp_dict:
            self.filename = resp_dict['filename']

    def compute(self, molecule, basis=None, scf_results=None, flag=None):
        """
        Computes RESP charges.

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
            flag is None or flag.lower() == 'resp',
            f'{type(self).__name__}.compute: Use of resp/esp flag is ' +
            'deprecated. Please use either RespChargesDriver or ' +
            'EspChargesDriver')

        if isinstance(molecule, list):
            # conformer-weighted resp charges
            if scf_results is not None:
                cur_str = 'Ignoring scf_results since multiple molecules were '
                cur_str += 'provided.'
                self.ostream.print_warning(cur_str)
            molecule_list = molecule
            return self.compute_multiple_molecules(molecule_list, basis)
        else:
            # single molecule resp charges
            return self.compute_single_molecule(molecule, basis, scf_results)

    def compute_single_molecule(self, molecule, basis=None, scf_results=None):
        """
        Computes RSP charges for a single molecule.

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

        # sanity check for RESP
        assert_msg_critical(
            self.grid_type == 'mk',
            f'{type(self).__name__}.compute: For RESP charges, ' +
            'grid_type must be \'mk\'')

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
        else:
            use_631gs_basis = (basis.get_label() in ['6-31G*', '6-31G_D_'])
            if not use_631gs_basis:
                cur_str = 'Recommended basis set 6-31G* is not used!'
                self.ostream.print_warning(cur_str)
                self.ostream.print_blank()

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
            q = self.compute_resp_charges([molecule], [grid_m], [esp_m], [1.0])

        return q

    def compute_multiple_molecules(self, molecule_list, basis_set_list=None):
        """
        Computes conformer-weighted RESP charges.

        :param molecule_list:
            The list of molecules.
        :param basis_set_list:
            The list of molecular basis set.

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
                f'{type(self).__name__}.compute: Invalid equal_charges')

        # sanity check for RESP
        assert_msg_critical(
            self.grid_type == 'mk',
            f'{type(self).__name__}.compute: For RESP charges, ' +
            'grid_type must be \'mk\'')

        molecules = molecule_list
        basis_sets = basis_set_list

        if basis_sets is None:
            if self.rank == mpi_master():
                basis_sets = [
                    MolecularBasis.read(mol, '6-31g*', verbose=False)
                    for mol in molecules
                ]
            basis_sets = self.comm.bcast(basis_sets, root=mpi_master())
        else:
            use_631gs_basis = all([
                bas.get_label() in ['6-31G*', '6-31G_D_'] for bas in basis_sets
            ])
            if not use_631gs_basis:
                cur_str = 'Recommended basis set 6-31G* is not used!'
                self.ostream.print_warning(cur_str)
                self.ostream.print_blank()

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
            q = self.compute_resp_charges(molecules, grids, esp, weights)
        else:
            q = None

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
        self.print_references(len(molecules))
        self.ostream.flush()

        return q

    def optimize_charges(self, molecules, grids, esp, weights, q0, constr,
                         restraint_strength):
        """
        Iterative procedure to optimize charges for a single stage of the ESP fit.

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
            The strength of the hyperbolic penalty function (used in RESP).

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

        # iterative ESP fitting procedure
        for iteration in range(self.fitting_max_iter):
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
            if dq_norm < self.fitting_threshold:
                fitting_converged = True
                break

        errmsg = f'{type(self).__name__}: Charge fitting is not converged!'
        assert_msg_critical(fitting_converged, errmsg)

        cur_str = '*** Charge fitting converged in '
        cur_str += f'{current_iteration + 1} iterations.'
        self.ostream.print_header(cur_str.ljust(40))
        self.ostream.print_blank()

        return q_new[:n_atoms]

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
        cur_str = 'Max. Number of Iterations    :  ' + str(
            self.fitting_max_iter)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Convergence Threshold (a.u.) :  ' + str(
            self.fitting_threshold)
        self.ostream.print_header(cur_str.ljust(str_width))
        self.ostream.print_blank()
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

    def print_references(self, n_conf):
        """
        Prints references.

        :param n_conf:
            The number of conformers.
        """

        str_width = 44

        title = 'Reference: '
        self.ostream.print_header(title.ljust(str_width))

        title = 'J. Phys. Chem. 1993, 97, 10269-10280.'
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
