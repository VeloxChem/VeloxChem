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
from contextlib import redirect_stdout, redirect_stderr
from pathlib import Path
from io import StringIO
from copy import deepcopy
import numpy as np
import time as tm
import tempfile
import math

from .veloxchemlib import mpi_master, hartree_in_kjpermol, bohr_in_angstrom
from .molecule import Molecule
from .molecularbasis import MolecularBasis
from .outputstream import OutputStream
from .optimizationengine import OptimizationEngine
from .scfrestdriver import ScfRestrictedDriver
from .scfunrestdriver import ScfUnrestrictedDriver
from .scfrestopendriver import ScfRestrictedOpenDriver
from .scfgradientdriver import ScfGradientDriver
from .scfhessiandriver import ScfHessianDriver
from .vibrationalanalysis import VibrationalAnalysis
from .tddftgradientdriver import TddftGradientDriver
from .xtbdriver import XtbDriver
from .xtbgradientdriver import XtbGradientDriver
from .openmmdriver import OpenMMDriver
from .openmmgradientdriver import OpenMMGradientDriver
from .mmdriver import MMDriver
from .mmgradientdriver import MMGradientDriver
from .inputparser import (parse_input, print_keywords,
                          get_random_string_parallel, unparse_input,
                          read_unparsed_input_from_hdf5)
from .errorhandler import assert_msg_critical
from .resultsio import read_molecule_and_basis, write_opt_results_to_hdf5

with redirect_stderr(StringIO()) as fg_err:
    import geometric


class OptimizationDriver:
    """
    Implements optimization driver.

    :param drv:
        The energy or gradient driver.

    Instance variables
        - rank: The rank of MPI process.
        - coordsys: The coordinate system.
        - constraints: The constraints.
        - check_interval: The interval (number of steps) for checking
          coordinate system.
        - max_iter: The maximum number of optimization steps.
        - filename: The filename that will be used by geomeTRIC.
        - grad_drv: The gradient driver.
        - transition: The flag for transition state searching.
        - hessian: The flag for computing Hessian.
    """

    def __init__(self, drv):
        """
        Initializes optimization driver.
        """

        grad_drv = self._pick_driver(drv)

        self.comm = grad_drv.comm
        self.rank = grad_drv.comm.Get_rank()
        self.ostream = grad_drv.ostream

        self.grad_drv = grad_drv

        self.coordsys = 'tric'
        self.constraints = None
        self.check_interval = 0
        self.trust = None
        self.tmax = None
        self.max_iter = 300
        self.conv_maxiter = False

        self.conv_energy = None
        self.conv_grms = None
        self.conv_gmax = None
        self.conv_drms = None
        self.conv_dmax = None

        self.transition = False
        self.irc = False
        self.hessian = 'never'

        self.ref_xyz = None

        self.keep_files = False

        self.filename = None

        self.restart = True

        self._debug = False

        # input keywords
        self.input_keywords = {
            'optimize': {
                'coordsys': ('str_lower', 'coordinate system'),
                'constraints': ('list', 'constraints'),
                'check_interval':
                    ('int', 'interval for checking coordinate system'),
                'trust': ('float', 'trust radius to begin with'),
                'tmax': ('float', 'maximum value of trust radius'),
                'max_iter': ('int', 'maximum number of optimization steps'),
                'transition': ('bool', 'transition state search'),
                'irc': ('bool', 'flag for intrinsic reaction coordinate'),
                'hessian': ('str_lower', 'hessian flag'),
                'ref_xyz': ('str', 'reference geometry'),
                'keep_files': ('bool', 'flag to keep output files'),
                'restart': ('bool', 'flag to restart from checkpoint'),
                'conv_maxiter':
                    ('bool', 'consider converged if max_iter is reached'),
                'conv_energy': ('float', ''),
                'conv_grms': ('float', ''),
                'conv_gmax': ('float', ''),
                'conv_drms': ('float', ''),
                'conv_dmax': ('float', ''),
                '_debug': ('bool', 'print debug info'),
            },
        }

    @property
    def is_scf(self):
        """
        Checks if optimization uses SCF driver.

        :return:
            True if optimization uses SCF driver.
        """

        if hasattr(self.grad_drv, 'scf_driver'):
            return isinstance(self.grad_drv, ScfGradientDriver)
        else:
            return False

    def print_keywords(self):
        """
        Prints input keywords in optimization driver.
        """

        print_keywords(self.input_keywords, self.ostream)

    def update_settings(self, opt_dict):
        """
        Updates settings in optimization driver.

        :param opt_dict:
            The dictionary of optimize input.
        """

        opt_keywords = {
            key: val[0] for key, val in self.input_keywords['optimize'].items()
        }

        parse_input(self, opt_keywords, opt_dict)

        if 'filename' in opt_dict:
            self.filename = opt_dict['filename']

        # update hessian option for transition state search
        if ('hessian' not in opt_dict) and (self.transition or self.irc):
            self.hessian = 'first'

    def read_settings(self, checkpoint_file):
        """
        Reads opt settings from checkpoint file.

        :param checkpoint_file:
            The checkpoint file to read settings from.
        """

        if self.rank == mpi_master():
            checkpoint_opt_input = read_unparsed_input_from_hdf5(
                checkpoint_file, group_name='opt_settings')
        else:
            checkpoint_opt_input = None

        checkpoint_opt_input = self.comm.bcast(checkpoint_opt_input,
                                               root=mpi_master())

        self.update_settings(checkpoint_opt_input)

        self.grad_drv.read_settings(checkpoint_file)

    def _pick_driver(self, drv):
        """
        Chooses the gradient driver.

        :param drv:
            The energy or gradient driver.
        """

        if isinstance(drv, (ScfRestrictedDriver, ScfUnrestrictedDriver,
                            ScfRestrictedOpenDriver)):
            grad_drv = ScfGradientDriver(drv)

        elif isinstance(drv, XtbDriver):
            grad_drv = XtbGradientDriver(drv)

        elif isinstance(drv, OpenMMDriver):
            grad_drv = OpenMMGradientDriver(drv)

        elif isinstance(drv, MMDriver):
            grad_drv = MMGradientDriver(drv)

        elif (isinstance(drv, ScfGradientDriver) or
              isinstance(drv, XtbGradientDriver) or
              isinstance(drv, OpenMMGradientDriver) or
              isinstance(drv, TddftGradientDriver) or
              isinstance(drv, MMGradientDriver)):
            grad_drv = drv

        else:
            assert_msg_critical(
                False,
                'OptimizationDriver: Invalid argument for initialization')

        return grad_drv

    def compute(self, molecule, *args):
        """
        Performs geometry optimization.

        :param molecule:
            The molecule.
        :param args:
            The same arguments as the "compute" function of the gradient driver.

        :return:
            A dictionary containing the results of the geometry optimization.
        """

        # update hessian option for transition state search
        if self.hessian == 'never' and (self.transition or self.irc):
            self.hessian = 'first'
        elif self.hessian == 'last' and (self.transition or self.irc):
            self.hessian = 'first+last'

        if (self.hessian != 'never') or (self.transition or self.irc):
            err_msg = (
                'The installed geometric package does not support\n' +
                '  Hessian or transition state search. Please install\n' +
                '  the latest geometric via pip or conda.\n')
            assert_msg_critical(hasattr(geometric, 'normal_modes'), err_msg)

        self.print_header()

        valstr = 'Using geomeTRIC for geometry optimization.'
        self.ostream.print_info(valstr)
        self.ostream.print_blank()
        self.ostream.print_reference(self.get_geometric_reference())
        self.ostream.print_blank()
        self.ostream.flush()

        start_time = tm.time()

        if self.is_scf:
            basis = args[0]
            if self.restart:
                valid_chkpnt = self.grad_drv.scf_driver.validate_checkpoint(
                    molecule.get_element_ids(), basis.get_label(),
                    self.grad_drv.scf_driver.scf_type)
                if valid_chkpnt:
                    if self.rank == mpi_master():
                        molecule, basis = read_molecule_and_basis(
                            self.grad_drv.scf_driver.get_checkpoint_file())
                    molecule = self.comm.bcast(molecule, root=mpi_master())
                    self.ostream.print_info(
                        'Reading molecular geometry from checkpoint file...')
                    self.ostream.print_blank()
                    self.ostream.flush()

        opt_engine = OptimizationEngine(self.grad_drv, molecule, *args)

        # save unparsed opt_dict in opt_engine for later writing to checkpoint
        opt_keywords = {
            key: val[0] for key, val in self.input_keywords['optimize'].items()
        }
        opt_engine.opt_unparsed_input = unparse_input(self, opt_keywords)

        if self._debug:
            opt_engine._debug = True

        hessian_exit = False
        final_mol = None

        # inherit filename from scf results

        # check that the args contain molecular basis and scf_results
        # note that we only check the type of scf_results on the master rank
        if (isinstance(self.grad_drv, ScfGradientDriver) and
                isinstance(args[0], MolecularBasis) and (len(args) >= 2)):
            args_filename = None
            if self.rank == mpi_master():
                # read filename from scf_results
                if isinstance(args[1], dict) and ('filename' in args[1]) and (
                        args[1]['filename'] is not None):
                    args_filename = args[1]['filename']
            args_filename = self.comm.bcast(args_filename, root=mpi_master())
            # update filename
            if args_filename is not None:
                self.filename = args_filename

        # run within temp_dir since geomeTRIC will generate intermediate files

        try:
            temp_dir = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        except TypeError:
            temp_dir = tempfile.TemporaryDirectory()
        temp_path = Path(temp_dir.name)

        if self.rank == mpi_master() and self.filename is not None:
            self.clean_up_file(Path(self.filename + '.log'))
            self.clean_up_file(
                Path(self.filename + '.tmp', 'hessian', 'hessian.txt'))
            self.clean_up_file(
                Path(self.filename + '.tmp', 'hessian', 'coords.xyz'))

        # filename is used by geomeTRIC to create .log and other files

        if self.filename is not None:
            base_fname = self.filename
        else:
            name_string = get_random_string_parallel(self.comm)
            base_fname = 'vlx_' + name_string

        if self.is_scf and self.grad_drv.scf_driver.filename is None:
            # make sure that the scfdriver has filename
            self.grad_drv.scf_driver.filename = base_fname

        if self.rank == mpi_master() and self.keep_files:
            filename = base_fname
        else:
            filename = Path(base_fname).name
            filename = str(temp_path / f'{filename}_{self.rank}')

        if self.constraints:
            constr_file = Path(filename + '.constr.txt')
            constr_dict = {'freeze': [], 'set': [], 'scan': []}
            for line in self.constraints:
                content = line.strip().split()
                key, val = content[0], ' '.join(content[1:])
                assert_msg_critical(
                    key in ['freeze', 'set', 'scan'],
                    'OptimizationDriver: Invalid constraint {:s}'.format(key))
                constr_dict[key].append(val)
            with constr_file.open('w') as fh:
                for key in ['freeze', 'set', 'scan']:
                    if constr_dict[key]:
                        print(f'${key}', file=fh)
                        for line in constr_dict[key]:
                            print(line, file=fh)
            constr_filename = constr_file.as_posix()

            # self.ostream.print_info('The following constraints are passed to geomeTRIC:')
            # self.ostream.print_blank()
            # with constr_file.open('r') as fh:
            #     for line in fh:
            #         self.ostream.print_header(line.rstrip().ljust(104))
            # self.ostream.print_blank()
            # self.ostream.flush()

        else:
            constr_filename = None

        optinp_filename = Path(filename + '.optinp').as_posix()
        basis = args[0] if args and isinstance(args[0], MolecularBasis) else None

        # prepare for post-opt Hessian

        need_scf_postopt_hessian = False
        if self.is_scf and self.hessian == 'last':
            need_scf_postopt_hessian = True
            self.hessian = 'never'
        elif self.is_scf and self.hessian == 'first+last':
            need_scf_postopt_hessian = True
            self.hessian = 'first'

        # pre-compute Hessian

        if self.is_scf and self.hessian == 'first':
            hessian_drv = ScfHessianDriver(self.grad_drv.scf_driver)
            # conservative choice of disabling restart
            # since geometry is likely changed during an opt calculation
            hessian_drv.cphf_dict = {'restart': False}
            hessian_drv.compute(molecule, args[0])
            if self.rank == mpi_master():
                hess_data = hessian_drv.hessian.copy()
            else:
                hess_data = None
            hess_data = self.comm.bcast(hess_data, root=mpi_master())

            hessian_dir = temp_path / f'rank_{self.rank}'
            hessian_dir.mkdir(parents=True, exist_ok=True)
            hessian_filename = (hessian_dir / 'hessian.txt').as_posix()
            np.savetxt(hessian_filename, hess_data)
            self.hessian = f'file:{hessian_filename}'

        # determine trust radius

        if self.trust is None:
            # from geomeTRIC params.py
            default_trust = 0.01 if self.transition else 0.1
        else:
            default_trust = self.trust

        if self.tmax is None:
            # from geomeTRIC params.py
            default_tmax = (0.03 if self.transition else
                            (default_trust if self.irc else 0.3))
        else:
            default_tmax = self.tmax

        # redirect geomeTRIC stdout/stderr

        with redirect_stdout(StringIO()) as fg_out, redirect_stderr(
                StringIO()) as fg_err:
            try:
                m = geometric.optimize.run_optimizer(
                    customengine=opt_engine,
                    coordsys=self.coordsys,
                    check=self.check_interval,
                    trust=default_trust,
                    tmax=default_tmax,
                    maxiter=self.max_iter,
                    converge=self.conv_flags(),
                    constraints=constr_filename,
                    transition=self.transition,
                    irc=self.irc,
                    hessian=self.hessian,
                    input=optinp_filename)
            except geometric.errors.HessianExit:
                hessian_exit = True

        # post-process and print results while temp_dir is still available

        if hessian_exit:
            if self.rank == mpi_master() and self.hessian == 'stop':
                self.print_vib_analysis(filename, 'vdata_first')
            final_mol = molecule

            opt_results = {'final_geometry': final_mol.get_xyz_string()}

        else:
            coords = m.xyzs[-1] / geometric.nifty.bohr2ang
            labels = molecule.get_labels()
            atom_basis_labels = molecule.get_atom_basis_labels()

            if self.rank == mpi_master():
                final_mol = Molecule(labels, coords.reshape(-1, 3), 'au',
                                     atom_basis_labels)
                final_mol.set_charge(molecule.get_charge())
                final_mol.set_multiplicity(molecule.get_multiplicity())
            else:
                final_mol = None
            final_mol = self.comm.bcast(final_mol, root=mpi_master())

            opt_results = {'final_geometry': final_mol.get_xyz_string()}
            opt_results['final_molecule'] = final_mol

            if self.rank == mpi_master():
                self.grad_drv.ostream.print_info(
                    'Geometry optimization completed.')
                self.ostream.print_blank()

                self.ostream.print_block(
                    final_mol.get_string(title='Final Geometry'))

                is_scan_job = False
                if self.constraints:
                    for line in self.constraints:
                        key = line.strip().split()[0]
                        is_scan_job = (key == 'scan')

                if is_scan_job:
                    self.print_scan_result(m)

                    all_energies = []
                    all_coords_au = []
                    for step, energy, xyz in zip(m.comms, m.qm_energies,
                                                 m.xyzs):
                        if step.split()[:2] == ['Iteration', '0']:
                            all_energies.append([])
                            all_coords_au.append([])
                        all_energies[-1].append(energy)
                        all_coords_au[-1].append(xyz / geometric.nifty.bohr2ang)

                    opt_results['scan_energies'] = [
                        opt_energies[-1] for opt_energies in all_energies
                    ]
                    opt_results['scan_coordinates_au'] = np.array(
                        [opt_coords_au[-1] for opt_coords_au in all_coords_au])

                    opt_results['scan_geometries'] = [
                        self._get_xyz_string(labels, opt_coords_au[-1])
                        for opt_coords_au in all_coords_au
                    ]

                elif self.irc:
                    self.print_irc_result(m)

                    opt_results['irc_energies'] = list(m.qm_energies)
                    opt_results['ts_index'] = int(np.argmax(m.qm_energies))

                    opt_coordinates_au = [
                        xyz / geometric.nifty.bohr2ang for xyz in m.xyzs
                    ]
                    opt_results['irc_geometries'] = [
                        self._get_xyz_string(labels, coords_au)
                        for coords_au in opt_coordinates_au
                    ]

                else:
                    self.print_opt_result(m)

                    opt_results['opt_energies'] = list(m.qm_energies)

                    opt_coordinates_au = [
                        xyz / geometric.nifty.bohr2ang for xyz in m.xyzs
                    ]
                    opt_results['opt_geometries'] = [
                        self._get_xyz_string(labels, coords_au)
                        for coords_au in opt_coordinates_au
                    ]
                    opt_results['opt_coordinates_au'] = np.array(
                        opt_coordinates_au)

                    if self.ref_xyz:
                        self.print_ic_rmsd(final_mol, self.ref_xyz)
                    else:
                        self.print_ic_rmsd(final_mol, molecule)

                if self.hessian in ['last', 'first+last', 'each']:
                    self.print_vib_analysis(filename, 'vdata_last')

                valstr = '*** Time spent in Optimization Driver: '
                valstr += '{:.2f} sec'.format(tm.time() - start_time)
                self.ostream.print_header(valstr)
                self.ostream.print_blank()
                self.ostream.flush()

                # Write opt results to final hdf5 file
                # Note: use base_fname so that the final h5 file is kept even
                # when keep_files is False
                final_h5_fname = base_fname + ".h5"
                # pass in the basis set object if it is in args
                if len(args) >= 1 and isinstance(args[0], MolecularBasis):
                    self._write_final_hdf5(final_h5_fname,
                                           final_mol,
                                           opt_results,
                                           basis=args[0])
                else:
                    self._write_final_hdf5(final_h5_fname, final_mol,
                                           opt_results)

            opt_results = self.comm.bcast(opt_results, root=mpi_master())

        # post-opt Hessian

        if self.is_scf and need_scf_postopt_hessian:
            vib_drv = VibrationalAnalysis(self.grad_drv.scf_driver)
            vib_drv.filename = self.grad_drv.scf_driver.filename
            # conservative choice of disabling restart
            # since geometry is likely changed during an opt calculation
            vib_drv.cphf_dict = {'restart': False}
            vib_results_not_used = vib_drv.compute(molecule, args[0])
            # restore Hessian option
            if self.hessian == 'never':
                need_scf_postopt_hessian = False
                self.hessian == 'last'
            elif self.hessian == 'first':
                need_scf_postopt_hessian = False
                self.hessian == 'first+last'

        try:
            temp_dir.cleanup()
        except (NotADirectoryError, PermissionError):
            pass

        return opt_results

    def conv_flags(self):
        """
        Generates convergence keywords for GeomTRIC.

        :return:
            The convergence keywords for GeomTRIC.
        """

        opt_flags = []
        if self.conv_energy is not None:
            opt_flags.append('energy')
            opt_flags.append(str(self.conv_energy))
        if self.conv_grms is not None:
            opt_flags.append('grms')
            opt_flags.append(str(self.conv_grms))
        if self.conv_gmax is not None:
            opt_flags.append('gmax')
            opt_flags.append(str(self.conv_gmax))
        if self.conv_drms is not None:
            opt_flags.append('drms')
            opt_flags.append(str(self.conv_drms))
        if self.conv_dmax is not None:
            opt_flags.append('dmax')
            opt_flags.append(str(self.conv_dmax))
        if self.conv_maxiter:
            opt_flags.append('maxiter')
        return opt_flags

    @staticmethod
    def clean_up_file(extfile):
        """
        Cleans up existing file.

        :param path_list:
            A list contains the path to the file.
        """

        if extfile.is_file():
            try:
                extfile.unlink()
            except PermissionError:
                pass

    @staticmethod
    def _get_xyz_string(labels, coords_au, precision=12):
        """
        Formats an xyz string from labels and Bohr coordinates.

        :param labels:
            The atomic labels.
        :param coords_au:
            Atomic coordinates in Bohr.
        :param precision:
            Decimal precision for coordinates in Angstrom.

        :return:
            The xyz string.
        """

        coords_angstrom = coords_au * bohr_in_angstrom()
        xyz = f'{len(labels)}\n\n'

        for label, (xa, ya, za) in zip(labels, coords_angstrom):
            xyz += f'{label:<6s}'
            xyz += f' {xa:{precision + 10}.{precision}f}'
            xyz += f' {ya:{precision + 10}.{precision}f}'
            xyz += f' {za:{precision + 10}.{precision}f}\n'

        return xyz

    @staticmethod
    def get_ic_rmsd(opt_mol, ref_mol):
        """
        Gets statistical deviation of bonds, angles and dihedral angles between
        optimized and reference geometry.

        :param opt_mol:
            The optimized molecule.
        :param ref_mol:
            The reference molecule (or xyz filename).

        :return:
            The statistical deviation of bonds, angles and dihedral angles.
        """

        if isinstance(ref_mol, str):
            errmsg = '*** Note: invalid reference xyz file!'
        else:
            errmsg = '*** Note: invalid reference molecule!'

        if isinstance(ref_mol, str):
            if Path(ref_mol).is_file():
                ref_mol = Molecule.read_xyz_file(ref_mol)
            else:
                return errmsg

        if ref_mol.get_labels() != opt_mol.get_labels():
            return errmsg

        g_mol = geometric.molecule.Molecule()
        g_mol.elem = opt_mol.get_labels()
        g_mol.xyzs = [
            opt_mol.get_coordinates_in_bohr() * geometric.nifty.bohr2ang
        ]

        ic = geometric.internal.DelocalizedInternalCoordinates(g_mol,
                                                               build=True)

        ref_geom = ref_mol.get_coordinates_in_bohr() * geometric.nifty.bohr2ang
        opt_geom = opt_mol.get_coordinates_in_bohr() * geometric.nifty.bohr2ang

        bonds = []
        angles = []
        dihedrals = []

        for internal in ic.Prims.Internals:
            if isinstance(internal, geometric.internal.Distance):
                v1 = internal.value(ref_geom)
                v2 = internal.value(opt_geom)
                bonds.append(abs(v1 - v2))
            elif isinstance(internal, geometric.internal.Angle):
                v1 = internal.value(ref_geom)
                v2 = internal.value(opt_geom)
                angles.append(abs(v1 - v2) * 180.0 / np.pi)
            elif isinstance(internal, geometric.internal.Dihedral):
                v1 = internal.value(ref_geom)
                v2 = internal.value(opt_geom)
                diff_in_deg = (v1 - v2) * 180.0 / np.pi
                if diff_in_deg > 180.0:
                    diff_in_deg -= 360.0
                elif diff_in_deg < -180.0:
                    diff_in_deg += 360.0
                dihedrals.append(abs(diff_in_deg))

        ic_rmsd = {'bonds': None, 'angles': None, 'dihedrals': None}

        if bonds:
            np_bonds = np.array(bonds)
            rms_bonds = np.sqrt(np.mean(np_bonds**2))
            max_bonds = np.max(np_bonds)
            ic_rmsd['bonds'] = {
                'rms': rms_bonds,
                'max': max_bonds,
                'unit': 'Angstrom'
            }

        if angles:
            np_angles = np.array(angles)
            rms_angles = np.sqrt(np.mean(np_angles**2))
            max_angles = np.max(np_angles)
            ic_rmsd['angles'] = {
                'rms': rms_angles,
                'max': max_angles,
                'unit': 'degree'
            }

        if dihedrals:
            np_dihedrals = np.array(dihedrals)
            rms_dihedrals = np.sqrt(np.mean(np_dihedrals**2))
            max_dihedrals = np.max(np_dihedrals)
            ic_rmsd['dihedrals'] = {
                'rms': rms_dihedrals,
                'max': max_dihedrals,
                'unit': 'degree'
            }

        return ic_rmsd

    def print_opt_result(self, progress):
        """
        Prints summary of geometry optimization.

        :param progress:
            The geomeTRIC progress of geometry optimization.
        """

        self.ostream.print_blank()
        self.ostream.print_header('Summary of Geometry Optimization')
        self.ostream.print_header('=' * 34)
        self.ostream.print_blank()

        energies = list(progress.qm_energies)
        coords = [xyz / geometric.nifty.bohr2ang for xyz in progress.xyzs]

        line = '{:>8s}{:>20s}  {:>25s}{:>30}'.format('Opt.Step',
                                                     'Energy (a.u.)',
                                                     'Energy Change (a.u.)',
                                                     'Displacement (RMS, Max)')
        self.ostream.print_header(line)
        self.ostream.print_header('-' * len(line))

        for i in range(len(energies)):
            if i > 0:
                delta_e = energies[i] - energies[i - 1]
                rmsd, maxd = geometric.optimize.calc_drms_dmax(
                    coords[i], coords[i - 1])
            else:
                delta_e = 0.0
                rmsd, maxd = 0.0, 0.0
            line = '{:>5d}   {:22.12f}{:22.12f}   {:15.3e}{:15.3e}'.format(
                i, energies[i], delta_e, rmsd, maxd)
            self.ostream.print_header(line)

        self.ostream.print_blank()
        self.ostream.flush()

    def print_scan_result(self, progress):
        """
        Prints summary of geometry scan.

        :param progress:
            The geomeTRIC progress of geometry scan.
        """

        self.ostream.print_blank()
        self.ostream.print_header('Summary of Geometry Scan')
        self.ostream.print_header('=' * 26)
        self.ostream.print_blank()

        scans = []

        for step, energy in zip(progress.comms, progress.qm_energies):
            if step.split()[:2] == ['Iteration', '0']:
                scans.append([])
            scans[-1].append(energy)

        energies = [s[-1] for s in scans]

        e_min = min(energies)
        relative_energies = [e - e_min for e in energies]

        line = '{:>5s}{:>20s}  {:>25s}{:>30s}'.format(
            'Scan', 'Energy (a.u.)', 'Relative Energy (a.u.)',
            'Relative Energy (kJ/mol)')
        self.ostream.print_header(line)
        self.ostream.print_header('-' * len(line))
        for i, (e, rel_e) in enumerate(zip(energies, relative_energies)):
            line = '{:>5d}{:22.12f}{:22.12f}   {:25.10f}     '.format(
                i + 1, e, rel_e, rel_e * hartree_in_kjpermol())
            self.ostream.print_header(line)

        self.ostream.print_blank()
        self.ostream.flush()

    def print_irc_result(self, progress):
        """
        Prints summary of IRC calculation.

        :param progress:
            The geomeTRIC progress of IRC calculation.
        """

        self.ostream.print_blank()
        self.ostream.print_header('Summary of IRC Calculation')
        self.ostream.print_header('=' * 28)
        self.ostream.print_blank()

        energies = list(progress.qm_energies)
        coords = [xyz / geometric.nifty.bohr2ang for xyz in progress.xyzs]

        ts_index = np.argmax(energies)
        e_ts = energies[ts_index]

        line = '{:>8s}{:>21s}{:>25s}{:>23s}{:>28s}{:>7s}'.format(
            'IRC Pt.', 'Energy (a.u.)', 'Energy Change (a.u.)',
            'E - E(TS) (kJ/mol)', 'Displacement (RMS, Max)', '')
        self.ostream.print_header(line)
        self.ostream.print_header('-' * len(line))

        for i in range(len(energies)):
            if i > 0:
                delta_e = energies[i] - energies[i - 1]
                rmsd, maxd = geometric.optimize.calc_drms_dmax(
                    coords[i], coords[i - 1])
            else:
                delta_e = 0.0
                rmsd, maxd = 0.0, 0.0
            
            rel_e_ts = (energies[i] - e_ts) * hartree_in_kjpermol()
            status = '<- TS' if i == ts_index else ''

            line = ('{:>5d}   {:22.12f}{:22.12f}   {:13.3f}'
                    '     {:15.3e}{:15.3e}{:>7s}').format(
                        i, energies[i], delta_e, rel_e_ts, rmsd, maxd,
                        status)
            self.ostream.print_header(line)

        self.ostream.print_blank()

        # Print summary statistics
        e_reactant = energies[0]
        e_product = energies[-1]
        barrier_fwd = (e_ts - e_reactant) * hartree_in_kjpermol()
        barrier_bwd = (e_ts - e_product) * hartree_in_kjpermol()
        reaction_energy = (e_product - e_reactant) * hartree_in_kjpermol()

        title = 'IRC Path Summary'
        self.ostream.print_header(title)
        self.ostream.print_header('-' * (len(title) + 2))

        summary_items = [
            ('Transition State Point', f'{ts_index:12d}'),
            ('Forward Barrier (TS - Reactant)',
             f'{barrier_fwd:12.3f} kJ/mol'),
            ('Backward Barrier (TS - Product)',
             f'{barrier_bwd:12.3f} kJ/mol'),
            ('Reaction Energy (Product - Reactant)',
             f'{reaction_energy:12.3f} kJ/mol'),
        ]
        lines = [
            f'{label:<38s} :    {value}' for label, value in summary_items
        ]
        maxlen = max(len(line) for line in lines)
        for line in lines:
            self.ostream.print_header(line.ljust(maxlen))

        self.ostream.print_blank()
        self.ostream.flush()

    def print_vib_analysis(self, filename, vdata_label):
        """
        Prints summary of vibrational analysis.

        :param vdata_label:
            The label of vdata filename.
        """

        self.ostream.print_blank()
        self.ostream.print_header('Summary of Vibrational Analysis')
        self.ostream.print_header('=' * 33)

        vdata_file = Path(filename + '.' + vdata_label)
        assert_msg_critical(
            vdata_file.is_file(),
            'OptimizationDriver: cannot find vdata file {:s}'.format(
                str(vdata_file)))

        text = []
        with vdata_file.open() as fh:
            for line in fh:
                if line[:2] == '# ':
                    text.append(line[2:].strip())

        maxlen = max([len(line) for line in text])
        for line in text:
            self.ostream.print_header(line.ljust(maxlen))

        self.ostream.print_blank()
        self.ostream.flush()

    def print_ic_rmsd(self, opt_mol, ref_mol):
        """
        Prints statistical deviation of bonds, angles and dihedral angles
        between the optimized geometry and reference geometry.

        :param opt_mol:
            The optimized molecule.
        :param ref_mol:
            The reference molecule (or xyz filename).
        """

        self.ostream.print_blank()

        xyz_filename = ref_mol if isinstance(ref_mol, str) else None

        ic_rmsd = self.get_ic_rmsd(opt_mol, ref_mol)

        if isinstance(ic_rmsd, str):
            self.ostream.print_header(ic_rmsd)
            self.ostream.print_blank()
            self.ostream.flush()
            return

        valstr = 'Statistical Deviation between'
        self.ostream.print_header(valstr)
        if xyz_filename:
            valstr = 'Optimized Geometry and Reference Geometry'
            valstr += ' ({})'.format(xyz_filename)
        else:
            valstr = 'Optimized Geometry and Initial Geometry'
        self.ostream.print_header(valstr)
        self.ostream.print_header((len(valstr) + 2) * '=')
        self.ostream.print_blank()

        valstr = '{:>15s} {:>20s}  {:>21s}'.format('Internal Coord.',
                                                   'RMS deviation',
                                                   'Max. deviation')
        self.ostream.print_header(valstr)
        self.ostream.print_header(len(valstr) * '-')

        if ic_rmsd['bonds']:
            valstr = '{:>12s}    {:12.3f} {:<8s} {:12.3f} {:<8s}'.format(
                'Bonds    ',
                ic_rmsd['bonds']['rms'],
                ic_rmsd['bonds']['unit'],
                ic_rmsd['bonds']['max'],
                ic_rmsd['bonds']['unit'],
            )
            self.ostream.print_header(valstr)
        if ic_rmsd['angles']:
            valstr = '{:>12s}    {:12.3f} {:<8s} {:12.3f} {:<8s}'.format(
                'Angles   ',
                ic_rmsd['angles']['rms'],
                ic_rmsd['angles']['unit'],
                ic_rmsd['angles']['max'],
                ic_rmsd['angles']['unit'],
            )
            self.ostream.print_header(valstr)
        if ic_rmsd['dihedrals']:
            valstr = '{:>12s}    {:12.3f} {:<8s} {:12.3f} {:<8s}'.format(
                'Dihedrals',
                ic_rmsd['dihedrals']['rms'],
                ic_rmsd['dihedrals']['unit'],
                ic_rmsd['dihedrals']['max'],
                ic_rmsd['dihedrals']['unit'],
            )
            self.ostream.print_header(valstr)

        self.ostream.print_blank()
        self.ostream.flush()

    def print_header(self):
        """
        Prints header for optimization driver.
        """

        self.ostream.print_blank()
        self.ostream.print_header('Optimization Driver Setup')
        self.ostream.print_header(27 * '=')
        self.ostream.print_blank()

        lines = []
        lines.append('Coordinate System       :    ' + self.coordsys.upper())
        lines.append('Constraints             :    ' +
                     ('Yes' if self.constraints else 'No'))
        lines.append('Max. Number of Steps    :    ' + str(self.max_iter))
        lines.append('Transition State        :    ' +
                     ('Yes' if self.transition else 'No'))
        lines.append('IRC                     :    ' +
                     ('Yes' if self.irc else 'No'))
        lines.append('Hessian                 :    ' + self.hessian)

        maxlen = max([len(line.split(':')[0]) * 2 for line in lines])
        for line in lines:
            self.ostream.print_header(line.ljust(maxlen))

        self.ostream.print_blank()
        self.ostream.flush()

    def get_geometric_reference(self):
        """
        Gets reference string for geomeTRIC.
        """

        return 'L.-P. Wang and C.C. Song, J. Chem. Phys. 2016, 144, 214108'

    def plot_convergence(self, opt_results):
        """
        Plot the convergence of the optimization (energy plot only).

        :param opt_results:
            The dictionary of optimize results.
        """

        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError('matplotlib is required for this functionality.')

        energies = opt_results['opt_energies']
        total_steps = len(energies) - 1

        min_energy = np.min(energies)
        rel_energies = energies - min_energy
        rel_energies_kJ = rel_energies * hartree_in_kjpermol()

        steps = range(len(rel_energies_kJ))
        x = np.linspace(0, total_steps, 100)
        y = np.interp(x, steps, rel_energies_kJ)

        plt.figure(figsize=(6.5, 4))
        plt.plot(x,
                 y,
                 color='black',
                 alpha=0.9,
                 linewidth=2.5,
                 ls='-',
                 zorder=0)
        plt.scatter(steps,
                    rel_energies_kJ,
                    color='black',
                    alpha=1.0,
                    s=120 / math.log(total_steps, 10),
                    facecolors="darkcyan",
                    edgecolor="darkcyan",
                    zorder=1)
        plt.xlabel('Iteration')
        plt.ylabel('Relative energy [kJ/mol]')
        plt.title("Geometry optimization")
        # Ensure x-axis displays as integers
        plt.xticks(np.arange(0, total_steps + 1, max(1, total_steps // 10)))
        plt.tight_layout()
        plt.show()

    def show_convergence(self, opt_results, atom_indices=False):
        """
        Plot the convergence of the optimization

        :param opt_results:
            The dictionary of optimize results.
        """

        try:
            import ipywidgets
        except ImportError:
            raise ImportError('ipywidgets is required for this functionality.')

        energies = opt_results['opt_energies']
        geometries = opt_results['opt_geometries']
        total_steps = len(energies) - 1
        ipywidgets.interact(self.show_iteration,
                            energies=ipywidgets.fixed(energies),
                            geometries=ipywidgets.fixed(geometries),
                            step=ipywidgets.IntSlider(min=0,
                                                      max=total_steps,
                                                      step=1,
                                                      value=total_steps),
                            atom_indices=ipywidgets.fixed(atom_indices))

    def show_iteration(self, energies, geometries, step=0, atom_indices=False):
        """
        Show the geometry at a specific iteration.
        """

        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError('matplotlib is required for this functionality.')

        min_energy = np.min(energies)
        rel_energies = energies - min_energy
        rel_energies_kJ = rel_energies * hartree_in_kjpermol()

        xyz_data_i = geometries[step]
        steps = range(len(rel_energies_kJ))
        total_steps = len(rel_energies_kJ) - 1
        x = np.linspace(0, total_steps, 100)
        y = np.interp(x, steps, rel_energies_kJ)
        plt.figure(figsize=(6.5, 4))
        plt.plot(x,
                 y,
                 color='black',
                 alpha=0.9,
                 linewidth=2.5,
                 ls='-',
                 zorder=0)
        plt.scatter(steps,
                    rel_energies_kJ,
                    color='black',
                    alpha=0.7,
                    s=120 / math.log(total_steps, 10),
                    facecolors="none",
                    edgecolor="darkcyan",
                    zorder=1)
        plt.scatter(step,
                    rel_energies_kJ[step],
                    marker='o',
                    color='darkcyan',
                    alpha=1.0,
                    s=120 / math.log(total_steps, 10),
                    zorder=2)
        plt.xlabel('Iteration')
        plt.ylabel('Relative energy [kJ/mol]')
        plt.title("Geometry optimization")
        # Ensure x-axis displays as integers
        plt.xticks(np.arange(0, total_steps + 1, max(1, total_steps // 10)))
        plt.tight_layout()
        plt.show()

        mol = Molecule.read_xyz_string(xyz_data_i)
        mol.show(atom_indices=atom_indices, width=640, height=360)

    def plot_scan(self, opt_results):
        """
        Plot the scan energies (energy plot only).

        :param opt_results:
            The dictionary of optimize results.
        """

        try:
            import matplotlib.pyplot as plt
            from scipy.interpolate import CubicSpline
        except ImportError:
            raise ImportError('matplotlib and scipy are required for this functionality.')

        energies = opt_results['scan_energies']
        total_points = len(energies)

        min_energy = np.min(energies)
        rel_energies = energies - min_energy
        rel_energies_kJ = rel_energies * hartree_in_kjpermol()

        scan_points = np.array(range(total_points))

        plt.figure(figsize=(6.5, 4))
        if len(scan_points) > 3:
            cs = CubicSpline(scan_points, rel_energies_kJ)
            x = np.linspace(0, total_points - 1, 300)
            y = cs(x)
            plt.plot(x, y, color='black', alpha=0.9, linewidth=2.5, ls='-', zorder=0)
        else:
            plt.plot(scan_points, rel_energies_kJ, color='black', alpha=0.9, linewidth=2.5, ls='-', zorder=0)
        plt.scatter(scan_points,
                    rel_energies_kJ,
                    color='black',
                    alpha=1.0,
                    s=120,
                    facecolors='darkcyan',
                    edgecolor='darkcyan',
                    zorder=1)
        plt.xlabel('Scan point')
        plt.ylabel('Relative energy [kJ/mol]')
        plt.title("Geometry scan")
        # Ensure x-axis displays as integers
        plt.xticks(np.arange(0, total_points, max(1, total_points // 10)))
        plt.tight_layout()
        plt.show()

    def show_scan(self, opt_results, atom_indices=False):
        """
        Plot the scan energies with interactive molecular visualization

        :param opt_results:
            The dictionary of optimize results.
        :param atom_indices:
            Flag for displaying atom indices.
        """

        try:
            import ipywidgets
        except ImportError:
            raise ImportError('ipywidgets is required for this functionality.')

        energies = opt_results['scan_energies']
        geometries = opt_results['scan_geometries']
        total_points = len(energies)
        ipywidgets.interact(self.show_scan_point,
                            energies=ipywidgets.fixed(energies),
                            geometries=ipywidgets.fixed(geometries),
                            point=ipywidgets.IntSlider(min=0,
                                                       max=total_points - 1,
                                                       step=1,
                                                       value=0),
                            atom_indices=ipywidgets.fixed(atom_indices))

    def show_scan_point(self, energies, geometries, point=0, atom_indices=False):
        """
        Show the geometry at a specific scan point.

        :param energies:
            The energies at each scan point.
        :param geometries:
            The geometries at each scan point.
        :param point:
            The scan point to display.
        :param atom_indices:
            Flag for displaying atom indices.
        """

        try:
            import matplotlib.pyplot as plt
            from scipy.interpolate import CubicSpline
        except ImportError:
            raise ImportError('matplotlib and scipy are required for this functionality.')

        min_energy = np.min(energies)
        rel_energies = energies - min_energy
        rel_energies_kJ = rel_energies * hartree_in_kjpermol()

        xyz_data_i = geometries[point]
        scan_points = np.array(range(len(rel_energies_kJ)))
        total_points = len(rel_energies_kJ)

        plt.figure(figsize=(6.5, 4))
        if len(scan_points) > 3:
            cs = CubicSpline(scan_points, rel_energies_kJ)
            x = np.linspace(0, total_points - 1, 300)
            y = cs(x)
            plt.plot(x, y, color='black', alpha=0.9, linewidth=2.5, ls='-', zorder=0)
        else:
            plt.plot(scan_points, rel_energies_kJ, color='black', alpha=0.9, linewidth=2.5, ls='-', zorder=0)
        plt.scatter(scan_points,
                    rel_energies_kJ,
                    color='black',
                    alpha=0.7,
                    s=120,
                    facecolors="none",
                    edgecolor="darkcyan",
                    zorder=1)
        plt.scatter(point,
                    rel_energies_kJ[point],
                    marker='o',
                    color='darkcyan',
                    alpha=1.0,
                    s=120,
                    zorder=2)
        plt.xlabel('Scan point')
        plt.ylabel('Relative energy [kJ/mol]')
        plt.title("Geometry scan")
        # Ensure x-axis displays as integers
        plt.xticks(np.arange(0, total_points, max(1, total_points // 10)))
        plt.tight_layout()
        plt.show()

        mol = Molecule.read_xyz_string(xyz_data_i)
        mol.show(atom_indices=atom_indices, width=640, height=360)

    def plot_irc(self, opt_results):
        """
        Plot the IRC energies (energy plot only).

        :param opt_results:
            The dictionary of optimize results.
        """

        try:
            import matplotlib.pyplot as plt
            from scipy.interpolate import CubicSpline
        except ImportError:
            raise ImportError('matplotlib and scipy are required for this functionality.')

        energies = opt_results['irc_energies']
        ts_index = opt_results['ts_index']
        total_points = len(energies)

        e_ts = energies[ts_index]
        rel_energies = np.array(energies) - e_ts
        rel_energies_kJ = rel_energies * hartree_in_kjpermol()

        irc_points = np.array(range(total_points))

        plt.figure(figsize=(6.5, 4))
        if len(irc_points) > 3:
            cs = CubicSpline(irc_points, rel_energies_kJ)
            x = np.linspace(0, total_points - 1, 300)
            y = cs(x)
            plt.plot(x, y, color='black', alpha=0.9, linewidth=2.5, ls='-', zorder=0)
        else:
            plt.plot(irc_points, rel_energies_kJ, color='black', alpha=0.9, linewidth=2.5, ls='-', zorder=0)
        
        # Plot all IRC points
        plt.scatter(irc_points,
                    rel_energies_kJ,
                    color='black',
                    alpha=0.7,
                    s=80,
                    facecolors='darkcyan',
                    edgecolor='darkcyan',
                    zorder=1)
        
        # Highlight transition state
        plt.scatter(ts_index,
                rel_energies_kJ[ts_index],
                color='#993399',
                alpha=0.7,
                s=80,
                edgecolor='#993399',
                linewidth=1.5,
                zorder=2,
                label='TS')
        
        plt.xlabel('IRC point')
        plt.ylabel('Energy - E(TS) [kJ/mol]')
        plt.title("Intrinsic Reaction Coordinate")
        plt.legend(frameon=False)
        plt.axhline(y=0, color='gray', linestyle='--', alpha=0.3)
        # Ensure x-axis displays as integers
        plt.xticks(np.arange(0, total_points, max(1, total_points // 10)))
        plt.tight_layout()
        plt.show()

    def show_irc(self, opt_results, atom_indices=False):
        """
        Plot the IRC energies with interactive molecular visualization

        :param opt_results:
            The dictionary of optimize results.
        :param atom_indices:
            Flag for displaying atom indices.
        """

        try:
            import ipywidgets
        except ImportError:
            raise ImportError('ipywidgets is required for this functionality.')

        energies = opt_results['irc_energies']
        geometries = opt_results['irc_geometries']
        ts_index = opt_results['ts_index']
        total_points = len(energies)
        ipywidgets.interact(self.show_irc_point,
                            energies=ipywidgets.fixed(energies),
                            geometries=ipywidgets.fixed(geometries),
                            ts_index=ipywidgets.fixed(ts_index),
                            point=ipywidgets.IntSlider(min=0,
                                                       max=total_points - 1,
                                                       step=1,
                                                       value=ts_index),
                            atom_indices=ipywidgets.fixed(atom_indices))

    def show_irc_point(self, energies, geometries, ts_index, point=0, atom_indices=False):
        """
        Show the geometry at a specific IRC point.

        :param energies:
            The energies at each IRC point.
        :param geometries:
            The geometries at each IRC point.
        :param ts_index:
            The index of the transition state.
        :param point:
            The IRC point to display.
        :param atom_indices:
            Flag for displaying atom indices.
        """

        try:
            import matplotlib.pyplot as plt
            from scipy.interpolate import CubicSpline
        except ImportError:
            raise ImportError('matplotlib and scipy are required for this functionality.')

        e_ts = energies[ts_index]
        rel_energies = np.array(energies) - e_ts
        rel_energies_kJ = rel_energies * hartree_in_kjpermol()

        xyz_data_i = geometries[point]
        irc_points = np.array(range(len(rel_energies_kJ)))
        total_points = len(rel_energies_kJ)

        plt.figure(figsize=(6.5, 4))
        if len(irc_points) > 3:
            cs = CubicSpline(irc_points, rel_energies_kJ)
            x = np.linspace(0, total_points - 1, 300)
            y = cs(x)
            plt.plot(x, y, color='black', alpha=0.9, linewidth=2.5, ls='-', zorder=0)
        else:
            plt.plot(irc_points, rel_energies_kJ, color='black', alpha=0.9, linewidth=2.5, ls='-', zorder=0)
        
        # Plot all IRC points
        plt.scatter(irc_points,
                    rel_energies_kJ,
                    color='black',
                    alpha=0.5,
                    s=80,
                    facecolors="none",
                    edgecolor="darkcyan",
                    zorder=1)
        
        # Highlight transition state
        plt.scatter(ts_index,
                rel_energies_kJ[ts_index],
                color='#993399',
                alpha=0.7,
                s=80,
                edgecolor='#993399',
                linewidth=1.5,
                zorder=2,
                label='TS')
        
        # Highlight current point
        plt.scatter(point,
                    rel_energies_kJ[point],
                    marker='o',
                    color='darkcyan',
                    alpha=1.0,
                    s=120,
                    zorder=3)
        
        plt.xlabel('IRC point')
        plt.ylabel('Energy - E(TS) [kJ/mol]')
        plt.title("Intrinsic Reaction Coordinate")
        plt.axhline(y=0, color='gray', linestyle='--', alpha=0.3)
        # Ensure x-axis displays as integers
        plt.xticks(np.arange(0, total_points, max(1, total_points // 10)))
        plt.tight_layout()
        plt.show()

        mol = Molecule.read_xyz_string(xyz_data_i)
        mol.show(atom_indices=atom_indices, width=640, height=360)

    def _write_final_hdf5(self, fname, molecule, opt_results, basis=None):
        """
        Creates a HDF5 file entry and saves the optimization results.

        :param fname:
            Name of the HDF5 file.
        :param molecule:
            The molecule.
        :param opt_results:
            The dictionary of optimzation results.
        :param basis:
            Optional AO basis set object (for taking care of ECP core electrons).
        """

        if (fname and isinstance(fname, str) and Path(fname).is_file()):
            # Keep Molecule as an in-memory convenience object, and serialize
            # only data-backed optimization results.
            h5_opt_results = {
                key: value
                for key, value in opt_results.items() if key != 'final_molecule'
            }
            write_opt_results_to_hdf5(fname, h5_opt_results)

            valstr = 'Optimization results written to file: '
            valstr += fname
            self.ostream.print_info(valstr)
            self.ostream.print_blank()
            self.ostream.flush()

    def __deepcopy__(self, memo):
        """
        Implements deepcopy.

        :param memo:
            The memo dictionary for deepcopy.

        :return:
            A deepcopy of self.
        """

        new_opt_drv = OptimizationDriver(deepcopy(self.grad_drv))

        for key, val in vars(self).items():
            if isinstance(val, (MPI.Intracomm, OutputStream)):
                continue
            setattr(new_opt_drv, key, deepcopy(val))

        return new_opt_drv
