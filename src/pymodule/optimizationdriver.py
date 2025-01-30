#
#                              VELOXCHEM
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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

from contextlib import redirect_stdout, redirect_stderr
from pathlib import Path
from io import StringIO
import numpy as np
import time as tm
import tempfile
import math

from .veloxchemlib import mpi_master, hartree_in_kcalpermol
from .molecule import Molecule
from .optimizationengine import OptimizationEngine
from .scfrestdriver import ScfRestrictedDriver
from .scfunrestdriver import ScfUnrestrictedDriver
from .scfrestopendriver import ScfRestrictedOpenDriver
from .scfgradientdriver import ScfGradientDriver
from .scfhessiandriver import ScfHessianDriver
from .xtbdriver import XtbDriver
from .xtbgradientdriver import XtbGradientDriver
from .openmmdriver import OpenMMDriver
from .openmmgradientdriver import OpenMMGradientDriver
from .mmdriver import MMDriver
from .mmgradientdriver import MMGradientDriver
from .inputparser import parse_input, print_keywords, get_random_string_parallel
from .errorhandler import assert_msg_critical

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

        self.coordsys = 'tric'
        self.constraints = None
        self.check_interval = 0
        self.trust = None
        self.tmax = None
        self.max_iter = 300

        self.conv_energy = None
        self.conv_grms = None
        self.conv_gmax = None
        self.conv_drms = None
        self.conv_dmax = None

        self.transition = False
        self.hessian = 'never'

        self.ref_xyz = None

        self.keep_files = False

        self.filename = None
        self.grad_drv = grad_drv

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
                'hessian': ('str_lower', 'hessian flag'),
                'ref_xyz': ('str', 'reference geometry'),
                'keep_files': ('bool', 'flag to keep output files'),
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
            return isinstance(self.grad_drv.scf_driver,
                              (ScfRestrictedDriver, ScfUnrestrictedDriver,
                               ScfRestrictedOpenDriver))
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
        if ('hessian' not in opt_dict) and self.transition:
            self.hessian = 'first'

    def _pick_driver(self, drv):
        """
        Chooses the gradient driver

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
            The tuple with final geometry, and energy of molecule.
        """

        # update hessian option for transition state search
        if self.hessian == 'never' and self.transition:
            self.hessian = 'first'

        if self.hessian or self.transition:
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

        opt_engine = OptimizationEngine(self.grad_drv, molecule, *args)

        if self._debug:
            opt_engine._debug = True

        hessian_exit = False
        final_mol = None

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

        if self.is_scf and self.grad_drv.scf_driver.checkpoint_file is None:
            # make sure that the scfdriver has checkpoint_file
            fpath = Path(base_fname)
            fpath = fpath.with_name(f'{fpath.stem}_scf{fpath.suffix}')
            fpath = fpath.with_suffix('.h5')
            self.grad_drv.scf_driver.checkpoint_file = str(fpath)

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
        else:
            constr_filename = None

        optinp_filename = Path(filename + '.optinp').as_posix()

        # pre-compute Hessian

        if self.is_scf and self.hessian == 'first':
            hessian_drv = ScfHessianDriver(self.grad_drv.scf_driver)
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
            default_tmax = 0.03 if self.transition else 0.3
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

            if self.rank == mpi_master():
                self.grad_drv.ostream.print_info(
                    'Geometry optimization completed.')
                self.ostream.print_blank()

                self.ostream.print_block(final_mol.get_string())

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

                    opt_results['scan_geometries'] = []
                    labels = molecule.get_labels()
                    for opt_coords_au in all_coords_au:
                        mol = Molecule(labels, opt_coords_au[-1], 'au',
                                       atom_basis_labels)
                        opt_results['scan_geometries'].append(
                            mol.get_xyz_string())

                else:
                    self.print_opt_result(m)

                    opt_results['opt_energies'] = list(m.qm_energies)

                    opt_results['opt_geometries'] = []
                    labels = molecule.get_labels()
                    for xyz in m.xyzs:
                        mol = Molecule(labels, xyz / geometric.nifty.bohr2ang,
                                       'au', atom_basis_labels)
                        opt_results['opt_geometries'].append(
                            mol.get_xyz_string())

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

            opt_results = self.comm.bcast(opt_results, root=mpi_master())

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
        return opt_flags

    @staticmethod
    def clean_up_file(extfile):
        """
        Cleans up existing file.

        :param path_list:
            A list contains the path to the file.
        """

        if extfile.is_file():
            extfile.unlink()

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
            'Relative Energy (kcal/mol)')
        self.ostream.print_header(line)
        self.ostream.print_header('-' * len(line))
        for i, (e, rel_e) in enumerate(zip(energies, relative_energies)):
            line = '{:>5d}{:22.12f}{:22.12f}   {:25.10f}     '.format(
                i + 1, e, rel_e, rel_e * hartree_in_kcalpermol())
            self.ostream.print_header(line)

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
        rel_energies_kcal = rel_energies * hartree_in_kcalpermol()

        xyz_data_i = geometries[step]
        steps = range(len(rel_energies_kcal))
        total_steps = len(rel_energies_kcal) - 1
        x = np.linspace(0, total_steps, 100)
        y = np.interp(x, steps, rel_energies_kcal)
        plt.figure(figsize=(6.5, 4))
        plt.plot(x,
                 y,
                 color='black',
                 alpha=0.9,
                 linewidth=2.5,
                 ls='-',
                 zorder=0)
        plt.scatter(steps,
                    rel_energies_kcal,
                    color='black',
                    alpha=0.7,
                    s=120 / math.log(total_steps, 10),
                    facecolors="none",
                    edgecolor="darkcyan",
                    zorder=1)
        plt.scatter(step,
                    rel_energies_kcal[step],
                    marker='o',
                    color='darkcyan',
                    alpha=1.0,
                    s=120 / math.log(total_steps, 10),
                    zorder=2)
        plt.xlabel('Iteration')
        plt.ylabel('Relative energy [kcal/mol]')
        plt.title("Geometry optimization")
        # Ensure x-axis displays as integers
        plt.xticks(np.arange(0, total_steps + 1, max(1, total_steps // 10)))
        plt.tight_layout()
        plt.show()

        mol = Molecule.read_xyz_string(xyz_data_i)
        mol.show(atom_indices=atom_indices, width=640, height=360)
