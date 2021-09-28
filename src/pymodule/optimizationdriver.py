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

from pathlib import PurePath
from pathlib import Path
from os import devnull
import time as tm
import sys
import tempfile
import contextlib
import geometric

from .veloxchemlib import mpi_master
from .veloxchemlib import hartree_in_kcalpermol
from .veloxchemlib import bohr_in_angstroms
from .molecule import Molecule
from .optimizationengine import OptimizationEngine
from .errorhandler import assert_msg_critical
from .inputparser import parse_input
from .veloxchemlib import CommonNeighbors


class OptimizationDriver:
    """
    Implements optimization driver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - rank: The rank of MPI process.
        - coordsys: The coordinate system.
        - constraints: The constraints.
        - check_interval: The interval (number of steps) for checking
          coordinate system.
        - max_iter: The maximum number of optimization steps.
        - filename: The filename that will be used by geomeTRIC.
        - grad_drv: The gradient driver.
        - flag: The type of the optimization driver.
        - cna: The flag for computation of Jaccard similarity index.
        - cna_bond: The cut-off radius for chemical bond in CNA analysis.
        - cna_rcut: The cut-off radius for chemical bonds environment in
          CNA analysis.
    """

    def __init__(self, filename, grad_drv, flag):
        """
        Initializes optimization driver.
        """

        self.comm = grad_drv.comm
        self.rank = grad_drv.comm.Get_rank()
        self.ostream = grad_drv.ostream

        self.coordsys = 'tric'
        self.constraints = None
        self.check_interval = 0
        self.max_iter = 300
 
        self.conv_energy = None 
        self.conv_grms = None 
        self.conv_gmax = None 
        self.conv_drms = None 
        self.conv_dmax = None

        self.ref_xyz = None

        self.filename = filename
        self.grad_drv = grad_drv
        self.flag = flag
        
        self.cna = False
        self.cna_bond = None
        self.cna_rcut = None
        

    def update_settings(self, opt_dict):
        """
        Updates settings in optimization driver.

        :param opt_dict:
            The input dictionary of optimize group.
        """

        opt_keywords = {
            'coordsys': 'str_lower',
            'constraints': 'list',
            'check_interval': 'int',
            'max_iter': 'int',
            'conv_energy': 'float', 
            'conv_grms': 'float', 
            'conv_gmax': 'float', 
            'conv_drms': 'float', 
            'conv_dmax': 'float',
            'ref_xyz': 'str',
            'cna': 'bool',
            'cna_bond': 'float',
            'cna_rcut': 'float',
        }

        parse_input(self, opt_keywords, opt_dict)
        
        # update CNA bond cut-off radius
        if self.cna_bond is None:
            self.cna_bond = 3.0
        else:
            self.cna_bond /= bohr_in_angstroms()
            
        # update CNA bond environment cut-off radius
        if self.cna_rcut is None:
            self.cna_rcut = 4.5
        else:
            self.cna_rcut /= bohr_in_angstroms() 

    def compute(self, molecule, ao_basis, min_basis=None):
        """
        Performs geometry optimization.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param min_basis:
            The minimal AO basis set.

        :return:
            The tuple with final geometry, and energy of molecule.
        """

        self.print_header()
        start_time = tm.time()

        opt_engine = OptimizationEngine(molecule, ao_basis, min_basis,
                                        self.grad_drv, self.flag)

        # filename is used by geomeTRIC to create .log and other files. On
        # master node filename is determined based on the input/output file.
        # On other nodes filename points to file in a temporary directory.

        with tempfile.TemporaryDirectory() as temp_dir:
            if self.rank == mpi_master():
                filename = self.filename
                self.clean_up_file(filename + '.log')
                self.clean_up_file(filename + '.tmp', 'hessian', 'hessian.txt')
                self.clean_up_file(filename + '.tmp', 'hessian', 'coords.xyz')
            else:
                filename = PurePath(self.filename).name
                filename = str(PurePath(temp_dir, f'{filename}_{self.rank}'))

            if self.constraints:
                constr_filename = Path(filename).with_suffix('.constr.txt')
                with open(str(constr_filename), 'w') as fh:
                    for line in self.constraints:
                        print(line, file=fh)
            else:
                constr_filename = None

            log_ini = Path(temp_dir, f'log.ini_{self.rank}')
            self.write_log_ini(log_ini)

            # geomeTRIC prints information to stdout and stderr. On master node
            # this is redirected to the output stream. On other nodes this is
            # redirected to devnull.

            with open(devnull, 'w') as f_devnull:

                if self.rank == mpi_master():
                    f_out = sys.stdout
                else:
                    f_out = f_devnull

                with contextlib.redirect_stdout(f_out):
                    m = geometric.optimize.run_optimizer(
                        customengine=opt_engine,
                        coordsys=self.coordsys,
                        check=self.check_interval,
                        maxiter=self.max_iter,
                        converge=self.conv_flags(), 
                        constraints=constr_filename,
                        input=filename,
                        logIni=str(log_ini))

        coords = m.xyzs[-1] / geometric.nifty.bohr2ang
        labels = molecule.get_labels()

        if self.rank == mpi_master():
            final_mol = Molecule(labels, coords.reshape(-1, 3), units='au') 
            final_ene = m.qm_energies[-1]
        else:
            final_mol = Molecule()
            final_ene = None
        final_mol.broadcast(self.rank, self.comm)

        if self.rank == mpi_master():
            self.grad_drv.ostream.print_info('Geometry optimization completed.')
            self.ostream.print_blank()

            self.ostream.print_block(final_mol.get_string())

            if self.constraints and '$scan' in self.constraints:
                self.print_scan_result(m)
            else:
                self.print_opt_result(m)
                if self.ref_xyz:
                    self.print_ic_rmsd(final_mol, self.ref_xyz)
                else:
                    self.print_ic_rmsd(final_mol, molecule)
                    
            if self.cna:
                self.cna_analysis(final_mol, molecule, self.ref_xyz)
    
            valstr = '*** Time spent in Optimization Driver: '
            valstr += '{:.2f} sec'.format(tm.time() - start_time)
            self.ostream.print_header(valstr)
            self.ostream.print_blank()
            self.ostream.flush()
            
        return (final_mol,final_ene)
    
    def conv_flags(self): 
        """
        Generates convergence keywords for GeomTRIC. 

        :return:
            The convergence keywords for GeomTRIC. 
        """

        opt_flags = []
        if self.conv_energy is not None: 
            opt_flags.append('energy')
            opt_flags.append(self.conv_energy) 
        if self.conv_grms is not None: 
            opt_flags.append('grms') 
            opt_flags.append(self.conv_grms) 
        if self.conv_gmax is not None:
            opt_flags.append('gmax')
            opt_flags.append(self.conv_gmax)
        if self.conv_drms is not None:
            opt_flags.append('drms')
            opt_flags.append(self.conv_drms)
        if self.conv_gmax is not None:
            opt_flags.append('dmax')
            opt_flags.append(self.conv_dmax)
        return opt_flags

    def clean_up_file(self, *path_list):
        """
        Cleans up existing geomeTRIC file.

        :param path_list:
            A list contains the path to the file.
        """

        extfile = Path(*path_list)

        if extfile.is_file():
            extfile.unlink()

    def write_log_ini(self, fname):

        lines = [
            '[loggers]',
            'keys=root',
            '[handlers]',
            'keys=stream_handler, file_handler',
            '[formatters]',
            'keys=formatter',
            '[logger_root]',
            'level=INFO',
            'handlers=stream_handler, file_handler',
            '[handler_stream_handler]',
            'class=geometric.nifty.RawStreamHandler',
            'level=INFO',
            'formatter=formatter',
            'args=(sys.stderr,)',
            '[handler_file_handler]',
            'class=geometric.nifty.RawFileHandler',
            'level=INFO',
            'formatter=formatter',
            'args=(r\'%(logfilename)s\',)',
            '[formatter_formatter]',
            'format=%(message)s',
        ]

        with fname.open('w') as f_ini:
            for line in lines:
                print(line, file=f_ini)

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

    def print_vib_analysis(self, vdata_label):
        """
        Prints summary of vibrational analysis.

        :param vdata_label:
            The label of vdata filename.
        """

        self.ostream.print_blank()
        self.ostream.print_header('Summary of Vibrational Analysis')
        self.ostream.print_header('=' * 33)

        hessian_file = Path(self.filename + '.tmp', 'hessian', 'hessian.txt')
        vdata_file = Path(self.filename + '.' + vdata_label)

        assert_msg_critical(
            hessian_file.is_file(),
            'OptimizationDriver: cannot find hessian file {:s}'.format(
                str(hessian_file)))
        assert_msg_critical(
            vdata_file.is_file(),
            'OptimizationDriver: cannot find vdata file {:s}'.format(
                str(vdata_file)))

        text = []
        with open(str(vdata_file)) as fh:
            for line in fh:
                if line[:2] == '# ':
                    text.append(line[2:].strip())

        maxlen = max([len(line) for line in text])
        for line in text:
            self.ostream.print_header(line.ljust(maxlen))

        valstr = 'Hessian saved in {:s}'.format(str(hessian_file))
        self.ostream.print_header(valstr.ljust(maxlen))
        valstr = 'Frequencies and normal modes saved in {:s}'.format(
            str(vdata_file))
        self.ostream.print_header(valstr.ljust(maxlen))

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

        ic_rmsd = opt_mol.get_ic_rmsd(ref_mol)

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

        maxlen = max([len(line.split(':')[0]) * 2 for line in lines])
        for line in lines:
            self.ostream.print_header(line.ljust(maxlen))

        self.ostream.print_blank()
        self.ostream.flush()

    def cna_analysis(self, last_mol, start_mol, ref_mol):
        """
        Performs common neighbor analysis between the initial, optimized,
        and reference geometries.

        :param last_mol:
            The last i.e. optimized molecule.
        :param start_mol:
            The first i.e. initial molecule.
        :param ref_mol:
            The reference molecule (or xyz filename).
        """

        self.ostream.print_blank()

        xyz_filename = ref_mol if isinstance(ref_mol, str) else None

        self.ostream.print_blank()
        self.ostream.print_header('Summary of Common Neighbor Analysis')
        self.ostream.print_header('=' * 37)
        self.ostream.print_blank()
        
        # compute CNA signatures for first/last molecules
        
        cna_first = CommonNeighbors(start_mol, self.cna_bond)
        cna_first.generate(self.cna_rcut)
        
        cna_last = CommonNeighbors(last_mol, self.cna_bond)
        cna_last.generate(self.cna_rcut)
        
        # compute CNA signatures for reference molecule
        
        cna_ref = None
        if xyz_filename is not None:
            if Path(xyz_filename).is_file():
                ref_mol = Molecule.read_xyz(xyz_filename)
            else:
                return '*** Note: invalid reference xyz file!'
            cna_ref = CommonNeighbors(ref_mol, self.cna_bond)
            cna_ref.generate(self.cna_rcut)
        
        # print signatures
        
        self.print_bonding_pattern(cna_first, 'Initial')
        self.print_bonding_pattern(cna_last, 'Optimized')
        if cna_ref is not None:
            self.print_bonding_pattern(cna_ref, 'Reference')
        
        # print Jaccard similarity indexes
        
        self.ostream.print_header('Jaccard Similarity Indexes')
        self.ostream.print_header('-' * 28)
        fact = cna_last.comp_cna(cna_first) * 100.0
        valstr = 'Initial/Optimized   {:3.2f}'.format(fact)
        self.ostream.print_header(valstr)
        if cna_ref is not None:
            fact = cna_ref.comp_cna(cna_first) * 100.0
            valstr = 'Reference/Initial   {:3.2f}'.format(fact)
            self.ostream.print_header(valstr)
            fact = cna_ref.comp_cna(cna_last) * 100.0
            valstr = 'Reference/Optimized {:3.2f}'.format(fact)
            self.ostream.print_header(valstr)
        self.ostream.print_blank()
        
        # print radius
        valstr = 'Bond cut-off radius: {:3.1f}'.format(self.cna_bond)
        self.ostream.print_header(valstr)
        valstr = 'Environment cut-off radius: {:3.1f}'.format(self.cna_rcut)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        
    def print_bonding_pattern(self, cna_data, label):
        """
        Prints bonding pattern for the given CNA data.

        :param cna_data:
            The CNA data with bonding pattern information.
        :param label:
            The label of bonding patterm header.
        """
        valstr = 'Bonding Pattern of ' + label + ' Geometry'
        self.ostream.print_header(valstr)
        self.ostream.print_header('-' * 41)
        valstr = '{:>9s} {:>15s} {:>15s}'.format('Bond Type',
                                                 'Signature',
                                                 'Repetitions')
        self.ostream.print_header(valstr)
        self.ostream.print_header(len(valstr) * '-')
        self.ostream.print_block(str(cna_data))
        self.ostream.print_blank()
