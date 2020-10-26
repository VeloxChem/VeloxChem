from pathlib import PurePath
from pathlib import Path
from os import devnull
import sys
import tempfile
import contextlib
import geometric

from .veloxchemlib import mpi_master
from .veloxchemlib import hartree_in_kcalpermol
from .molecule import Molecule
from .optimizationengine import OptimizationEngine
from .errorhandler import assert_msg_critical


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
        - max_iter: The maximum number of optimization steps
        - transition: The flag for transition state searching.
        - hessian: The flag for computing Hessian.
        - filename: The filename that will be used by geomeTRIC.
        - grad_drv: The gradient driver.
        - flag: The type of the optimization driver.
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

        self.transition = False
        self.hessian = 'never'

        self.filename = filename
        self.grad_drv = grad_drv
        self.flag = flag

    def update_settings(self, opt_dict):
        """
        Updates settings in optimization driver.

        :param opt_dict:
            The input dictionary of optimize group.
        """

        if 'coordsys' in opt_dict:
            self.coordsys = opt_dict['coordsys'].lower()
        if 'constraints' in opt_dict:
            self.constraints = list(opt_dict['constraints'])

        if 'check_interval' in opt_dict:
            self.check_interval = int(opt_dict['check_interval'])
        elif 'check' in opt_dict:
            self.check_interval = int(opt_dict['check'])

        if 'max_iter' in opt_dict:
            self.max_iter = int(opt_dict['max_iter'])
        elif 'maxiter' in opt_dict:
            self.max_iter = int(opt_dict['maxiter'])

        if 'transition' in opt_dict:
            key = opt_dict['transition'].lower()
            self.transition = True if key == 'yes' else False

        if 'hessian' in opt_dict:
            self.hessian = opt_dict['hessian'].lower()
            if self.hessian == 'only':
                self.hessian = 'stop'
        elif self.transition:
            self.hessian = 'first'

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
            The molecule with final geometry.
        """

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
                filename = str(
                    PurePath(temp_dir, '{:s}_{:d}'.format(filename, self.rank)))

            if self.constraints:
                constr_filename = filename + '_constr.txt'
                with open(constr_filename, 'w') as fh:
                    for line in self.constraints:
                        print(line, file=fh)
            else:
                constr_filename = None

            # geomeTRIC prints information to stdout and stderr. On master node
            # this is redirected to the output stream. On other nodes this is
            # redirected to devnull.

            with open(devnull, 'w') as f_devnull:

                if self.rank == mpi_master():
                    f_out = sys.stdout
                else:
                    f_out = f_devnull

                with contextlib.redirect_stdout(f_out):
                    with contextlib.redirect_stderr(f_out):
                        try:
                            m = geometric.optimize.run_optimizer(
                                customengine=opt_engine,
                                coordsys=self.coordsys,
                                check=self.check_interval,
                                maxiter=self.max_iter,
                                constraints=constr_filename,
                                transition=self.transition,
                                hessian=self.hessian,
                                input=filename)
                        except geometric.errors.HessianExit:
                            if (self.rank == mpi_master() and
                                    self.hessian == 'stop'):
                                self.print_vib_analysis('vdata_first')
                            return molecule

        coords = m.xyzs[-1] / geometric.nifty.bohr2ang
        labels = molecule.get_labels()

        if self.rank == mpi_master():
            final_mol = Molecule(labels, coords.reshape(-1, 3), units='au')
        else:
            final_mol = Molecule()
        final_mol.broadcast(self.rank, self.comm)

        if self.rank == mpi_master():
            if self.constraints and '$scan' in self.constraints:
                self.print_scan_result(m)
            else:
                self.print_opt_result(m)
            if self.hessian in ['last', 'first+last', 'each']:
                self.print_vib_analysis('vdata_last')

        return final_mol

    def clean_up_file(self, *path_list):
        """
        Cleans up existing geomeTRIC file.

        :param path_list:
            A list contains the path to the file.
        """

        extfile = Path(*path_list)

        if extfile.is_file():
            extfile.unlink()

    def print_opt_result(self, progress):
        """
        Prints summary of geometry optimization.

        :param progress:
            The geomeTRIC progress of geometry scan.
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
                rmsd, maxd = geometric.step.calc_drms_dmax(
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
