import os
import sys
import tempfile
import contextlib
import geometric

from .veloxchemlib import mpi_master
from .veloxchemlib import hartree_in_kcalpermol
from .molecule import Molecule
from .optimizationengine import OptimizationEngine


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
            else:
                filename = os.path.join(temp_dir, 'tmp_{:d}'.format(self.rank))

            if self.constraints:
                constr_filename = filename + '_constr.txt'
                with open(constr_filename, 'w') as fh:
                    fh.write(os.linesep.join(self.constraints))
            else:
                constr_filename = None

            # geomeTRIC prints information to stdout and stderr. On master node
            # this is redirected to the output stream. On other nodes this is
            # redirected to os.devnull.

            with open(os.devnull, 'w') as devnull:

                if self.rank == mpi_master():
                    fh = sys.stdout
                else:
                    fh = devnull

                with contextlib.redirect_stdout(fh):
                    with contextlib.redirect_stderr(fh):
                        m = geometric.optimize.run_optimizer(
                            customengine=opt_engine,
                            coordsys=self.coordsys,
                            check=self.check_interval,
                            maxiter=self.max_iter,
                            constraints=constr_filename,
                            transition=self.transition,
                            hessian=self.hessian,
                            input=filename)

        coords = m.xyzs[-1] / geometric.nifty.bohr2ang
        labels = molecule.get_labels()

        if self.rank == mpi_master():
            final_mol = Molecule(labels, coords.reshape(-1, 3), units='au')
        else:
            final_mol = Molecule()
        final_mol.broadcast(self.rank, self.comm)

        if self.constraints and '$scan' in self.constraints:
            if self.rank == mpi_master():
                self.print_scan_result(m)

        return final_mol

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
