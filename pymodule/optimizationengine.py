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
from copy import deepcopy
from contextlib import redirect_stderr
from io import StringIO
import numpy as np
import time as tm

from .veloxchemlib import mpi_master
from .veloxchemlib import XCFunctional, MolecularGrid
from .outputstream import OutputStream
from .molecule import Molecule
from .profiler import Profiler

with redirect_stderr(StringIO()) as fg_err:
    import geometric


class OptimizationEngine(geometric.engine.Engine):
    """
    Implements optimization engine for geomeTRIC.

    :param grad_drv:
        The gradient driver.
    :param molecule:
        The molecule.
    :params **:
        The input parameters of the main energy driver

    Instance variables
        - molecule: The molecule.
        - grad_drv: The gradient driver.
        - comm: The MPI communicator.
        - rank: The rank of MPI process.
    """

    def __init__(self, grad_drv, molecule, *args):
        """
        Initializes optimization engine for geomeTRIC.
        """

        g_molecule = geometric.molecule.Molecule()
        g_molecule.elem = molecule.get_labels()
        g_molecule.xyzs = [
            molecule.get_coordinates_in_bohr() * geometric.nifty.bohr2ang
        ]

        super().__init__(g_molecule)

        self.molecule = molecule
        self.args = args
        self.grad_drv = grad_drv

        self.comm = grad_drv.comm
        self.rank = grad_drv.comm.Get_rank()

        self._debug = False

    def lower(self):
        """
        Required in order to get MECI optimization working in geomeTRIC
        """

        return 'custom engine'

    def calc_new(self, coords, dirname):
        """
        Implements calc_new method for the engine.

        :param coords:
            The coordinates.
        :param dirname:
            The relative path.

        :return:
            A dictionary containing energy and gradient.
        """

        start_time = tm.time()

        labels = self.molecule.get_labels()
        atom_basis_labels = self.molecule.get_atom_basis_labels()

        if self.rank == mpi_master():
            new_mol = Molecule(labels, coords.reshape(-1, 3), 'au',
                               atom_basis_labels)
            new_mol.set_charge(self.molecule.get_charge())
            new_mol.set_multiplicity(self.molecule.get_multiplicity())
        else:
            new_mol = Molecule()
        new_mol = self.comm.bcast(new_mol, root=mpi_master())

        self.grad_drv.ostream.print_info('Computing energy and gradient...')
        self.grad_drv.ostream.flush()

        if self._debug:
            profiler = Profiler()
            self.grad_drv.ostream.print_blank()
            self.grad_drv.ostream.print_info(
                '==DEBUG==   available memory before gradient: ' +
                profiler.get_available_memory())
            self.grad_drv.ostream.print_blank()
            self.grad_drv.ostream.flush()

        if not self._debug:
            #self.grad_drv.ostream.mute()
            pass

        energy = self.grad_drv.compute_energy(new_mol, *self.args)
        self.grad_drv.compute(new_mol, *self.args)
        gradient = self.grad_drv.get_gradient()

        if not self._debug:
            #self.grad_drv.ostream.unmute()
            pass

        energy = self.comm.bcast(energy, root=mpi_master())
        gradient = self.comm.bcast(gradient, root=mpi_master())

        if self.rank == mpi_master():
            grad2 = np.sum(gradient**2, axis=1)
            rms_grad = np.sqrt(np.mean(grad2))
            max_grad = np.max(np.sqrt(grad2))
            valstr = '  Energy   : {:.10f} a.u.'.format(energy)
            self.grad_drv.ostream.print_info(valstr)
            valstr = '  Gradient : {:.6e} a.u. (RMS)'.format(rms_grad)
            self.grad_drv.ostream.print_info(valstr)
            valstr = '             {:.6e} a.u. (Max)'.format(max_grad)
            self.grad_drv.ostream.print_info(valstr)
            valstr = '  Time     : {:.2f} sec'.format(tm.time() - start_time)
            self.grad_drv.ostream.print_info(valstr)
            self.grad_drv.ostream.print_blank()
            self.grad_drv.ostream.flush()

        if self._debug:
            profiler = Profiler()
            self.grad_drv.ostream.print_info(
                '==DEBUG==   available memory after  gradient: ' +
                profiler.get_available_memory())
            self.grad_drv.ostream.print_blank()
            self.grad_drv.ostream.flush()

        return {
            'energy': energy,
            'gradient': gradient.flatten(),
        }

    def copy_scratch(self, src, dest):
        """
        Implements copy_scratch method for the engine.

        :param src:
            The source.
        :param dest:
            The destination.
        """

        return

    def __deepcopy__(self, memo):
        """
        Implements deepcopy.

        :param memo:
            The memo dictionary for deepcopy.

        :return:
            A deepcopy of self.
        """

        new_engine = OptimizationEngine(deepcopy(self.grad_drv),
                                        deepcopy(self.molecule),
                                        deepcopy(self.args))

        for key, val in vars(self).items():
            if isinstance(val, (MPI.Intracomm, OutputStream)):
                pass
            elif isinstance(val, XCFunctional):
                new_engine.key = XCFunctional(val)
            elif isinstance(val, MolecularGrid):
                new_engine.key = MolecularGrid(val)
            else:
                new_engine.key = deepcopy(val)

        return new_engine
