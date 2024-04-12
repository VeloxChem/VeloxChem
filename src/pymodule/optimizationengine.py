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

    def lower(self):
        """ Required in order to get MECI optimization working in geomeTRIC"""
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

        if self.rank == mpi_master():
            new_mol = Molecule(labels, coords.reshape(-1, 3), units='au')
            new_mol.set_charge(self.molecule.get_charge())
            new_mol.set_multiplicity(self.molecule.get_multiplicity())
        else:
            new_mol = Molecule()
        new_mol.broadcast(self.rank, self.comm)

        self.grad_drv.ostream.print_info('Computing energy and gradient...')
        self.grad_drv.ostream.flush()

        self.grad_drv.ostream.mute()

        energy = self.grad_drv.compute_energy(new_mol, *self.args)
        self.grad_drv.compute(new_mol, *self.args)
        gradient = self.grad_drv.get_gradient()

        self.grad_drv.ostream.unmute()

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
