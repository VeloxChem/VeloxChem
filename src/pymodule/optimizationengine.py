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
from .outputstream import OutputStream
from .scfgradientdriver import ScfGradientDriver
from .molecule import Molecule
from .inputparser import write_unparsed_input_to_hdf5
from .errorhandler import assert_msg_critical

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

    def __init__(self, grad_drv, molecule, *args, frozen_idx=None,
                 frozen_method='zero'):
        """
        Initializes optimization engine for geomeTRIC.

        Frozen atoms are hidden from geomeTRIC. With the 'reduce' method only
        the active (non-frozen) atoms are exposed; with the 'zero' method the
        full molecule is exposed but the frozen-atom gradient is zeroed.
        """

        labels = molecule.get_labels()
        natoms = len(labels)

        if frozen_idx is None:
            frozen_idx = np.array([], dtype=int)
        self.frozen_idx = np.asarray(frozen_idx, dtype=int)
        frozen_set = set(self.frozen_idx.tolist())
        self.active_idx = np.array(
            [i for i in range(natoms) if i not in frozen_set], dtype=int)
        self.frozen_method = frozen_method

        # fixed positions of the frozen atoms (Bohr), held constant throughout
        coords_in_bohr = molecule.get_coordinates_in_bohr()
        self.frozen_coords = coords_in_bohr[self.frozen_idx]

        g_molecule = geometric.molecule.Molecule()
        if self.frozen_idx.size > 0 and frozen_method == 'reduce':
            # expose only the active atoms to geomeTRIC
            g_molecule.elem = [labels[i] for i in self.active_idx]
            g_molecule.xyzs = [
                coords_in_bohr[self.active_idx] * geometric.nifty.bohr2ang
            ]
        else:
            # full molecule ('zero' method, or no frozen atoms)
            g_molecule.elem = labels
            g_molecule.xyzs = [coords_in_bohr * geometric.nifty.bohr2ang]
        g_molecule.build_topology()

        super().__init__(g_molecule)

        self.molecule = molecule
        self.args = args
        self.grad_drv = grad_drv

        self.comm = grad_drv.comm
        self.rank = grad_drv.comm.Get_rank()

        self.opt_unparsed_input = None
        self.opt_current_step = 0

    def lower(self):
        """
        Required in order to get MECI optimization working in geomeTRIC
        """

        return 'custom engine'

    def _validate_step_results(self, energy, gradient, natoms):
        """Validates optimizer step results before returning to geomeTRIC.

        :param energy:
            The energy returned by the driver stack.
        :param gradient:
            The gradient returned by the driver stack.
        :param natoms:
            Number of atoms in the current molecule.
        """

        assert_msg_critical(
            np.isscalar(energy),
            'OptimizationEngine.calc_new: expected scalar energy after MPI '
            'synchronization')

        try:
            gradient = np.asarray(gradient, dtype=float)
        except (TypeError, ValueError):
            gradient = None

        assert_msg_critical(
            gradient is not None,
            'OptimizationEngine.calc_new: expected array-like gradient after '
            'MPI synchronization')

        return float(energy), gradient

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

        # geomeTRIC passes only the active atoms with the 'reduce' method;
        # splice them back into the full geometry with the frozen atoms held
        # at their fixed positions
        if self.frozen_idx.size > 0 and self.frozen_method == 'reduce':
            full_coords = self.molecule.get_coordinates_in_bohr()
            full_coords[self.active_idx] = coords.reshape(-1, 3)
            full_coords[self.frozen_idx] = self.frozen_coords
        else:
            full_coords = coords.reshape(-1, 3)
            # 'zero' method: geomeTRIC sees the full molecule and may drift the
            # frozen atoms (their gradient is zeroed); always evaluate energy
            # and gradient with the frozen atoms pinned to their fixed positions
            if self.frozen_idx.size > 0:
                full_coords = full_coords.copy()
                full_coords[self.frozen_idx] = self.frozen_coords

        if self.rank == mpi_master():
            new_mol = Molecule(labels, full_coords, 'au', atom_basis_labels)
            new_mol.set_charge(self.molecule.get_charge())
            new_mol.set_multiplicity(self.molecule.get_multiplicity())
        else:
            new_mol = Molecule()
        new_mol = self.comm.bcast(new_mol, root=mpi_master())

        title_txt = f'Optimization Step {self.opt_current_step}'
        self.grad_drv.ostream.print_header(title_txt)
        self.grad_drv.ostream.print_header('=' * (len(title_txt) + 2))
        self.grad_drv.ostream.print_blank()
        self.grad_drv.ostream.print_info('Computing energy and gradient...')
        self.grad_drv.ostream.flush()

        energy = self.grad_drv.compute_energy(new_mol, *self.args)
        self.grad_drv.compute(new_mol, *self.args)
        gradient = self.grad_drv.get_gradient()

        if hasattr(self.grad_drv, "scf_driver") and isinstance(
                self.grad_drv, ScfGradientDriver):
            checkpoint_file = self.grad_drv.scf_driver.get_checkpoint_file()
            if checkpoint_file is not None and self.opt_unparsed_input is not None:
                if self.rank == mpi_master():
                    write_unparsed_input_to_hdf5(checkpoint_file,
                                                 self.opt_unparsed_input,
                                                 group_name='opt_settings')

        energy = self.comm.bcast(energy, root=mpi_master())
        gradient = self.comm.bcast(gradient, root=mpi_master())
        energy, gradient = self._validate_step_results(
            energy, gradient, new_mol.number_of_atoms())

        self.opt_current_step += 1

        # build the gradient seen by geomeTRIC: active rows only for 'reduce',
        # full gradient with zeroed frozen rows for 'zero'
        if self.frozen_idx.size > 0 and self.frozen_method == 'reduce':
            eff_gradient = gradient[self.active_idx]
        elif self.frozen_idx.size > 0:
            eff_gradient = gradient.copy()
            eff_gradient[self.frozen_idx] = 0.0
        else:
            eff_gradient = gradient

        if self.rank == mpi_master():
            grad2 = np.sum(eff_gradient**2, axis=1)
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
            'gradient': eff_gradient.flatten(),
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
                                        *deepcopy(self.args),
                                        frozen_idx=deepcopy(self.frozen_idx),
                                        frozen_method=self.frozen_method)

        for key, val in vars(self).items():
            if isinstance(val, (MPI.Intracomm, OutputStream)):
                continue
            setattr(new_engine, key, deepcopy(val))

        return new_engine
