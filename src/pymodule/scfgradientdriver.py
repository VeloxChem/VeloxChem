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

from mpi4py import MPI
from copy import deepcopy
import numpy as np
import time as tm

from .veloxchemlib import (OverlapGeom100Driver, KineticEnergyGeom100Driver,
                           NuclearPotentialGeom100Driver,
                           NuclearPotentialGeom010Driver, FockGeom1000Driver)
from .veloxchemlib import XCFunctional, MolecularGrid, XCMolecularGradient
from .veloxchemlib import mpi_master, mat_t
from .veloxchemlib import partition_atoms, make_matrix
from .veloxchemlib import parse_xc_func
from .griddriver import GridDriver
from .outputstream import OutputStream
from .gradientdriver import GradientDriver
from .dftutils import get_default_grid_level
from .errorhandler import assert_msg_critical


class ScfGradientDriver(GradientDriver):
    """
    Implements SCF gradient driver.

    :param scf_drv:
        The SCF driver.

    Instance variables
        - flag: The driver flag.
        - delta_h: The displacement for finite difference.
    """

    def __init__(self, scf_drv):
        """
        Initializes gradient driver.
        """

        super().__init__(scf_drv.comm, scf_drv.ostream)

        self.scf_driver = scf_drv
        self.flag = 'SCF Gradient Driver'

        self.numerical = False
        self.delta_h = 0.001

    def compute(self, molecule, basis, scf_results):
        """
        Performs calculation of gradient.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_results:
            The dictionary containing converged SCF results.
        """

        start_time = tm.time()
        self.print_header()

        if self.rank == mpi_master():
            nocc = molecule.number_of_alpha_electrons()

            D = scf_results['D_alpha']

            ene_occ = scf_results['E_alpha'][:nocc]
            mo_occ = scf_results['C_alpha'][:, :nocc].copy()
            W = np.linalg.multi_dot([mo_occ, np.diag(ene_occ), mo_occ.T])
        else:
            D = None
            W = None

        D = self.comm.bcast(D, root=mpi_master())
        W = self.comm.bcast(W, root=mpi_master())

        natoms = molecule.number_of_atoms()

        self.gradient = np.zeros((natoms, 3))

        local_atoms = partition_atoms(natoms, self.rank, self.nodes)

        # kinetic energy contribution to gradient

        kin_grad_drv = KineticEnergyGeom100Driver()

        for iatom in local_atoms:
            gmats = kin_grad_drv.compute(molecule, basis, iatom)

            for i, label in enumerate(['X', 'Y', 'Z']):
                gmat = gmats.matrix(label).full_matrix().to_numpy()
                self.gradient[iatom, i] += 2.0 * np.sum((gmat + gmat.T) * D)

        # nuclear potential contribution to gradient

        npot_grad_100_drv = NuclearPotentialGeom100Driver()
        npot_grad_010_drv = NuclearPotentialGeom010Driver()

        for iatom in local_atoms:
            gmats_100 = npot_grad_100_drv.compute(molecule, basis, iatom)
            gmats_010 = npot_grad_010_drv.compute(molecule, basis, iatom)

            for i, label in enumerate(['X', 'Y', 'Z']):
                gmat_100 = gmats_100.matrix(label).full_matrix().to_numpy()
                gmat_010 = gmats_010.matrix(label).full_matrix().to_numpy()

                # TODO: move minus sign into function call (such as in oneints)
                self.gradient[iatom, i] -= 2.0 * np.sum(
                    (gmat_100 + gmat_100.T) * D)
                self.gradient[iatom, i] -= 2.0 * np.sum(gmat_010 * D)

        # orbital contribution to gradient

        ovl_grad_drv = OverlapGeom100Driver()

        for iatom in local_atoms:
            gmats = ovl_grad_drv.compute(molecule, basis, iatom)

            for i, label in enumerate(['X', 'Y', 'Z']):
                gmat = gmats.matrix(label).full_matrix().to_numpy()
                # Note: minus sign for energy weighted density
                self.gradient[iatom, i] -= 2.0 * np.sum((gmat + gmat.T) * W)

        # ERI contribution to gradient

        if self.rank == mpi_master():
            use_dft = 'xcfun' in scf_results
        else:
            use_dft = False
        use_dft = self.comm.bcast(use_dft, root=mpi_master())

        # determine fock_type and exchange_scaling_factor
        if use_dft:
            if self.rank == mpi_master():
                xcfun = scf_results['xcfun']
            else:
                xcfun = None
            xcfun = self.comm.bcast(xcfun, root=mpi_master())
            xcfun = parse_xc_func(xcfun)

            if xcfun.is_hybrid():
                fock_type = '2jkx'
                exchange_scaling_factor = xcfun.get_frac_exact_exchange()
            else:
                fock_type = 'j'
                exchange_scaling_factor = 0.0
        else:
            fock_type = '2jk'
            exchange_scaling_factor = 1.0

        # TODO: range-separated Fock

        den_mat_for_fock = make_matrix(basis, mat_t.symmetric)
        den_mat_for_fock.set_values(D)

        fock_grad_drv = FockGeom1000Driver()

        for iatom in local_atoms:
            gmats = fock_grad_drv.compute(basis, molecule, den_mat_for_fock,
                                          iatom, fock_type,
                                          exchange_scaling_factor, 0.0)

            factor = 2.0 if fock_type == 'j' else 1.0

            for i, label in enumerate(['X', 'Y', 'Z']):
                gmat = gmats.matrix(label).full_matrix().to_numpy()
                self.gradient[iatom, i] += np.sum(gmat * D) * factor

        # XC contribution to gradient

        if use_dft:

            if self.rank == mpi_master():
                xcfun = scf_results['xcfun']
                grid_level = scf_results.get('grid_level', None)
            else:
                xcfun = None
                grid_level = None

            xcfun, grid_level = self.comm.bcast((xcfun, grid_level),
                                                root=mpi_master())

            grid_level = (get_default_grid_level(xcfun)
                          if grid_level is None else grid_level)

            # TODO: take molecular grid from scf
            grid_drv = GridDriver(self.comm)
            grid_drv.set_level(grid_level)
            mol_grid = grid_drv.generate(molecule)

            grad_drv = XCMolecularGradient()
            self.gradient += grad_drv.integrate_vxc_gradient(
                molecule, basis, D, mol_grid, xcfun)

        # nuclear contribution to gradient (only added on master rank)

        if self.rank == mpi_master():
            self.gradient += self.grad_nuc_contrib(molecule)

        # collect gradient

        self.gradient = self.comm.allreduce(self.gradient, op=MPI.SUM)

        # print gradient
        self.print_geometry(molecule)
        self.print_gradient(molecule)

        valstr = '*** Time spent in gradient calculation: '
        valstr += '{:.2f} sec ***'.format(tm.time() - start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

    def compute_numerical_gradient(self, molecule, ao_basis, scf_results):
        """
        Performs calculation of gradient.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param scf_results:
            The dictionary containing converged SCF results.
        """

        start_time = tm.time()
        self.print_header()

        self.ostream.mute()
        # Currently, only numerical gradients activated
        self.compute_numerical(molecule, ao_basis, scf_results)
        self.ostream.unmute()

        # print gradient
        self.print_geometry(molecule)
        self.print_gradient(molecule)

        valstr = '*** Time spent in gradient calculation: '
        valstr += '{:.2f} sec ***'.format(tm.time() - start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

    def compute_energy(self, molecule, ao_basis, scf_results):
        """
        Computes the energy at current geometry.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param scf_results:
            The dictionary containing converged SCF results.

        :return:
            The energy.
        """

        self.ostream.mute()
        self.scf_driver.restart = False
        new_scf_results = self.scf_driver.compute(molecule, ao_basis)
        assert_msg_critical(self.scf_driver.is_converged,
                            'ScfGradientDriver: SCF did not converge')
        self.ostream.unmute()

        if self.rank == mpi_master():
            scf_results.update(new_scf_results)

        return self.scf_driver.get_scf_energy()

    def __deepcopy__(self, memo):
        """
        Implements deepcopy.

        :param memo:
            The memo dictionary for deepcopy.

        :return:
            A deepcopy of self.
        """

        new_grad_drv = ScfGradientDriver(self.comm, self.ostream)

        for key, val in vars(self).items():
            if isinstance(val, (MPI.Intracomm, OutputStream)):
                pass
            elif isinstance(val, XCFunctional):
                new_grad_drv.key = XCFunctional(val)
            elif isinstance(val, MolecularGrid):
                new_grad_drv.key = MolecularGrid(val)
            else:
                new_grad_drv.key = deepcopy(val)

        return new_grad_drv
