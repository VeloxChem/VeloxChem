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

import time as tm
from mpi4py import MPI
import numpy as np
import sys

from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .veloxchemlib import XCFunctional
from .veloxchemlib import XCIntegrator
from veloxchem.veloxchemlib import GridDriver
from veloxchem.veloxchemlib import XCMolecularGradient
from .veloxchemlib import MolecularGrid
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import denmat
from .veloxchemlib import parse_xc_func
from .distributedarray import DistributedArray
from .errorhandler import assert_msg_critical
from .molecule import Molecule
from .outputstream import OutputStream


class GradientDriver:
    """
    Implements gradient driver.

    :param energy_drv:
        The energy driver.
    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - gradient: The gradient.
        - flag: The type of gradient driver.
        - numerical: Perform numerical gradient calculation.
        - do_four_point: Perform numerical gradient calculation using the five-point stencil method?
        - dft: The flag for running DFT.
        - grid_level: The accuracy level of DFT grid.
        - xcfun: The XC functional.

    """

    def __init__(self, energy_drv, comm=None, ostream=None):
        """
        Initializes gradient driver.
        """

        if comm is None:
            if hasattr(energy_drv, 'comm'):
                comm = energy_drv.comm
            else:
                comm = MPI.COMM_WORLD

        if ostream is None:
            if hasattr(energy_drv, 'ostream'):
                ostream = energy_drv.ostream
            else:
                if comm.Get_rank() == mpi_master():
                    ostream = OutputStream(sys.stdout)
                else:
                    ostream = OutputStream(None)

        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        self.ostream = ostream

        self.gradient = None
        self.flag = None
        self.numerical = False
        # flag for two-point or four-point approximation
        self.do_four_point = False

        # DFT information
        self.dft = False
        self.grid_level = 4
        self.xcfun = XCFunctional()

    def update_settings(self, grad_dict, method_dict):
        """
        Updates settings in GradientDriver.

        :param grad_dict:
            The input dictionary of gradient settings group.
        :param method_dict:
            The input dicitonary of method settings group.
        """


        # if this is True, numerical must also be True
        if 'do_four_point' in grad_dict:
            key = grad_dict['do_four_point'].lower()
            self.do_four_point = True if key in ['yes', 'y'] else False
            # if four-point is desired, numerical is also set to True
            if self.do_four_point:
                self.numerical = True

        # Numerical gradient?
        if 'numerical' in grad_dict:
            key = grad_dict['numerical'].lower()
            self.numerical = True if key in ['yes', 'y'] else False

        # DFT
        if 'dft' in method_dict:
            key = method_dict['dft'].lower()
            self.dft = True if key in ['yes', 'y'] else False
        if 'grid_level' in method_dict:
            self.grid_level = int(method_dict['grid_level'])
        if 'xcfun' in method_dict:
            if 'dft' not in method_dict:
                self.dft = True
            self.xcfun = parse_xc_func(method_dict['xcfun'].upper())
            assert_msg_critical(not self.xcfun.is_undefined(),
                                'Gradient driver: Undefined XC functional')

        # TODO: Analytical DFT gradient is not implemented yet
        #if self.dft and not self.numerical:
        #    self.numerical = True
        #    warn_msg = '*** Warning: Analytical DFT gradient is not yet implemented.'
        #    self.ostream.print_blank()
        #    self.ostream.print_header(warn_msg.ljust(56))
        #    warn_msg = '              Gradient will be calculated numerically instead.'
        #    self.ostream.print_header(warn_msg.ljust(56))
        #    self.ostream.flush()

        # step size for finite differences
        if 'delta_h' in grad_dict:
            self.delta_h = float(grad_dict['delta_h'])

        # Numerical derivative of dipole moment
        if 'dipole_deriv' in grad_dict:
            key = grad_dict['dipole_deriv'].lower()
            self.dipole_deriv = True if key in ['yes', 'y'] else False
            if self.dipole_deriv and not self.numerical:
                self.ostream.print_blank()
                warn_msg = '*** Warning: Dipole moment derivatives requested.'
                self.ostream.print_header(warn_msg.ljust(56))
                warn_msg = '             Gradient will be calculated numerically.'
                self.ostream.print_header(warn_msg.ljust(56))
                self.numerical = True



    def compute(self, molecule, *args):
        """
        Performs calculation of numerical gradient.

        :param molecule:
            The molecule.
        :param args:
            The arguments.
        """

        return

    def compute_numerical(self, molecule, *args):
        """
        Performs calculation of numerical gradient.

        :param molecule:
            The molecule.
        :param args:
            The same arguments as the "compute" function.
        """

        # atom labels
        labels = molecule.get_labels()

        # atom coordinates (nx3)
        coords = molecule.get_coordinates()

        # numerical gradient
        self.gradient = np.zeros((molecule.number_of_atoms(), 3))

        for i in range(molecule.number_of_atoms()):
            for d in range(3):
                coords[i, d] += self.delta_h
                new_mol = Molecule(labels, coords, units='au')
                e_plus = self.compute_energy(new_mol, *args)

                coords[i, d] -= 2.0 * self.delta_h
                new_mol = Molecule(labels, coords, units='au')
                e_minus = self.compute_energy(new_mol, *args)

                if self.do_four_point:
                    coords[i, d] -= self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    e_minus2 = self.compute_energy(new_mol, *args)
                    coords[i, d] += 4.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    e_plus2 = self.compute_energy(new_mol, *args)
                    coords[i, d] -= 2.0 * self.delta_h
                    self.gradient[i, d] = ( (e_minus2 - 8 * e_minus 
                                           + 8 * e_plus - e_plus2) 
                                           / (12.0 * self.delta_h) )
                else:
                    coords[i, d] += self.delta_h
                    self.gradient[i, d] = (e_plus - e_minus) / (2.0 * self.delta_h)
                    
        # restore energy driver
        self.compute_energy(molecule, *args)

    def compute_energy(self, molecule, *args):
        """
        Computes the energy at current geometry.

        :param molecule:
            The molecule.
        :param args:
            The same arguments as the "compute" function.

        :return:
            The energy.
        """

        return

    def init_dft(self, molecule, scf_tensors):
        """
        Initializes DFT.

        :param molecule:
            The molecule.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.

        :return:
            The dictionary of DFT information.
        """

        # generate integration grid
        if self.dft:
            grid_drv = GridDriver(self.comm)
            grid_drv.set_level(self.grid_level)

            grid_t0 = tm.time()
            molgrid = grid_drv.generate(molecule)

            n_grid_points = molgrid.number_of_points()
            self.ostream.print_info(
                'Molecular grid with {0:d} points generated in {1:.2f} sec.'.
                format(n_grid_points,
                       tm.time() - grid_t0))
            self.ostream.print_blank()

            if self.rank == mpi_master():
                gs_density = AODensityMatrix([scf_tensors['D_alpha']],
                                             denmat.rest)
            else:
                gs_density = AODensityMatrix()
            gs_density.broadcast(self.rank, self.comm)
            molgrid.broadcast(self.rank, self.comm) # TODO duble check

            dft_func_label = self.xcfun.get_func_label().upper()
        else:
            molgrid = MolecularGrid()
            gs_density = AODensityMatrix()
            dft_func_label = 'HF'

        return {
            'molgrid': molgrid,
            'gs_density': gs_density,
            'dft_func_label': dft_func_label,
        }

    def grad_vxc_contrib(self, molecule, ao_basis, rhow_density, gs_density, xcfun_label):
        """
        Calculates the ground-state contribution of the exchange-correlation
        energy to the analytical gradient.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param rhow_density:
            The perturbed density.
        :param gs_density:
            The ground state density.
        :param xcfun_label:
            The label of the xc functional.

        :return:
            The ground-state exchange-correlation contribution to the gradient.
        """

        grid_drv = GridDriver(self.comm)
        grid_drv.set_level(self.grid_level)
        mol_grid = grid_drv.generate(molecule)
        mol_grid.distribute(self.rank, self.nodes, self.comm)

        xc_molgrad_drv = XCMolecularGradient(self.comm)
        atom_ids = list(range(molecule.number_of_atoms()))
        vxc_contrib = xc_molgrad_drv.integrate_vxc_gradient(
            atom_ids, rhow_density, gs_density, molecule, ao_basis, mol_grid,
            xcfun_label)
        vxc_contrib = self.comm.reduce(vxc_contrib, root=mpi_master())

        return vxc_contrib

    def grad_nuc_contrib(self, molecule):
        """
        Calculates the contribution of the nuclear-nuclear repulsion
        to the analytical nuclear gradient.

        :param molecule:
            The molecule.

        :return:
            The nuclear contribution to the gradient.
        """

        # number of atoms
        natm = molecule.number_of_atoms()

        # nuclear repulsion energy contribution to nuclear gradient
        nuc_contrib = np.zeros((natm, 3))

        # atom coordinates (nx3)
        coords = molecule.get_coordinates()

        # atomic charges
        nuclear_charges = molecule.elem_ids_to_numpy()

        # loop over all distinct atom pairs and add energy contribution
        for i in range(natm):
            z_a = nuclear_charges[i]
            r_a = coords[i]
            for j in range(i+1, natm):
                z_b = nuclear_charges[j]
                r_b = coords[j]
                r = np.sqrt(np.dot(r_a - r_b, r_a - r_b))
                f_ij = z_a * z_b * (r_b - r_a) / r**3
                nuc_contrib[i] += f_ij
                nuc_contrib[j] -= f_ij #z_a * z_b * (r_a - r_b) / r**3

        return nuc_contrib


    def get_gradient(self):
        """
        Gets the gradient.

        :return:
            The gradient.
        """

        return self.gradient

    def print_geometry(self, molecule):
        """
        Prints the geometry.

        :param molecule:
            The molecule.
        """

        self.ostream.print_block(molecule.get_string())


    def print_gradient(self, molecule):
        """
        Prints the gradient.

        :param molecule:
            The molecule.
        """

        # atom labels
        labels = molecule.get_labels()

        if self.numerical:
            title = 'Numerical '
        else:
            title = 'Analytical '

        title += 'Gradient (Hartree/Bohr)'
        self.ostream.print_header(title)
        self.ostream.print_header('-' * (len(title) + 2))
        self.ostream.print_blank()

        valstr = '  Atom '
        valstr += '{:>20s}  '.format('Gradient X')
        valstr += '{:>20s}  '.format('Gradient Y')
        valstr += '{:>20s}  '.format('Gradient Z')
        self.ostream.print_header(valstr)
        self.ostream.print_blank()

        for i in range(molecule.number_of_atoms()):
            valstr = '  {:<4s}'.format(labels[i])
            for d in range(3):
                valstr += '{:22.12f}'.format(self.gradient[i, d])
            self.ostream.print_header(valstr)

        self.ostream.print_blank()
        self.ostream.flush()


    def print_header(self, state_deriv_index=None):
        """
        Prints gradient calculation setup details to output stream.
        """

        str_width = 60

        self.ostream.print_blank()
        self.ostream.print_header(self.flag)
        self.ostream.print_header((len(self.flag) + 2) * '=')
        self.ostream.flush()

        cur_str = 'Gradient Type               : '

        if self.numerical:
            cur_str += 'Numerical'
            cur_str2 = 'Numerical Method            : '
            if self.do_four_point:
                cur_str2 += 'Five-Point Stencil'
            else:
                cur_str2 += 'Symmetric Difference Quotient'
            cur_str3 = 'Finite Difference Step Size : '
            cur_str3 += str(self.delta_h) + ' a.u.'
        else:
            cur_str += 'Analytical'

        if self.numerical or state_deriv_index is not None:
            self.ostream.print_blank()
            self.ostream.print_header(cur_str.ljust(str_width))
        if self.numerical:
            self.ostream.print_header(cur_str2.ljust(str_width))
            self.ostream.print_header(cur_str3.ljust(str_width))

        # print solver-specific info

        if state_deriv_index is not None:
            cur_str = 'Excited State of Interest   : ' + str(state_deriv_index +
                                                                 1)
            self.ostream.print_header(cur_str.ljust(str_width))

        self.ostream.print_blank()
        self.ostream.flush()

