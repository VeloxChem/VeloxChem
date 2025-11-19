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

import time as tm

from .dftutils import get_default_grid_level
from .errorhandler import assert_msg_critical
from .hessiandriver import HessianDriver
from .tdaeigensolver import TdaEigenSolver
from .tddftgradientdriver import TddftGradientDriver
from .veloxchemlib import mpi_master


class TddftHessianDriver(HessianDriver):
    """
    Implements the TDDFT Hessian driver which is used
    to determine the Hessian for a specific excited state.

    :param scf_driver:
        The SCF driver.
    :param rsp_driver:
        The response driver.

    Instance variables
        - hessian: The Hessian in Hartree per Bohr**2.
        - state_deriv_index: The index of the excited state of interest.
    """

    def __init__(self, scf_drv, rsp_drv):
        """
        Initializes the TDDFT Hessian driver.

        :param scf_drv:
            The SCF driver.
        :param rsp_drv:
            The linear response driver.
        """
        super().__init__(scf_drv.comm, scf_drv.ostream)

        self._xcfun_ldstaging = scf_drv._xcfun_ldstaging

        if scf_drv._dft:
            self.flag = "TDDFT Hessian Driver"
        else:
            self.flag = "TDHF Hessian Driver"
        self.scf_driver = scf_drv
        self.rsp_driver = rsp_drv
        self.do_print_hessian = False

        # excited state information; if set to None,
        # the first excited state will be calculated.
        self.state_deriv_index = None

        # excited-state relaxed dipole moment for IR
        self.relaxed_dipole_moment = None

        # option dictionaries from input
        # TODO: cleanup
        self.cphf_dict = {}
        self.grad_dict = {}

        if isinstance(rsp_drv, TdaEigenSolver):
            self.tamm_dancoff = True
        else:
            self.tamm_dancoff = False

        self._input_keywords["hessian"].update(
            {
                "state_deriv_index": ("int", "excited state index"),
            }
        )

    def update_settings(
        self, method_dict, hessian_dict=None, grad_dict=None, cphf_dict=None
    ):
        """
        Updates settings in TddftHessianDriver.

        :param method_dict:
            The input dictionary of method settings group.
        :param hessian_dict:
            The input dictionary of Hessian settings group.
        :param grad_dict:
            The input dictionary for TDDFT gradient settings.
        :param cphf_dict:
            The input dictionary of orbital response settings.
        """
        if hessian_dict is None:
            hessian_dict = {}

        super().update_settings(method_dict, hessian_dict)

        if grad_dict is None:
            grad_dict = {}

        if cphf_dict is None:
            cphf_dict = {}

        self.grad_dict = dict(grad_dict)
        self.cphf_dict = dict(cphf_dict)

    def compute(self, molecule, basis):
        """
        Performs the computation of the TDDFT Hessian.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        """

        if self.rank == mpi_master():
            self.print_header()
            if self.state_deriv_index is None:
                self.state_deriv_index = 1

        self.state_deriv_index = self.comm.bcast(
            self.state_deriv_index, root=mpi_master()
        )

        start_time = tm.time()

        # compute hessian
        self.compute_numerical(molecule, basis)

        if self.rank == mpi_master():
            if self.do_print_hessian:
                self.print_geometry(molecule)
                self.ostream.print_blank()
                title = "Numerical Hessian based on analytical "
                title += "gradient (Hartree/Bohr**2)"
                self.print_hessian(molecule, title=title)

            valstr = "*** Time spent in Hessian calculation: "
            valstr += "{:.2f} sec ***".format(tm.time() - start_time)
            self.ostream.print_header(valstr)
            self.ostream.print_blank()
            self.ostream.print_blank()
            self.ostream.flush()

    def compute_energy(self, molecule, basis):
        """
        Computes the TDDFT energy.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        """
        self.scf_driver.restart = False
        self.scf_driver.ostream.mute()
        scf_results = self.scf_driver.compute(molecule, basis)
        self.scf_driver.ostream.unmute()
        assert_msg_critical(
            self.scf_driver.is_converged, "TddftHessianDriver: SCF did not converge"
        )

        self.rsp_driver.restart = False
        self.rsp_driver.ostream.mute()
        rsp_results = self.rsp_driver.compute(molecule, basis, scf_results)
        self.rsp_driver.ostream.unmute()
        assert_msg_critical(
            self.rsp_driver.is_converged,
            "TddftHessianDriver: response did not converge",
        )

        if self.rank == mpi_master():
            index = self.state_deriv_index - 1
            scf_energy = self.scf_driver.get_scf_energy()
            excitation_energy = rsp_results["eigenvalues"][index]

            energy = scf_energy + excitation_energy

            return energy
        else:
            return None

    def compute_gradient(self, molecule, basis):
        """
        Computes the TDDFT gradient at the current molecular geometry.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        """
        # Recalculate the ground state energy for the new geometry
        self.scf_driver.restart = False
        self.scf_driver.ostream.mute()
        scf_results = self.scf_driver.compute(molecule, basis)
        self.scf_driver.ostream.unmute()
        assert_msg_critical(
            self.scf_driver.is_converged, "TddftHessianDriver: SCF did not converge"
        )

        # Recalculate the excited state energies for the new geometry
        self.rsp_driver.restart = False
        self.rsp_driver.ostream.mute()
        rsp_results = self.rsp_driver.compute(molecule, basis, scf_results)
        self.rsp_driver.ostream.unmute()
        assert_msg_critical(
            self.rsp_driver.is_converged,
            "TddftHessianDriver: response did not converge",
        )

        # Calculate the excited state gradient.
        tddft_gradient_driver = TddftGradientDriver(self.scf_driver)
        tddft_gradient_driver.update_settings(
            grad_dict=self.grad_dict, orbrsp_dict=self.cphf_dict
        )

        tddft_gradient_driver.state_deriv_index = self.state_deriv_index

        if self.rank == mpi_master():
            # Calculate excited-state relaxed dipole moment
            if self.do_dipole_gradient:
                tddft_gradient_driver.do_first_order_prop = True

        tddft_gradient_driver.ostream.mute()
        tddft_gradient_driver.compute(
            molecule, basis, self.scf_driver, self.rsp_driver, rsp_results
        )
        tddft_gradient_driver.ostream.unmute()

        if self.rank == mpi_master():
            if self.do_dipole_gradient:
                self.relaxed_dipole_moment = (
                    tddft_gradient_driver.relaxed_dipole_moment[0].copy()
                )

        return tddft_gradient_driver.gradient[0].copy()

    # The relaxed dipole moment is calculated at the same time as
    # the gradient so, to avoid re-calculating everything, here the value
    # is simply returned. Used by compute_numerical in HessianDriver.
    def compute_electric_dipole_moment(self, molecule, basis):
        """
        Returns the electric dipole moment calculated with the
        TDDFT gradient in compute_gradient.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        """

        if self.rank == mpi_master():
            return self.relaxed_dipole_moment
        else:
            return None

    def print_header(self):
        """
        Prints Hessian calculation setup details to output stream.
        """

        str_width = 70

        self.ostream.print_blank()
        self.ostream.print_header(self.flag)
        self.ostream.print_header((len(self.flag) + 2) * "=")
        self.ostream.flush()

        cur_str = "Hessian Type                    : "
        cur_str += "Numerical using analytical gradient"
        cur_str2 = "Numerical Method                : "
        if self.do_four_point:
            cur_str2 += "Five-Point Stencil"
        else:
            cur_str2 += "Symmetric Difference Quotient"
        cur_str3 = "Finite Difference Step Size     : "
        cur_str3 += str(self.delta_h) + " a.u."
        if self.state_deriv_index is None:
            s = 1
        else:
            s = self.state_deriv_index

        cur_str4 = "Exited State                    : %d" % s

        self.ostream.print_blank()
        self.ostream.print_header(cur_str.ljust(str_width))
        self.ostream.print_header(cur_str2.ljust(str_width))
        self.ostream.print_header(cur_str3.ljust(str_width))
        self.ostream.print_header(cur_str4.ljust(str_width))

        if self.scf_driver._dft:
            cur_str = "Exchange-Correlation Functional : "
            cur_str += self.scf_driver.xcfun.get_func_label().upper()
            self.ostream.print_header(cur_str.ljust(str_width))
            grid_level = (
                get_default_grid_level(self.scf_driver.xcfun)
                if self.scf_driver.grid_level is None
                else self.scf_driver.grid_level
            )
            cur_str = "Molecular Grid Level            : " + str(grid_level)
            self.ostream.print_header(cur_str.ljust(str_width))

        self.ostream.print_blank()
        self.ostream.flush()
