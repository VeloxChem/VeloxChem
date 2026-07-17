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

from .veloxchemlib import mpi_master
from .gradientdriver import GradientDriver


class OpenQPGradientDriver(GradientDriver):
    """
    Implements the OpenQP ground-state (SCF) gradient driver.

    :param openqp_drv:
        The OpenQP SCF driver.

    Instance variables
        - flag: The driver flag.
    """

    def __init__(self, openqp_drv):
        """
        Initializes the OpenQP gradient driver.
        """

        super().__init__(openqp_drv.comm, openqp_drv.ostream)

        self.openqp_driver = openqp_drv
        self.flag = 'OpenQP Gradient Driver'

    def compute(self, molecule, basis=None, scf_results=None):
        """
        Performs the calculation of the OpenQP ground-state gradient.

        :param molecule:
            The molecule.
        :param basis:
            Unused; present for interface compatibility.
        :param scf_results:
            Unused; present for interface compatibility.
        """

        if self.numerical:
            if self.rank == mpi_master():
                self.compute_numerical(molecule)
            else:
                self.gradient = None
        else:
            self.openqp_driver.compute_gradient(molecule)
            self.gradient = self.openqp_driver.get_gradient()

        self.gradient = self.comm.bcast(self.gradient, root=mpi_master())

        self.print_geometry(molecule)
        self.print_gradient(molecule)

        self.openqp_driver._invalidate_cache()

        self.ostream.print_blank()
        self.ostream.flush()

    def compute_energy(self, molecule, *args):
        """
        Computes the OpenQP ground-state (SCF) energy at the current geometry.

        :param molecule:
            The molecule.

        :return:
            The SCF energy.
        """

        self.openqp_driver.compute(molecule)

        return self.openqp_driver.get_energy()
