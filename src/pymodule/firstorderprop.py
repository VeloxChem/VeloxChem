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
import numpy as np
import sys

from .veloxchemlib import mpi_master
from .veloxchemlib import dipole_in_debye
from .outputstream import OutputStream
from .firstorderpropdriver import FirstOrderPropertyDriver


class FirstOrderProperties:
    """
    Implements first-order properties.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables:
        - properties: The dictionary of properties.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes SCF first-order properties.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        self.comm = comm
        self.rank = comm.Get_rank()

        self.ostream = ostream

        self.properties = {}

        self._prop_drv = FirstOrderPropertyDriver(self.comm, self.ostream)

    def _set_dipole_moment(self, dipole_moment):
        """
        Stores dipole moment under the canonical property key.

        :param dipole_moment:
            The dipole moment array for one or multiple states.
        """

        self.properties['dipole moment'] = deepcopy(dipole_moment)

    def _get_dipole_moment(self):
        """
        Gets dipole moment from the canonical property key.

        :return:
            The dipole moment array.
        """

        return deepcopy(self.properties['dipole moment'])

    def compute_scf_prop(self, molecule, basis, scf_results):
        """
        Computes first-order properties for SCF.

        :param molecule:
            The molecule
        :param basis:
            The AO basis set.
        :param scf_results:
            The dictionary containing SCF results.
        """

        if self.rank == mpi_master():
            total_density = scf_results['D_alpha'] + scf_results['D_beta']
        else:
            total_density = None

        self.compute(molecule, basis, total_density)

    def compute(self, molecule, basis, total_density):
        """
        Computes first-order properties.

        :param molecule:
            The molecule
        :param basis:
            The AO basis set.
        :param total_density:
            The total electron density.
        """

        dipole_dict = self._prop_drv.compute_electric_dipole_moment(
            molecule, basis, total_density)

        self._set_dipole_moment(dipole_dict['total'])

    def get_property(self, key):
        """
        Gets first-order property.

        :param key:
            The name of the property.

        :return:
            The property.
        """

        dipole_keys = [
            'dipole moment',
            'dipole_moment',
            'electric dipole moment',
            'electric_dipole_moment',
        ]

        if key in dipole_keys:
            return self._get_dipole_moment()

        return self.properties[key]

    def print_properties(self, molecule, title=None, states=None):
        """
        Prints first-order properties.

        :param molecule:
            The molecule.
        :param title:
            The title of the printout, giving information about the method,
            relaxed/unrelaxed density, which excited state.
        :param states:
            The excited states for which first order properties ar calculated.
        """

        if title is None:
            title = "Ground State Dipole Moment"
        self.ostream.print_blank()

        self.ostream.print_header(title)
        self.ostream.print_header('-' * (len(title) + 2))

        # prints warning if the molecule is charged
        if molecule.get_charge() != 0:
            warn_msg = '*** Warning: Molecule has non-zero charge. Dipole'
            self.ostream.print_header(warn_msg.ljust(56))
            warn_msg = '    moment will be dependent on the choice of origin.'
            self.ostream.print_header(warn_msg.ljust(56))
            warn_msg = '    Center of nuclear charge is chosen as the origin.'
            self.ostream.print_header(warn_msg.ljust(56))

        self.ostream.print_blank()

        dip = self._get_dipole_moment()
        if states is None:
            dip_au = list(dip) + [np.linalg.norm(dip)]
            dip_debye = [m * dipole_in_debye() for m in dip_au]

            for i, a in enumerate(['  X', '  Y', '  Z', 'Total']):
                valstr = '{:<5s} :'.format(a)
                valstr += '{:17.6f} a.u.'.format(dip_au[i])
                valstr += '{:17.6f} Debye   '.format(dip_debye[i])
                self.ostream.print_header(valstr)
        else:
            for index, s in enumerate(states):
                self.ostream.print_blank()
                state_text = 'Excited State %d' % s
                self.ostream.print_header(state_text)
                self.ostream.print_blank()
                dip_au = list(dip[index]) + [np.linalg.norm(dip[index])]
                dip_debye = [m * dipole_in_debye() for m in dip_au]

                for i, a in enumerate(['  X', '  Y', '  Z', 'Total']):
                    valstr = '{:<5s} :'.format(a)
                    valstr += '{:17.6f} a.u.'.format(dip_au[i])
                    valstr += '{:17.6f} Debye   '.format(dip_debye[i])
                    self.ostream.print_header(valstr)

        self.ostream.print_blank()
        self.ostream.flush()

    def __deepcopy__(self, memo):
        """
        Implements deepcopy.

        :param memo:
            The memo dictionary for deepcopy.

        :return:
            A deepcopy of self.
        """

        new_scf_prop = FirstOrderProperties(self.comm, self.ostream)

        for key, val in vars(self).items():
            if isinstance(val, (MPI.Intracomm, OutputStream)):
                continue
            setattr(new_scf_prop, key, deepcopy(val))

        return new_scf_prop
