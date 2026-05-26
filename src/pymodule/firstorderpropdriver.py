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
from .oneeints import compute_electric_dipole_integrals
from .errorhandler import assert_msg_critical


class FirstOrderPropertyDriver:
    """
    Implements first-order property driver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes SCF first-order property driver.
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

        self.property = 'electric dipole moment'

        self.origin = None

    def compute(self, molecule, basis, scf_results):
        """
        Computes first-order property.

        :param molecule:
            The molecule
        :param basis:
            The AO basis set.
        :param scf_results:
            The dictionary containing SCF results.
        """

        assert_msg_critical(
            self.property.lower() in [
                'electric dipole moment',
            ],
            f'{type(self).__name__}: Property {self.property} not yet supported'
        )

        if self.rank == mpi_master():
            total_density = scf_results['D_alpha'] + scf_results['D_beta']
        else:
            total_density = None

        if self.property.lower() in [
                'electric dipole moment',
        ]:
            dipole_dict = self.compute_electric_dipole_moment(
                molecule, basis, total_density)
        else:
            dipole_dict = None

        return {
            'electric dipole moment': dipole_dict,
        }

    def compute_electric_dipole_moment(self, molecule, basis, total_density):
        """
        Computes electric dipole moment.

        :param molecule:
            The molecule
        :param basis:
            The AO basis set.
        :param total_density:
            The total electron density.
        :return:
            The electric dipole moment.
        """

        # choose center of nuclear charges as origin
        coords = molecule.get_coordinates_in_bohr()
        nuclear_charges = molecule.get_effective_nuclear_charges(basis)

        if self.origin is None:
            origin = np.sum(coords.T * nuclear_charges,
                            axis=1) / np.sum(nuclear_charges)
        else:
            origin = np.array(self.origin)

        dipole_moment = None
        if self.rank == mpi_master():
            dipole_ints = compute_electric_dipole_integrals(
                molecule, basis, origin)

            # electronic contribution

            # multiple states:
            if len(total_density.shape) > 2:
                dof = total_density.shape[0]
                electronic_dipole = np.zeros((dof, 3))
                for i in range(dof):
                    electronic_dipole[i] = np.array([
                        np.sum(dipole_ints[d] * total_density[i])
                        for d in range(3)
                    ])

            # or single state:
            else:
                electronic_dipole = np.array(
                    [np.sum(dipole_ints[d] * total_density) for d in range(3)])

            # nuclear contribution
            coords = molecule.get_coordinates_in_bohr()
            nuclear_charges = molecule.get_effective_nuclear_charges(basis)
            nuclear_dipole = np.sum((coords - origin).T * nuclear_charges,
                                    axis=1)

            dipole_moment = nuclear_dipole + electronic_dipole

            if molecule.get_charge() != 0:
                warning = [
                    '*** Warning: Molecule has non-zero charge. Dipole',
                    '    moment will be dependent on the choice of origin.',
                    '    Center of nuclear charge is chosen as the origin.',
                ]
            else:
                warning = None

            dipole_dict = {
                'total': dipole_moment,
                'nuclear': nuclear_dipole,
                'electronic': electronic_dipole,
                'origin': origin,
                'units': 'a.u.',
                'warning': warning,
            }
        else:
            dipole_dict = None

        dipole_dict = self.comm.bcast(dipole_dict, root=mpi_master())

        return dipole_dict

    def print_property(self, prop_dict, title=None):
        """
        Prints first-order property.

        :param prop_dict:
            The dictionary from the compute method.
        :param title:
            The title of the printout.
        """

        assert_msg_critical(
            self.property.lower() in [
                'electric dipole moment',
            ],
            f'{type(self).__name__}: Property {self.property} not yet supported'
        )

        total = prop_dict['electric dipole moment']['total']

        if title is None:
            title = "Ground State Dipole Moment"

        self.ostream.print_blank()
        self.ostream.print_header(title)
        self.ostream.print_header('-' * (len(title) + 2))
        self.ostream.print_blank()

        warning = prop_dict['electric dipole moment'].get('warning')
        if warning is not None:
            for warn_msg in warning:
                self.ostream.print_header(warn_msg.ljust(56))
            self.ostream.print_blank()

        if total.ndim == 1:
            dip_au = list(total) + [np.linalg.norm(total)]
            dip_debye = [m * dipole_in_debye() for m in dip_au]

            for i, a in enumerate(['  X', '  Y', '  Z', 'Total']):
                valstr = '{:<5s} :'.format(a)
                valstr += '{:17.6f} a.u.'.format(dip_au[i])
                valstr += '{:17.6f} Debye   '.format(dip_debye[i])
                self.ostream.print_header(valstr)
        else:
            for index in range(total.shape[0]):
                self.ostream.print_blank()
                state_text = 'Excited State %d' % (index + 1)
                self.ostream.print_header(state_text)
                self.ostream.print_blank()
                dip_au = list(total[index]) + [np.linalg.norm(total[index])]
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

        new_prop_drv = FirstOrderPropertyDriver(self.comm, self.ostream)

        for key, val in vars(self).items():
            if isinstance(val, (MPI.Intracomm, OutputStream)):
                continue
            setattr(new_prop_drv, key, deepcopy(val))

        return new_prop_drv
