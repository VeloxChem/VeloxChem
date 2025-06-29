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

from .veloxchemlib import AODensityMatrix
from .veloxchemlib import denmat
from .veloxchemlib import mpi_master


def _AODensityMatrix_broadcast(self, comm, root=mpi_master()):
    """
    Broadcasts AODensityMatrix object.

    :param comm:
        The MPI communicator.
    :param root:
        The root rank to broadcast from.

    :return:
        The AODensityMatrix object.
    """

    densities = []

    if comm.Get_rank() == root:
        density_type = ('rest'
                        if self.get_density_type() == denmat.rest else 'unrest')
        num_dens = self.number_of_density_matrices()
    else:
        density_type = None
        num_dens = None

    density_type = comm.bcast(density_type, root=root)
    num_dens = comm.bcast(num_dens, root=root)

    if comm.Get_rank() == root:
        for i in range(num_dens):
            densities.append(self.alpha_to_numpy(i))
            if density_type == 'unrest':
                densities.append(self.beta_to_numpy(i))

    else:
        for i in range(num_dens):
            densities.append(None)
            if density_type == 'unrest':
                densities.append(None)

    if density_type == 'rest':
        for i in range(num_dens):
            densities[i] = comm.bcast(densities[i], root=root)
    elif density_type == 'unrest':
        for i in range(num_dens * 2):
            densities[i] = comm.bcast(densities[i], root=root)

    density_type = (denmat.rest if density_type == 'rest' else denmat.unrest)

    return AODensityMatrix(densities, density_type)


AODensityMatrix.broadcast = _AODensityMatrix_broadcast
