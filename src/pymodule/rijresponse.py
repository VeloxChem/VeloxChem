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


"""The resolution of the identity for the Coulomb matrices of a response
calculation (RI-J), shared by the linear solver and the nonlinear one.

The integrals are held in memory, cut into parts which the ranks divide. The
setup predicts the most a rank reaches and falls back to the four-center
integrals where that is more than the budget. A build takes its densities in
chunks, so that what it holds for them beside the integrals stays within the
room the prediction left.
"""

from mpi4py import MPI
import numpy as np

from .veloxchemlib import mpi_master
from .veloxchemlib import PackedMatrix
from .veloxchemlib import SimdRIJFockDriver
from .veloxchemlib import mat_t
from .veloxchemlib import rimode
from .molecularbasis import MolecularBasis
from .profiler import Profiler
from . import rijkresponse


def initialize(solver, molecule, basis):
    """
    Prepares the simd Coulomb only driver the solver builds with, or falls back
    to the four-center integrals where it does not fit. Every rank must call
    this.

    :param solver:
        The solver, which is given the driver, or _ri_j_fallback where it does
        not fit.
    :param molecule:
        The molecule.
    :param basis:
        The AO basis set.
    """

    solver._ri_j_fallback = False
    solver._ri_j_memory = None
    solver._ri_drv = None

    rijkresponse.mark_setup(solver)

    basis_ri = MolecularBasis.read(molecule,
                                   solver.ri_auxiliary_basis,
                                   ostream=None,
                                   auxiliary_basis=True)

    ri_drv = SimdRIJFockDriver()

    # NOTE: the parts are cut from the budget, and every rank has to cut the
    # same ones or a rank would sweep a part no other rank knows of. The budget
    # is therefore the smallest of the ranks' and not each one's own.
    budget = min(
        solver.comm.allgather(
            rijkresponse.memory_budget(solver, molecule, basis)))

    # NOTE: the integrals are held in memory, the way a response build wants:
    # it sweeps them for many densities a call. A calculation whose peak in
    # that way exceeds the budget falls back to the four-center integrals, as
    # the SCF does, rather than taking the direct way the automatic choice
    # would pick without a check of the peak. ri, ri_coulomb and ri_jk are left
    # as they were set.
    peak = ri_drv.peak_memory(molecule, basis, basis_ri, solver.eri_thresh,
                              budget, solver.rank, solver.nodes)
    peaks = solver.comm.allgather(peak)

    where = ('available' if solver.nodes == 1 else 'on the largest rank')

    # NOTE: an SCF asked for force_ri hands it on through its results.
    if (max(peaks) > budget) and getattr(solver, '_force_ri', False):
        solver.ostream.print_warning(
            'RI-J would need roughly ' +
            f'{Profiler.memory_to_string(max(peaks))} of ' +
            f'{Profiler.memory_to_string(budget)} {where}; force_ri is set ' +
            'on the SCF, so the response uses RI-J.')
        solver.ostream.print_blank()
        solver.ostream.flush()

    elif max(peaks) > budget:
        solver._ri_j_fallback = True

        solver.ostream.print_warning(
            'RI-J would need roughly ' +
            f'{Profiler.memory_to_string(max(peaks))} of ' +
            f'{Profiler.memory_to_string(budget)} {where}, so the response ' +
            'Coulomb matrices are formed from the four-center integrals.')
        solver.ostream.print_blank()
        solver.ostream.flush()
        return

    # NOTE: the metric is inverted once on the master and handed to the ranks,
    # as the Fock build does it: the inversion picks its own fallback from the
    # matrix, and two ranks choosing differently would build with metrics which
    # are not the same.
    metric = PackedMatrix()

    if solver.rank == mpi_master():
        metric = ri_drv.make_metric(molecule, basis_ri,
                                    solver.ri_metric_threshold)

    metric = metric.broadcast(solver.comm, root=mpi_master())

    # NOTE: the parts are what the ranks divide, so there have to be at least
    # as many of them as there are ranks.
    ri_drv.prepare(molecule, basis, basis_ri, solver.eri_thresh, budget,
                   solver.ri_metric_threshold, rimode.in_memory, metric,
                   solver.rank, solver.nodes)

    del metric

    solver._ri_drv = ri_drv

    # NOTE: the same on every rank, as the chunks of a build meet in a
    # reduction each and every rank has to cut the same ones.
    solver._ri_j_memory = (budget, max(peaks),
                           basis_ri.get_dimensions_of_basis())

    solver.ostream.print_info(
        'Using the SIMD resolution of the identity (RI-J) for response.')
    solver.ostream.print_info(
        'RI-J needs roughly ' + f'{Profiler.memory_to_string(max(peaks))} of ' +
        f'{Profiler.memory_to_string(budget)} {where}; three-center ' +
        f'integrals are held in memory, in {ri_drv.number_of_parts()} parts.')
    solver.ostream.print_blank()
    solver.ostream.flush()


def _max_densities(solver, nao, nspins):
    """
    The most densities one chunk of a build may take.

    :param solver:
        The solver holding the driver and the memory it was checked with.
    :param nao:
        The number of basis functions.
    :param nspins:
        One for a restricted reference, two for an unrestricted one.

    :return:
        The number of densities, at least one, the same on every rank.
    """

    if solver._ri_j_memory is None:
        return None

    budget, peak, naux = solver._ri_j_memory

    # NOTE: what a chunk holds for each density beside what the driver's peak
    # counts for one: the density packed for the fitting and its fitting
    # coefficients, local and reduced, and for an open shell the two spins'
    # densities added. The Fock matrices returned are the build's output, which
    # the four-center build holds as well.
    triangle = nao * (nao + 1) // 2 * 8
    per_density = triangle + 2 * naux * 8 + (nspins - 1) * nao * nao * 8

    return max((budget - peak) // per_density, 1)


def coulomb_matrices(solver, dens, comm):
    """
    Computes this rank's share of the Coulomb matrices of a batch of densities.

    :param solver:
        The solver holding the driver.
    :param dens:
        The densities, as numpy arrays.
    :param comm:
        The communicator the build is divided over.

    :return:
        This rank's share of each Coulomb matrix, doubled as every matrix of
        the `j` kind is, as numpy arrays. The shares are added by the reduction
        the caller makes at the end of the build.
    """

    if not dens:
        return []

    nao = dens[0].shape[0]

    count = _max_densities(solver, nao, 1)

    fock_arrays = []

    for start, end in rijkresponse._chunks(len(dens), count):
        # NOTE: twice, as every other way of building a matrix of the `j`
        # kind is doubled by its caller; this path does not go through the
        # loop which doubles what the conventional driver answers.
        fock_arrays += [
            2.0 * coulomb
            for coulomb in _coulomb_of_chunk(solver, dens[start:end], comm)
        ]

    return fock_arrays


def coulomb_matrices_unrestricted(solver, dens_a, dens_b, comm):
    """
    Computes this rank's share of the Coulomb matrices of an open shell batch.

    :param solver:
        The solver holding the driver.
    :param dens_a:
        The alpha densities, as numpy arrays.
    :param dens_b:
        The beta densities, as numpy arrays.
    :param comm:
        The communicator the build is divided over.

    :return:
        This rank's share of the Coulomb matrices, an alpha and a beta for each
        density in the order the caller appends them.
    """

    if not dens_a:
        return []

    nao = dens_a[0].shape[0]

    count = _max_densities(solver, nao, 2)

    fock_arrays = []

    for start, end in rijkresponse._chunks(len(dens_a), count):
        # NOTE: the Coulomb of an open shell is of the two densities added and
        # is the same matrix for both spins, so one is formed and handed to
        # each. There is no factor of two here, where the closed shell carries
        # one: the total density already counts both spins.
        totals = [
            a + b for a, b in zip(dens_a[start:end], dens_b[start:end])
        ]

        for coulomb in _coulomb_of_chunk(solver, totals, comm):
            # NOTE: a copy each, not the same array twice. What the caller
            # does to one spin's matrix afterwards must not reach the other's.
            fock_arrays.append(coulomb.copy())
            fock_arrays.append(coulomb)

        del totals

    return fock_arrays


def _coulomb_of_chunk(solver, dens, comm):
    """
    Computes this rank's share of the Coulomb matrices of one chunk.

    :param solver:
        The solver holding the driver.
    :param dens:
        The densities of the chunk, as numpy arrays.
    :param comm:
        The communicator the build is divided over.

    :return:
        This rank's share of each Coulomb matrix, not doubled.
    """

    nao = dens[0].shape[0]

    ri_drv = solver._ri_drv

    mine = ri_drv.owned_parts()

    # NOTE: the density is symmetrized, as it is for the conventional driver of
    # this approximation: a response density is not symmetric and the fitting
    # sees only the symmetric part of it.
    packed = []

    for density in dens:
        matrix = PackedMatrix(nao, nao, mat_t.symmetric)
        matrix.from_numpy(np.ascontiguousarray(0.5 * (density + density.T)))
        packed.append(matrix)

    # NOTE: the right hand side of every fitting of the chunk is gathered into
    # one array and reduced once, a single collective for the chunk where a
    # density at a time would be one each.
    local = np.array([ri_drv.compute_gamma(matrix, mine) for matrix in packed],
                     dtype=np.float64)

    del packed

    total = np.zeros_like(local)
    comm.Allreduce(local, total, op=MPI.SUM)

    del local

    coulombs = []

    for idx in range(total.shape[0]):
        gamma = ri_drv.solve_fitting(list(total[idx]))

        fock = PackedMatrix(nao, nao, mat_t.symmetric)
        fock.zero()

        ri_drv.compute_coulomb(gamma, mine, fock)

        coulombs.append(fock.to_numpy())

    return coulombs
