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

"""The resolution of the identity for the Fock matrices of a response calculation.

The linear solver and the nonlinear one are separate classes which do not share a
base, and both of them need the same three things: the B vectors formed once, the
memory they are allowed, and a Fock matrix built from densities given as their
factors. They are written here as functions of the solver rather than twice as
methods of it.

A density is given as the factors it was made from, and three shapes are taken:

    (left, rights)                     D = left rights[k]^T
    (left, rights, transposed_rights)  D = left rights[k]^T + transposed[k] left^T
    (left_a, rights_a, left_b, rights_b)
                                       D = left_a rights_a[k]^T + left_b rights_b[k]^T

The last of them carries **any** density and not only a block diagonal one. Split the
left index of the molecular orbital density into its occupied and virtual halves and

    r_a = C(occ) M(oo)^T + C(vir) M(ov)^T
    r_b = C(occ) M(vo)^T + C(vir) M(vv)^T

give it exactly, whatever blocks of M are nonzero. A batch which mixes the orders,
as a three-time perturbed calculation does, therefore needs no telling apart of its
densities: one set of factors with C(occ) and C(vir) as the shared halves carries all
of them. It costs the basis times the orbitals where a density which is only between
the occupied orbitals and the virtual ones would cost the basis times twice the
occupied, so it is the general form and not the cheapest one.

The first is a trial vector of the Tamm-Dancoff approximation, the second one of
linear response or a density of third order in the perturbation, and the third a
density of second order, which is block diagonal in the orbitals: the occupied block
carries the occupied orbitals on both sides and the virtual block the virtual ones.
The third is two calls of the first added together, the Fock matrix being linear in
the density.

A three-time perturbed calculation has densities of both the second and the third
order in one batch, cut from two arrays with strides of their own. That is not a
shape of factors but two of them, and the caller holds them apart and calls this
twice rather than this module learning about batches it cannot see.
"""

from contextlib import contextmanager

import numpy as np
import psutil

from .veloxchemlib import SimdRIJKFockDriver
from .veloxchemlib import SimdRIJKResponseDriver
from .veloxchemlib import PackedMatrix
from .veloxchemlib import mat_t
from .veloxchemlib import rimode
from .molecularbasis import MolecularBasis
from .environment import get_ri_free_memory
from .errorhandler import assert_msg_critical
from .profiler import Profiler


def memory_budget(solver, molecule, basis):
    """
    Gets the memory the simd RI driver of a response solver may hold, per rank,
    in bytes. Every rank must call this. RI-J and RI-JK both take it.

    :param solver:
        The solver.
    :param molecule:
        The molecule, whose open shell doubles the rank reserve.
    :param basis:
        The AO basis set.

    :return:
        The memory budget in bytes: ri_memory_budget where the solver sets it
        (in GB, taken as it is), and otherwise the memory free to this rank less
        the reserves of the node and of the rank. The rank's reserve is kept on
        the solver, for memory profiling to report.
    """

    ri_memory_budget = getattr(solver, 'ri_memory_budget', None)

    if ri_memory_budget is not None:
        solver._ri_memory_reserve = 0
        return int(ri_memory_budget * 1.0e9)

    # NOTE: the rank's reserve is taken from what is free to the rank after the
    # node's reserve, so it is worked out here and not in get_ri_free_memory,
    # which the SCF and its gradient share with a fixed reserve.
    free = get_ri_free_memory(solver.comm, 0)

    reserve = rank_reserve(molecule, basis, free)

    solver._ri_memory_reserve = reserve

    return max(free - reserve, 0)


def rank_reserve(molecule, basis, free):
    """
    Gets the memory a rank keeps back from the simd RI driver of a response
    solver, for what the response allocates after the budget is taken, in bytes.

    :param molecule:
        The molecule.
    :param basis:
        The AO basis set.
    :param free:
        The memory free to this rank after the node's reserve, in bytes.

    :return:
        The reserve in bytes.
    """

    # NOTE: part of what the response allocates after the budget is taken goes
    # as the square of the basis: the quadrature of the functional and its
    # ground state density, the densities of a build broadcast to every rank
    # and the Fock matrices formed from them, and their reduction. Each spin of
    # an open shell carries its own.
    nao = basis.get_dimensions_of_basis()
    nspins = (1 if molecule.number_of_alpha_electrons()
              == molecule.number_of_beta_electrons() else 2)

    matrices = _RANK_RESERVE_MATRICES * nspins * nao * nao * 8

    # NOTE: the larger part does not go as the square of the basis: the
    # solver's reduced space (trial and sigma vectors, each of the size of the
    # occupied times the virtual orbitals, growing with the states or
    # frequencies and with the iterations) and a nonlinear driver's batch of
    # second-order densities. Neither is known when the budget is taken, so a
    # share of the free memory is kept for them, which grows with the machine
    # as a calculation asking for more of them does. For noradrenaline in
    # def2-TZVPD (552 basis functions), B3LYP, on one rank with 12 threads,
    # this part was about 0.4 GB for 5 RPA states, 1.2 GB for 40 complex
    # polarizability frequencies and 2.1 GB for 50 RPA states.
    share = free // _RANK_RESERVE_SHARE_DIVISOR

    return max(matrices, share, 512 * 1024**2)


# NOTE: the matrices of the basis squared the rank reserve keeps for each spin;
# see rank_reserve.
_RANK_RESERVE_MATRICES = 12

# NOTE: the rank reserve keeps at least this fraction (its inverse) of what is
# free to the rank; see rank_reserve.
_RANK_RESERVE_SHARE_DIVISOR = 4


def mark_setup(solver):
    """
    Records what a rank holds before the simd RI driver of a response solver is
    set up, against which memory profiling holds its peak after the response.
    With memory profiling, the peak is reset here, so that what follows is not
    read against the SCF's peak.

    :param solver:
        The solver.
    """

    # NOTE: what the solver holds for the densities of a batch, recorded by
    # the builds after this setup (LinearSolver._note_ri_batch_memory).
    solver._ri_batch_memory = 0

    if not getattr(solver, 'memory_profiling', False):
        return

    Profiler.reset_peak_memory()

    solver._ri_rss_at_setup = psutil.Process().memory_info().rss


def print_memory_check(solver):
    """
    Prints, with memory profiling, the most the simd RI driver of a response
    solver was predicted to reach against what each rank's peak grew by after
    it was set up. Every rank must call this.

    :param solver:
        The solver.
    """

    if (not getattr(solver, 'memory_profiling', False) or
            getattr(solver, '_ri_rss_at_setup', None) is None):
        return

    if solver._uses_ri_jk() and solver._ri_jk_memory is not None:
        label = 'RI-JK'
        predicted = getattr(solver, '_ri_jk_predicted_peak', 0)
    elif solver._uses_ri_coulomb() and solver._ri_j_memory is not None:
        label = 'RI-J'
        predicted = solver._ri_j_memory[1]
    else:
        return

    predicted = max(solver.comm.allgather(predicted))

    # NOTE: recorded on the master rank only, which holds the most for each
    # density of a batch.
    batch_memory = max(
        solver.comm.allgather(getattr(solver, '_ri_batch_memory', 0)))

    # NOTE: the growth of the peak and not the peak itself, as in the SCF: what
    # a rank held before the setup was not free when the budget was taken.
    grown = Profiler.get_peak_memory() - solver._ri_rss_at_setup
    grown_all = solver.comm.allgather(grown)
    largest = max(grown_all)

    where = ('' if solver.nodes == 1 else
             f' on rank {grown_all.index(largest):d}, the largest of ' +
             f'{solver.nodes:d}')

    solver.ostream.print_info(
        f'{label} was predicted to need roughly ' +
        f'{Profiler.memory_to_string(predicted)} at its largest build, and ' +
        f'the solver {Profiler.memory_to_string(batch_memory)} for the ' +
        'densities of its largest batch, with ' +
        f'{Profiler.memory_to_string(solver._ri_memory_reserve)} reserved ' +
        'for the rest of the response.')
    solver.ostream.print_info('Peak memory grew by ' +
                              f'{Profiler.memory_to_string(largest)} after ' +
                              f'the {label} setup{where}.')
    solver.ostream.print_blank()
    solver.ostream.flush()


def range_separation(solver):
    """The range separation parameter the solver's functional asks for.

    :param solver:
        The solver.

    :return:
        The parameter, or zero where the functional is not range separated and the
        plain B vectors are the whole of what is needed.

    .. note::
        This is read from the functional and not from the parameters of a build.
        Those are settled per call and say what that call wants; what the B vectors
        are is settled once, before any of them.
    """

    if not getattr(solver, '_dft', False):
        return 0.0

    xcfun = getattr(solver, 'xcfun', None)

    if xcfun is None or not xcfun.is_range_separated():
        return 0.0

    return xcfun.get_rs_omega()


def initialize(solver, molecule, basis, reuse=False):
    """
    Forms the B vectors the response driver contracts, on the solver, or checks
    that the solver holds the ones it is to use again.

    :param solver:
        The solver, which is given the drivers and the auxiliary basis.
    :param molecule:
        The molecule.
    :param basis:
        The AO basis set.
    :param reuse:
        True where the solver formed the B vectors, or fell back, earlier in
        the same calculation with the same molecule and basis, and uses that
        again: nothing is formed, and an error is raised if there is nothing to
        use. Callers pass it for the second setup of one calculation, which is
        a nonlinear driver's after it formed the B vectors to lend them (see
        lend) and an orbital response's omega after its compute.
    """

    _check_not_lent(solver)

    # NOTE: a solver which borrows the B vectors of the driver that drives it
    # (see lend) uses them for the call they were lent for and forms none of
    # its own. The driver formed them for its own functional; every driver
    # copies its xcfun to a solver among the keywords it sets before lending,
    # and the check below fails if one ever does not.
    if solver._ri_jk_borrowing:
        assert_msg_critical(
            solver._ri_jk_fallback or
            solver._ri_jk_drv.get_omega() == range_separation(solver),
            f'{type(solver).__name__}: the RI-JK B vectors it borrows are ' +
            'not formed for the range separation of its functional')
        return

    if reuse:
        assert_msg_critical(
            solver._ri_jk_fallback or solver._ri_jk_drv is not None,
            f'{type(solver).__name__}: no RI-JK B vectors were formed ' +
            'earlier in this calculation')
        return

    # NOTE: B vectors held from an earlier calculation are let go before new
    # ones are formed. They belong to the molecule and basis of that
    # calculation, which may not be this one (a solver computed again at
    # another geometry), and let go here they are freed before the new ones
    # are allocated rather than after. The fallback is decided again as well.
    _release(solver)

    mark_setup(solver)

    assert_msg_critical(
        'jkfit' in solver.ri_auxiliary_basis.lower(),
        f'{type(solver).__name__}: RI-JK needs a jkfit auxiliary basis, and ' +
        f'{solver.ri_auxiliary_basis} fits the Coulomb alone')

    solver._ri_jk_aux_basis = MolecularBasis.read(molecule,
                                                  solver.ri_auxiliary_basis,
                                                  ostream=None,
                                                  auxiliary_basis=True)

    solver._ri_jk_drv = SimdRIJKFockDriver()

    # NOTE: a hybrid range-separated functional needs the B vectors of the
    # attenuated operator beside the plain ones, each fitted in its own metric, so
    # it holds twice this and the check below is made against twice this.
    omega = range_separation(solver)

    # NOTE: the atoms of the auxiliary basis this rank holds. Every term of the
    # Coulomb and of the exchange is a sum over the auxiliary basis, so a rank given
    # a share of its atoms forms a share of each Fock matrix and the reduction the
    # solver already makes adds the shares. One rank is given the whole of it.
    #
    # NOTE: divided by the work the atoms carry and not by their number, as the Fock
    # build divides them: dealing them round robin gave the ranks equal counts of
    # auxiliary functions and unequal counts of values, and a build ends when its
    # last rank ends.
    #
    # NOTE: a list of atoms is those atoms. One rank is given every atom by name,
    # and a rank dealt none of them has an empty list, which holds nothing and needs
    # no memory for it.
    if solver.nodes == 1:
        aux_atoms = list(range(molecule.number_of_atoms()))
    else:
        weights = solver._ri_jk_drv.aux_atom_weights(molecule, basis,
                                                     solver._ri_jk_aux_basis,
                                                     solver.eri_thresh)
        aux_atoms = molecule.partition_atoms_by_weight(solver.comm, weights)

    # NOTE: the response driver contracts the B vectors and cannot form them
    # again, so the mode which holds them is the only one it can use.
    #
    # NOTE: the budget is the smallest of the ranks', so that every rank takes
    # the same decision against the largest of their peaks.
    budget = min(solver.comm.allgather(memory_budget(solver, molecule, basis)))

    # NOTE: the exchange work arrays are a transient on top of the B vectors and
    # the rest of the live set, so they take a fraction of the rank's share and
    # not all of it. This keeps the peak from growing with the size of the node
    # and avoids the driver's fixed 4 GiB default, which is neither derived from
    # the machine nor from the rank's share. The bounds keep the batches sane on
    # very small and very large nodes.
    exchange_budget = min(
        budget, max(1 * 1024**3, min(8 * 1024**3, budget // 8)))

    response_drv = SimdRIJKResponseDriver(solver.eri_thresh, 0,
                                          exchange_budget)

    # NOTE: the peak is not the B vectors alone. Forming them holds the
    # integrals of a batch beside them, and a build holds the transformed
    # factors, the Coulomb buffers and its densities. The densities of a build
    # are taken in chunks small enough to fit (see _max_densities), so the
    # check here is that forming fits and that a build of one density of the
    # largest shape the solvers hand over does.
    forming, build, nheld = response_drv.fixed_memory(molecule, basis,
                                                      solver._ri_jk_aux_basis,
                                                      solver.eri_thresh,
                                                      aux_atoms, omega > 0.0)

    peak = max(
        response_drv.peak_memory(forming, build, *shape, 1, nheld)
        for shape in _worst_shapes(molecule, basis))

    peaks = solver.comm.allgather(peak)

    if max(peaks) > budget:
        # NOTE: the solver falls back to the four-center integrals, and so does
        # every solver this one lends its B vectors to (see lend). ri, ri_jk
        # and ri_coulomb are left as they were set.
        solver._ri_jk_fallback = True
        solver._ri_jk_drv = None
        solver._ri_jk_response_drv = None
        solver._ri_jk_aux_basis = None
        solver._ri_jk_memory = None

        where = ('available' if solver.nodes == 1 else 'on the largest rank')
        solver.ostream.print_warning(
            'RI-JK would need roughly ' +
            f'{Profiler.memory_to_string(max(peaks))} of ' +
            f'{Profiler.memory_to_string(budget)} {where}, so the response ' +
            'Fock matrices are formed from the four-center integrals.')
        solver.ostream.print_blank()
        solver.ostream.flush()
        return

    # NOTE: the metric is inverted through its Cholesky factor (the False passed
    # as use_inverse_square_root here and to prepare) with the solver's own
    # ri_metric_threshold. The SCF's ri_metric_route and ri_metric_threshold are
    # not handed to response, so a response after an SCF which inverted its
    # metric through its eigenvalues still takes the Cholesky route. With the
    # SCF's defaults (cholesky, 1e-12) the two invert the metric the same way.
    metric_erf = PackedMatrix()

    if omega > 0.0:
        # NOTE: the two come out of one call, which forms both operators over one
        # set of primitive pairs. There is no way to answer here as make_metric
        # answers one: the way which holds the B vectors is the only way a range
        # separated build has, and the response path has no other way either.
        metric, metric_erf = solver._ri_jk_drv.make_metric_rs(
            molecule, solver._ri_jk_aux_basis, solver.ri_metric_threshold, False,
            rimode.in_memory, omega)
        mode = rimode.in_memory
    else:
        metric, mode = solver._ri_jk_drv.make_metric(molecule,
                                                     solver._ri_jk_aux_basis,
                                                     solver.ri_metric_threshold,
                                                     False, rimode.in_memory)

    # NOTE: an empty share means this rank was dealt no auxiliary atoms, and not
    # that it should take all of them, as in the Fock build's own call.
    solver._ri_jk_drv.prepare(molecule, basis, solver._ri_jk_aux_basis,
                              solver.eri_thresh, budget,
                              solver.ri_metric_threshold, False, mode,
                              aux_atoms, metric, solver.nodes, omega,
                              metric_erf)

    # NOTE: the ranks' shares have to cover the auxiliary basis once between them,
    # as in the SCF, or the exchange the ranks add up is not the whole one and
    # nothing says so.
    naux = solver._ri_jk_aux_basis.get_dimensions_of_basis()

    held = solver.comm.allgather(solver._ri_jk_drv.number_of_aux_functions())

    assert_msg_critical(
        sum(held) == naux,
        f'{type(solver).__name__}: the ranks hold {held} of {naux} ' +
        'auxiliary functions, which is not the auxiliary basis once over ' +
        'the ranks')

    solver._ri_jk_response_drv = response_drv

    # NOTE: what a build may take, which _max_densities cuts its densities by.
    solver._ri_jk_memory = (budget, forming, build, nheld)

    # NOTE: the most this rank was predicted to reach so far, which the builds
    # raise as they take more densities (see _note_peak).
    solver._ri_jk_predicted_peak = peak

    # NOTE: the B vectors are a share of the auxiliary basis on more than one rank,
    # and the response driver is told which functions that share carries. Left to
    # sweep the whole basis it would form, zero and multiply a matrix for every
    # function another rank holds, which is that rank's work done again as zeros and
    # is why the response build did not fall with the ranks while the SCF build did.
    # A rank dealt no atoms is told an empty share, which is not the same as the
    # whole basis and is why the functions are named rather than counted.
    solver._ri_jk_response_drv.set_aux_functions(
        solver._ri_jk_drv.get_aux_functions())

    solver.ostream.print_info(
        'Using the SIMD resolution of the identity (RI-JK) for response.')
    if omega > 0.0:
        solver.ostream.print_info(
            'Range-separated functional: two sets of B vectors at ' +
            f'omega = {omega:.3f}.')

    if solver.nodes > 1:
        solver.ostream.print_info(
            f'The auxiliary basis is divided over {solver.nodes} ranks.')

    where = ('available' if solver.nodes == 1 else 'on the largest rank')
    solver.ostream.print_info(
        'RI-JK needs roughly ' +
        f'{Profiler.memory_to_string(max(peaks))} of ' +
        f'{Profiler.memory_to_string(budget)} {where} for a build of one ' +
        'density; larger batches are cut to fit.')
    solver.ostream.print_blank()
    solver.ostream.flush()


def _worst_shapes(molecule, basis):
    """The shapes of a build of one density which the setup checks.

    :param molecule:
        The molecule.
    :param basis:
        The AO basis set.

    :return:
        The arguments of the response driver's peak_memory after forming and
        build, as (nao, nvec, nterms, ncalls, nspins), for each shape a solver
        may hand over.

    .. note::
        A linear response trial vector has two terms over the occupied orbitals;
        a block diagonal density is two calls, the second over the virtual
        orbitals, which may be nearly as many as the basis; an open shell has
        two spins of two terms. The virtual orbitals are bounded by the basis
        less the occupied ones, as the orbitals may be fewer than the basis.
    """

    nao = basis.get_dimensions_of_basis()
    nalpha = molecule.number_of_alpha_occupied_orbitals(basis)
    nbeta = molecule.number_of_beta_occupied_orbitals(basis)

    shapes = [
        (nao, nalpha, 2, 1, 1),
        (nao, max(nalpha, nao - nbeta), 1, 2, 1),
    ]

    if nalpha != nbeta:
        shapes.append((nao, max(nalpha, nbeta), 2, 1, 2))

    return shapes


def _max_densities(solver, nao, nvec, nterms, ncalls, nspins):
    """The most densities one call of a build may take within the budget.

    :param solver:
        The solver holding the B vectors and the memory they were checked with.
    :param nao:
        The number of basis functions.
    :param nvec:
        The columns of the largest left factor of the call.
    :param nterms:
        One for densities of one term, two for densities of two terms.
    :param ncalls:
        One, or two for a block diagonal density.
    :param nspins:
        One for a restricted reference, two for an unrestricted one.

    :return:
        The number of densities, at least one, or None where the solver was
        handed B vectors without the memory they were checked with and nothing
        is cut.
    """

    memory = getattr(solver, '_ri_jk_memory', None)

    if memory is None:
        return None

    budget, forming, build, nheld = memory

    drv = solver._ri_jk_response_drv

    # NOTE: at least one. The setup checked that one density of every shape a
    # solver hands over fits, so this is a bound that holds and not a guess.
    count = max(
        drv.max_densities(budget, forming, build, nao, nvec, nterms, ncalls,
                          nspins, nheld), 1)

    return count


def _note_peak(solver, nao, nvec, nterms, ncalls, nspins, ndens):
    """Records the most a build was predicted to reach, for memory profiling.

    :param solver:
        The solver holding the B vectors.
    :param ndens:
        The densities of the call.

    The other parameters are those of _max_densities.
    """

    memory = getattr(solver, '_ri_jk_memory', None)

    if memory is None:
        return

    _, forming, build, nheld = memory

    peak = solver._ri_jk_response_drv.peak_memory(forming, build, nao, nvec,
                                                  nterms, ncalls, nspins,
                                                  ndens, nheld)

    solver._ri_jk_predicted_peak = max(
        getattr(solver, '_ri_jk_predicted_peak', 0), peak)


def _chunks(ndens, count):
    """The ranges of the densities a batch is cut into.

    :param ndens:
        The densities of the batch.
    :param count:
        The most densities of a chunk, or None for the whole batch at once.

    :return:
        The (start, end) of each chunk.
    """

    if count is None or count >= ndens:
        return [(0, ndens)]

    return [(start, min(start + count, ndens))
            for start in range(0, ndens, count)]


def _spin_factors(dens_factors, spin):
    """The left factor and the one or two right factors of one spin.

    Two shapes are taken, as in the restricted case: a trial vector of the
    Tamm-Dancoff approximation carries one right factor, and one of linear response
    carries two, the second standing to the left of the shared factor transposed.
    """
    if len(dens_factors) == 2:
        left, rights = dens_factors
        transposed_rights = []
    elif len(dens_factors) == 3:
        left, rights, transposed_rights = dens_factors
    else:
        assert_msg_critical(
            False, 'rijkresponse: the factors of the ' + spin + ' spin are of a ' +
            'shape an unrestricted build does not take')

    return left, rights, transposed_rights


def fock_matrices_unrestricted(solver,
                               basis,
                               dens_factors_alpha,
                               dens_factors_beta,
                               exchange_scaling_factor,
                               erf_exchange_scaling_factor=0.0):
    """
    Computes the two-electron part for a batch of factorised densities of an
    unrestricted reference.

    :param solver:
        The solver holding the B vectors and the response driver.
    :param basis:
        The AO basis set.
    :param dens_factors_alpha:
        The factors of the alpha densities: the left factor the batch shares, which
        is that spin's occupied orbitals, and the right factor of each density.
    :param dens_factors_beta:
        The same for the beta spin, whose left factor has a number of columns of its
        own.
    :param exchange_scaling_factor:
        The fraction of exact exchange.
    :param erf_exchange_scaling_factor:
        The coefficient of the exchange of the attenuated operator, zero for a
        functional which is not range separated.

    :return:
        The Fock matrices as numpy arrays, **interleaved alpha and beta**, which is
        the order the dense path appends them in and the order the callers unpack.

    .. note::
        The driver has one routine for a range separated functional and another for
        the rest, and this picks between them. They are separate there so that a
        plain calculation cannot reach the attenuated code at all.
    """
    _check_not_lent(solver)

    drv = solver._ri_jk_response_drv
    bq = solver._ri_jk_drv.get_bq_vectors()
    aux = solver._ri_jk_aux_basis

    left_a, rights_a, transposed_a = _spin_factors(dens_factors_alpha, 'alpha')
    left_b, rights_b, transposed_b = _spin_factors(dens_factors_beta, 'beta')

    assert_msg_critical(
        len(rights_a) == len(rights_b),
        'rijkresponse: the two spins do not carry the same number of densities')

    if not rights_a:
        return []

    bq_erf = (solver._ri_jk_drv.get_bq_vectors_erf()
              if erf_exchange_scaling_factor != 0.0 else None)

    # NOTE: the densities are taken in chunks small enough that a call stays
    # within the budget the B vectors were checked with. The Fock matrix is
    # linear in the density, so the chunks are independent.
    nao = basis.get_dimensions_of_basis()
    nvec = max(left_a.shape[1], left_b.shape[1])
    nterms = 2 if (transposed_a or transposed_b) else 1
    ndens = len(rights_a)

    count = _max_densities(solver, nao, nvec, nterms, 1, 2)

    # NOTE: alpha and beta of one density side by side, which is what the four
    # centre path appends and what every caller of it unpacks.
    fock_arrays = []

    for start, end in _chunks(ndens, count):
        _note_peak(solver, nao, nvec, nterms, 1, 2, end - start)

        shared = (basis, aux, _packed(left_a),
                  [_packed(r) for r in rights_a[start:end]],
                  [_packed(r) for r in transposed_a[start:end]],
                  _packed(left_b), [_packed(r) for r in rights_b[start:end]],
                  [_packed(r) for r in transposed_b[start:end]])

        if bq_erf is not None:
            focks_a, focks_b = drv.compute_unrestricted_rs(
                bq, bq_erf, *shared, exchange_scaling_factor,
                erf_exchange_scaling_factor)
        else:
            focks_a, focks_b = drv.compute_unrestricted(
                bq, *shared, exchange_scaling_factor)

        for fock_a, fock_b in zip(focks_a, focks_b):
            fock_arrays.append(fock_a.to_numpy())
            fock_arrays.append(fock_b.to_numpy())

        # NOTE: freed before the next chunk is formed, which the count was
        # taken for, rather than when the names are bound again after it.
        del shared, focks_a, focks_b

    return fock_arrays


def general_factors(mo, nocc, mo_density):
    """The two right factors any density is carried by.

    Splitting the left index of a density in the molecular orbitals into its
    occupied and its virtual half gives

        D = C(occ) r_a^T + C(vir) r_b^T

    with r_a and r_b as below, whatever blocks of the density are nonzero. It is
    the general form: a density which lives only between the occupied orbitals and
    the virtual ones is carried more cheaply by the two-term shape, and a block
    diagonal one by its own blocks, but a batch which mixes them needs one shape
    for all of them and this is it.

    :param mo:
        The molecular orbital coefficients.
    :param nocc:
        The number of occupied orbitals.
    :param mo_density:
        The density in the molecular orbitals, real.

    :return:
        The right factor of the occupied half and of the virtual half.
    """

    mo_occ = mo[:, :nocc]
    mo_vir = mo[:, nocc:]

    right_a = (np.matmul(mo_occ, mo_density[:nocc, :nocc].T) +
               np.matmul(mo_vir, mo_density[:nocc, nocc:].T))

    right_b = (np.matmul(mo_occ, mo_density[nocc:, :nocc].T) +
               np.matmul(mo_vir, mo_density[nocc:, nocc:].T))

    return right_a, right_b


def _release(solver):
    """Lets go of the B vectors the solver holds and of its fallback.

    :param solver:
        The solver.
    """

    solver._ri_jk_drv = None
    solver._ri_jk_response_drv = None
    solver._ri_jk_aux_basis = None
    solver._ri_jk_memory = None
    solver._ri_jk_predicted_peak = 0
    solver._ri_jk_fallback = False


def _check_not_lent(solver):
    """Raises an error where the solver's B vectors are lent out (see lend).

    :param solver:
        The solver.
    """

    assert_msg_critical(
        not solver._ri_jk_lent,
        f'{type(solver).__name__}: its RI-JK B vectors are lent to a solver ' +
        'it drives and cannot be used or formed again until they are returned')


@contextmanager
def lend(driver, solver):
    """Lends the B vectors of a driver to a solver it drives, for one call.

    A response driver drives linear solvers of its own, and every one of them
    needs the same B vectors of the same molecule in the same basis. Formed once
    and lent they cost one transformation of the integrals; formed by each
    solver they cost one for each.

    Used as ``with lend(driver, solver): solver.compute(...)``. Inside the block
    the solver borrows the driver's B vectors and forms none of its own (see
    initialize), and the driver, which stays their owner, raises an error if
    it uses or forms them. On leaving the block, also on an error, the solver
    holds no reference to them and the driver can use them again.

    A driver which fell back to the four-center integrals lends that, so its
    solvers do not form the B vectors it could not hold or warn again. A
    driver which holds neither lends nothing, and the solver sets itself up.

    :param driver:
        The driver which formed the B vectors.
    :param solver:
        The solver the driver is about to compute.
    """

    _check_not_lent(driver)

    assert_msg_critical(
        not driver._ri_jk_borrowing,
        f'{type(driver).__name__}: it lends only RI-JK B vectors it formed')

    assert_msg_critical(
        not solver._ri_jk_lent and not solver._ri_jk_borrowing and
        solver._ri_jk_drv is None and not solver._ri_jk_fallback,
        f'{type(solver).__name__}: it holds RI-JK B vectors or a fallback of ' +
        'its own and cannot borrow those of the driver')

    if not driver._ri_jk_fallback and driver._ri_jk_drv is None:
        yield
        return

    solver._ri_jk_fallback = driver._ri_jk_fallback
    solver._ri_jk_drv = driver._ri_jk_drv
    solver._ri_jk_response_drv = driver._ri_jk_response_drv
    solver._ri_jk_aux_basis = driver._ri_jk_aux_basis
    solver._ri_jk_memory = driver._ri_jk_memory
    solver._ri_jk_predicted_peak = getattr(driver, '_ri_jk_predicted_peak', 0)
    solver._ri_rss_at_setup = getattr(driver, '_ri_rss_at_setup', None)
    solver._ri_memory_reserve = getattr(driver, '_ri_memory_reserve', 0)

    solver._ri_jk_borrowing = True
    driver._ri_jk_lent = True

    try:
        yield
    finally:
        _release(solver)
        solver._ri_jk_borrowing = False
        driver._ri_jk_lent = False


def slice_factors(dens_factors, start, end):
    """The factors of the densities of one batch.

    The Fock build takes the densities in batches, so the factors are cut the same
    way. The shared left factors are the same for every batch and only the lists of
    right factors are cut.

    :param dens_factors:
        The factors of the whole set.
    :param start:
        The first density of the batch.
    :param end:
        The density past the last of the batch.

    :return:
        The factors of that batch, in the shape they were given in.
    """

    if dens_factors is None:
        return None

    if len(dens_factors) == 4:
        left_a, rights_a, left_b, rights_b = dens_factors
        return (left_a, rights_a[start:end], left_b, rights_b[start:end])

    if len(dens_factors) == 2:
        left, rights = dens_factors
        return (left, rights[start:end])

    left, rights, transposed_rights = dens_factors
    return (left, rights[start:end], transposed_rights[start:end])


def _packed(array):
    """The array as a general packed matrix."""

    array = np.ascontiguousarray(array)
    matrix = PackedMatrix(array.shape[0], array.shape[1], mat_t.general)
    matrix.from_numpy(array)
    return matrix


def fock_matrices(solver, basis, dens_factors, exchange_scaling_factor,
                  erf_exchange_scaling_factor=0.0):
    """
    Computes the two-electron part for a batch of factorised densities.

    :param solver:
        The solver holding the B vectors and the response driver.
    :param basis:
        The AO basis set.
    :param dens_factors:
        The factors, in one of the three shapes this module describes.
    :param exchange_scaling_factor:
        The fraction of exact exchange.
    :param erf_exchange_scaling_factor:
        The coefficient of the exchange of the attenuated operator, which is the
        erf coefficient of a hybrid range-separated functional and zero for every
        other calculation. It is subtracted as the plain exchange is, so this is
        the same number the four-centre way is passed.

    :return:
        The Fock matrices as numpy arrays, one for each density.
    """

    _check_not_lent(solver)

    drv = solver._ri_jk_response_drv
    bq = solver._ri_jk_drv.get_bq_vectors()
    aux = solver._ri_jk_aux_basis

    # NOTE: the attenuated B vectors, or nothing where the functional does not ask
    # for them. The driver refuses a coefficient without them rather than leaving
    # the long-range term out in silence.
    bq_erf = (solver._ri_jk_drv.get_bq_vectors_erf()
              if erf_exchange_scaling_factor != 0.0 else None)

    # NOTE: what comes back is twice the Coulomb less the scaled exchange
    # already, which is what the builders are asked for, so nothing is scaled
    # here. A pure functional asks for a scaling of zero and is given twice the
    # Coulomb, where the dense path forms the Coulomb and doubles it.

    # NOTE: the densities are taken in chunks small enough that a call stays
    # within the budget the B vectors were checked with. The Fock matrix is
    # linear in the density, so the chunks are independent; each one transforms
    # the shared left factor again, which is the price of the bound.
    nao = basis.get_dimensions_of_basis()

    fock_arrays = []

    if len(dens_factors) == 4:
        # NOTE: a block diagonal density, taken as the sum of its blocks. The
        # Fock matrix is linear in the density, so two batches added is the Fock
        # matrix of the sum, and neither block needs a factor the driver does not
        # already take.
        left_a, rights_a, left_b, rights_b = dens_factors

        nvec = max(left_a.shape[1], left_b.shape[1])

        count = _max_densities(solver, nao, nvec, 1, 2, 1)

        for start, end in _chunks(len(rights_a), count):
            _note_peak(solver, nao, nvec, 1, 2, 1, end - start)

            first = drv.compute(bq, basis, aux, _packed(left_a),
                                [_packed(r) for r in rights_a[start:end]],
                                exchange_scaling_factor, bq_erf,
                                erf_exchange_scaling_factor)

            second = drv.compute(bq, basis, aux, _packed(left_b),
                                 [_packed(r) for r in rights_b[start:end]],
                                 exchange_scaling_factor, bq_erf,
                                 erf_exchange_scaling_factor)

            fock_arrays += [
                a.to_numpy() + b.to_numpy() for a, b in zip(first, second)
            ]

            # NOTE: freed before the next chunk is formed.
            del first, second

        return fock_arrays

    if len(dens_factors) == 2:
        left, rights = dens_factors
        transposed_rights = None
    else:
        left, rights, transposed_rights = dens_factors

    nvec = left.shape[1]
    nterms = 1 if transposed_rights is None else 2

    count = _max_densities(solver, nao, nvec, nterms, 1, 1)

    for start, end in _chunks(len(rights), count):
        _note_peak(solver, nao, nvec, nterms, 1, 1, end - start)

        args = [
            bq, basis, aux,
            _packed(left), [_packed(r) for r in rights[start:end]]
        ]

        if transposed_rights is not None:
            args.append(
                [_packed(r) for r in transposed_rights[start:end]])

        focks = drv.compute(*args, exchange_scaling_factor, bq_erf,
                            erf_exchange_scaling_factor)

        fock_arrays += [fock.to_numpy() for fock in focks]

        # NOTE: freed before the next chunk is formed.
        del args, focks

    return fock_arrays
