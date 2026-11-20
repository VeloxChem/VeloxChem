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

"""What the resolution of the identity in a response calculation shares between
RI-J (rijresponse) and RI-JK (rijkresponse): the budget of the simd RI driver
of a solver and the reserve it keeps back, the memory check printed with memory
profiling, the chunks a build cuts its densities into, and the lending of a
driver's B vectors to the solvers it drives (lend_ri).
"""

from contextlib import contextmanager

import psutil

from .veloxchemlib import release_free_memory
from .environment import get_ri_free_memory
from .errorhandler import assert_msg_critical
from .profiler import Profiler


def get_ri_memory_budget(solver, molecule, basis, reserve_floor=0):
    """
    Gets the memory the simd RI driver of a response solver may hold, per rank,
    in bytes. Every rank must call this. RI-J and RI-JK both take it.

    :param solver:
        The solver.
    :param molecule:
        The molecule, whose open shell doubles the rank reserve.
    :param basis:
        The AO basis set.
    :param reserve_floor:
        The least the rank's reserve may be, in bytes: what a nonlinear
        driver's smallest Fock batch holds (see
        NonlinearSolver._set_ri_reserve_floor), and zero for the others.

    :return:
        The memory budget in bytes: ri_memory_budget where the solver sets it
        (in GB, taken as it is), and otherwise the memory free to this rank less
        the reserves of the node and of the rank. The rank's reserve and the
        memory free to it are kept on the solver, for memory profiling and the
        fallback warning to report.
    """

    ri_memory_budget = getattr(solver, 'ri_memory_budget', None)

    if ri_memory_budget is not None:
        solver._ri_memory_reserve = 0
        solver._ri_memory_free = None
        return int(ri_memory_budget * 1.0e9)

    # NOTE: the rank's reserve is taken from what is free to the rank after the
    # node's reserve, so it is worked out here and not in get_ri_free_memory,
    # which the SCF and its gradient share with a fixed reserve.
    free = get_ri_free_memory(solver.comm, 0)

    # NOTE: the larger of the two and not their sum. The share rank_reserve
    # keeps is for the reduced space of a solver, and a nonlinear driver's Fock
    # build runs after its solvers have returned.
    reserve = max(rank_reserve(molecule, basis, free), reserve_floor)

    solver._ri_memory_reserve = reserve
    solver._ri_memory_free = free

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


def get_ri_reserve_floor_text(reserve_floor, budget_free, budget_reserve):
    """
    The part of a fallback warning that says the floor of a nonlinear driver
    set the reserve of the rank with the smallest budget.

    :param reserve_floor:
        The floor of the rank reserve, in bytes.
    :param budget_free:
        The memory free to the rank with the smallest budget, in bytes, or
        None where ri_memory_budget set the budget.
    :param budget_reserve:
        The reserve of that rank, in bytes.

    :return:
        The text, with a leading space, or an empty string where the floor
        did not set that reserve.
    """

    # NOTE: the budget is the difference of the free memory and the floor, or
    # zero where the floor exceeds the free memory.
    if (reserve_floor > 0 and budget_free is not None and
            budget_reserve == reserve_floor):
        return (f' ({Profiler.memory_to_string(budget_free)} free less ' +
                f'{Profiler.memory_to_string(reserve_floor)} for the ' +
                'smallest batch of the nonlinear Fock build)')

    return ''


def print_ri_reserve_floor(solver, reserve_floor):
    """
    Prints the memory each rank keeps for the smallest batch of a nonlinear
    driver's Fock builds, where that and not the share of rank_reserve set the
    reserve.

    :param solver:
        The solver.
    :param reserve_floor:
        The floor of the rank reserve, in bytes.
    """

    # NOTE: a budget the solver sets (ri_memory_budget) is taken as it is,
    # with no reserve.
    if (reserve_floor > 0 and
            getattr(solver, 'ri_memory_budget', None) is None and
            reserve_floor >= getattr(solver, '_ri_memory_reserve', 0)):
        solver.ostream.print_info(
            f'{Profiler.memory_to_string(reserve_floor)} of each rank is ' +
            'kept for the smallest batch of the nonlinear Fock build.')


def mark_ri_setup(solver):
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

    # NOTE: the largest growth of the peak reached before a solver this one
    # drives reset the peak (see lend_ri).
    solver._ri_peak_grown = 0


def print_ri_memory_check(solver):
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
    #
    # NOTE: at least zero. The start is the resident memory read at the setup
    # and the end the kernel's high-water mark, and the two counters can differ
    # by a few hundred kB, which turns a growth near zero negative.
    grown = max(Profiler.get_peak_memory() - solver._ri_rss_at_setup, 0,
                getattr(solver, '_ri_peak_grown', 0))
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


def get_ri_chunks(ndens, count):
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


def release_ri_jk(solver):
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


def release_ri_j(solver):
    """Lets go of the simd RI-J driver the solver set up for itself
    (rijresponse.initialize_ri_j) and of the memory it was checked with.

    :param solver:
        The solver.
    """

    if getattr(solver, '_ri_drv', None) is None:
        return

    solver._ri_drv = None
    solver._ri_j_memory = None

    # NOTE: the three-center integrals are freed in blocks the allocator may
    # keep resident; they are returned here.
    release_free_memory()


def check_ri_jk_not_lent(solver):
    """Raises an error where the solver's B vectors are lent out (see lend_ri).

    :param solver:
        The solver.
    """

    assert_msg_critical(
        not solver._ri_jk_lent,
        f'{type(solver).__name__}: its RI-JK B vectors are lent to a solver ' +
        'it drives and cannot be used or formed again until they are returned')


def _measures_own_window(driver):
    """
    Checks whether a driver which lends no B vectors holds the memory check of
    a simd RI-J driver it set up, which is what memory profiling prints at the
    end of its calculation (print_ri_memory_check).

    :param driver:
        The driver.

    :return:
        True if the driver measures the growth of the peak after its RI-J
        setup.
    """

    return (getattr(driver, 'memory_profiling', False) and
            getattr(driver, '_ri_rss_at_setup', None) is not None and
            driver._uses_ri_coulomb() and driver._ri_j_memory is not None)


def _keep_peak_growth(driver):
    """
    Keeps the largest growth of the peak after the driver's RI setup reached
    so far, before a solver it drives resets the peak.

    :param driver:
        The driver.
    """

    grown = max(Profiler.get_peak_memory() - driver._ri_rss_at_setup, 0)

    driver._ri_peak_grown = max(getattr(driver, '_ri_peak_grown', 0), grown)


@contextmanager
def lend_ri(driver, solver):
    """Lends the B vectors of a driver to a solver it drives, for one call.

    A response driver drives linear solvers of its own, and every one of them
    needs the same B vectors of the same molecule in the same basis. Formed once
    and lent they cost one transformation of the integrals; formed by each
    solver they cost one for each.

    Used as ``with lend_ri(driver, solver): solver.compute(...)``. Inside the
    block the solver borrows the driver's B vectors and forms none of its own
    (see rijkresponse.initialize_ri_jk), and the driver, which stays their
    owner, raises an error if it uses or forms them. On leaving the block, also
    on an error, the solver holds no reference to them and the driver can use
    them again.

    A driver which fell back to the four-center integrals lends that, so its
    solvers do not form the B vectors it could not hold or warn again. A
    driver which holds neither lends nothing, and the solver sets itself up.
    Where that driver set up the simd RI-J driver and profiles its memory, its
    memory check leaves out the solver, which checks its own builds. On
    leaving the block, also on an error, the solver lets go of the simd RI-J
    driver it set up.

    :param driver:
        The driver which formed the B vectors.
    :param solver:
        The solver the driver is about to compute.
    """

    check_ri_jk_not_lent(driver)

    assert_msg_critical(
        not driver._ri_jk_borrowing,
        f'{type(driver).__name__}: it lends only RI-JK B vectors it formed')

    assert_msg_critical(
        not solver._ri_jk_lent and not solver._ri_jk_borrowing and
        solver._ri_jk_drv is None and not solver._ri_jk_fallback,
        f'{type(solver).__name__}: it holds RI-JK B vectors or a fallback of ' +
        'its own and cannot borrow those of the driver')

    if not driver._ri_jk_fallback and driver._ri_jk_drv is None:
        measures = _measures_own_window(driver)
        if measures:
            _keep_peak_growth(driver)
        try:
            yield
        finally:
            # NOTE: the callers keep the solver until their own calculation
            # returns and read only is_converged from it after the block. The
            # solver's own RI-J integrals are not used after its compute, and
            # kept, they stay resident: the driver's RI-J budget and batches,
            # fixed at its setup, do not count them, and the next solver the
            # driver computes finds them taken from its free memory.
            release_ri_j(solver)
        # NOTE: a solver which sets itself up resets the peak (mark_ri_setup)
        # and checks its own builds. The driver's check then reads, beside the
        # growth it kept before the solver, the peak of what follows the
        # solver, and not the solver's peak against the driver's setup and
        # prediction.
        if measures:
            Profiler.reset_peak_memory()
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
        # NOTE: the driver's memory check measures from its own setup, which
        # comes before the solver's builds, so it takes the largest RI build
        # the solver was predicted to reach and the largest batch it held.
        driver._ri_jk_predicted_peak = max(
            getattr(driver, '_ri_jk_predicted_peak', 0),
            getattr(solver, '_ri_jk_predicted_peak', 0))
        driver._ri_batch_memory = max(getattr(driver, '_ri_batch_memory', 0),
                                      getattr(solver, '_ri_batch_memory', 0))

        release_ri_jk(solver)
        solver._ri_jk_borrowing = False
        driver._ri_jk_lent = False
