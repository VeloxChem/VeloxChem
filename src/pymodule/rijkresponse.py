#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#

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

import numpy as np

from .veloxchemlib import SimdRIJKFockDriver
from .veloxchemlib import SimdRIJKResponseDriver
from .veloxchemlib import PackedMatrix
from .veloxchemlib import mat_t
from .veloxchemlib import rimode
from .molecularbasis import MolecularBasis
from .errorhandler import assert_msg_critical


def memory_budget(solver):
    """
    Gets the memory the simd RI-JK driver may hold, per rank, in bytes.

    :param solver:
        The solver, whose communicator says how many ranks share this machine.

    :return:
        The memory budget in bytes.
    """

    try:
        import psutil
        available = psutil.virtual_memory().available
    except ImportError:
        available = 8 * 1024**3

    reserve = 4 * 1024**3

    # NOTE: the ranks of one node see the same free memory and would each claim the
    # whole of it, so the node would be promised several times over. The ranks
    # sharing a node are counted by their host name, as the Fock build counts them,
    # and the memory is divided between them. This asked for the whole of it on
    # every rank, which on eight ranks of a node promised that node its memory eight
    # times.
    import platform
    here = platform.node()
    on_this_node = max(solver.comm.allgather(here).count(here), 1)

    return int(max(available - reserve, 0.25 * available) / on_this_node)


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


def initialize(solver, molecule, basis):
    """
    Forms the B vectors the response driver contracts, on the solver.

    :param solver:
        The solver, which is given the drivers and the auxiliary basis.
    :param molecule:
        The molecule.
    :param basis:
        The AO basis set.
    """

    # NOTE: the B vectors are the dearest thing this path forms and they depend on
    # the molecule, the basis and the fitting set alone. A driver which drives
    # solvers of its own hands them down with share below, and each of those
    # solvers initializes itself in turn, so this is reached again with the work
    # already done. Doing it again would cost what it cost the first time.
    if is_prepared(solver):
        return

    assert_msg_critical(
        'jkfit' in solver.ri_auxiliary_basis.lower(),
        f'{type(solver).__name__}: RI-JK needs a jkfit auxiliary basis, and ' +
        f'{solver.ri_auxiliary_basis} fits the Coulomb alone')

    solver._ri_jk_aux_basis = MolecularBasis.read(molecule,
                                                  solver.ri_auxiliary_basis,
                                                  ostream=None)

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
    if solver.nodes == 1:
        aux_atoms = []
    else:
        weights = solver._ri_jk_drv.aux_atom_weights(molecule, basis,
                                                     solver._ri_jk_aux_basis,
                                                     solver.eri_thresh)
        aux_atoms = molecule.partition_atoms_by_weight(solver.comm, weights)


    needed = solver._ri_jk_drv.required_memory(molecule, basis,
                                               solver._ri_jk_aux_basis,
                                               solver.eri_thresh, aux_atoms,
                                               omega > 0.0)

    # NOTE: the response driver contracts the B vectors and cannot form them
    # again, so the mode which holds them is the only one it can use. The memory
    # is checked here rather than left to the allocator.
    budget = memory_budget(solver)

    assert_msg_critical(
        needed <= budget,
        f'{type(solver).__name__}: the B vectors need ' +
        f'{needed / 1024**3:.2f} GB of {budget / 1024**3:.2f} GB available')

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

    solver._ri_jk_drv.prepare(molecule, basis, solver._ri_jk_aux_basis,
                              solver.eri_thresh, budget,
                              solver.ri_metric_threshold, False, mode,
                              aux_atoms, metric, solver.nodes, omega,
                              metric_erf)

    solver._ri_jk_response_drv = SimdRIJKResponseDriver(solver.eri_thresh)

    solver.ostream.print_info(
        'Using the SIMD resolution of the identity (RI-JK) for response.')
    if omega > 0.0:
        solver.ostream.print_info(
            'Range-separated functional: two sets of B vectors at ' +
            f'omega = {omega:.3f}.')

    if solver.nodes > 1:
        solver.ostream.print_info(
            f'The auxiliary basis is divided over {solver.nodes} ranks.')

    solver.ostream.print_info(
        f'B vectors need {needed / 1024**3:.2f} GB of ' +
        f'{budget / 1024**3:.2f} GB available, a rank.')
    solver.ostream.print_blank()
    solver.ostream.flush()


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

    shared = (basis, aux, _packed(left_a), [_packed(r) for r in rights_a],
              [_packed(r) for r in transposed_a], _packed(left_b),
              [_packed(r) for r in rights_b], [_packed(r) for r in transposed_b])

    if erf_exchange_scaling_factor != 0.0:
        bq_erf = solver._ri_jk_drv.get_bq_vectors_erf()

        focks_a, focks_b = drv.compute_unrestricted_rs(
            bq, bq_erf, *shared, exchange_scaling_factor,
            erf_exchange_scaling_factor)
    else:
        focks_a, focks_b = drv.compute_unrestricted(bq, *shared,
                                                    exchange_scaling_factor)

    # NOTE: alpha and beta of one density side by side, which is what the four
    # centre path appends and what every caller of it unpacks.
    fock_arrays = []

    for fock_a, fock_b in zip(focks_a, focks_b):
        fock_arrays.append(fock_a.to_numpy())
        fock_arrays.append(fock_b.to_numpy())

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


def is_prepared(solver):
    """Whether the solver already holds B vectors it can contract.

    :param solver:
        The solver.

    :return:
        True where the B vectors are formed and ready.
    """

    drv = getattr(solver, '_ri_jk_drv', None)

    if (drv is None or not drv.is_prepared() or
            getattr(solver, '_ri_jk_response_drv', None) is None or
            getattr(solver, '_ri_jk_aux_basis', None) is None):
        return False

    # NOTE: and prepared for the functional this solver has. A driver handed down
    # by another one carries the operators that solver needed, which is the same
    # set in every case that arises today -- a driver drives solvers of its own
    # functional -- but a plain set answering a range separated build would drop
    # the long-range term, and one which asks here is told no and forms its own.
    return drv.get_omega() == range_separation(solver)


def share(source, target):
    """Hands the B vectors of one solver to another.

    A response driver drives linear solvers of its own, and every one of them
    needs the same B vectors of the same molecule in the same basis. Formed once
    and handed down they cost one transformation of the integrals; left to each
    solver they cost one for each, which on a molecule of seventy atoms was eight
    seconds apiece.

    :param source:
        The solver which holds them, or which does not, in which case nothing is
        handed over and the target forms its own.
    :param target:
        The solver which is given them.
    """

    # NOTE: is_prepared asks whether the source holds vectors for its own
    # functional. A source and a target of different functionals is not a case
    # which arises -- a driver drives solvers of its own -- and if it ever did, the
    # target's own is_prepared would reject what it was handed and it would form
    # the set it needs.
    if not is_prepared(source):
        return

    target._ri_jk_drv = source._ri_jk_drv
    target._ri_jk_response_drv = source._ri_jk_response_drv
    target._ri_jk_aux_basis = source._ri_jk_aux_basis


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

    if len(dens_factors) == 4:
        # NOTE: a block diagonal density, taken as the sum of its blocks. The
        # Fock matrix is linear in the density, so two batches added is the Fock
        # matrix of the sum, and neither block needs a factor the driver does not
        # already take.
        left_a, rights_a, left_b, rights_b = dens_factors

        first = drv.compute(bq, basis, aux, _packed(left_a),
                            [_packed(r) for r in rights_a],
                            exchange_scaling_factor, bq_erf,
                            erf_exchange_scaling_factor)

        second = drv.compute(bq, basis, aux, _packed(left_b),
                             [_packed(r) for r in rights_b],
                             exchange_scaling_factor, bq_erf,
                             erf_exchange_scaling_factor)

        return [a.to_numpy() + b.to_numpy() for a, b in zip(first, second)]

    if len(dens_factors) == 2:
        left, rights = dens_factors
        transposed_rights = None
    else:
        left, rights, transposed_rights = dens_factors

    args = [bq, basis, aux, _packed(left), [_packed(r) for r in rights]]

    if transposed_rights is not None:
        args.append([_packed(r) for r in transposed_rights])

    focks = drv.compute(*args, exchange_scaling_factor, bq_erf,
                        erf_exchange_scaling_factor)

    return [fock.to_numpy() for fock in focks]
