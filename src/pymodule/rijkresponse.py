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

The first is a trial vector of the Tamm-Dancoff approximation, the second one of
linear response, and the third a density of second order in the perturbation, which
is block diagonal in the orbitals: the occupied block carries the occupied orbitals
on both sides and the virtual block the virtual ones. The third is two calls of the
first added together, the Fock matrix being linear in the density.
"""

import numpy as np

from .veloxchemlib import SimdRIJKFockDriver
from .veloxchemlib import SimdRIJKResponseDriver
from .veloxchemlib import PackedMatrix
from .veloxchemlib import mat_t
from .veloxchemlib import rimode
from .molecularbasis import MolecularBasis
from .errorhandler import assert_msg_critical


def memory_budget():
    """
    Gets the memory the simd RI-JK driver may hold, in bytes.

    :return:
        The memory budget in bytes.
    """

    try:
        import psutil
        available = psutil.virtual_memory().available
    except ImportError:
        available = 8 * 1024**3

    reserve = 4 * 1024**3

    return int(max(available - reserve, 0.25 * available))


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

    assert_msg_critical(
        'jkfit' in solver.ri_auxiliary_basis.lower(),
        f'{type(solver).__name__}: RI-JK needs a jkfit auxiliary basis, and ' +
        f'{solver.ri_auxiliary_basis} fits the Coulomb alone')

    solver._ri_jk_aux_basis = MolecularBasis.read(molecule,
                                                  solver.ri_auxiliary_basis,
                                                  ostream=None)

    solver._ri_jk_drv = SimdRIJKFockDriver()

    needed = solver._ri_jk_drv.required_memory(molecule, basis,
                                               solver._ri_jk_aux_basis,
                                               solver.eri_thresh, [])

    # NOTE: the response driver contracts the B vectors and cannot form them
    # again, so the mode which holds them is the only one it can use. The memory
    # is checked here rather than left to the allocator.
    budget = memory_budget()

    assert_msg_critical(
        needed <= budget,
        f'{type(solver).__name__}: the B vectors need ' +
        f'{needed / 1024**3:.2f} GB of {budget / 1024**3:.2f} GB available')

    metric, mode = solver._ri_jk_drv.make_metric(molecule,
                                                 solver._ri_jk_aux_basis,
                                                 solver.ri_metric_threshold,
                                                 False, rimode.in_memory)

    solver._ri_jk_drv.prepare(molecule, basis, solver._ri_jk_aux_basis,
                              solver.eri_thresh, budget,
                              solver.ri_metric_threshold, False, mode, [],
                              metric, 1)

    solver._ri_jk_response_drv = SimdRIJKResponseDriver(solver.eri_thresh)

    solver.ostream.print_info(
        'Using the SIMD resolution of the identity (RI-JK) for response.')
    solver.ostream.print_info(
        f'B vectors need {needed / 1024**3:.2f} GB of ' +
        f'{budget / 1024**3:.2f} GB available.')
    solver.ostream.print_blank()
    solver.ostream.flush()


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


def fock_matrices(solver, basis, dens_factors, exchange_scaling_factor):
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

    :return:
        The Fock matrices as numpy arrays, one for each density.
    """

    drv = solver._ri_jk_response_drv
    bq = solver._ri_jk_drv.get_bq_vectors()
    aux = solver._ri_jk_aux_basis

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
                            exchange_scaling_factor)

        second = drv.compute(bq, basis, aux, _packed(left_b),
                             [_packed(r) for r in rights_b],
                             exchange_scaling_factor)

        return [a.to_numpy() + b.to_numpy() for a, b in zip(first, second)]

    if len(dens_factors) == 2:
        left, rights = dens_factors
        transposed_rights = None
    else:
        left, rights, transposed_rights = dens_factors

    args = [bq, basis, aux, _packed(left), [_packed(r) for r in rights]]

    if transposed_rights is not None:
        args.append([_packed(r) for r in transposed_rights])

    focks = drv.compute(*args, exchange_scaling_factor)

    return [fock.to_numpy() for fock in focks]
