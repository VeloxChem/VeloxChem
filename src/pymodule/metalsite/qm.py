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
"""
The QM steps of the metal site force field builder: the constrained
optimization, the Hessian and the RESP charges.

These are the three collective steps -- every MPI rank enters them with the
same molecule -- and the only functions of the package that build a
VeloxChem driver. Everything else is in core.
"""

from mpi4py import MPI
from contextlib import contextmanager
import numpy as np

from ..veloxchemlib import mpi_master
from ..molecularbasis import MolecularBasis
from ..respchargesdriver import RespChargesDriver
from ..scfrestdriver import ScfRestrictedDriver
from ..scfunrestdriver import ScfUnrestrictedDriver
from ..scfhessiandriver import ScfHessianDriver
from ..optimizationdriver import OptimizationDriver
from ..outputstream import OutputStream
from ..errorhandler import assert_msg_critical
from ..molecule import Molecule
from . import core
from . import printing


@contextmanager
def _muted(driver, mute_scf=True):
    """
    Mutes the output stream of a driver for the duration of a run.

    OutputStream.mute is reference counted, so a mute that is not balanced
    by its unmute swallows every later warning of that stream outright.
    Doing it here means the balance holds even when the driver raises.

    :param driver:
        The driver whose stream to mute.
    :param mute_scf:
        Whether to mute at all.
    """

    if mute_scf:
        driver.ostream.mute()
    try:
        yield
    finally:
        if mute_scf:
            driver.ostream.unmute()


def _get_scf_driver(molecule,
                    basis_set_label='def2-svp',
                    scf_drv=None,
                    xcfun='PBE0',
                    comm=None,
                    ostream=None):
    """
    Returns the SCF driver and basis set for the active site

    The driver assigned to scf_drv is used exactly as given, so it carries
    every QM setting beyond the functional and the basis set.

    :param molecule:
        The active site molecule.

    :return:
        The tuple of the SCF driver and the basis set.
    """

    comm = MPI.COMM_WORLD if comm is None else comm
    ostream = OutputStream(None) if ostream is None else ostream

    if scf_drv is None:
        if molecule.get_multiplicity() != 1:
            scf_drv = ScfUnrestrictedDriver(comm, ostream)
        else:
            scf_drv = ScfRestrictedDriver(comm, ostream)
        if xcfun is not None:
            scf_drv.xcfun = xcfun

    basis = MolecularBasis.read(molecule, basis_set_label)

    return scf_drv, basis


def _run_scf(scf_drv, molecule, basis, mute_scf=True):
    """
    Runs the SCF for a geometry, whatever the driver already holds.

    Both the gradient and the Hessian driver reuse whatever is already in
    scf_drv.scf_results and only fall back to running an SCF when it is
    empty. They do not check that those results belong to the geometry they
    were handed, so the SCF is run explicitly whenever the geometry may have
    moved since the driver was last used.

    :param scf_drv:
        The SCF driver.
    :param molecule:
        The active site molecule.
    :param basis:
        The basis set.
    :param mute_scf:
        Whether to mute the driver while it runs.
    """

    with _muted(scf_drv, mute_scf):
        scf_drv.compute(molecule, basis)


def optimize_active_site(active_site,
                         frozen_indices=None,
                         mute_scf=True,
                         ostream=None,
                         basis_set_label='def2-svp',
                         scf_drv=None,
                         xcfun='PBE0',
                         constrain_capping_hydrogens=False,
                         comm=None):
    """
    Optimizes the active site with the beta carbons frozen.

    Freezing them keeps the spatial arrangement imposed by the protein
    backbone. Without it the site relaxes to a gas-phase geometry, and
    since the Seminario method takes the equilibrium values straight from
    the geometry, every fitted bond length and angle would then describe
    the wrong structure.

    :param frozen_indices:
        The zero-based active site indices to freeze. Defaults to
        constrained_indices().

    :return:
        The tuple of the optimized molecule and the results of the
        optimization driver. The caller decides whether to put the molecule
        back into the active site.
    """

    if frozen_indices is None:
        frozen_indices = core.constrained_indices(
            active_site,
            constrain_capping_hydrogens=constrain_capping_hydrogens)

    molecule = active_site['molecule']

    constraint = core.freeze_constraints(active_site,
                                         frozen_indices=frozen_indices)
    constraints = None if constraint is None else [constraint]

    printing.print_muted_notice(
        f'the constrained optimization with {len(frozen_indices)} '
        'atom(s) frozen',
        mute_scf=mute_scf,
        ostream=ostream)

    scf_drv, basis = _get_scf_driver(molecule,
                                     basis_set_label=basis_set_label,
                                     scf_drv=scf_drv,
                                     xcfun=xcfun,
                                     comm=comm,
                                     ostream=ostream)

    opt_drv = OptimizationDriver(scf_drv)
    opt_drv.constraints = constraints

    with _muted(opt_drv, mute_scf):
        opt_results = opt_drv.compute(molecule, basis)

    optimized = Molecule.read_xyz_string(opt_results['final_geometry'])
    optimized.set_charge(molecule.get_charge())
    optimized.set_multiplicity(molecule.get_multiplicity())

    return optimized, opt_results


def hessian_pairs(active_site,
                  bond_count=2,
                  partial_hessian_cutoff=core.PARTIAL_HESSIAN_CUTOFF):
    """
    Finds the atom pairs a partial Hessian has to hold blocks for.

    The pairs extract_pairs walks out of the connectivity cover exactly the
    metal terms the fit makes on the geometry as it stands. That is one
    geometry's worth of perception, and a coordination the QM optimization
    opened up past metal_bond_cutoff is precisely what add_metal_bond is
    reached for afterwards -- at which point the Hessian holds an all-zero
    block for the new bond and Seminario gives it no force constant, with
    the whole Hessian to pay for again to repair it. The walk therefore
    starts from a connectivity that additionally treats every donor atom
    within partial_hessian_cutoff of a metal as bonded to it, so a bond
    added later is already covered.

    Only what is computed is widened. Nothing reads this connectivity back:
    the active site keeps the bonding it was perceived with, and the fit
    keeps fitting that.

    :param active_site:
        The active site whose Hessian is being computed. Not modified.
    :param bond_count:
        The number of bonds to walk, as extract_pairs takes it.
    :param partial_hessian_cutoff:
        The distance in Angstrom within which an unbonded donor atom is
        covered along with the metal anyway. None walks the connectivity as
        it stands.

    :return:
        The tuple of the sorted pair list and the sorted atom list.
    """

    matrix = np.array(active_site['connectivity_matrix'], dtype=bool)
    metals = list(active_site['metal_indices'])

    if partial_hessian_cutoff is not None:
        molecule = active_site['molecule']
        coordinates = molecule.get_coordinates_in_angstrom()
        labels = molecule.get_labels()

        donors = [
            index for index, label in enumerate(labels)
            if label in core.DONOR_ELEMENTS
        ]

        for metal in metals:
            for donor in donors:
                if matrix[metal, donor]:
                    continue
                distance = np.linalg.norm(coordinates[donor] -
                                          coordinates[metal])
                if distance <= partial_hessian_cutoff:
                    matrix[metal, donor] = True
                    matrix[donor, metal] = True

    return core.extract_pairs(matrix, metals, bond_count=bond_count)


def compute_hessian(active_site,
                    atom_pairs=None,
                    mute_scf=True,
                    basis_set_label='def2-svp',
                    scf_drv=None,
                    xcfun='PBE0',
                    ostream=None,
                    comm=None):
    """
    Computes the nuclear Hessian of the active site.

    With atom_pairs the analytical Hessian is restricted to those blocks
    and everything else is left at zero, which is all the Seminario method
    needs when only the metal terms are being fitted. The diagonal blocks
    are added by the Hessian driver itself, so only the off-diagonal pairs
    need to be given.

    compute() passes the pairs of hessian_pairs unless
    calculate_partial_hessian is off, in which case it asks for the whole
    Hessian instead.

    :param atom_pairs:
        The list of zero-based (i, j) tuples, typically from
        hessian_pairs. None computes the full Hessian.

    :return:
        The Hessian as a (3N, 3N) numpy array in Hartree per Bohr squared.
    """

    molecule = active_site['molecule']

    scf_drv, basis = _get_scf_driver(molecule,
                                     basis_set_label=basis_set_label,
                                     scf_drv=scf_drv,
                                     xcfun=xcfun,
                                     comm=comm,
                                     ostream=ostream)

    assert_msg_critical(
        scf_drv.solvation_model is None, 'compute_hessian: ScfHessianDriver '
        'does not support a solvation model')

    # The Hessian driver reuses scf_drv.scf_results without checking
    # which geometry they belong to, so the SCF is run here for the
    # current one. This costs nothing: the driver would otherwise run the
    # same SCF itself.
    printing.print_muted_notice('the SCF for the Hessian',
                                mute_scf=mute_scf,
                                ostream=ostream)
    _run_scf(scf_drv, molecule, basis, mute_scf=mute_scf)

    printing.print_muted_notice('the Hessian',
                                mute_scf=mute_scf,
                                ostream=ostream)

    hessian_drv = ScfHessianDriver(scf_drv)
    # the numerical path ignores atom_pairs entirely and would silently
    # compute the full Hessian instead
    hessian_drv.numerical = False
    if atom_pairs is None:
        hessian_drv.atom_pairs = None
    else:
        hessian_drv.atom_pairs = [tuple(pair) for pair in atom_pairs]

    with _muted(hessian_drv, mute_scf):
        hessian_drv.compute(molecule, basis)

    hessian = np.copy(hessian_drv.hessian)

    return hessian


def compute_resp_charges(active_site, mute_scf=True, comm=None, ostream=None):
    """
    Computes RESP charges for the active site.

    :return:
        The partial charges as an (N,) numpy array.
    """

    molecule = active_site['molecule']

    printing.print_muted_notice('the RESP charge fit at Hartree-Fock/6-31G*',
                                mute_scf=mute_scf,
                                ostream=ostream)

    resp_drv = RespChargesDriver(comm, ostream)

    # Neither a basis nor SCF results are passed: the driver then defaults
    # to Hartree-Fock with 6-31G*, which is what RESP charges are meant to
    # be fitted to, and runs its own SCF. Handing it the active site's own
    # functional and basis would silently fit the charges at a level the
    # RESP parameters were never derived for.
    with _muted(resp_drv, mute_scf):
        charges = resp_drv.compute(molecule)

    charges = comm.bcast(charges, root=mpi_master())
    charges = np.array(charges)

    return charges
