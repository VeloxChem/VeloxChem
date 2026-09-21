"""One RI-J calculation, timed and recorded.

The Coulomb only approximation serves a **pure** functional: there is no exact
exchange to form, so the fitting never has to be closed for an orbital and the
metric never multiplies the tensor. It multiplies a vector of one value per
auxiliary function, twice a build. That is the whole of the difference from RI-JK,
and it is why this suite is separate rather than another column of `scfbench`.

**The fitting set is not the same one.** RI-JK is measured against
`def2-universal-jkfit`, which has to describe the products of orbitals an exchange
needs; a Coulomb only fitting wants `def2-universal-jfit`, which is a third the
size. Reading a row of this file against a row of the RI-JK tables therefore
compares two approximations **and** two auxiliary bases, and the honest comparison
within a row is the one against four centres.

What is timed is the driver's own profiler rather than a wrapper: `FockERI` is the
two-electron build and `FockXC` the quadrature, summed over the iterations. For a
pure functional those two are the whole of a Fock build, and the point of the suite
is what happens to their ratio.
"""
import time

from mpi4py import MPI

from veloxchem.outputstream import OutputStream
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.veloxchemlib import mpi_master

import veloxchem as vlx

from scfbench import _fock_timings, geometry, provenance, write  # noqa: F401

# NOTE: the Coulomb fitting set and not the exchange one. It is fixed here rather
# than in a runner so that two runners cannot come to disagree about which pair was
# measured, which is why `scfbench.fitting_set` exists for the other suite.
AUX = "def2-universal-jfit"

# method -> the attributes it sets on the driver
METHODS = {
    "full": {},
    "ri_j_conventional": {
        "ri_coulomb": True
    },
    "ri_j_simd": {
        "ri_coulomb": True,
        "ri_coulomb_simd": True,
        "ri_mode": "in_memory"
    },
    "ri_j_simd_direct": {
        "ri_coulomb": True,
        "ri_coulomb_simd": True,
        "ri_mode": "direct"
    },
}


def run(molecule_name,
        basis_name,
        method,
        functional="BLYP",
        aux_name=AUX,
        conv_thresh=1.0e-8,
        max_iter=100,
        charge=0,
        multiplicity=1,
        comm=None):
    """Runs one calculation and returns its record.

    :param functional:
        A **pure** functional. A hybrid is refused by the driver rather than
        answered without its exchange, and the refusal is the right answer: this
        approximation has no exchange to give it.
    :param multiplicity:
        Anything but one is run unrestricted, which is two Fock matrices an
        iteration and not one, so an open shell row is not comparable with a closed
        shell one however alike the two read.
    :param max_iter:
        A hundred and not the driver's fifty. A fitted calculation converges along a
        slightly different path from a four-centre one -- the surface it is on is
        not quite the same with respect to the rotations of the orbitals -- and it
        takes a few more iterations for it. A calculation which did not converge has
        no timing worth recording, so it raises rather than returning a row.
    """
    molecule = vlx.Molecule.read_xyz_file(str(geometry(molecule_name)))
    molecule.set_charge(charge)
    molecule.set_multiplicity(multiplicity)

    basis = vlx.MolecularBasis.read(molecule, basis_name.upper(), ostream=None)

    settings = METHODS[method]

    scf_class = (ScfRestrictedDriver
                 if multiplicity == 1 else ScfUnrestrictedDriver)

    driver = scf_class(comm=comm, ostream=OutputStream(None))
    driver.timing = True
    driver.conv_thresh = conv_thresh
    driver.max_iter = max_iter
    driver.xcfun = functional

    for key, value in settings.items():
        setattr(driver, key, value)

    if settings:
        driver.ri_auxiliary_basis = aux_name.upper()

    t0 = time.time()
    results = driver.compute(molecule, basis)
    wall = time.time() - t0

    if driver.rank != mpi_master():
        return None

    if not driver.scf_results:
        raise SystemExit(f"{molecule_name} {basis_name} {method}: the SCF did "
                         f"not converge in {max_iter} iterations")

    eri, xc, builds = _fock_timings(driver)

    aux = None
    if settings:
        aux = vlx.MolecularBasis.read(molecule, aux_name.upper(), ostream=None)

    # NOTE: what the driver settled on rather than what was asked for. `automatic`
    # answers one of the two, and a named way is honoured rather than checked
    # against the budget, so what ran is worth recording either way.
    mode = None
    if settings.get("ri_coulomb_simd") and getattr(driver, "_ri_drv", None):
        mode = str(driver._ri_drv.get_mode()).replace("rimode.", "")

    return {
        "molecule": molecule_name,
        "atoms": molecule.number_of_atoms(),
        "charge": charge,
        "multiplicity": multiplicity,
        "basis": basis_name,
        "nao": basis.get_dimensions_of_basis(),
        "aux_basis": aux_name if settings else None,
        "naux": aux.get_dimensions_of_basis() if aux is not None else None,
        "method": method,
        "functional": functional,
        "ri_mode": mode,
        "iterations": driver.num_iter,
        "builds": builds,
        "energy": results["scf_energy"],
        "wall": wall,
        "coulomb": eri,
        "xc": xc,
        "rest": wall - eri - xc,
    }
