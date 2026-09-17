"""One molecular gradient, one record.

The record is the calculation, its provenance, and the wall time of the gradient
alone -- the SCF that precedes it is timed too, but separately: what this suite
compares is the derivative, and an RI-JK calculation has already won on the energy
before the gradient starts.

Speedups are not recorded. They are computed when a table is rendered, against the
`full` row of the same molecule, basis and functional, so that re-measuring one
method cannot leave a stale ratio behind.

The gradient is repeated and the best taken, and both ways run in one process, so
that the comparison is never made across runs.
"""
import json
import time
from pathlib import Path

import numpy as np
import veloxchem as vlx
from veloxchem.outputstream import OutputStream
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.scfgradientdriver import ScfGradientDriver

from scfbench import GEOMETRIES, geometry, provenance, write

# method -> (ri_jk, ri_jk_simd, ri_mode)
METHODS = {
    "full": (False, False, None),
    "ri_jk_simd": (True, True, "in_memory"),
}


def run(molecule_name, basis_name, aux_name, method, functional,
        repeats=2, conv_thresh=1.0e-8, charge=0, multiplicity=1,
        max_iter=100):
    """Runs one SCF and times its gradient, and returns the record.

    :param charge:
        The charge of the molecule.
    :param multiplicity:
        The spin multiplicity. Anything but one is run unrestricted, whose
        gradient contracts an exchange for each spin where a closed shell
        contracts one and doubles it.
    :param max_iter:
        The iterations the SCF is allowed. A hundred and not the driver's fifty:
        a radical cation of caffeine takes 45 by the four centre way and 51 by
        either resolution of the identity, and a gradient of a calculation which
        did not converge is not a gradient.
    """
    molecule = vlx.Molecule.read_xyz_file(str(geometry(molecule_name)))
    molecule.set_charge(charge)
    molecule.set_multiplicity(multiplicity)

    basis = vlx.MolecularBasis.read(molecule, basis_name.upper(), ostream=None)

    ri_jk, simd, ri_mode = METHODS[method]

    scf_class = ScfRestrictedDriver if multiplicity == 1 else ScfUnrestrictedDriver

    driver = scf_class(ostream=OutputStream(None))
    driver.conv_thresh = conv_thresh
    driver.max_iter = max_iter
    if functional.upper() != "HF":
        driver.xcfun = functional
    if ri_jk:
        driver.ri_jk = True
        driver.ri_jk_simd = simd
        driver.ri_auxiliary_basis = aux_name.upper()
        if ri_mode is not None:
            driver.ri_mode = ri_mode

    t0 = time.time()
    results = driver.compute(molecule, basis)
    scf_wall = time.time() - t0

    # NOTE: a run which did not converge leaves no results behind, and the
    # gradient driver indexes them and raises a KeyError from the middle of
    # itself. What went wrong is said here instead.
    if not driver.scf_results:
        raise SystemExit(
            f"{molecule_name} {basis_name} {functional} {method}: the SCF did "
            f"not converge in {max_iter} iterations, so it has no gradient")

    grad_driver = ScfGradientDriver(driver)
    grad_driver.ostream = OutputStream(None)

    walls = []
    for _ in range(repeats):
        t0 = time.time()
        grad_driver.compute(molecule, basis)
        walls.append(time.time() - t0)

    gradient = grad_driver.gradient.copy()

    aux = None
    if ri_jk:
        aux = vlx.MolecularBasis.read(molecule, aux_name.upper(), ostream=None)

    mode = None
    if simd and getattr(driver, "_ri_drv", None) is not None:
        mode = str(driver._ri_drv.get_mode()).replace("rimode.", "")

    return {
        "molecule": molecule_name,
        "atoms": molecule.number_of_atoms(),
        "basis": basis_name,
        "nao": basis.get_dimensions_of_basis(),
        "aux_basis": aux_name if ri_jk else None,
        "naux": aux.get_dimensions_of_basis() if aux is not None else None,
        "occupied": molecule.number_of_alpha_electrons(),
        "method": method,
        "ri_mode": mode,
        "functional": functional,
        "conv_thresh": conv_thresh,
        "energy": results["scf_energy"],
        "iterations": driver.num_iter,
        "scf_wall": round(scf_wall, 3),
        "charge": int(charge),
        "multiplicity": int(multiplicity),
        "scf_type": "restricted" if multiplicity == 1 else "unrestricted",
        "grad_wall": round(min(walls), 3),
        "grad_walls": [round(w, 3) for w in walls],
        "grad_max": round(float(np.abs(gradient).max()), 9),
        # NOTE: the sum over the atoms, which is zero for a gradient of every
        # atom of a molecule that is not in a field. It is kept as a check that
        # travels with the number rather than as a claim made once in a message.
        "translational": float(np.abs(gradient.sum(axis=0)).max()),
        "_gradient": gradient.tolist(),
    }
