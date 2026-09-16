"""One geometry optimization, one record.

The record is the optimization, its provenance, the steps it took and the wall
time of the whole of it. A step is an SCF and a gradient, so the time per step is
the quantity the methods differ in and the step count is the quantity they must
agree on -- two paths which converge in different numbers of steps are not being
compared on speed alone, and the record keeps both so a table cannot hide it.

Speedups are not recorded. They are computed when a table is rendered, against
the `full` row of the same molecule, basis and functional.
"""
import time
from pathlib import Path

import numpy as np
import veloxchem as vlx
from veloxchem.outputstream import OutputStream
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.optimizationdriver import OptimizationDriver

from scfbench import GEOMETRIES

# method -> (ri_jk, ri_jk_simd, ri_mode)
METHODS = {
    "full": (False, False, None),
    "ri_jk_simd": (True, True, "in_memory"),
}


def run(molecule_name, basis_name, aux_name, method, functional,
        conv_thresh=1.0e-8):
    """Runs one geometry optimization and returns its record."""
    molecule = vlx.Molecule.read_xyz_file(str(GEOMETRIES / f"{molecule_name}.xyz"))
    basis = vlx.MolecularBasis.read(molecule, basis_name.upper(), ostream=None)

    ri_jk, simd, ri_mode = METHODS[method]

    driver = ScfRestrictedDriver(ostream=OutputStream(None))
    driver.conv_thresh = conv_thresh
    if functional.upper() != "HF":
        driver.xcfun = functional
    if ri_jk:
        driver.ri_jk = True
        driver.ri_jk_simd = simd
        driver.ri_auxiliary_basis = aux_name.upper()
        if ri_mode is not None:
            driver.ri_mode = ri_mode

    opt = OptimizationDriver(driver)
    opt.ostream = OutputStream(None)

    t0 = time.time()
    results = opt.compute(molecule, basis)
    wall = time.time() - t0

    energies = results["opt_energies"]
    steps = len(energies)

    final = vlx.Molecule.read_xyz_string(results["final_geometry"])

    aux = None
    if ri_jk:
        aux = vlx.MolecularBasis.read(molecule, aux_name.upper(), ostream=None)

    return {
        "molecule": molecule_name,
        "atoms": molecule.number_of_atoms(),
        "basis": basis_name,
        "nao": basis.get_dimensions_of_basis(),
        "aux_basis": aux_name if ri_jk else None,
        "naux": aux.get_dimensions_of_basis() if aux is not None else None,
        "method": method,
        "ri_mode": ri_mode if ri_jk else None,
        "functional": functional,
        "conv_thresh": conv_thresh,
        "energy": energies[-1],
        "steps": steps,
        "wall": round(wall, 3),
        "wall_per_step": round(wall / steps, 3),
        "geometry": final.get_coordinates_in_bohr().tolist(),
        "labels": final.get_labels(),
    }
