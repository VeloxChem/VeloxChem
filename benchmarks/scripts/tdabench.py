"""One excited state calculation, one record.

The record is the calculation, its provenance, the wall time of the excited state
part alone, and the excitation energies it found. The ground state is timed too but
separately: what this suite compares is the response, and an RI-JK calculation has
already won on the ground state before the first trial vector is formed.

The iterations are kept beside the time. A calculation which converges in a
different number of them is not being compared on speed alone, and the record keeps
both so a table cannot hide it.
"""
import time
from pathlib import Path

import numpy as np
import veloxchem as vlx
from veloxchem.outputstream import OutputStream
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tdaeigensolver import TdaEigenSolver
from veloxchem.lreigensolver import LinearResponseEigenSolver

from scfbench import GEOMETRIES

# NOTE: the two ways of solving for the same states. The Tamm-Dancoff
# approximation drops the de-excitation block, so its trial vector gives a density
# of one term where a full linear response vector gives two.
SOLVERS = {
    "tda": TdaEigenSolver,
    "rpa": LinearResponseEigenSolver,
}

# method -> (ri_jk, ri_jk_simd, ri_mode)
METHODS = {
    "full": (False, False, None),
    "ri_jk_simd": (True, True, "in_memory"),
}


def run(molecule_name, basis_name, aux_name, method, functional, nstates=5,
        conv_thresh=1.0e-8, tda_thresh=1.0e-5, solver="tda"):
    """Runs one ground state and its TDA, and returns the record."""
    molecule = vlx.Molecule.read_xyz_file(str(GEOMETRIES / f"{molecule_name}.xyz"))
    basis = vlx.MolecularBasis.read(molecule, basis_name.upper(), ostream=None)

    ri_jk, simd, ri_mode = METHODS[method]

    scf = ScfRestrictedDriver(ostream=OutputStream(None))
    scf.conv_thresh = conv_thresh
    if functional.upper() != "HF":
        scf.xcfun = functional
    if ri_jk:
        scf.ri_jk = True
        scf.ri_jk_simd = simd
        scf.ri_auxiliary_basis = aux_name.upper()
        if ri_mode is not None:
            scf.ri_mode = ri_mode

    t0 = time.time()
    scf.compute(molecule, basis)
    scf_wall = time.time() - t0

    tda = SOLVERS[solver](ostream=OutputStream(None))
    tda.nstates = nstates
    tda.conv_thresh = tda_thresh
    if functional.upper() != "HF":
        tda.xcfun = functional
    if ri_jk:
        tda.ri_jk = True
        tda.ri_jk_simd = simd
        tda.ri_auxiliary_basis = aux_name.upper()

    t0 = time.time()
    results = tda.compute(molecule, basis, scf.scf_tensors)
    tda_wall = time.time() - t0

    aux = None
    if ri_jk:
        aux = vlx.MolecularBasis.read(molecule, aux_name.upper(), ostream=None)

    return {
        "molecule": molecule_name,
        "solver": solver,
        "atoms": molecule.number_of_atoms(),
        "basis": basis_name,
        "nao": basis.get_dimensions_of_basis(),
        "aux_basis": aux_name if ri_jk else None,
        "naux": aux.get_dimensions_of_basis() if aux is not None else None,
        "occupied": molecule.number_of_alpha_electrons(),
        "method": method,
        "ri_mode": ri_mode if ri_jk else None,
        "functional": functional,
        "nstates": nstates,
        "conv_thresh": conv_thresh,
        "tda_thresh": tda_thresh,
        "scf_wall": round(scf_wall, 3),
        "tda_wall": round(tda_wall, 3),
        # NOTE: the iterations the excited state part took, which is the other
        # half of what its time is made of.
        "iterations": int(getattr(tda, "_cur_iter", -1)) + 1,
        "converged": bool(getattr(tda, "_is_converged", False)),
        "excitation_energies": [float(e) for e in results["eigenvalues"]],
    }
