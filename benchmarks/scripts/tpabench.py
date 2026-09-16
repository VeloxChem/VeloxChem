"""One two-photon absorption calculation, one record.

The record keeps the cross sections and the strengths beside the timings, because
what a two-photon calculation is for is the spectrum and not the speed, and a table
of times which does not carry the numbers cannot be checked afterwards.

The excited state part is timed apart from the ground state: the resolution of the
identity has already won on the ground state before the first trial vector.
"""
import time
from pathlib import Path

import numpy as np
import veloxchem as vlx
from veloxchem.outputstream import OutputStream
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tpatransitiondriver import TpaTransitionDriver

from scfbench import GEOMETRIES

# method -> (ri_jk, ri_jk_simd, ri_mode)
METHODS = {
    "full": (False, False, None),
    "ri_jk_simd": (True, True, "in_memory"),
}


def _listed(values):
    """A numpy array of reals as a list, or None."""
    if values is None:
        return None
    return [float(np.real(v)) for v in np.asarray(values).ravel()]


def run(molecule_name, basis_name, aux_name, method, functional, nstates=5,
        conv_thresh=1.0e-8, rsp_thresh=1.0e-5):
    """Runs one ground state and its two-photon absorption."""
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

    drv = TpaTransitionDriver(ostream=OutputStream(None))
    drv.nstates = nstates
    drv.conv_thresh = rsp_thresh
    if functional.upper() != "HF":
        drv.xcfun = functional
    if ri_jk:
        drv.ri_jk = True
        drv.ri_jk_simd = simd
        drv.ri_auxiliary_basis = aux_name.upper()

    t0 = time.time()
    results = drv.compute(molecule, basis, scf.scf_tensors)
    tpa_wall = time.time() - t0

    aux = None
    if ri_jk:
        aux = vlx.MolecularBasis.read(molecule, aux_name.upper(), ostream=None)

    strengths = results.get("tpa_strengths", {})

    return {
        "molecule": molecule_name,
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
        "rsp_thresh": rsp_thresh,
        "scf_wall": round(scf_wall, 3),
        "tpa_wall": round(tpa_wall, 3),
        # NOTE: the spectrum itself, kept so the table can be checked against the
        # numbers it was made from rather than only against another table.
        "cross_sections": _listed(results.get("cross_sections")),
        "photon_energies": _listed(results.get("photon_energies")),
        "oscillator_strengths": _listed(results.get("oscillator_strengths")),
        "tpa_strengths_linear": _listed(strengths.get("linear")),
        "tpa_strengths_circular": _listed(strengths.get("circular")),
    }
