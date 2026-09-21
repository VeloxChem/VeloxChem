"""One reduced two-photon absorption calculation, one record.

The record keeps the gamma of every frequency beside the timings. A two-photon
calculation is for the response function and not for the clock, and a table of times
which does not carry the numbers cannot be checked against anything afterwards.

Gamma is complex and is kept as its two parts, so the record is plain json and the
comparison of two runs is a comparison of numbers rather than of formatting.
"""
import time
from pathlib import Path

import numpy as np
import veloxchem as vlx
from veloxchem.outputstream import OutputStream
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tpareddriver import TpaReducedDriver

from scfbench import GEOMETRIES

# method -> (ri_jk, ri_jk_simd, ri_mode)
METHODS = {
    "full": (False, False, None),
    "ri_jk_simd": (True, True, "in_memory"),
}


def run(molecule_name, basis_name, aux_name, method, functional, frequencies,
        conv_thresh=1.0e-8, rsp_thresh=1.0e-6):
    """Runs one ground state and its reduced two-photon absorption."""
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

    drv = TpaReducedDriver(ostream=OutputStream(None))
    drv.frequencies = list(frequencies)
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

    # NOTE: gamma is keyed by the triple of frequencies the driver used, which is
    # (w, -w, w). It is stored against the frequency itself, in order, so two runs
    # can be compared without matching keys.
    gamma = results.get("gamma", {})
    ordered = [gamma[key] for key in sorted(gamma, key=lambda k: k[0])]

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
        "frequencies": [float(w) for w in frequencies],
        "conv_thresh": conv_thresh,
        "rsp_thresh": rsp_thresh,
        "scf_wall": round(scf_wall, 3),
        "tpa_wall": round(tpa_wall, 3),
        "gamma_real": [float(np.real(g)) for g in ordered],
        "gamma_imag": [float(np.imag(g)) for g in ordered],
        "cross_sections": [float(np.real(v))
                           for v in np.asarray(
                               results.get("cross_sections", [])).ravel()],
    }
