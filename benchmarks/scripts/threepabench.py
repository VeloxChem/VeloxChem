"""One three-photon absorption calculation, one record.

The record keeps the strengths beside the timings, because what a three-photon
calculation is for is the spectrum and not the speed, and a table of times which does
not carry the numbers cannot be checked afterwards.

The driver reports strengths and not cross sections: it computes cross sections and
does not return them, the line being commented out above a note to check them. So
none are stored here rather than one being made up from the strengths.

Everything compared here is invariant to the phase of an excited state vector. The
transition moments are not -- the two runs pick the phase differently and half of
them come back with the opposite sign meaning the same thing -- so they are recorded
and not compared.
"""
import time
from pathlib import Path

import numpy as np
import veloxchem as vlx
from veloxchem.outputstream import OutputStream
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.threepatransitiondriver import ThreePATransitionDriver

from scfbench import GEOMETRIES

# method -> (ri_jk, ri_jk_simd, ri_mode)
METHODS = {
    "full": (False, False, None),
    "ri_jk_simd": (True, True, "in_memory"),
}


def _listed(values):
    """A sequence of reals as a list, or None."""
    if values is None:
        return None
    return [float(np.real(v)) for v in np.asarray(values).ravel()]


def _strengths(results, kind):
    """The strengths of one polarization, in the order of the photon energies.

    They are returned keyed by the frequency, which is the negative of the photon
    energy, and a dictionary has no order of its own. Matched by value so that a
    key which is a float cannot miss.
    """
    table = (results.get("3pa_strengths") or {}).get(kind)
    energies = results.get("photon_energies")
    if not table or energies is None:
        return None

    out = []
    for energy in energies:
        key = next((k for k in table if abs(float(k) + float(energy)) < 1.0e-10),
                   None)
        if key is None:
            return None
        out.append(float(np.real(table[key])))
    return out


def run(molecule_name, basis_name, aux_name, method, functional, nstates=5,
        conv_thresh=1.0e-8, rsp_thresh=1.0e-5):
    """Runs one ground state and its three-photon absorption."""
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

    drv = ThreePATransitionDriver(ostream=OutputStream(None))
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
    rsp_wall = time.time() - t0

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
        "occupied": molecule.number_of_alpha_electrons(),
        "method": method,
        "ri_mode": ri_mode if ri_jk else None,
        "functional": functional,
        "nstates": nstates,
        "conv_thresh": conv_thresh,
        "rsp_thresh": rsp_thresh,
        "scf_wall": round(scf_wall, 3),
        "rsp_wall": round(rsp_wall, 3),
        # NOTE: the spectrum itself, kept so the table can be checked against the
        # numbers it was made from rather than only against another table.
        "photon_energies": _listed(results.get("photon_energies")),
        "oscillator_strengths": _listed(results.get("oscillator_strengths")),
        "elec_trans_dipoles": _listed(results.get("elec_trans_dipoles")),
        "strengths_linear": _strengths(results, "linear"),
        "strengths_circular": _strengths(results, "circular"),
    }
