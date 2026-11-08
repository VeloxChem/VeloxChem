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
"""Pure-NumPy geometry and hydration-shell analysis.

Minimum-image displacements (box vectors as rows; one cell or a cell per
frame), Kabsch superposition and RMSD, direction assignment, angle sets, and
framework atom selection and measurements. Atoms are selected by 0-based
index or by ``<element><ordinal>`` labels such as ``AL1`` (first Al) or ``O3``
(third O), counted in input order. Positions are in nm.

The analysis section's ``analyze_trajectory`` turns aqueous production frames
into IOD, CN and the geometry RMSD properties; IOD and distance/shape RMSDs
are returned in angstrom, angles in degrees.
"""

from __future__ import annotations

import itertools
import math
from typing import Sequence, Any
import re

import numpy as np

from .models import ANGSTROM_TO_NM, GEOMETRY_RMSD_PROPERTIES

# ==========================================================================
# geometry
# ==========================================================================

# ------------------------------------------------------------------- cells


def minimum_image(delta_nm: np.ndarray,
                  box_nm: np.ndarray | None) -> np.ndarray:
    """Apply the minimum-image convention to displacement vectors (rows).

    ``box_nm`` is one (3, 3) cell or a stack of cells (frames, 3, 3) matching
    ``delta_nm`` of shape (frames, n, 3); None leaves the vectors unchanged.
    """
    if box_nm is None:
        return delta_nm
    fractional = np.asarray(delta_nm) @ np.linalg.inv(box_nm)
    fractional -= np.round(fractional)
    return fractional @ box_nm


def box_widths(box_nm: np.ndarray) -> np.ndarray:
    """Perpendicular widths of a (possibly triclinic) cell."""
    volume = abs(np.linalg.det(box_nm))
    a, b, c = box_nm
    return np.array([
        volume / np.linalg.norm(np.cross(b, c)),
        volume / np.linalg.norm(np.cross(c, a)),
        volume / np.linalg.norm(np.cross(a, b))
    ])


# --------------------------------------------------------------- superposition


def kabsch_rotation(mobile_centered: np.ndarray,
                    reference_centered: np.ndarray) -> np.ndarray:
    """Proper rotation R minimizing |mobile @ R.T - reference| for centred point sets."""
    u, _, vt = np.linalg.svd(mobile_centered.T @ reference_centered)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    return vt.T @ np.diag([1.0, 1.0, d]) @ u.T


def kabsch_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    """Minimum RMSD between matched point sets P and Q over rigid rotation.

    Both sets are centred first, so the result ignores translation and
    rotation but keeps the chirality/torsion that distances and angles miss.
    """
    Pc = np.asarray(P, float) - np.mean(P, axis=0)
    Qc = np.asarray(Q, float) - np.mean(Q, axis=0)
    R = kabsch_rotation(Pc, Qc)
    return float(np.sqrt(np.mean(np.sum((Pc @ R.T - Qc)**2, axis=1))))


def superpose_all(points: np.ndarray, mobile: np.ndarray,
                  reference: np.ndarray) -> np.ndarray:
    """Apply the Kabsch transform taking ``mobile`` onto ``reference`` to ``points``."""
    center = mobile.mean(axis=0)
    rotation = kabsch_rotation(mobile - center,
                               reference - reference.mean(axis=0))
    return (points - center) @ rotation.T + reference.mean(axis=0)


def element_matched_rmsd(local, elements, reference, ref_elements) -> float:
    """Kabsch RMSD after matching atoms of equal element by optimal assignment.

    Row 0 (the metal) is fixed in both sets; the rest are reassigned within
    each element while alternating with superposition, so
    symmetry-equivalent ligands in any order are handled. Angstrom.
    """
    from scipy.optimize import linear_sum_assignment

    n = len(elements)
    order = [0] + sorted(range(1, n), key=lambda i: elements[i])
    ref_sorted = [0] + sorted(range(1, n), key=lambda j: ref_elements[j])
    # element-consistent start
    order = [order[ref_sorted.index(j)] for j in range(n)]
    for _ in range(5):  # alternate assignment and superposition
        mobile = superpose_all(local, local[order], reference)
        order = [0] + [None] * (n - 1)
        for element in set(elements[1:]):
            mine = [i for i in range(1, n) if elements[i] == element]
            theirs = [j for j in range(1, n) if ref_elements[j] == element]
            cost = np.linalg.norm(mobile[mine][:, None] -
                                  reference[theirs][None],
                                  axis=2)
            rows, cols = linear_sum_assignment(cost)
            for row, col in zip(rows, cols):
                order[theirs[col]] = mine[row]
    return kabsch_rmsd(local[order], reference)


# --------------------------------------------------------------- directions


def assign_to_reference(obs_dirs: np.ndarray,
                        ref_dirs: np.ndarray) -> np.ndarray:
    """For each observed direction, the ideal vertex it matches (Hungarian, max overlap)."""
    from scipy.optimize import linear_sum_assignment
    rows, cols = linear_sum_assignment(-(obs_dirs @ ref_dirs.T))
    perm = np.arange(len(obs_dirs))
    perm[rows] = cols
    return perm


def align_template(template: np.ndarray, observed: np.ndarray) -> np.ndarray:
    """Rotate unit ``template`` directions onto ``observed`` unit vectors.

    The rotation minimizing the RMSD (with the best direction assignment) is
    found from 64 deterministic starting rotations.
    """
    from scipy.optimize import linear_sum_assignment
    from scipy.spatial.transform import Rotation

    best, best_rmsd = None, np.inf
    for start in Rotation.random(64, random_state=0):
        rotation = start
        for _ in range(5):
            rotated = rotation.apply(template)
            rows, cols = linear_sum_assignment(-rotated @ observed.T)
            rotation, _ = Rotation.align_vectors(observed[cols], template[rows])
        rotated = rotation.apply(template)
        rows, cols = linear_sum_assignment(-rotated @ observed.T)
        rmsd = np.sqrt(
            np.mean(np.sum((rotated[rows] - observed[cols])**2, axis=1)))
        if rmsd < best_rmsd - 1e-9:
            best, best_rmsd = rotated, rmsd
    return best


# ------------------------------------------------------------------- angles


def pairwise_angles_deg(unit_dirs: np.ndarray) -> np.ndarray:
    """All C(n,2) angles (degrees) between unit direction vectors."""
    return np.asarray([
        math.degrees(
            math.acos(
                float(np.clip(np.dot(unit_dirs[i], unit_dirs[j]), -1.0, 1.0))))
        for i, j in itertools.combinations(range(len(unit_dirs)), 2)
    ])


def angle_rmsd_deg(unit_dirs: np.ndarray,
                   ideal_sorted_deg: np.ndarray) -> float:
    """RMSD between the sorted observed and ideal angle multisets."""
    observed = np.sort(pairwise_angles_deg(unit_dirs))
    n = min(len(observed), len(ideal_sorted_deg))
    if n == 0:
        return 0.0
    return float(np.sqrt(np.mean((observed[:n] - ideal_sorted_deg[:n])**2)))


# ---------------------------------------------------------------- selection


def atom_labels(elements: Sequence[str]) -> list[str]:
    """Return ``<ELEMENT><ordinal>`` labels in input order (AL1, O1, O2, ...)."""
    counts: dict[str, int] = {}
    labels = []
    for element in elements:
        key = element.upper()
        counts[key] = counts.get(key, 0) + 1
        labels.append(f"{key}{counts[key]}")
    return labels


def resolve_atom(selector: int | str, labels: Sequence[str]) -> int:
    """Map an index or label to a 0-based atom index, with a clear error."""
    if isinstance(selector,
                  (int, np.integer)) and not isinstance(selector, bool):
        if not 0 <= int(selector) < len(labels):
            raise ValueError(
                f"atom index {selector} is outside the structure (0..{len(labels) - 1})"
            )
        return int(selector)
    if isinstance(selector, str):
        key = selector.strip().upper()
        if key in labels:
            return list(labels).index(key)
        raise ValueError(
            f"unknown atom label {selector!r}; labels are element + ordinal, "
            f"e.g. {', '.join(labels[:4])}")
    raise ValueError(
        f"atom selector must be an index or a label, got {selector!r}")


def nearest_heavy_atoms(positions_nm, center: int, elements, excluded, box_nm,
                        size: int) -> list[int]:
    """The ``size`` nearest non-hydrogen atoms to ``center``, skipping ``excluded``."""
    delta = minimum_image(positions_nm - positions_nm[center], box_nm)
    distances = np.linalg.norm(delta, axis=1)
    candidates = [
        i for i in np.argsort(distances)
        if i != center and elements[i] != "H" and i not in excluded
    ]
    if len(candidates) < size:
        raise ValueError(
            f"the target has fewer than {size} coordinating heavy atoms")
    return [int(i) for i in candidates[:size]]


# ----------------------------------------------------------- measurements


def distance_angstrom(xyz_nm, atoms, box_nm) -> float:
    """Minimum-image distance between ``atoms[0]`` and ``atoms[1]``, in angstrom."""
    return float(
        np.linalg.norm(
            minimum_image(xyz_nm[atoms[1]] - xyz_nm[atoms[0]], box_nm)) /
        ANGSTROM_TO_NM)


def angle_degrees(xyz_nm, atoms, box_nm) -> float:
    """Minimum-image angle ``atoms[0]-atoms[1]-atoms[2]`` (vertex in the middle), in degrees."""
    left = minimum_image(xyz_nm[atoms[0]] - xyz_nm[atoms[1]], box_nm)
    right = minimum_image(xyz_nm[atoms[2]] - xyz_nm[atoms[1]], box_nm)
    cosine = left @ right / (np.linalg.norm(left) * np.linalg.norm(right))
    return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))


def count_within(xyz_nm, center: int, partners, cutoff_nm: float,
                 box_nm) -> float:
    """Number of ``partners`` closer than ``cutoff_nm`` to ``center`` (minimum image), as a float."""
    delta = minimum_image(xyz_nm[partners] - xyz_nm[center], box_nm)
    return float(np.sum(np.linalg.norm(delta, axis=1) < cutoff_nm))


# ==========================================================================
# analysis
# ==========================================================================

WATER_RESIDUES = frozenset(
    {"S01", "HOH", "SOL", "WAT", "TIP3", "TIP4", "TIP4PEW"})

# Metal-O RDF: 0.01 angstrom bins out to 5 angstrom, Gaussian-smoothed over 3 bins.
RDF_EDGES_NM = np.arange(0.0, 0.501, 0.001)
RDF_SMOOTHING_BINS = 3.0
SHELL_ONSET_RDF = 0.05  # g(r) at which the first shell starts (ion radius R1)

# ------------------------------------------------------------------ topology


def cdam_indices(topology) -> tuple[int | None, list[int]]:
    """Index of the CDAM metal and of its dummy sites (residue ``CDM``).

    Dummy sites are named ``D1`` ... ``Dn``; the metal is named by its element,
    which may itself start with D (Dy, Db, Ds).
    """
    metal, dummies = None, []
    atoms = topology.atoms() if callable(topology.atoms) else topology.atoms
    for atom in atoms:
        if atom.residue.name == "CDM":
            if re.fullmatch(r"D\d+", atom.name):
                dummies.append(atom.index)
            else:
                metal = atom.index
    return metal, dummies


def water_oxygens(topology) -> list[int]:
    """Indices of water oxygens (atom ``O`` in a water residue); raises ``ValueError`` if none."""
    oxygens = [
        a.index
        for a in topology.atoms()
        if a.name == "O" and a.residue.name in WATER_RESIDUES
    ]
    if not oxygens:
        raise ValueError("Could not isolate water oxygens. Found residues: "
                         f"{set(r.name for r in topology.residues())}")
    return oxygens


def count_waters(topology) -> int:
    """Number of water residues in the topology."""
    return sum(1 for r in topology.residues() if r.name in WATER_RESIDUES)


# ------------------------------------------------------------------ analysis


def analyze_trajectory(
        frames: np.ndarray,
        boxes: np.ndarray,
        topology,
        targets,
        reference_positions: np.ndarray | None = None) -> dict[str, Any]:
    """Measure the requested shell properties from a periodic metal-O RDF.

    No coordination number is assumed: the RDF is smoothed, its first peak
    gives IOD, and waters inside the first resolved low-density minimum give
    CN. Returns property values, ``R1`` (shell onset, nm), ``IOD_std`` and
    ``CN_std`` spreads, and the RDF itself under ``_analysis``.
    """
    from scipy.ndimage import gaussian_filter1d
    from scipy.signal import find_peaks

    metal, dummies = cdam_indices(topology)
    metal = 0 if metal is None else metal
    oxygens = water_oxygens(topology)
    displacement = minimum_image(
        frames[:, oxygens, :] - frames[:, metal, None, :], boxes)
    distances = np.linalg.norm(displacement, axis=-1)

    edges = RDF_EDGES_NM
    r = 0.5 * (edges[:-1] + edges[1:])
    counts, _ = np.histogram(distances, bins=edges)
    shell_volumes = (4.0 * math.pi / 3.0) * (edges[1:]**3 - edges[:-1]**3)
    density = len(oxygens) / np.linalg.det(boxes)
    g_r = counts / (shell_volumes * density.sum())

    results: dict[str, Any] = {}
    onset = np.where(g_r > SHELL_ONSET_RDF)[0]
    if not len(onset):
        raise ValueError(
            "No solvent molecules detected within 0.5 nm of the metal ion. "
            "The proposed parameters caused complete solvent depletion "
            "(unphysical state).")
    results["R1"] = float(r[onset[0]])

    smooth = gaussian_filter1d(g_r, sigma=RDF_SMOOTHING_BINS)
    peaks, _ = find_peaks(smooth, prominence=0.1 * smooth.max())
    if not len(peaks):
        raise ValueError("No resolved first hydration-shell RDF peak")
    peak = int(peaks[0])
    iod_nm = float(r[peak])
    left, middle, right = smooth[peak - 1:peak + 2]
    curvature = left - 2 * middle + right
    if curvature < 0:  # parabolic interpolation of the peak position
        iod_nm += float(np.clip(0.5 * (left - right) / curvature, -0.5,
                                0.5)) * (r[1] - r[0])

    analysis = {
        "method": "smoothed_rdf_first_peak_and_resolved_minimum_v1",
        "bin_width_angstrom": 0.01,
        "smoothing_sigma_angstrom": 0.03,
        "n_frames": len(distances),
        "r_angstrom": (r * 10).tolist(),
        "rdf": g_r.tolist(),
        "smoothed_rdf": smooth.tolist(),
    }

    if "IOD" in targets:
        results["IOD"] = iod_nm * 10.0
        shell = distances[distances < iod_nm + 0.05]
        results["IOD_std"] = float(shell.std() * 10.0) if len(shell) else 0.0

    if "CN" in targets:
        valleys, _ = find_peaks(-smooth)
        valleys = [
            i for i in valleys if r[i] > iod_nm +
            0.02 and np.max(smooth[max(0, i - 5):i + 6]) < 0.2 * smooth[peak]
        ]
        if not valleys:
            raise ValueError(
                "No resolved first-shell RDF minimum; CN cannot be measured reliably"
            )
        cutoff_nm = float(r[valleys[0]])
        coordination = (distances <= cutoff_nm).sum(axis=1)
        results["CN"] = float(coordination.mean())
        results["CN_std"] = float(coordination.std())
        analysis["cn_cutoff_angstrom"] = cutoff_nm * 10
        analysis["cn_per_frame"] = coordination.tolist()

    results.update(
        geometry_rmsd(frames, metal, dummies, displacement, distances,
                      reference_positions, targets))
    results["_analysis"] = analysis
    return results


def geometry_rmsd(frames: np.ndarray, metal: int, dummies: list[int],
                  displacement: np.ndarray, distances: np.ndarray,
                  reference_positions: np.ndarray | None,
                  targets) -> dict[str, float]:
    """Geometry-fidelity RMSDs against the ideal polyhedron, averaged over frames.

    SHELL_*: the N nearest water oxygens (``displacement``/``distances`` from
    the metal) vs the ideal vertices. DUMMY_*: the CDAM dummy star vs its
    pre-MD geometry in ``reference_positions``. N is the number of dummies.
    """
    requested = {t for t in targets if t in GEOMETRY_RMSD_PROPERTIES}
    if not requested:
        return {}
    if reference_positions is None or not dummies:
        raise ValueError(
            "geometry RMSD targets require the CDAM dummy sites and a reference "
            "geometry; none were available for this system")

    ref = np.asarray(reference_positions, float)
    ref_dummy = ref[dummies] - ref[metal]
    ref_radii = np.linalg.norm(ref_dummy, axis=1)
    ref_dirs = ref_dummy / ref_radii[:, None]
    ideal_angles = np.sort(pairwise_angles_deg(ref_dirs))
    results = {}

    if requested & {"DUMMY_DIST_RMSD", "DUMMY_ANGLE_RMSD", "DUMMY_SHAPE_RMSD"}:
        vectors = frames[:, dummies, :] - frames[:, metal, None, :]
        radii = np.linalg.norm(vectors, axis=-1)
        if "DUMMY_DIST_RMSD" in requested:
            results["DUMMY_DIST_RMSD"] = float(
                np.sqrt(np.mean((radii - ref_radii[None, :])**2)) * 10.0)
        if "DUMMY_ANGLE_RMSD" in requested:
            results["DUMMY_ANGLE_RMSD"] = float(
                np.mean([
                    angle_rmsd_deg(v / rad[:, None], ideal_angles)
                    for v, rad in zip(vectors, radii)
                ]))
        if "DUMMY_SHAPE_RMSD" in requested:
            results["DUMMY_SHAPE_RMSD"] = float(
                np.mean([kabsch_rmsd(v, ref_dummy) for v in vectors]) * 10.0)

    if requested & {"SHELL_DIST_RMSD", "SHELL_ANGLE_RMSD", "SHELL_SHAPE_RMSD"}:
        shell_dist, shell_angle, shell_shape = [], [], []
        for frame_disp, frame_dist in zip(displacement, distances):
            nearest = np.argsort(frame_dist)[:len(dummies)]
            points, radii = frame_disp[nearest], frame_dist[nearest]
            if np.any(radii <= 0):
                continue
            dirs = points / radii[:, None]
            shell_dist.append(np.sqrt(np.mean((radii - radii.mean())**2)))
            if "SHELL_ANGLE_RMSD" in requested:
                shell_angle.append(angle_rmsd_deg(dirs, ideal_angles))
            if "SHELL_SHAPE_RMSD" in requested:
                ideal = ref_dirs[assign_to_reference(dirs,
                                                     ref_dirs)] * radii.mean()
                shell_shape.append(kabsch_rmsd(points, ideal))
        if "SHELL_DIST_RMSD" in requested and shell_dist:
            results["SHELL_DIST_RMSD"] = float(np.mean(shell_dist) * 10.0)
        if "SHELL_ANGLE_RMSD" in requested and shell_angle:
            results["SHELL_ANGLE_RMSD"] = float(np.mean(shell_angle))
        if "SHELL_SHAPE_RMSD" in requested and shell_shape:
            results["SHELL_SHAPE_RMSD"] = float(np.mean(shell_shape) * 10.0)
    return results
