"""
molecule_toolkit.py
====================
A unified toolkit for building and modifying molecules in VeloxChem.

Merges the metal-complex assembly machinery from MoleculeModifier with the
organic substitution and fragment-coupling tools from MoleculeBuilder into one
coherent class.

Capabilities
------------
Metal complex assembly
    add_monodentate_ligand      — place any single-donor ligand at a target site
    add_universal_bidentate     — build en-type N-C-C-N chelates from scratch
    add_bidentate_rigid         — place rigid/aromatic bidentates (bipy, ppy, acac…)
    prepare_bidentate_conformation — rotate backbone to align donor lone-pairs
    find_chelate_donors         — auto-locate donor indices (N-N, N-C, C-N…)
    remove_H_from_donor         — deprotonate a cyclometallating carbon
    remove_ligand               — strip a ligand by cutting the M-donor bond
    transmetalate               — swap the metal, rescaling all M-L distances

Organic chemistry / fragment coupling
    substitute                  — swap a functional group matched by SMILES
    remove_all_fragment_matches — strip all occurrences of a SMILES fragment
    join_molecules              — couple two VeloxChem molecules at given indices
    link_at_vector              — attach a SMILES fragment at an explicit vector
    run_radical_addition        — add a SMILES radical to every aromatic C

Geometry helpers
    get_rotation_matrix         — Rodrigues rotation A→B (antiparallel-safe)
    get_geometry_vectors        — octahedral / tetrahedral site vectors
    get_lone_pair_vector        — LP direction for any donor atom
    rotate_vector               — Rodrigues single-vector rotation
"""

from __future__ import annotations

import numpy as np
import networkx as nx
import veloxchem as vlx
from rdkit import Chem
from rdkit.Chem import AllChem
from networkx.algorithms import isomorphism


class MoleculeToolkit:
    """
    Unified molecule-building and modification toolkit for VeloxChem.

    **All atom indices are 1-based**, following VeloxChem convention.
    The first atom in a molecule is atom 1, not atom 0.

    Instantiate once and use all methods on any VeloxChem Molecule objects:

        tk = MoleculeToolkit()
        mol = vlx.Molecule.read_molecule_string("Pt  0.0  0.0  0.0")
        mol = tk.add_monodentate_ligand(mol, 1, nh3, 1, [1,0,0], 2.05)
        ...
    """

    # ------------------------------------------------------------------
    # Radius tables
    # ------------------------------------------------------------------

    # Covalent radii (Å) for bond perception — extended with common metals
    cov_radii: dict[str, float] = {
        'H':  0.31, 'C':  0.76, 'N':  0.71, 'O':  0.66,
        'F':  0.57, 'P':  1.07, 'S':  1.05, 'Cl': 1.02,
        'Br': 1.20, 'I':  1.39,
        # Transition metals
        'Fe': 1.32, 'Co': 1.26, 'Ni': 1.24, 'Cu': 1.32,
        'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39,
        'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36,
    }
    _default_cov_radius: float = 0.76

    # Van der Waals radii (Å) for steric clash scoring (organic)
    vdw_radii: dict[str, float] = {
        'H':  1.20, 'C':  1.70, 'N':  1.55, 'O':  1.52,
        'F':  1.47, 'S':  1.80, 'Cl': 1.75, 'Br': 1.85,
    }
    _default_vdw_radius: float = 1.50

    # Default single-bond lengths when attaching organic fragments
    _bond_defaults: dict[str, float] = {
        'C': 1.54, 'N': 1.47, 'O': 1.41, 'S': 1.81,
    }

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    def __init__(self) -> None:
        pass   # all data lives at the class level; nothing to initialise

    # ------------------------------------------------------------------
    # Index convention
    # ------------------------------------------------------------------

    @staticmethod
    def _i(idx: int) -> int:
        """
        Convert a 1-based atom index (as used in VeloxChem) to the 0-based
        index used internally by NumPy arrays.  Every public method converts
        its user-facing index parameters through this helper immediately on
        entry, so all internal logic remains 0-based throughout.
        """
        if idx < 1:
            raise ValueError(
                f"Atom index {idx} is invalid — MoleculeToolkit uses "
                "1-based indices (VeloxChem convention). "
                "The first atom is index 1."
            )
        return idx - 1

    # ==================================================================
    # ── Section 1: Low-level geometry helpers ─────────────────────────
    # ==================================================================

    def _cov_radius(self, element: str) -> float:
        """Covalent radius for *element*, with a sensible fallback."""
        return self.cov_radii.get(element, self._default_cov_radius)

    def _vdw_radius(self, element: str) -> float:
        """Van der Waals radius for *element*, with a sensible fallback."""
        return self.vdw_radii.get(element, self._default_vdw_radius)

    def _get_molecule_graph(self, vlx_mol) -> nx.Graph:
        """
        Build a NetworkX connectivity graph from a VeloxChem Molecule using
        covalent radii with a 1.3× tolerance.  Each node carries an ``element``
        attribute and a ``mol_idx`` attribute (its original 0-based index).
        """
        coords = np.array(vlx_mol.get_coordinates_in_angstrom(), dtype=float)
        labels = vlx_mol.get_labels()
        G = nx.Graph()
        for i, lbl in enumerate(labels):
            G.add_node(i, element=lbl, mol_idx=i)
        for i in range(len(coords)):
            for j in range(i + 1, len(coords)):
                r_sum = self._cov_radius(labels[i]) + self._cov_radius(labels[j])
                if np.linalg.norm(coords[i] - coords[j]) < 1.3 * r_sum:
                    G.add_edge(i, j)
        return G

    def _get_neighbors_by_radii(
        self, coords: np.ndarray, labels: list, atom_idx: int
    ) -> np.ndarray:
        """Return indices of atoms covalently bonded to *atom_idx*."""
        r_i = self._cov_radius(labels[atom_idx])
        return np.array(
            [j for j, lbl in enumerate(labels)
             if j != atom_idx
             and np.linalg.norm(coords[j] - coords[atom_idx])
                 < 1.3 * (r_i + self._cov_radius(lbl))],
            dtype=int,
        )

    # ------------------------------------------------------------------
    # Core rotation helpers (used throughout)
    # ------------------------------------------------------------------

    def get_rotation_matrix(self, a: np.ndarray, b: np.ndarray) -> np.ndarray:
        """
        Return the 3×3 rotation matrix **R** such that ``R @ â = b̂``.

        Uses Rodrigues' formula.  The antiparallel edge case (a ≈ -b) is handled
        with a proper 180° rotation (det = +1), not a reflection.
        """
        a = np.array(a, dtype=float) / np.linalg.norm(a)
        b = np.array(b, dtype=float) / np.linalg.norm(b)
        v = np.cross(a, b)
        c = np.dot(a, b)

        if np.linalg.norm(v) < 1e-9:
            if c > 0:
                return np.eye(3)
            # Antiparallel: 180° around an arbitrary perpendicular axis
            perp = np.array([1., 0., 0.]) if abs(a[0]) < 0.9 else np.array([0., 1., 0.])
            ax = np.cross(a, perp); ax /= np.linalg.norm(ax)
            K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
            return np.eye(3) + 2 * K @ K

        s = np.linalg.norm(v)
        K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        return np.eye(3) + K + K @ K * ((1 - c) / s ** 2)

    def rotate_vector(
        self, v: np.ndarray, axis: np.ndarray, angle_rad: float
    ) -> np.ndarray:
        """
        Rotate vector *v* around *axis* by *angle_rad* (Rodrigues formula).
        Returns a new vector; does not modify *v*.
        """
        axis = axis / np.linalg.norm(axis)
        return (v * np.cos(angle_rad)
                + np.cross(axis, v) * np.sin(angle_rad)
                + axis * np.dot(axis, v) * (1.0 - np.cos(angle_rad)))

    def get_geometry_vectors(self, geometry: str) -> np.ndarray:
        """
        Unit vectors for standard coordination geometries.

        Parameters
        ----------
        geometry : ``'octahedral'`` or ``'tetrahedral'``

        Returns
        -------
        ndarray of shape (n_sites, 3)
        """
        geoms = {
            'tetrahedral': np.array(
                [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=float
            ) / np.sqrt(3),
            'octahedral': np.array(
                [[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]], dtype=float
            ),
        }
        key = geometry.lower()
        if key not in geoms:
            raise ValueError(f"Unknown geometry '{geometry}'. Choose from: {list(geoms)}")
        return geoms[key]

    # ------------------------------------------------------------------
    # Lone-pair vector
    # ------------------------------------------------------------------

    def _lone_pair_and_nbonds(
        self, coords: np.ndarray, labels: list, atom_idx: int
    ) -> tuple[np.ndarray, int]:
        """
        Core LP implementation.  Returns ``(lp_unit_vec, n_bonds)``.

        The LP direction is the *negative* of the sum of outgoing unit bond
        vectors — correct for sp3, sp2, and sp donors alike.
        ``n_bonds == 0`` signals an isolated atom (caller should handle).
        """
        pos_i  = coords[atom_idx]
        r_i    = self._cov_radius(labels[atom_idx])
        bv_sum = np.zeros(3)
        n      = 0
        last_j = 0
        for j, lbl in enumerate(labels):
            if j == atom_idx:
                continue
            r_sum = r_i + self._cov_radius(lbl)
            dist  = np.linalg.norm(coords[j] - pos_i)
            if dist < 1.3 * r_sum:
                bv_sum += (coords[j] - pos_i) / dist
                n += 1
                last_j = j

        if n == 0:
            return np.array([0., 0., 1.]), 0

        norm = np.linalg.norm(bv_sum)
        if norm < 1e-9:
            v0   = (coords[last_j] - pos_i) / np.linalg.norm(coords[last_j] - pos_i)
            perp = np.array([1., 0., 0.]) if abs(v0[0]) < 0.9 else np.array([0., 1., 0.])
            lp   = np.cross(v0, perp)
            return lp / np.linalg.norm(lp), n

        return -bv_sum / norm, n

    def get_lone_pair_vector(
        self, coords: np.ndarray, labels: list, atom_idx: int
    ) -> np.ndarray:
        """
        Public lone-pair interface.  Returns a unit vector pointing toward
        where a metal (or bond partner) would coordinate to *atom_idx*.

        *atom_idx* is **1-based** (VeloxChem convention).

        Works for any hybridisation:
        * sp³ N (NH₃)      — from 3 N–H vectors
        * sp² N (pyridine) — in-plane, away from ring
        * sp  C (CO)       — sigma LP on C
        """
        lp, _ = self._lone_pair_and_nbonds(coords, labels, self._i(atom_idx))
        return lp

    # ==================================================================
    # ── Section 2: Metal complex assembly ─────────────────────────────
    # ==================================================================

    def add_monodentate_ligand(
        self,
        complex_mol,
        metal_idx: int,
        ligand_mol,
        coord_atom_idx: int,
        target_vec,
        bond_length: float = 1.9,
    ):
        """
        Place a monodentate ligand so its donor atom sits at
        ``metal_pos + target_unit × bond_length``.

        The ligand body is oriented away from the metal using the centroid of
        all non-donor atoms.  A 360-step torsion scan minimises steric clash
        with the existing complex.  Single-atom ligands (Cl⁻, Br⁻, …) skip
        both the body rotation and the torsion scan.

        Parameters
        ----------
        complex_mol     : current complex
        metal_idx       : 1-based index of the metal atom in *complex_mol*
        ligand_mol      : ligand to attach
        coord_atom_idx  : 1-based donor atom index in *ligand_mol*
        target_vec      : direction from metal to donor site (need not be normalised)
        bond_length     : M–donor distance in Å
        """
        metal_idx      = self._i(metal_idx)
        coord_atom_idx = self._i(coord_atom_idx)
        assert metal_idx      < complex_mol.number_of_atoms()
        assert coord_atom_idx < ligand_mol.number_of_atoms()

        c_labels = list(complex_mol.get_labels())
        c_coords = np.array(complex_mol.get_coordinates_in_angstrom(), dtype=float)
        l_labels = list(ligand_mol.get_labels())
        l_coords = np.array(ligand_mol.get_coordinates_in_angstrom(), dtype=float)
        t_unit   = np.array(target_vec, dtype=float) / np.linalg.norm(target_vec)

        print(f"[add_monodentate] ligand={l_labels}, donor={coord_atom_idx+1}, "
              f"target={np.round(t_unit, 2)}, bl={bond_length}")

        # Centre on donor
        l_coords -= l_coords[coord_atom_idx].copy()

        # Body rotation (skip for single-atom ligands)
        other_idx = [i for i in range(len(l_labels)) if i != coord_atom_idx]
        if other_idx:
            body_vec  = np.mean(l_coords[other_idx], axis=0)
            body_unit = body_vec / np.linalg.norm(body_vec)
            rot       = self.get_rotation_matrix(body_unit, t_unit)
            l_coords  = l_coords @ rot.T

        # Translate donor to target position
        l_coords += c_coords[metal_idx] + t_unit * bond_length

        # Torsion scan to minimise steric clash (skip for single-atom ligands)
        if len(l_labels) > 1 and len(c_coords) > 1:
            axis  = t_unit
            pivot = c_coords[metal_idx] + t_unit * bond_length
            best_angle, best_score = 0.0, -1.0
            for deg in np.linspace(0, 360, 360):
                angle = np.radians(deg)
                K     = np.array([[0, -axis[2], axis[1]],
                                  [axis[2], 0, -axis[0]],
                                  [-axis[1], axis[0], 0]])
                R     = np.eye(3) + np.sin(angle)*K + (1-np.cos(angle))*(K@K)
                trial = pivot + (l_coords - pivot) @ R.T
                dists = [np.linalg.norm(trial[i] - c_coords[j])
                         for i in range(len(trial))
                         for j in range(len(c_coords)) if j != metal_idx]
                score = min(dists) if dists else 0.0
                if score > best_score:
                    best_score = score
                    best_angle = angle
            K = np.array([[0, -axis[2], axis[1]],
                          [axis[2], 0, -axis[0]],
                          [-axis[1], axis[0], 0]])
            R = np.eye(3) + np.sin(best_angle)*K + (1-np.cos(best_angle))*(K@K)
            l_coords = pivot + (l_coords - pivot) @ R.T

        return vlx.Molecule(c_labels + l_labels,
                            np.vstack([c_coords, l_coords]), 'angstrom')

    # ------------------------------------------------------------------

    def _find_ring_carbons_2d(
        self, n1_2d, n2_2d, bis_2d, n_c_bond, c_c_bond
    ):
        """
        Solve for the C1/C2 backbone positions in the chelate-ring plane.

        Scans 72 000 angles; returns ``(c1_2d, c2_2d)`` or ``None``.
        """
        def mirror_about(p, axis):
            a = axis / np.linalg.norm(axis)
            return 2 * np.dot(p, a) * a - p

        best, best_err = None, 1e9
        for theta in np.linspace(0, 2*np.pi, 72000):
            v  = np.array([np.cos(theta), np.sin(theta)])
            c1 = n1_2d + n_c_bond * v
            c2 = mirror_about(c1, bis_2d)
            err = (abs(np.linalg.norm(c2 - n2_2d) - n_c_bond)
                 + abs(np.linalg.norm(c2 - c1)    - c_c_bond))
            if err < best_err and np.dot(c1, bis_2d) > np.dot(n1_2d, bis_2d) * 0.5:
                best_err = err
                best = (c1.copy(), c2.copy())
        return best if best_err < 0.02 else None

    def add_universal_bidentate(
        self,
        complex_mol,
        metal_idx: int,
        ligand_mol,
        idx1: int,
        idx2: int,
        target_vec1,
        target_vec2,
        bond_length: float = 2.0,
        bond_length2: float = None,
    ):
        """
        Build an N-C-C-N flexible bidentate chelate (e.g. ethylenediamine) from
        scratch.  The backbone geometry is solved analytically to close a
        5-membered ring; sp³ H atoms are placed at correct M-N-H angles (~110°).

        Use this for flexible, non-aromatic bidentates where the backbone can
        freely adopt the required conformation.  For rigid/aromatic bidentates
        (bipy, phen, ppy…) use :meth:`add_bidentate_rigid` instead.

        Parameters
        ----------
        idx1, idx2      : 1-based donor atom indices in *ligand_mol*
        target_vec1/2   : direction vectors from metal to each donor site
        bond_length     : M–donor₁ distance
        bond_length2    : M–donor₂ distance (defaults to *bond_length*)
        """
        metal_idx = self._i(metal_idx)
        idx1      = self._i(idx1)
        idx2      = self._i(idx2)
        assert metal_idx < complex_mol.number_of_atoms()
        if bond_length2 is None:
            bond_length2 = bond_length

        c_labels = list(complex_mol.get_labels())
        c_coords = np.array(complex_mol.get_coordinates_in_angstrom(), dtype=float)
        l_labels = list(ligand_mol.get_labels())
        l_coords = np.array(ligand_mol.get_coordinates_in_angstrom(), dtype=float)

        m_pos = c_coords[metal_idx]
        t1    = np.array(target_vec1, dtype=float); t1 /= np.linalg.norm(t1)
        t2    = np.array(target_vec2, dtype=float); t2 /= np.linalg.norm(t2)
        N1    = m_pos + t1 * bond_length
        N2    = m_pos + t2 * bond_length2

        ring_normal = np.cross(t1, t2)
        if np.linalg.norm(ring_normal) < 1e-9:
            perp = np.array([1.,0.,0.]) if abs(t1[0]) < 0.9 else np.array([0.,1.,0.])
            ring_normal = np.cross(t1, perp)
        ring_normal /= np.linalg.norm(ring_normal)

        e1 = t1.copy()
        e2 = t2 - np.dot(t2, e1)*e1; e2 /= np.linalg.norm(e2)

        bl = bond_length
        n1_2d = np.array([bl, 0.])
        n2_2d = np.array([bl*np.dot(t2, e1), bl*np.dot(t2, e2)])
        bis_2d = n1_2d + n2_2d
        bis_2d = bis_2d / np.linalg.norm(bis_2d) if np.linalg.norm(bis_2d) > 1e-9 \
                 else np.array([0.707, 0.707])

        cr = self._cov_radius
        def bonded_to(idx):
            return [j for j in range(len(l_labels)) if j != idx and
                    np.linalg.norm(l_coords[j]-l_coords[idx])
                    < 1.3*(cr(l_labels[idx])+cr(l_labels[j]))]

        c1_cands = [j for j in bonded_to(idx1) if l_labels[j] == 'C']
        c2_cands = [j for j in bonded_to(idx2) if l_labels[j] == 'C']
        if c1_cands and c2_cands:
            c1i, c2i = c1_cands[0], c2_cands[0]
            n_c_bond = (np.linalg.norm(l_coords[idx1]-l_coords[c1i])
                      + np.linalg.norm(l_coords[idx2]-l_coords[c2i])) / 2.0
            c_c_bond = np.linalg.norm(l_coords[c1i]-l_coords[c2i])
        else:
            n_c_bond, c_c_bond = 1.47, 1.52

        result2d = self._find_ring_carbons_2d(n1_2d, n2_2d, bis_2d, n_c_bond, c_c_bond)
        if result2d is None:
            mid = (n1_2d + n2_2d) / 2.0
            off = np.array([-mid[1], mid[0]]) / np.linalg.norm(mid) * n_c_bond
            c1_2d, c2_2d = mid + off*0.5, mid - off*0.5
        else:
            c1_2d, c2_2d = result2d

        to3d = lambda p: m_pos + p[0]*e1 + p[1]*e2
        C1, C2 = to3d(c1_2d), to3d(c2_2d)

        def nh2_h(n, metal, c_back, rnorm, h=1.01):
            vnm = (metal-n)/np.linalg.norm(metal-n)
            vnc = (c_back-n)/np.linalg.norm(c_back-n)
            bis = vnm+vnc; bis = bis/np.linalg.norm(bis) if np.linalg.norm(bis)>1e-9 \
                  else np.cross(vnm, rnorm)
            bis /= np.linalg.norm(bis)
            ang = np.radians(55)
            return (n + h*(-np.cos(ang)*bis + np.sin(ang)*rnorm),
                    n + h*(-np.cos(ang)*bis - np.sin(ang)*rnorm))

        def ch2_h(c, na, nb, rnorm, h=1.09):
            v1 = (na-c)/np.linalg.norm(na-c); v2 = (nb-c)/np.linalg.norm(nb-c)
            bis = v1+v2; bis = bis/np.linalg.norm(bis) if np.linalg.norm(bis)>1e-9 \
                  else np.cross(v1, rnorm)
            bis /= np.linalg.norm(bis)
            ang = np.radians(55)
            return (c + h*(-np.cos(ang)*bis + np.sin(ang)*rnorm),
                    c + h*(-np.cos(ang)*bis - np.sin(ang)*rnorm))

        hN1a, hN1b = nh2_h(N1, m_pos, C1, ring_normal)
        hN2a, hN2b = nh2_h(N2, m_pos, C2, ring_normal)
        hC1a, hC1b = ch2_h(C1, N1,   C2, ring_normal)
        hC2a, hC2b = ch2_h(C2, N2,   C1, ring_normal)

        new_ll = ['N','C','C','N','H','H','H','H','H','H','H','H']
        new_lc = np.array([N1,C1,C2,N2, hN1a,hN1b,hC1a,hC1b,hC2a,hC2b,hN2a,hN2b])
        return vlx.Molecule(c_labels+new_ll, np.vstack([c_coords, new_lc]), 'angstrom')

    # ------------------------------------------------------------------

    def find_chelate_donors(
        self, ligand_mol, el1: str, el2: str
    ) -> tuple[int, int]:
        """
        Automatically identify the two donor atom indices for a bidentate
        ligand, without relying on fragile SMILES-order assumptions.

        * Symmetric case (``el1 == el2``, e.g. ``'N','N'`` for bipy):
          returns the indices of the first two atoms of that element.

        * Asymmetric case (e.g. ``'N','C'`` for ppy):
          finds the N atom, then identifies the ortho phenyl C whose lone-pair
          vector is on the *same face* as the N lone-pair — this is the
          cyclometallating carbon.

        Returns
        -------
        ``(idx1, idx2)`` — **1-based** indices of the el1 and el2 donor atoms.
        """
        labels = list(ligand_mol.get_labels())
        coords = np.array(ligand_mol.get_coordinates_in_angstrom(), dtype=float)
        G = self._get_molecule_graph(ligand_mol)

        el1_atoms = [i for i, l in enumerate(labels) if l == el1]

        if el1 == el2:
            if len(el1_atoms) < 2:
                raise ValueError(f"Need ≥2 '{el1}' atoms, found {len(el1_atoms)}")
            return el1_atoms[0] + 1, el1_atoms[1] + 1

        if not el1_atoms:
            raise ValueError(f"No '{el1}' atom found in ligand")
        idx_N = el1_atoms[0]

        rings    = list(nx.cycle_basis(G))
        N_ring   = next((r for r in rings if idx_N in r), None)
        if N_ring is None:
            raise ValueError("el1 atom is not in any ring")
        other_rings = [r for r in rings if idx_N not in r]
        if not other_rings:
            raise ValueError("Could not find a second ring for C^N donor detection")
        other_ring = other_rings[0]

        bridging = next(
            ((u, v) for u, v in G.edges()
             if (u in N_ring and v in other_ring)
             or (v in N_ring and u in other_ring)),
            None,
        )
        if bridging is None:
            raise ValueError("No inter-ring bond found")

        C_bridge_ph = bridging[1] if bridging[0] in N_ring else bridging[0]
        ortho = [n for n in G.neighbors(C_bridge_ph)
                 if n in other_ring and labels[n] == el2]
        if not ortho:
            raise ValueError("No ortho el2 candidates found")

        lp_N = self._lone_pair_and_nbonds(coords, labels, idx_N)[0]
        best, best_dot = ortho[0], -2.0
        for c in ortho:
            dot = np.dot(lp_N, self._lone_pair_and_nbonds(coords, labels, c)[0])
            if dot > best_dot:
                best_dot = dot; best = c

        return idx_N + 1, best + 1

    # ------------------------------------------------------------------

    def add_bidentate_ligand(
        self,
        complex_mol,
        metal_idx: int,
        ligand_mol,
        idx1: int,
        idx2: int,
        target_vec1,
        target_vec2,
        bond_length: float = 2.0,
        bond_length2: float = None,
    ):
        """
        Place a bidentate ligand whose geometry is already correct for chelation.

        The ligand is treated as a rigid body — no bonds are broken.  Use this
        when the input ligand is already in the right conformation (e.g. a
        pre-built or pre-rotated flat cis-bipy).  For ligands fresh from
        ``read_smiles`` that need conformation fixing, use
        :meth:`add_bidentate_rigid` instead.

        Algorithm
        ---------
        1. Two-point rigid superposition: align the donor–donor axis to
           ``(target1 − target2)`` and translate the donor midpoint to the
           target midpoint.  Ligand internal geometry is fully preserved.
        2. Rotate around the donor–donor axis (3600 steps) to maximise the
           sum of LP·(donor→M) for both donors, placing the coordinating face
           of the ligand toward the metal.

        Parameters
        ----------
        complex_mol     : current complex
        metal_idx       : 1-based index of the metal atom in *complex_mol*
        ligand_mol      : bidentate ligand in chelating conformation
        idx1, idx2      : 1-based donor atom indices in *ligand_mol*
        target_vec1/2   : direction vectors from metal to each donor site
        bond_length     : M–donor₁ distance in Å
        bond_length2    : M–donor₂ distance (defaults to *bond_length*)
        """
        metal_idx = self._i(metal_idx)
        idx1      = self._i(idx1)
        idx2      = self._i(idx2)
        assert metal_idx < complex_mol.number_of_atoms()
        if bond_length2 is None:
            bond_length2 = bond_length

        c_labels = list(complex_mol.get_labels())
        c_coords = np.array(complex_mol.get_coordinates_in_angstrom(), dtype=float)
        l_labels = list(ligand_mol.get_labels())
        l_coords = np.array(ligand_mol.get_coordinates_in_angstrom(), dtype=float)
        m_pos    = c_coords[metal_idx]

        t1 = np.array(target_vec1, dtype=float); t1 /= np.linalg.norm(t1)
        t2 = np.array(target_vec2, dtype=float); t2 /= np.linalg.norm(t2)
        target1 = m_pos + t1 * bond_length
        target2 = m_pos + t2 * bond_length2

        # Step 1: two-point rigid superposition
        donor_mid  = (l_coords[idx1] + l_coords[idx2]) / 2.0
        l_c        = l_coords - donor_mid
        R_align    = self.get_rotation_matrix(l_c[idx1] - l_c[idx2], target1 - target2)
        l_c        = l_c @ R_align.T
        l_c       += (target1 + target2) / 2.0

        # Step 2: rotate around donor–donor axis to face LP toward metal
        nn_axis = (target1 - target2) / np.linalg.norm(target1 - target2)
        pivot   = (target1 + target2) / 2.0
        best_angle, best_score = 0.0, -1e9
        for deg in np.linspace(0.0, 360.0, 3600, endpoint=False):
            angle = np.radians(deg)
            K     = np.array([[0, -nn_axis[2], nn_axis[1]],
                               [nn_axis[2], 0, -nn_axis[0]],
                               [-nn_axis[1], nn_axis[0], 0]])
            R_try = np.eye(3) + np.sin(angle)*K + (1-np.cos(angle))*(K@K)
            trial = (l_c - pivot) @ R_try.T + pivot
            lp1   = self._lone_pair_and_nbonds(trial, l_labels, idx1)[0]
            lp2   = self._lone_pair_and_nbonds(trial, l_labels, idx2)[0]
            v1    = (m_pos - trial[idx1]) / np.linalg.norm(m_pos - trial[idx1])
            v2    = (m_pos - trial[idx2]) / np.linalg.norm(m_pos - trial[idx2])
            score = np.dot(lp1, v1) + np.dot(lp2, v2)
            if score > best_score:
                best_score = score; best_angle = angle

        K      = np.array([[0, -nn_axis[2], nn_axis[1]],
                            [nn_axis[2], 0, -nn_axis[0]],
                            [-nn_axis[1], nn_axis[0], 0]])
        R_best = np.eye(3) + np.sin(best_angle)*K + (1-np.cos(best_angle))*(K@K)
        l_c    = (l_c - pivot) @ R_best.T + pivot

        return vlx.Molecule(c_labels + l_labels,
                            np.vstack([c_coords, l_c]), 'angstrom')

    # ------------------------------------------------------------------

    def add_bidentate_rigid(
        self,
        complex_mol,
        metal_idx: int,
        ligand_mol,
        idx1: int,
        idx2: int,
        target_vec1,
        target_vec2,
        bond_length: float = 2.0,
        bond_length2: float = None,
    ):
        """
        Place a rigid/aromatic bidentate ligand (bipy, phen, ppy, acac…).

        The ligand is split at the bond between the two donor-containing
        fragments (e.g. the inter-ring C–C bond in bipy).  Each fragment is
        placed independently via body-centroid rotation — identical to
        :meth:`add_monodentate_ligand` — so each donor lands exactly at its
        target position with no backbone atoms unphysically close to the metal.
        The inter-fragment bond will be slightly stretched; the QM optimiser
        restores it.

        For flexible N-C-C-N donors (en) use :meth:`add_universal_bidentate`.

        Parameters
        ----------
        idx1, idx2      : 1-based donor atom indices in *ligand_mol*
        target_vec1/2   : direction vectors from metal to each donor site
        bond_length     : M–donor₁ distance
        bond_length2    : M–donor₂ distance (defaults to *bond_length*)
        """
        metal_idx = self._i(metal_idx)
        idx1      = self._i(idx1)
        idx2      = self._i(idx2)
        assert metal_idx < complex_mol.number_of_atoms()
        if bond_length2 is None:
            bond_length2 = bond_length

        c_labels = list(complex_mol.get_labels())
        c_coords = np.array(complex_mol.get_coordinates_in_angstrom(), dtype=float)
        l_labels = list(ligand_mol.get_labels())
        l_coords = np.array(ligand_mol.get_coordinates_in_angstrom(), dtype=float)
        n_lig    = len(l_labels)

        m_pos   = c_coords[metal_idx]
        t1      = np.array(target_vec1, dtype=float); t1 /= np.linalg.norm(t1)
        t2      = np.array(target_vec2, dtype=float); t2 /= np.linalg.norm(t2)
        target1 = m_pos + t1 * bond_length
        target2 = m_pos + t2 * bond_length2

        # Build ligand connectivity and split at the inter-donor bond
        G = nx.Graph()
        for i in range(n_lig): G.add_node(i)
        for i in range(n_lig):
            for j in range(i+1, n_lig):
                r = self._cov_radius(l_labels[i]) + self._cov_radius(l_labels[j])
                if np.linalg.norm(l_coords[i]-l_coords[j]) < 1.3*r:
                    G.add_edge(i, j)

        path = nx.shortest_path(G, idx1, idx2)
        mid  = len(path) // 2
        G.remove_edge(path[mid-1], path[mid])
        frag1 = set(nx.node_connected_component(G, idx1))
        frag2 = set(nx.node_connected_component(G, idx2))

        def place_fragment(frag_indices, donor_idx, target_pos, target_unit):
            frag      = sorted(frag_indices)
            f_coords  = l_coords[frag].copy()
            loc_donor = frag.index(donor_idx)
            f_coords -= f_coords[loc_donor]
            others    = [i for i in range(len(frag)) if i != loc_donor]
            if others:
                body = np.mean(f_coords[others], axis=0)
                body /= np.linalg.norm(body)
                f_coords = f_coords @ self.get_rotation_matrix(body, target_unit).T
            f_coords += target_pos
            return frag, f_coords

        frag1_idx, frag1_coords = place_fragment(frag1, idx1, target1, t1)
        frag2_idx, frag2_coords = place_fragment(frag2, idx2, target2, t2)

        l_final = l_coords.copy()
        for li, gi in enumerate(frag1_idx): l_final[gi] = frag1_coords[li]
        for li, gi in enumerate(frag2_idx): l_final[gi] = frag2_coords[li]

        return vlx.Molecule(c_labels + l_labels,
                            np.vstack([c_coords, l_final]), 'angstrom')

    # ------------------------------------------------------------------

    def prepare_bidentate_conformation(
        self,
        ligand_mol,
        idx1: int,
        idx2: int,
        n_steps: int = 120,
    ):
        """
        Rotate the backbone between the two donor atoms to maximise LP
        alignment (both lone-pairs pointing toward the same side of space),
        as required for chelation.

        Uses a fully-vectorised grid scan over all rotatable bonds along the
        shortest donor–donor path; no Python loops over angle combinations.

        Parameters
        ----------
        n_steps : angular resolution per bond (default 120 → 3°)
        """
        labels = list(ligand_mol.get_labels())
        coords = np.array(ligand_mol.get_coordinates_in_angstrom(), dtype=float)
        G      = self._get_molecule_graph(ligand_mol)
        idx1   = self._i(idx1)
        idx2   = self._i(idx2)

        path = nx.shortest_path(G, idx1, idx2)
        rotatable = [(a, b) for a, b in zip(path[:-1], path[1:])
                     if not nx.has_path(G.copy().remove_edge(a, b) or G, a, b)]
        # Simpler check: bond is rotatable if removing it disconnects the graph
        rotatable = []
        for a, b in zip(path[:-1], path[1:]):
            Gt = G.copy(); Gt.remove_edge(a, b)
            if not nx.has_path(Gt, a, b):
                rotatable.append((a, b))

        if not rotatable:
            return ligand_mol

        angles     = np.linspace(0.0, 2*np.pi, n_steps, endpoint=False)
        lp1        = self._lone_pair_and_nbonds(coords, labels, idx1)[0]
        lp2_init   = self._lone_pair_and_nbonds(coords, labels, idx2)[0]
        axes_init  = [(coords[b]-coords[a])/np.linalg.norm(coords[b]-coords[a])
                      for a, b in rotatable]

        def batch_R(axes, angs):
            N, M = axes.shape[0], angs.shape[0]
            K = np.zeros((N, 3, 3))
            K[:,0,1]=-axes[:,2]; K[:,0,2]= axes[:,1]
            K[:,1,0]= axes[:,2]; K[:,1,2]=-axes[:,0]
            K[:,2,0]=-axes[:,1]; K[:,2,1]= axes[:,0]
            K2 = np.einsum('nij,njk->nik', K, K)
            s, c = np.sin(angs), np.cos(angs)
            return (np.eye(3)[None,None]
                    + s[None,:,None,None]*K[:,None]
                    + (1-c)[None,:,None,None]*K2[:,None])

        n_combos = 1
        cumR     = np.eye(3)[None]
        for k in range(len(rotatable)):
            axis_k = np.einsum('nij,j->ni', cumR, axes_init[k])
            R_k    = batch_R(axis_k, angles)
            new_cumR = np.einsum('nmij,nmjk->nmik',
                                 R_k, cumR[:,None].repeat(len(angles), axis=1))
            n_combos *= len(angles)
            cumR = new_cumR.reshape(n_combos, 3, 3)

        lp2_all = np.einsum('nij,j->ni', cumR, lp2_init)
        lp2_all /= np.maximum(np.linalg.norm(lp2_all, axis=1, keepdims=True), 1e-10)
        scores   = lp2_all @ lp1
        best_i   = int(np.argmax(scores))

        best_indices = []
        rem = best_i
        for _ in rotatable:
            best_indices.append(rem % len(angles)); rem //= len(angles)
        best_indices.reverse()

        result = coords.copy()
        for k, (a, b) in enumerate(rotatable):
            axis  = result[b]-result[a]; axis /= np.linalg.norm(axis)
            pivot = result[a].copy()
            angle = angles[best_indices[k]]
            K     = np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]])
            R     = np.eye(3) + np.sin(angle)*K + (1-np.cos(angle))*(K@K)
            Gt    = G.copy(); Gt.remove_edge(a, b)
            moving = list(nx.node_connected_component(Gt, b))
            result[moving] = (result[moving]-pivot) @ R.T + pivot

        return vlx.Molecule(labels, result, 'angstrom')

    # ------------------------------------------------------------------

    def remove_H_from_donor(self, ligand_mol, donor_idx: int):
        """
        Remove the H bonded to *donor_idx* (deprotonate a cyclometallating C
        or any other donor that loses a proton on coordination).

        Raises ``ValueError`` if no H is found on that atom.
        """
        labels    = list(ligand_mol.get_labels())
        coords    = np.array(ligand_mol.get_coordinates_in_angstrom(), dtype=float)
        donor_idx = self._i(donor_idx)
        r_d       = self._cov_radius(labels[donor_idx])
        h_idx  = next(
            (i for i, lbl in enumerate(labels)
             if i != donor_idx and lbl == 'H'
             and np.linalg.norm(coords[i]-coords[donor_idx])
                 < 1.3*(r_d + self._cov_radius('H'))),
            None,
        )
        if h_idx is None:
            raise ValueError(f"No H found on atom {donor_idx} ({labels[donor_idx]})")
        keep = [i for i in range(len(labels)) if i != h_idx]
        return vlx.Molecule([labels[i] for i in keep], coords[keep], 'angstrom')

    def remove_ligand(self, complex_mol, metal_idx: int, donor_idx: int):
        """
        Remove a ligand by severing the M–donor bond and deleting every atom
        in the ligand fragment.
        """
        metal_idx  = self._i(metal_idx)
        donor_idx  = self._i(donor_idx)
        assert metal_idx  < complex_mol.number_of_atoms()
        assert donor_idx  < complex_mol.number_of_atoms()
        labels = list(complex_mol.get_labels())
        coords = np.array(complex_mol.get_coordinates_in_angstrom(), dtype=float)
        G      = self._get_molecule_graph(complex_mol)
        if not G.has_edge(metal_idx, donor_idx):
            raise ValueError(f"Atoms {metal_idx} and {donor_idx} are not bonded")
        Gc = G.copy(); Gc.remove_edge(metal_idx, donor_idx)
        lig = set(nx.node_connected_component(Gc, donor_idx))
        keep = [i for i in range(len(labels)) if i not in lig]
        return vlx.Molecule([labels[i] for i in keep], coords[keep], 'angstrom')

    def transmetalate(self, vlx_mol, metal_idx: int, new_metal_symbol: str):
        """
        Replace the metal at *metal_idx* with *new_metal_symbol*, scaling all
        M–L bond lengths to match the new covalent radius.
        """
        metal_idx = self._i(metal_idx)
        assert metal_idx < vlx_mol.number_of_atoms()
        labels = list(vlx_mol.get_labels())
        coords = np.array(vlx_mol.get_coordinates_in_angstrom(), dtype=float)
        G      = self._get_molecule_graph(vlx_mol)
        delta  = self._cov_radius(new_metal_symbol) - self._cov_radius(labels[metal_idx])
        labels[metal_idx] = new_metal_symbol
        metal_pos = coords[metal_idx]
        shifts  = np.zeros_like(coords)
        shifted = set()
        for nb in G.neighbors(metal_idx):
            Gc = G.copy(); Gc.remove_edge(metal_idx, nb)
            u  = coords[nb]-metal_pos; u /= np.linalg.norm(u)
            for i in nx.node_connected_component(Gc, nb):
                if i not in shifted:
                    shifts[i] = u * delta; shifted.add(i)
        coords += shifts
        return vlx.Molecule(labels, coords, 'angstrom')

    # ==================================================================
    # ── Section 3: Organic chemistry / fragment coupling ──────────────
    # ==================================================================

    # ------------------------------------------------------------------
    # Internal helpers shared by organic methods
    # ------------------------------------------------------------------

    def _smiles_to_nx_with_hs(self, smiles: str) -> tuple[nx.Graph, int]:
        """
        Parse *smiles* (which may contain a ``*`` wildcard anchor) into a
        NetworkX graph with explicit H atoms.  Returns ``(G, star_idx)`` where
        ``star_idx`` is the index of the ``*`` node (or ``-1`` if absent).
        """
        mol      = Chem.MolFromSmiles(smiles)
        mol      = Chem.AddHs(mol)
        G        = nx.Graph()
        star_idx = -1
        for atom in mol.GetAtoms():
            sym = atom.GetSymbol(); idx = atom.GetIdx()
            G.add_node(idx, element=sym, is_star=(sym == '*'))
            if sym == '*':
                star_idx = idx
        for bond in mol.GetBonds():
            G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        return G, star_idx

    def _prepare_organic_fragment(self, obj, index=None):
        """
        Return ``(labels, coords, join_idx, exit_vec)`` for an organic fragment.

        *obj* may be:
        * a SMILES string with ``*`` marking the attachment point, or
        * a VeloxChem Molecule (``index`` is the bonding atom).

        One H is removed from the bonding atom; ``exit_vec`` is the vector
        that pointed from the bonding atom toward that removed H.
        """
        if isinstance(obj, str):
            full = obj.replace('*', '[H]')
            mol  = Chem.MolFromSmiles(full); mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            conf   = mol.GetConformer().GetPositions()
            labels = [a.GetSymbol() for a in mol.GetAtoms()]
            # First heavy atom is the join point
            join   = next(i for i, a in enumerate(mol.GetAtoms()) if a.GetSymbol() != 'H')
            h_nb   = next(n.GetIdx() for n in mol.GetAtoms()[join].GetNeighbors()
                          if n.GetSymbol() == 'H')
            exit_vec = conf[h_nb] - conf[join]
            f_lab = [l for i, l in enumerate(labels) if i != h_nb]
            f_coo = np.delete(conf, h_nb, axis=0)
            new_join = join if join < h_nb else join - 1
            return f_lab, f_coo, new_join, exit_vec

        # VeloxChem Molecule path
        labels = list(obj.get_labels())
        coords = np.array(obj.get_coordinates_in_angstrom())
        h_idx  = next(
            (i for i, (l, p) in enumerate(zip(labels, coords))
             if l == 'H' and np.linalg.norm(p-coords[index]) < 1.3),
            -1,
        )
        if h_idx != -1:
            exit_vec = coords[h_idx] - coords[index]
            labels.pop(h_idx); coords = np.delete(coords, h_idx, axis=0)
            if index > h_idx: index -= 1
        else:
            nbs      = [coords[i] for i in range(len(coords))
                        if i != index and np.linalg.norm(coords[i]-coords[index]) < 1.8]
            exit_vec = -np.mean([n-coords[index] for n in nbs], axis=0) if nbs \
                       else np.array([0., 0., 1.])
        return labels, coords, index, exit_vec

    def _clash_torsion_scan(
        self,
        coords_a: np.ndarray, labels_a: list,
        coords_b: np.ndarray, labels_b: list,
        pivot: np.ndarray, axis: np.ndarray,
        step_deg: int = 10,
    ) -> np.ndarray:
        """
        Scan *step_deg*-degree torsion steps around *axis* at *pivot* for
        *coords_b*, minimising steric clash against *coords_a*.
        Returns the best-rotated *coords_b*.
        """
        r_a = np.array([self._vdw_radius(l) for l in labels_a])
        r_b = np.array([self._vdw_radius(l) for l in labels_b])
        best, min_clash = coords_b.copy(), float('inf')
        for deg in range(0, 360, step_deg):
            trial = np.array([
                self.rotate_vector(p - pivot, axis, np.radians(deg)) + pivot
                for p in coords_b
            ])
            dists  = np.linalg.norm(coords_a[:,None,:] - trial[None,:,:], axis=2)
            clash  = np.sum(np.maximum(0, (r_a[:,None]+r_b[None,:]) - dists)**2)
            if clash < min_clash:
                min_clash = clash; best = trial
        return best

    # ------------------------------------------------------------------
    # Public organic API
    # ------------------------------------------------------------------

    def link_at_vector(
        self,
        mol_a,
        smiles_b: str,
        index_a: int,
        exit_vec_a: np.ndarray,
        bond_length: float,
    ):
        """
        Attach SMILES fragment *smiles_b* to *mol_a* at atom *index_a*, with
        the new bond pointing along *exit_vec_a*.

        One H is removed from the ``*`` atom of the fragment; the resulting
        radical is rotated to face *mol_a*, translated to *bond_length*, then
        optimised for steric clash by a 10° torsion scan.

        Returns a new VeloxChem Molecule.
        """
        l_b, c_b, j_b, v_b = self._prepare_organic_fragment(smiles_b)
        l_a = list(mol_a.get_labels())
        c_a = np.array(mol_a.get_coordinates_in_angstrom())
        index_a = self._i(index_a)

        # Centre mol_a at bonding atom; rotate so exit_vec points along +Z
        c_a -= c_a[index_a]
        c_a  = c_a @ self.get_rotation_matrix(exit_vec_a, [0,0,1.]).T

        # Centre fragment at its join; rotate so v_b points along -Z
        c_b -= c_b[j_b]
        c_b  = c_b @ self.get_rotation_matrix(v_b, [0,0,-1.]).T
        c_b += np.array([0., 0., bond_length])

        pivot = np.array([0., 0., bond_length])
        c_b   = self._clash_torsion_scan(c_a, l_a, c_b, l_b, pivot, np.array([0.,0.,1.]))

        return vlx.Molecule(l_a + l_b, np.concatenate([c_a, c_b], axis=0), 'angstrom')

    def join_molecules(
        self,
        mol_a,
        mol_b,
        idx_a: int,
        idx_b: int,
        bond_length: float = 1.54,
    ):
        """
        Couple two VeloxChem molecules at specific atom indices.

        One H is removed from each bonding atom, the molecules are aligned
        along Z, and a torsion scan minimises steric overlap.

        Parameters
        ----------
        idx_a, idx_b : 1-based bonding atom indices in *mol_a* and *mol_b*
        bond_length  : new bond length in Å (default 1.54 for C–C)
        """
        l_a, c_a, j_a, v_a = self._prepare_organic_fragment(mol_a, self._i(idx_a))
        l_b, c_b, j_b, v_b = self._prepare_organic_fragment(mol_b, self._i(idx_b))

        c_a -= c_a[j_a]
        c_a  = c_a @ self.get_rotation_matrix(v_a, [0,0,1.]).T
        c_b -= c_b[j_b]
        c_b  = c_b @ self.get_rotation_matrix(v_b, [0,0,-1.]).T
        c_b += np.array([0., 0., bond_length])

        pivot = np.array([0., 0., bond_length])
        c_b   = self._clash_torsion_scan(c_a, l_a, c_b, l_b, pivot,
                                          np.array([0.,0.,1.]), step_deg=15)

        return vlx.Molecule(l_a + l_b, np.concatenate([c_a, c_b], axis=0), 'angstrom')

    def substitute(
        self,
        vlx_mol,
        old_smiles: str,
        new_smiles: str,
        bond_length: float = None,
    ) -> dict:
        """
        Replace every occurrence of *old_smiles* in *vlx_mol* with *new_smiles*.

        *old_smiles* should contain a ``*`` atom marking the anchor (the atom
        that remains in the molecule after substitution).  The star atom is
        matched to each heavy atom in the molecule and subgraph isomorphism is
        used to locate the group to remove.

        Returns a dict ``{'molecules': [...], 'radical_indices': [...]}``.
        """
        coords = np.array(vlx_mol.get_coordinates_in_angstrom())
        labels = list(vlx_mol.get_labels())

        target_G = nx.Graph()
        for i, lbl in enumerate(labels):
            target_G.add_node(i, element=lbl)
        for i in range(len(coords)):
            for j in range(i+1, len(coords)):
                if np.linalg.norm(coords[i]-coords[j]) < 1.75:
                    target_G.add_edge(i, j)

        query_G, star_idx = self._smiles_to_nx_with_hs(old_smiles)
        GM = isomorphism.GraphMatcher(
            target_G, query_G,
            node_match=lambda n1, n2: True if n2['is_star'] else n1['element']==n2['element'],
        )

        unique, seen = [], set()
        for mapping in GM.subgraph_isomorphisms_iter():
            inv        = {v: k for k, v in mapping.items()}
            anchor_idx = inv[star_idx]
            to_remove  = set(mapping.keys()) - {anchor_idx}
            key        = tuple(sorted(to_remove))
            if key not in seen:
                unique.append((to_remove, anchor_idx)); seen.add(key)

        results = {'molecules': [], 'radical_indices': []}
        for to_remove, anchor_idx in unique:
            nbs = [n for n in target_G.neighbors(anchor_idx) if n in to_remove]
            if not nbs: continue
            exit_vec = coords[nbs[0]] - coords[anchor_idx]

            new_labels, new_coords, old_to_new = [], [], {}
            curr = 0
            for i in range(len(labels)):
                if i not in to_remove:
                    new_labels.append(labels[i]); new_coords.append(coords[i])
                    old_to_new[i] = curr; curr += 1

            pruned    = vlx.Molecule(new_labels, np.array(new_coords), 'angstrom')
            new_anch  = old_to_new[anchor_idx]
            tmp_lab, _, _, _ = self._prepare_organic_fragment(new_smiles)
            dist      = bond_length or self._bond_defaults.get(tmp_lab[0], 1.54)
            final_mol = self.link_at_vector(pruned, new_smiles, new_anch, exit_vec, dist)
            results['molecules'].append(final_mol)
            results['radical_indices'].append(new_anch)

        return results

    def remove_all_fragment_matches(
        self,
        vlx_mol,
        fragment_smiles: str,
    ) -> dict:
        """
        Strip every terminal occurrence of *fragment_smiles* from *vlx_mol*.

        ``*`` in the SMILES marks the anchor atom (kept in the molecule).
        Only fragments connected by exactly one bond are removed (terminal
        groups).  Returns ``{'molecules': [...], 'radical_indices': [...]}``.
        """
        coords = np.array(vlx_mol.get_coordinates_in_angstrom())
        labels = list(vlx_mol.get_labels())

        target_G = nx.Graph()
        for i, lbl in enumerate(labels):
            target_G.add_node(i, element=lbl, mol_idx=i)
        for i in range(len(coords)):
            for j in range(i+1, len(coords)):
                r1 = self._cov_radius(labels[i])
                r2 = self._cov_radius(labels[j])
                if np.linalg.norm(coords[i]-coords[j]) < 1.3*(r1+r2):
                    target_G.add_edge(i, j)

        query_G, star_idx = self._smiles_to_nx_with_hs(fragment_smiles)

        unique, seen = [], set()
        for candidate in range(len(labels)):
            if labels[candidate] == 'H': continue

            def node_match(n1, n2):
                if n2.get('is_star'): return n1.get('mol_idx') == candidate
                return n1['element'] == n2['element']

            GM = isomorphism.GraphMatcher(target_G, query_G, node_match=node_match)
            for mapping in GM.subgraph_isomorphisms_iter():
                inv       = {v: k for k, v in mapping.items()}
                anchor    = inv[star_idx]
                to_remove = set(mapping.keys()) - {anchor}
                key       = tuple(sorted(to_remove))
                if key in seen: break
                cut = sum(1 for u in to_remove
                          for v in target_G.neighbors(u) if v not in to_remove)
                if cut == 1:
                    unique.append((to_remove, anchor))
                    seen.add(key)
                    print(f"[remove_fragment] found {fragment_smiles} at anchor {anchor}")
                break  # one match per candidate

        output = {'molecules': [], 'radical_indices': []}
        for to_remove, rad_idx in unique:
            new_labels, new_coords, old_to_new = [], [], {}
            curr = 0
            for i in range(len(labels)):
                if i not in to_remove:
                    new_labels.append(labels[i]); new_coords.append(coords[i])
                    old_to_new[i] = curr; curr += 1
            output['molecules'].append(
                vlx.Molecule(new_labels, np.array(new_coords), 'angstrom'))
            output['radical_indices'].append(old_to_new[rad_idx])

        return output

    def identify_aromatic_carbons_gaff(self, vlx_mol) -> tuple[list, dict]:
        """
        Identify aromatic carbon sites in *vlx_mol* using GAFF atom types.

        Returns ``(ca_indices, gaff_dict)`` where *ca_indices* is a list of
        ``(idx, gaff_type)`` tuples for aromatic-carbon atoms.
        Indices are **1-based** (VeloxChem convention).
        """
        at_id = vlx.AtomTypeIdentifier()
        at_id.generate_gaff_atomtypes(vlx_mol)
        gaff_dict = at_id.atom_types_dict
        ca_types  = {'ca', 'cc', 'cd', 'cp', 'cq', 'ce', 'cf'}
        ca_indices = []
        for key, val in gaff_dict.items():
            if val.get('gaff') in ca_types:
                # GAFF keys are 1-based (e.g. "C1", "C2"); keep as 1-based for the user
                idx = int(''.join(filter(str.isdigit, key)))
                ca_indices.append((idx, val.get('gaff')))
        return ca_indices, gaff_dict

    def add_smiles_radical(
        self,
        vlx_mol,
        target_idx: int,
        smiles_with_star: str,
        gaff_dict: dict,
        bond_length: float = 1.45,
    ):
        """
        Add a SMILES radical fragment to the aromatic carbon *target_idx*.

        The radical is oriented perpendicular to the ring plane at *target_idx*
        and steric clash is minimised by a torsion scan.  Existing substituents
        on *target_idx* are deflected to accommodate the sp³ geometry.

        Returns a new VeloxChem Molecule, or ``None`` on failure.
        """
        coords = np.array(vlx_mol.get_coordinates_in_angstrom())
        labels = list(vlx_mol.get_labels())
        target_idx = self._i(target_idx)
        tpos   = coords[target_idx]

        nbs     = [i for i, p in enumerate(coords)
                   if i != target_idx and np.linalg.norm(p-tpos) < 1.7]
        ring_nbs, sub_idx = [], -1
        for ni in nbs:
            g = gaff_dict.get(f"{labels[ni]}{ni+1}", {}).get('gaff', '')
            if (g.startswith('c') or g.startswith('n')) and len(ring_nbs) < 2:
                ring_nbs.append(ni)
            else:
                sub_idx = ni
        if len(ring_nbs) < 2:
            return None

        v1 = coords[ring_nbs[0]]-tpos; v2 = coords[ring_nbs[1]]-tpos
        rnorm = np.cross(v1, v2); rnorm /= np.linalg.norm(rnorm)
        if sub_idx != -1:
            vsub = (coords[sub_idx]-tpos)/np.linalg.norm(coords[sub_idx]-tpos)
        else:
            vsub = -(v1/np.linalg.norm(v1)+v2/np.linalg.norm(v2))
            vsub /= np.linalg.norm(vsub)

        bend_axis  = np.cross(rnorm, vsub)
        v_bond_out = self.rotate_vector(vsub, bend_axis, np.radians(-54.75))

        coords_copy = coords.copy()
        if sub_idx != -1:
            new_vsub = self.rotate_vector(vsub, bend_axis, np.radians(54.75))
            rot_mat  = self.get_rotation_matrix(vsub, new_vsub)
            branch   = [sub_idx]; added = True
            while added:
                added = False
                for i, p in enumerate(coords):
                    if i in branch or i==target_idx or i in ring_nbs: continue
                    if any(np.linalg.norm(p-coords[b])<1.7 for b in branch):
                        branch.append(i); added = True
            for i in branch:
                coords_copy[i] = tpos + rot_mat @ (coords[i]-tpos)

        frag = Chem.MolFromSmiles(smiles_with_star.replace('*', '[Xe]'))
        frag = Chem.AddHs(frag); AllChem.EmbedMolecule(frag, AllChem.ETKDG())
        fc   = frag.GetConformer().GetPositions()
        xe   = [a.GetIdx() for a in frag.GetAtoms() if a.GetSymbol()=='Xe'][0]
        rad  = frag.GetAtoms()[xe].GetNeighbors()[0].GetIdx()
        fc  -= fc[rad]
        fc   = fc @ self.get_rotation_matrix(fc[xe]/np.linalg.norm(fc[xe]), -v_bond_out).T
        f_labels = [a.GetSymbol() for i, a in enumerate(frag.GetAtoms()) if i != xe]
        f_coords = np.delete(fc, xe, axis=0)

        f_coords = self._clash_torsion_scan(
            coords_copy, labels, f_coords, f_labels,
            pivot=np.zeros(3), axis=v_bond_out,
        )
        f_coords += tpos + v_bond_out * bond_length

        return vlx.Molecule(list(labels)+f_labels,
                            np.concatenate([coords_copy, f_coords], axis=0), 'angstrom')

    def run_radical_addition(
        self,
        molecule,
        smiles_fragment: str,
    ) -> list:
        """
        Add *smiles_fragment* to every aromatic carbon site in *molecule*.

        Returns a list of VeloxChem Molecules (one per site), skipping any
        sites where the addition fails.
        """
        ca_sites, gaff_info = self.identify_aromatic_carbons_gaff(molecule)
        return [iso for idx, gtype in ca_sites
                for iso in [self.add_smiles_radical(molecule, idx, smiles_fragment, gaff_info)]
                if iso is not None]
