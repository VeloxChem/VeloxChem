from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Tuple

import numpy as np
from mpi4py import MPI


@dataclass(frozen=True)
class RankLayout:
    counts: np.ndarray
    displs: np.ndarray

    @staticmethod
    def from_n_items(n_items: int, size: int) -> "RankLayout":
        ave, rem = divmod(int(n_items), int(size))
        counts = np.array([ave + 1 if p < rem else ave for p in range(size)], dtype=np.int64)
        displs = np.zeros(size, dtype=np.int64)
        if size > 1:
            displs[1:] = np.cumsum(counts[:-1], dtype=np.int64)
        return RankLayout(counts=counts, displs=displs)

    def local_slice(self, rank: int) -> slice:
        start = int(self.displs[rank])
        stop = start + int(self.counts[rank])
        return slice(start, stop)


@dataclass(frozen=True)
class StepPacket:
    org_int_coords: np.ndarray
    b_matrix: np.ndarray
    use_cosine_dihedral: bool


@dataclass(frozen=True)
class SymmetryTask:
    task_id: int
    sym_idx: int
    mask0: np.ndarray
    mask: np.ndarray


@dataclass
class StaticLabelPack:
    label: str
    natm: int
    n_ic: int
    energies: np.ndarray
    int_grad: np.ndarray
    int_hess: np.ndarray
    int_coords: np.ndarray
    tasks: List[SymmetryTask]
    layout: RankLayout
    bond_end: int
    angle_end: int
    dihedral_start: int
    dihedral_end: int
    sym3_rows: np.ndarray
    sym3_center_ids: np.ndarray
    n_sym3: int
    r_cut: float


@dataclass
class LocalLabelPack:
    global_pack: StaticLabelPack
    local_tasks: List[SymmetryTask]


@dataclass
class LabelReduction:
    S: float
    N: float
    dS: np.ndarray
    dN: np.ndarray
    rmsd_bond_sum: float
    rmsd_angle_sum: float
    rmsd_dihedral_sum: float
    rmsd_count: int

    @staticmethod
    def zeros(natm: int) -> "LabelReduction":
        return LabelReduction(
            S=0.0,
            N=0.0,
            dS=np.zeros((natm, 3), dtype=np.float64),
            dN=np.zeros((natm, 3), dtype=np.float64),
            rmsd_bond_sum=0.0,
            rmsd_angle_sum=0.0,
            rmsd_dihedral_sum=0.0,
            rmsd_count=0,
        )


def _allreduce_reduction(comm, local_red: LabelReduction) -> LabelReduction:
    return LabelReduction(
        S=comm.allreduce(local_red.S, op=MPI.SUM),
        N=comm.allreduce(local_red.N, op=MPI.SUM),
        dS=comm.allreduce(local_red.dS, op=MPI.SUM),
        dN=comm.allreduce(local_red.dN, op=MPI.SUM),
        rmsd_bond_sum=comm.allreduce(local_red.rmsd_bond_sum, op=MPI.SUM),
        rmsd_angle_sum=comm.allreduce(local_red.rmsd_angle_sum, op=MPI.SUM),
        rmsd_dihedral_sum=comm.allreduce(local_red.rmsd_dihedral_sum, op=MPI.SUM),
        rmsd_count=int(comm.allreduce(local_red.rmsd_count, op=MPI.SUM)),
    )


def _finalize_reduction(global_red: LabelReduction) -> Tuple[float, np.ndarray, Dict[str, float]]:
    if global_red.S <= 0.0:
        raise ValueError("MPI preload reduction produced non-positive total weight.")

    energy = global_red.N / global_red.S
    gradient = (global_red.dN * global_red.S - global_red.N * global_red.dS) / (global_red.S ** 2)

    if global_red.rmsd_count > 0:
        rmsd_diag = {
            'bond_mean': float(global_red.rmsd_bond_sum / global_red.rmsd_count),
            'angle_mean': float(global_red.rmsd_angle_sum / global_red.rmsd_count),
            'dihedral_mean': float(global_red.rmsd_dihedral_sum / global_red.rmsd_count),
            'count': float(global_red.rmsd_count),
        }
    else:
        rmsd_diag = {
            'bond_mean': 0.0,
            'angle_mean': 0.0,
            'dihedral_mean': 0.0,
            'count': 0.0,
        }

    return float(energy), gradient, rmsd_diag


def _grouped_signature_distance_and_gradient(
        grouped_current: List[np.ndarray],
        grouped_candidate: List[np.ndarray],
        grouped_brows: List[np.ndarray],
        ncart: int) -> Tuple[float, np.ndarray]:
    n_rotors = len(grouped_current)
    if n_rotors == 0:
        return 0.0, np.zeros((ncart,), dtype=np.float64)

    d2_total = 0.0
    grad_total = np.zeros((ncart,), dtype=np.float64)

    for rid in range(n_rotors):
        sig_c = grouped_current[rid]
        sig_r = grouped_candidate[rid]
        brows = grouped_brows[rid]

        delta = sig_c - sig_r
        terms = 2.0 * (1.0 - np.cos(delta))
        d2_r = float(np.mean(terms))

        grad_terms = 2.0 * np.sin(delta)[:, None] * brows
        grad_r = np.mean(grad_terms, axis=0)

        d2_total += d2_r
        grad_total += grad_r

    inv = 1.0 / max(float(n_rotors), 1.0)
    return d2_total * inv, grad_total * inv


def _signature_inverse_weight_and_gradient(
        d2: float,
        grad_d2: np.ndarray,
        confidence_radius: float = 1.0,
        power: int = 2,
        eps: float = 1.0e-12) -> Tuple[float, np.ndarray]:
    rho2 = max(confidence_radius * confidence_radius, eps)
    x = d2 / rho2 + eps

    w = 1.0 / (2.0 * (x ** power))
    dw_dd2 = -(power / (2.0 * rho2)) * (x ** (-(power + 1)))
    grad_w = dw_dd2 * grad_d2
    return float(w), np.asarray(grad_w, dtype=np.float64)


def _segment_rmsd(segment: np.ndarray) -> float:
    if segment.size == 0:
        return 0.0
    return float(np.sqrt(np.mean(np.sum(segment ** 2))))


def evaluate_candidate_taylor(
        *,
        energy: float,
        grad: np.ndarray,
        hess: np.ndarray,
        ref_coords: np.ndarray,
        mask0: np.ndarray,
        mask: np.ndarray,
        org_int_coords: np.ndarray,
        b_matrix: np.ndarray,
        use_cosine_dihedral: bool,
        dihedral_start: int,
        dihedral_end: int) -> Tuple[float, np.ndarray, np.ndarray, np.ndarray]:
    """
    Evaluates one masked Taylor candidate in a parity-safe way shared by
    serial and MPI preload paths.
    """
    grad_eff = np.asarray(grad, dtype=np.float64).copy()
    hess_eff = np.asarray(hess, dtype=np.float64).copy()
    ref = np.asarray(ref_coords, dtype=np.float64)
    mask0_arr = np.asarray(mask0, dtype=np.int64)
    mask_arr = np.asarray(mask, dtype=np.int64)
    org = np.asarray(org_int_coords, dtype=np.float64)
    bmat = np.asarray(b_matrix, dtype=np.float64)

    grad_eff[mask0_arr] = grad_eff[mask_arr]
    hess_eff[np.ix_(mask0_arr, mask0_arr)] = hess_eff[np.ix_(mask_arr, mask_arr)]

    dist_org = org - ref[mask_arr]
    dist_check = dist_org.copy()
    if not use_cosine_dihedral:
        dist_check[dihedral_start:dihedral_end] = np.sin(
            dist_org[dihedral_start:dihedral_end])

    dist_correlation = dist_check.copy()

    U_i = float(
        energy
        + dist_check @ grad_eff
        + 0.5 * np.linalg.multi_dot([dist_check, hess_eff, dist_check])
    )

    dist_hess = dist_check @ hess_eff
    if not use_cosine_dihedral:
        c = np.cos(dist_org[dihedral_start:dihedral_end])
        grad_eff[dihedral_start:dihedral_end] *= c
        dist_hess[dihedral_start:dihedral_end] *= c

    natm = bmat.shape[1] // 3
    dU_i = (bmat.T @ (grad_eff + dist_hess)).reshape(natm, 3)

    return U_i, dU_i, dist_org, dist_correlation


def _build_pack_for_label(driver: Any, dp_label: str, size: int) -> StaticLabelPack:
    entries = driver._build_candidates_for_label(dp_label)
    if len(entries) == 0:
        raise ValueError(f"No candidate entries found for label '{dp_label}'.")

    uniq = []
    idx_by_obj = {}
    for sym_dp, _, _ in entries:
        key = id(sym_dp)
        if key not in idx_by_obj:
            idx_by_obj[key] = len(uniq)
            uniq.append(sym_dp)

    n_sym = len(uniq)
    n_ic = len(uniq[0].internal_coordinates_values)
    natm = uniq[0].cartesian_coordinates.shape[0]

    energies = np.empty((n_sym,), dtype=np.float64)
    int_grad = np.empty((n_sym, n_ic), dtype=np.float64)
    int_hess = np.empty((n_sym, n_ic, n_ic), dtype=np.float64)
    int_coords = np.empty((n_sym, n_ic), dtype=np.float64)

    for i, dp in enumerate(uniq):
        energies[i] = float(dp.energy)
        int_grad[i] = np.asarray(dp.internal_gradient, dtype=np.float64)
        int_hess[i] = np.asarray(dp.internal_hessian, dtype=np.float64)
        int_coords[i] = np.asarray(dp.internal_coordinates_values, dtype=np.float64)

    tasks: List[SymmetryTask] = []
    for task_id, (sym_dp, mask0, mask) in enumerate(entries):
        tasks.append(
            SymmetryTask(
                task_id=task_id,
                sym_idx=idx_by_obj[id(sym_dp)],
                mask0=np.ascontiguousarray(mask0, dtype=np.int64),
                mask=np.ascontiguousarray(mask, dtype=np.int64),
            )
        )

    use_symmetry = bool(driver._use_symmetry_for_label(dp_label))
    if use_symmetry:
        torsion_meta = driver._get_symmetry_torsion_meta()
    else:
        torsion_meta = {
            'sym3_rows': np.array([], dtype=np.int64),
            'sym3_center_ids': np.array([], dtype=np.int64),
            'sym3_keys': (),
        }

    bounds = driver._get_internal_coordinate_partitions()

    layout = RankLayout.from_n_items(len(tasks), size)
    return StaticLabelPack(
        label=dp_label,
        natm=natm,
        n_ic=n_ic,
        energies=energies,
        int_grad=int_grad,
        int_hess=int_hess,
        int_coords=int_coords,
        tasks=tasks,
        layout=layout,
        bond_end=int(bounds['bond_end']),
        angle_end=int(bounds['angle_end']),
        dihedral_start=int(bounds['dihedral_start']),
        dihedral_end=int(bounds['dihedral_end']),
        sym3_rows=np.asarray(torsion_meta['sym3_rows'], dtype=np.int64),
        sym3_center_ids=np.asarray(torsion_meta['sym3_center_ids'], dtype=np.int64),
        n_sym3=int(len(torsion_meta['sym3_keys'])),
        r_cut=1.5,
    )


def _serialize_static_pack(pack: StaticLabelPack) -> Dict[str, Any]:
    return {
        'label': pack.label,
        'natm': int(pack.natm),
        'n_ic': int(pack.n_ic),
        'energies': pack.energies,
        'int_grad': pack.int_grad,
        'int_hess': pack.int_hess,
        'int_coords': pack.int_coords,
        'layout_counts': pack.layout.counts,
        'layout_displs': pack.layout.displs,
        'bond_end': int(pack.bond_end),
        'angle_end': int(pack.angle_end),
        'dihedral_start': int(pack.dihedral_start),
        'dihedral_end': int(pack.dihedral_end),
        'sym3_rows': pack.sym3_rows,
        'sym3_center_ids': pack.sym3_center_ids,
        'n_sym3': int(pack.n_sym3),
        'r_cut': float(pack.r_cut),
        'tasks': [
            {
                'task_id': int(t.task_id),
                'sym_idx': int(t.sym_idx),
                'mask0': t.mask0,
                'mask': t.mask,
            }
            for t in pack.tasks
        ],
    }


def _deserialize_static_pack(payload: Dict[str, Any]) -> StaticLabelPack:
    tasks = [
        SymmetryTask(
            task_id=int(item['task_id']),
            sym_idx=int(item['sym_idx']),
            mask0=np.ascontiguousarray(item['mask0'], dtype=np.int64),
            mask=np.ascontiguousarray(item['mask'], dtype=np.int64),
        )
        for item in payload['tasks']
    ]

    layout = RankLayout(
        counts=np.asarray(payload['layout_counts'], dtype=np.int64),
        displs=np.asarray(payload['layout_displs'], dtype=np.int64),
    )

    return StaticLabelPack(
        label=str(payload['label']),
        natm=int(payload['natm']),
        n_ic=int(payload['n_ic']),
        energies=np.asarray(payload['energies'], dtype=np.float64),
        int_grad=np.asarray(payload['int_grad'], dtype=np.float64),
        int_hess=np.asarray(payload['int_hess'], dtype=np.float64),
        int_coords=np.asarray(payload['int_coords'], dtype=np.float64),
        tasks=tasks,
        layout=layout,
        bond_end=int(payload['bond_end']),
        angle_end=int(payload['angle_end']),
        dihedral_start=int(payload['dihedral_start']),
        dihedral_end=int(payload['dihedral_end']),
        sym3_rows=np.asarray(payload['sym3_rows'], dtype=np.int64),
        sym3_center_ids=np.asarray(payload['sym3_center_ids'], dtype=np.int64),
        n_sym3=int(payload['n_sym3']),
        r_cut=float(payload['r_cut']),
    )


def _serialize_step_packet(packet: StepPacket) -> Dict[str, Any]:
    return {
        'org_int_coords': np.asarray(packet.org_int_coords, dtype=np.float64),
        'b_matrix': np.asarray(packet.b_matrix, dtype=np.float64),
        'use_cosine_dihedral': bool(packet.use_cosine_dihedral),
    }


def _deserialize_step_packet(payload: Dict[str, Any]) -> StepPacket:
    return StepPacket(
        org_int_coords=np.asarray(payload['org_int_coords'], dtype=np.float64),
        b_matrix=np.asarray(payload['b_matrix'], dtype=np.float64),
        use_cosine_dihedral=bool(payload['use_cosine_dihedral']),
    )


def _evaluate_local_tasks(local_pack: LocalLabelPack, packet: StepPacket) -> LabelReduction:
    gp = local_pack.global_pack
    red = LabelReduction.zeros(gp.natm)

    natm = gp.natm
    ncart = packet.b_matrix.shape[1]

    grouped_rows = [
        np.ascontiguousarray(gp.sym3_rows[gp.sym3_center_ids == rid], dtype=np.int64)
        for rid in range(gp.n_sym3)
    ]
    grouped_current = [
        np.asarray(packet.org_int_coords[rows], dtype=np.float64)
        for rows in grouped_rows
    ]
    grouped_brows = [
        np.asarray(packet.b_matrix[rows, :], dtype=np.float64)
        for rows in grouped_rows
    ]

    for task in local_pack.local_tasks:
        sym = task.sym_idx

        energy = gp.energies[sym]
        grad = gp.int_grad[sym]
        hess = gp.int_hess[sym]
        ref = gp.int_coords[sym]

        if grouped_rows:
            ref_masked = ref[task.mask]
            grouped_candidate = [
                np.asarray(ref_masked[rows], dtype=np.float64)
                for rows in grouped_rows
            ]
            d2_i, grad_d2_i = _grouped_signature_distance_and_gradient(
                grouped_current,
                grouped_candidate,
                grouped_brows,
                ncart,
            )
            W_i, grad_w_i = _signature_inverse_weight_and_gradient(
                d2_i,
                grad_d2_i,
                confidence_radius=1.0,
                power=2,
                eps=1.0e-12,
            )
            dW_i = grad_w_i.reshape(natm, 3)
        else:
            W_i = 1.0
            dW_i = np.zeros((natm, 3), dtype=np.float64)

        U_i, dU_i, dist_org, dist_correlation = evaluate_candidate_taylor(
            energy=energy,
            grad=grad,
            hess=hess,
            ref_coords=ref,
            mask0=task.mask0,
            mask=task.mask,
            org_int_coords=packet.org_int_coords,
            b_matrix=packet.b_matrix,
            use_cosine_dihedral=packet.use_cosine_dihedral,
            dihedral_start=gp.dihedral_start,
            dihedral_end=gp.dihedral_end,
        )

        red.S += W_i
        red.N += W_i * U_i
        red.dS += dW_i
        red.dN += dW_i * U_i + W_i * dU_i

        red.rmsd_bond_sum += _segment_rmsd(dist_org[:gp.bond_end])
        red.rmsd_angle_sum += _segment_rmsd(dist_org[gp.bond_end:gp.angle_end])
        red.rmsd_dihedral_sum += _segment_rmsd(dist_correlation[gp.dihedral_start:gp.dihedral_end])
        red.rmsd_count += 1

    return red


class InterpolationMPIPreloadEngine:
    def __init__(self, driver: Any, comm):
        self.driver = driver
        self.comm = comm
        self.rank = comm.Get_rank()
        self.size = comm.Get_size()
        self._dirty = True
        self._local_packs: Dict[str, LocalLabelPack] = {}

    def mark_dirty(self):
        self._dirty = True

    def preload_static_data(self, force=False):
        if (not force) and (not self._dirty):
            return

        if self.driver.qm_data_points is None:
            self.driver.qm_data_points = self.driver.read_qm_data_points()

        self.driver.prepare_runtime_data_cache(force=True)

        if self.rank == 0:
            labels = [dp.point_label for dp in self.driver.qm_data_points]
        else:
            labels = None
        labels = self.comm.bcast(labels, root=0)

        packs = {}
        for label in labels:
            if self.rank == 0:
                gp = _build_pack_for_label(self.driver, label, self.size)
                payload = _serialize_static_pack(gp)
            else:
                payload = None

            payload = self.comm.bcast(payload, root=0)
            gp = _deserialize_static_pack(payload)

            slc = gp.layout.local_slice(self.rank)
            packs[label] = LocalLabelPack(
                global_pack=gp,
                local_tasks=list(gp.tasks[slc]),
            )

        self._local_packs = packs
        self._dirty = False

    def broadcast_step_packet(self, packet: StepPacket | None, root: int = 0) -> StepPacket:
        if self.rank == root:
            if packet is None:
                raise ValueError("Root rank must provide StepPacket for broadcast.")
            payload = _serialize_step_packet(packet)
        else:
            payload = None

        payload = self.comm.bcast(payload, root=root)
        return _deserialize_step_packet(payload)

    def compute_label_bcast(
            self,
            dp_label: str,
            packet: StepPacket | None,
            root: int = 0) -> Tuple[float, np.ndarray, Dict[str, float]]:
        synced_packet = self.broadcast_step_packet(packet, root=root)
        return self.compute_label(dp_label, synced_packet)

    def compute_label(self, dp_label: str, packet: StepPacket) -> Tuple[float, np.ndarray, Dict[str, float]]:
        if self._dirty:
            self.preload_static_data(force=False)
        local_pack = self._local_packs[dp_label]
        local_red = _evaluate_local_tasks(local_pack, packet)
        global_red = _allreduce_reduction(self.comm, local_red)
        return _finalize_reduction(global_red)
