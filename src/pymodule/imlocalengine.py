from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from typing import Any, Dict, Tuple
import os

from .immpiengine import (
    StepPacket,
    RankLayout,
    LocalLabelPack,
    LabelReduction,
    _build_pack_for_label,
    _serialize_static_pack,
    _deserialize_static_pack,
    _serialize_step_packet,
    _deserialize_step_packet,
    _evaluate_local_tasks,
    _finalize_reduction,
)

_local_worker_storage = {}


def _set_worker_thread_env(omp_threads: int | None):

    if omp_threads is None:
        return

    n = str(max(1, int(omp_threads)))

    os.environ["OMP_NUM_THREADS"] = n
    os.environ["MKL_NUM_THREADS"] = n
    os.environ["OPENBLAS_NUM_THREADS"] = n
    os.environ["NUMEXPR_NUM_THREADS"] = n


def _init_local_worker(serialized_packs: Dict[str, dict], omp_threads: int | None):

    _set_worker_thread_env(omp_threads)
    packs = {}

    for label, payload in serialized_packs.items():

        packs[label] = _deserialize_static_pack(payload)

    _local_worker_storage["packs"] = packs


def _eval_chunk(payload: Tuple[str, dict, int, int]):

    dp_label, packet_payload, start, stop = payload
    packs = _local_worker_storage["packs"]
    gp = packs[dp_label]
    packet = _deserialize_step_packet(packet_payload)

    local_pack = LocalLabelPack(
        global_pack=gp,
        local_tasks=list(gp.tasks[int(start):int(stop)]),
    )

    red = _evaluate_local_tasks(local_pack, packet)

    return (
        float(red.S),
        float(red.N),
        red.dS,
        red.dN,
        float(red.rmsd_bond_sum),
        float(red.rmsd_angle_sum),
        float(red.rmsd_dihedral_sum),
        int(red.rmsd_count),
    )


class InterpolationLocalPreloadEngine:
    """
    Class to load all necessary information for a single compute potential call for a set of
    datapoints. Also handles the updating note when underlying database changed (mark dirty).
    """

    def __init__(self, driver: Any, n_workers: int = 0, min_tasks: int = 8, omp_threads: int = 1):

        self.driver = driver
        self.n_workers = n_workers
        self.min_tasks = min_tasks
        self.omp_threads = omp_threads

        self._dirty = True
        self._pool = None
        self._active_workers = 1
        self._global_packs: Dict[str, Any] = {}
        self._serialized_packs: Dict[str, dict] = {}

    def __del__(self):  # destructor of the current worker
        try:
            self.close()
        except Exception:
            pass

    def _resolve_workers(self) -> int:
        if self.n_workers > 0:
            return max(1, int(self.n_workers))
        return max(1, os.cpu_count() or 1)

    def mark_dirty(self):
        self._dirty = True

    def close(self):
        if self._pool is not None:
            self._pool.shutdown(wait=True)
            self._pool = None

    def preload_static_data(self, force: bool = False):
        if (not force) and (not self._dirty):
            return

        if self.driver.qm_data_points is None:
            self.driver.qm_data_points = self.driver.read_qm_data_points()

        self.driver.prepare_runtime_data_cache(force=True)
        resolved_workers = self._resolve_workers()
        payloads: Dict[str, dict] = {}
        packs: Dict[str, Any] = {}

        for dp in self.driver.qm_data_points:
            label = dp.point_label
            if label in payloads:
                continue
            gp = _build_pack_for_label(self.driver, label, resolved_workers)
            packs[label] = gp
            payloads[label] = _serialize_static_pack(gp)

        self._global_packs = packs
        self._serialized_packs = payloads

        self.close()
        self._active_workers = resolved_workers

        if self._active_workers > 1 and len(self._serialized_packs) > 0:
            self._pool = ProcessPoolExecutor(
                max_workers=self._active_workers,
                initializer=_init_local_worker,
                initargs=(self._serialized_packs, self.omp_threads),
            )

        self._dirty = False

    def compute_label(self, dp_label: str, packet: StepPacket):
        if self._dirty:
            self.preload_static_data(force=False)

        gp = self._global_packs.get(dp_label, None)
        if gp is None:
            raise KeyError(f"Local preload: label '{dp_label}' not found")

        n_tasks = len(gp.tasks)
        if n_tasks == 0:
            raise ValueError(f"Local preload: label '{dp_label}' has zero tasks.")

        run_serial = (
            self._pool is None
            or self._active_workers <= 1
            or n_tasks < self.min_tasks
        )

        if run_serial:
            local_pack = LocalLabelPack(global_pack=gp, local_tasks=list(gp.tasks))
            red = _evaluate_local_tasks(local_pack, packet)
            return _finalize_reduction(red)

        packet_payload = _serialize_step_packet(packet)
        layout = RankLayout.from_n_items(n_tasks, self._active_workers)

        chunks = []
        for w_idx in range(self._active_workers):
            count = int(layout.counts[w_idx])
            if count <= 0:
                continue
            start = int(layout.displs[w_idx])
            stop = start + count
            chunks.append((dp_label, packet_payload, start, stop))

        red_total = LabelReduction.zeros(gp.natm)

        try:
            parts = self._pool.map(_eval_chunk, chunks, chunksize=1)
            for part in parts:
                S_i, N_i, dS_i, dN_i, rb_i, ra_i, rd_i, rc_i = part
                red_total.S += S_i
                red_total.N += N_i
                red_total.dS += dS_i
                red_total.dN += dN_i
                red_total.rmsd_bond_sum += rb_i
                red_total.rmsd_angle_sum += ra_i
                red_total.rmsd_dihedral_sum += rd_i
                red_total.rmsd_count += rc_i
        except Exception:
            # Robust fallback
            local_pack = LocalLabelPack(global_pack=gp, local_tasks=list(gp.tasks))
            red = _evaluate_local_tasks(local_pack, packet)
            return _finalize_reduction(red)

        return _finalize_reduction(red_total)
