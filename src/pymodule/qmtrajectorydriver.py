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

from pathlib import Path
from mpi4py import MPI
import sys
import numpy as np
import h5py

from .outputstream import OutputStream
from .molecule import Molecule
from .molecularbasis import MolecularBasis
from .scfrestdriver import ScfRestrictedDriver
from .scfunrestdriver import ScfUnrestrictedDriver
from .scfrestopendriver import ScfRestrictedOpenDriver
from .lreigensolver import LinearResponseEigenSolver
from .cppsolver import ComplexResponse
from .veloxchemlib import mpi_master


class QMTrajectoryDriver:
    """
    Driver for QM calculations over a molecular dynamics trajectory.

    Designed for trajectories where the time interval between frames is small
    (e.g. Born–Oppenheimer or Car–Parrinello AIMD), so that the electronic
    structure changes smoothly from frame to frame.  The converged molecular
    orbitals from each frame are reused as the initial guess for the next
    frame via checkpoint-based restart, exploiting the small geometry change
    between consecutive AIMD steps to reduce the total number of SCF
    iterations significantly compared to independent SAD-guess calculations.

    Unlike the ensemble-based drivers, this driver:

    - Contains **no environment** (PE/NPE) infrastructure — every atom in each
      snapshot is treated at the QM level.
    - Processes frames **sequentially** and in order, preserving the
      temporal relationship between frames.
    - Propagates the SCF guess from frame *N* to frame *N + 1* via a
      a persistent internal checkpoint file.

    Snapshots are expected to come from :class:`QMTrajectoryParser`.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Usage::

        parser = QMTrajectoryParser()
        frames = parser.structures("traj.xtc", topology_file="top.tpr",
                                   qm_region="resname LIG")

        drv = QMTrajectoryDriver()
        results = drv.compute(
            frames,
            basis_set="def2-SVP",
            scf_options={"scf_type": "scf", "xcfun": "b3lyp"},
            property_options={"nstates": 5},
        )
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initialises the QMTrajectoryDriver.
        """
        if comm is None:
            comm = MPI.COMM_WORLD
        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()
        self.ostream = ostream

    # ── Static utilities ─────────────────────────────────────────────────────
    # These methods handle
    # NumPy array manipulation and HDF5 payload construction.  They are kept
    # self-contained here to avoid a circular dependency.

    @staticmethod
    def _fill_value_for_dtype(dtype):
        """Return a sensible padding value for a given dtype."""
        kind = np.dtype(dtype).kind
        if kind in ("f", "c"):
            return np.nan
        if kind in ("i", "u"):
            return -1
        if kind == "b":
            return False
        if kind in ("S", "U"):
            return b""
        return 0

    @staticmethod
    def _result_value_to_array(value):
        """
        Convert a result value to a NumPy array suitable for HDF5 writing.

        Returns ``None`` for unsupported object arrays.
        """
        if isinstance(value, np.ndarray):
            arr = value
        else:
            try:
                arr = np.asarray(value)
            except Exception:
                return None

        if arr.dtype.kind == "O":
            return None

        if arr.dtype.kind == "U":
            arr = arr.astype("S")

        return arr

    @staticmethod
    def _get_scf_hdf5_payload_local(scf_results, scf_history=None):
        """Build SCF payload dict for stacked HDF5 writing."""
        payload = {}

        keys = ["S"] + [
            f"{x}_{y}"
            for x in ["C", "E", "occ", "D", "F"]
            for y in ["alpha", "beta"]
        ]
        for key in keys:
            if key in scf_results:
                payload[key] = scf_results[key]

        payload["dipole_moment"] = scf_results["dipole_moment"]
        payload["scf_type"] = np.bytes_([scf_results["scf_type"]])
        payload["scf_energy"] = np.array([scf_results["scf_energy"]])

        if scf_history:
            keys = list(scf_history[0].keys())
            for key in keys:
                payload[f"scf_history_{key}"] = np.array(
                    [step[key] for step in scf_history]
                )

        return payload

    @staticmethod
    def _get_lr_rsp_hdf5_payload_local(rsp_results):
        """Build LR payload dict for stacked HDF5 writing."""
        payload = {}
        for key, value in rsp_results.items():
            if "vector" in key or "cube" in key or "file" in key or "details" in key:
                continue
            payload[key] = value
        return payload

    @staticmethod
    def _get_cpp_rsp_hdf5_payload_local(rsp_results, cpp_driver):
        """Build CPP payload dict for stacked HDF5 writing."""
        payload = {}

        freqs = rsp_results.get("frequencies", [])
        payload["frequencies"] = np.array(freqs)
        cpp_property = str(getattr(cpp_driver, "property", ""))
        payload["property"] = np.bytes_([cpp_property])

        spectrum = cpp_driver.get_spectrum(rsp_results, "au")
        if spectrum is not None:
            y_data = np.array(spectrum.get("y_data", []), dtype=float)
            payload["spectrum_y_data"] = y_data
            payload["spectrum_x_label"] = np.bytes_(
                [str(spectrum.get("x_label", ""))]
            )
            payload["spectrum_y_label"] = np.bytes_(
                [str(spectrum.get("y_label", ""))]
            )

            prop = cpp_property.lower()
            if prop == "absorption":
                payload["sigma"] = y_data
            elif prop == "ecd":
                payload["delta-epsilon"] = y_data

        return payload

    @staticmethod
    def _atomic_description_for_dataset(group_name, key):
        """Return metadata description for known atom-resolved datasets."""
        root_labels = {
            "nuclear_charges": "Nuclear Charges",
        }
        grouped_labels = {
            ("scf", "charges_resp"): "RESP Charges",
            ("scf", "population_analysis"): "Population Analysis",
            ("rsp", "detachment_charges"): "Detachment Charges",
            ("rsp", "attachment_charges"): "Attachment Charges",
        }

        if group_name is None:
            return root_labels.get(key)

        label = grouped_labels.get((str(group_name).lower(), str(key)))
        if label is not None:
            return label

        group = str(group_name).lower()
        key_s = str(key)
        if group == "rsp" and key_s.startswith("detach_charges_"):
            return "Detachment Charges"
        if group == "rsp" and key_s.startswith("attach_charges_"):
            return "Attachment Charges"

        return None

    @staticmethod
    def _apply_atomic_metadata(dset, group_name, key):
        """Attach standard metadata for atom-resolved scalar datasets."""
        description = QMTrajectoryDriver._atomic_description_for_dataset(
            group_name, key
        )
        if description is None:
            return
        dset.attrs["atomic_property"] = description

    @staticmethod
    def _temporal_description_for_dataset(group_name, key):
        """Return metadata description for known time-series datasets."""
        root_labels = {
            "nuclear_repulsion": "Nuclear Repulsion",
        }
        grouped_labels = {
            ("scf", "scf_energy"): "SCF Energy",
        }

        if group_name is None:
            return root_labels.get(key)
        return grouped_labels.get((str(group_name).lower(), str(key)))

    @staticmethod
    def _apply_temporal_metadata(dset, group_name, key):
        """Attach standard metadata for time-series (frame-stacked) datasets."""
        description = QMTrajectoryDriver._temporal_description_for_dataset(
            group_name, key
        )
        if description is None:
            return
        dset.attrs["temporal_property"] = description

    @staticmethod
    def _stack_result_arrays(values):
        """
        Stack per-frame arrays, padding variable shapes and returning lengths.

        :param values:
            Sequence of NumPy arrays, one per frame.
        :return:
            Tuple ``(stacked_array, lengths_or_none)``.
        """
        if not values:
            return None, None

        dtypes = [v.dtype for v in values]
        common_dtype = np.result_type(*dtypes)
        casted = [np.asarray(v, dtype=common_dtype) for v in values]

        shapes = [v.shape for v in casted]
        ndims = [v.ndim for v in casted]

        if len(set(shapes)) == 1:
            return np.stack(casted, axis=0), None

        if len(set(ndims)) != 1:
            flat = [v.reshape(-1) for v in casted]
            max_len = max(v.size for v in flat)
            out = np.full(
                (len(flat), max_len),
                QMTrajectoryDriver._fill_value_for_dtype(common_dtype),
                dtype=common_dtype,
            )
            lengths = np.zeros(len(flat), dtype=np.int32)
            for i, arr in enumerate(flat):
                out[i, : arr.size] = arr
                lengths[i] = arr.size
            return out, lengths

        max_shape = tuple(
            max(s[d] for s in shapes) for d in range(len(shapes[0]))
        )
        out_shape = (len(casted),) + max_shape
        out = np.full(
            out_shape,
            QMTrajectoryDriver._fill_value_for_dtype(common_dtype),
            dtype=common_dtype,
        )
        lengths = np.zeros(len(casted), dtype=np.int32)
        for i, arr in enumerate(casted):
            slices = (i,) + tuple(slice(0, s) for s in arr.shape)
            out[slices] = arr
            lengths[i] = arr.size

        return out, lengths

    # ── HDF5 output ──────────────────────────────────────────────────────────

    def _write_stacked_h5(
        self,
        snapshots,
        results,
        h5_file: str | Path,
        basis_source,
        dft_func_label,
        prop_driver=None,
    ):
        """
        Write trajectory SCF/property results to a single stacked HDF5 file.

        The layout mirrors the stacked HDF5 format so that
        the same post-processing tools can be applied to both ensemble and
        trajectory outputs.

        :param snapshots:
            List of snapshot dictionaries (or a single dict) as returned by
            :class:`QMTrajectoryParser`.
        :param results:
            Results dictionary returned by :meth:`compute`.
        :param h5_file:
            Output HDF5 path.
        :param basis_source:
            Basis label string or :class:`MolecularBasis` object.
        :param dft_func_label:
            DFT functional label written at root level.
        :param prop_driver:
            Optional property driver used to reconstruct LR extras.
        """
        if self.rank != mpi_master():
            return

        if isinstance(snapshots, dict):
            snapshots = [snapshots]

        h5_path = Path(h5_file)
        h5_path.parent.mkdir(parents=True, exist_ok=True)
        if h5_path.exists():
            h5_path.unlink()

        frames = np.array(
            [int(frame) for frame, _ in results.get("scf_all", [])],
            dtype=np.int32,
        )
        nframes = frames.size

        snapshot_by_frame = {int(s["frame"]): s for s in snapshots}
        ordered_snapshots = [snapshot_by_frame[int(fid)] for fid in frames]

        frame_indices = np.asarray(
            [
                int(s.get("frame_id", s.get("frame", -1)))
                for s in ordered_snapshots
            ],
            dtype=np.int32,
        )

        with h5py.File(h5_path, "w") as hf:
            hf.create_dataset("frame_id", data=frame_indices)
            total_frames = int(
                ordered_snapshots[0].get("total_frames", nframes)
            )
            hf.create_dataset(
                "total_frames", data=np.array([total_frames], dtype=np.int32)
            )

            traj_name = str(ordered_snapshots[0].get("trajectory_name", ""))
            hf.create_dataset("trajectory_name", data=np.bytes_([traj_name]))

            qm_idx_arrays = [
                np.asarray(s.get("qm_atom_indices", []), dtype=np.int32)
                for s in ordered_snapshots
            ]
            _stacked_idx, _lengths_idx = self._stack_result_arrays(qm_idx_arrays)
            if _stacked_idx is not None:
                hf.create_dataset("qm_atom_indices", data=_stacked_idx)
                if _lengths_idx is not None:
                    hf.create_dataset("qm_atom_indices_lengths", data=_lengths_idx)

            # ── Root-level metadata (mirrors create_hdf5 per frame, stacked) ──
            _basis_label = (
                basis_source if isinstance(basis_source, str)
                else basis_source.get_label()
            )
            _dft_label_str = str(dft_func_label) if dft_func_label else "HF"

            # Read basis once (mol0) to assign atom basis labels to the molecule
            # so molecule.get_atom_basis_labels() works correctly.
            _labels0 = [str(x) for x in ordered_snapshots[0]["qm_elements"]]
            _coords0 = np.asarray(ordered_snapshots[0]["qm_coords"], dtype=float)
            _mol0 = Molecule(_labels0, _coords0)
            _mol0.set_charge(int(ordered_snapshots[0].get("qm_charge", 0)))
            _mol0.set_multiplicity(int(ordered_snapshots[0].get("qm_multiplicity", 1)))
            MolecularBasis.read(_mol0, _basis_label)
            _abl_flat = []
            for _a, _e in _mol0.get_atom_basis_labels():
                _abl_flat += [_a, _e]
            _abl_bytes = np.bytes_(_abl_flat)

            per_nuclear_repulsion = []
            per_nuclear_charges = []
            per_atom_coordinates = []
            per_number_of_atoms = []
            per_number_of_alpha_electrons = []
            per_number_of_beta_electrons = []
            per_molecular_charge = []
            per_spin_multiplicity = []
            per_basis_set = []
            per_atom_basis_labels = []
            per_dft_func_label_list = []
            per_potfile_text = []

            for snap in ordered_snapshots:
                _lbl = [str(x) for x in snap["qm_elements"]]
                _crd = np.asarray(snap["qm_coords"], dtype=float)
                _mol = Molecule(_lbl, _crd)
                _mol.set_charge(int(snap.get("qm_charge", 0)))
                _mol.set_multiplicity(int(snap.get("qm_multiplicity", 1)))

                per_nuclear_repulsion.append(
                    np.array([_mol.nuclear_repulsion_energy(basis_source)])
                )
                per_nuclear_charges.append(
                    np.asarray(_mol.get_element_ids())
                )
                per_atom_coordinates.append(
                    np.asarray(_mol.get_coordinates_in_bohr())
                )
                per_number_of_atoms.append(
                    np.array([_mol.number_of_atoms()])
                )
                per_number_of_alpha_electrons.append(
                    np.array([_mol.number_of_alpha_electrons()])
                )
                per_number_of_beta_electrons.append(
                    np.array([_mol.number_of_beta_electrons()])
                )
                per_molecular_charge.append(
                    np.array([_mol.get_charge()])
                )
                per_spin_multiplicity.append(
                    np.array([_mol.get_multiplicity()])
                )
                per_basis_set.append(np.bytes_([_basis_label]))
                per_atom_basis_labels.append(_abl_bytes)
                per_dft_func_label_list.append(np.bytes_([_dft_label_str]))
                per_potfile_text.append(np.bytes_(['']))

            for _key, _arrays in [
                ("nuclear_repulsion", per_nuclear_repulsion),
                ("nuclear_charges", per_nuclear_charges),
                ("atom_coordinates", per_atom_coordinates),
                ("number_of_atoms", per_number_of_atoms),
                ("number_of_alpha_electrons", per_number_of_alpha_electrons),
                ("number_of_beta_electrons", per_number_of_beta_electrons),
                ("molecular_charge", per_molecular_charge),
                ("spin_multiplicity", per_spin_multiplicity),
                ("basis_set", per_basis_set),
                ("atom_basis_labels_flattened", per_atom_basis_labels),
                ("dft_func_label", per_dft_func_label_list),
                ("potfile_text", per_potfile_text),
            ]:
                _stacked, _lengths = self._stack_result_arrays(_arrays)
                if _stacked is not None:
                    _dset = hf.create_dataset(_key, data=_stacked)
                    if _key == "nuclear_charges":
                        self._apply_atomic_metadata(_dset, None, _key)
                    self._apply_temporal_metadata(_dset, None, _key)
                    if _lengths is not None:
                        hf.create_dataset(f"{_key}_lengths", data=_lengths)

            # Stacked SCF and RSP group payloads.
            scf_history_all = results.get("scf_history_all", [])
            if not scf_history_all:
                scf_history_all = [None] * len(results.get("scf_all", []))
            group_payloads = []
            if results.get("scf_all") is not None:
                group_payloads.append((
                    "scf",
                    [
                        self._get_scf_hdf5_payload_local(res, hist)
                        for (_, res), hist in zip(
                            results["scf_all"], scf_history_all
                        )
                        if res is not None
                    ],
                ))
            if results.get("prop_all") is not None:
                is_cpp = hasattr(prop_driver, "get_spectrum")
                if is_cpp:
                    rsp_payloads = [
                        self._get_cpp_rsp_hdf5_payload_local(res, prop_driver)
                        for _, res in results["prop_all"]
                    ]
                else:
                    rsp_payloads = [
                        self._get_lr_rsp_hdf5_payload_local(res)
                        for _, res in results["prop_all"]
                    ]
                group_payloads.append(("rsp", rsp_payloads))

            for group_name, frame_payloads in group_payloads:
                if not frame_payloads:
                    continue

                common_keys = sorted(
                    set.intersection(*(set(p.keys()) for p in frame_payloads))
                )
                grp = hf.require_group(group_name)

                for key in common_keys:
                    arrays = []
                    supported = True
                    for i, payload in enumerate(frame_payloads):
                        arr = self._result_value_to_array(payload[key])
                        if arr is None:
                            supported = False
                            break
                        arrays.append(arr)

                    if not supported:
                        continue

                    stacked, lengths = self._stack_result_arrays(arrays)
                    if stacked is None:
                        continue

                    dset = grp.create_dataset(key, data=stacked)
                    self._apply_atomic_metadata(dset, group_name, key)
                    self._apply_temporal_metadata(dset, group_name, key)
                    if lengths is not None:
                        grp.create_dataset(f"{key}_lengths", data=lengths)

            # Reconstruct LR-specific extras: state vectors (S1, S2, …) and
            # NTO datasets that are not present in the standard rsp payload.
            if (
                prop_driver is not None
                and results.get("prop_all") is not None
                and not hasattr(prop_driver, "get_spectrum")
                and hasattr(prop_driver, "get_full_solution_vector")
            ):
                rsp_grp = hf.require_group("rsp")
                state_vectors = {}
                nto_arrays = {}

                for i, (_, rsp_res) in enumerate(results["prop_all"]):
                    snap = ordered_snapshots[i]
                    scf_res = results["scf_all"][i][1]

                    labels = [str(x) for x in snap["qm_elements"]]
                    coords = np.asarray(snap["qm_coords"], dtype=float)
                    molecule = Molecule(labels, coords)
                    molecule.set_charge(int(snap.get("qm_charge", 0)))
                    molecule.set_multiplicity(
                        int(snap.get("qm_multiplicity", 1))
                    )

                    if isinstance(basis_source, str):
                        basis_for_snap = MolecularBasis.read(
                            molecule, basis_source
                        )
                    else:
                        basis_for_snap = basis_source

                    nocc = molecule.number_of_alpha_occupied_orbitals(
                        basis_for_snap
                    )

                    eigvecs = rsp_res.get("eigenvectors_distributed", {})
                    nstates = int(
                        rsp_res.get("number_of_states", len(eigvecs))
                    )

                    for s in range(nstates):
                        if s not in eigvecs:
                            continue
                        eigvec = prop_driver.get_full_solution_vector(
                            eigvecs[s]
                        )
                        if eigvec is None:
                            continue

                        s_key = f"S{s + 1}"
                        state_vectors.setdefault(s_key, []).append(
                            np.asarray(eigvec)
                        )

                        if getattr(prop_driver, "nto", False) and hasattr(prop_driver, "get_nto"):
                            if getattr(prop_driver, "core_excitation", False):
                                occ_count = int(
                                    getattr(prop_driver, "num_core_orbitals", 0)
                                )
                                mo_occ = scf_res["C_alpha"][
                                    :, :occ_count
                                ].copy()
                                mo_vir = scf_res["C_alpha"][:, nocc:].copy()
                                z_mat = eigvec[: eigvec.size // 2].reshape(
                                    occ_count, -1
                                )
                                y_mat = eigvec[eigvec.size // 2 :].reshape(
                                    occ_count, -1
                                )
                            elif getattr(
                                prop_driver, "restricted_subspace", False
                            ):
                                num_core = int(
                                    getattr(prop_driver, "num_core_orbitals", 0)
                                )
                                num_val = int(
                                    getattr(
                                        prop_driver, "num_valence_orbitals", 0
                                    )
                                )
                                num_vir = int(
                                    getattr(
                                        prop_driver, "num_virtual_orbitals", 0
                                    )
                                )
                                occ_count = num_core + num_val
                                mo_occ = np.hstack((
                                    scf_res["C_alpha"][:, :num_core].copy(),
                                    scf_res["C_alpha"][
                                        :, nocc - num_val : nocc
                                    ].copy(),
                                ))
                                mo_vir = scf_res["C_alpha"][
                                    :, nocc : nocc + num_vir
                                ].copy()
                                z_mat = eigvec[: eigvec.size // 2].reshape(
                                    occ_count, -1
                                )
                                y_mat = eigvec[eigvec.size // 2 :].reshape(
                                    occ_count, -1
                                )
                            else:
                                mo_occ = scf_res["C_alpha"][:, :nocc].copy()
                                mo_vir = scf_res["C_alpha"][:, nocc:].copy()
                                z_mat = eigvec[: eigvec.size // 2].reshape(
                                    nocc, -1
                                )
                                y_mat = eigvec[eigvec.size // 2 :].reshape(
                                    nocc, -1
                                )

                            nto_mo = prop_driver.get_nto(
                                z_mat - y_mat, mo_occ, mo_vir
                            )
                            nto_prefix = f"NTO_S{s + 1}_"
                            nto_alpha_orb = np.asarray(
                                nto_mo.alpha_to_numpy()
                            )
                            nto_arrays.setdefault(
                                nto_prefix + "alpha_orbitals", []
                            ).append(nto_alpha_orb)
                            nto_arrays.setdefault(
                                nto_prefix + "alpha_energies", []
                            ).append(np.asarray(nto_mo.ea_to_numpy()))
                            nto_arrays.setdefault(
                                nto_prefix + "alpha_occupations", []
                            ).append(np.asarray(nto_mo.occa_to_numpy()))

                for key, arrays in state_vectors.items():
                    if len(arrays) != nframes:
                        continue
                    stacked, lengths = self._stack_result_arrays(arrays)
                    if stacked is None:
                        continue
                    rsp_grp.create_dataset(key, data=stacked)
                    if lengths is not None:
                        rsp_grp.create_dataset(f"{key}_lengths", data=lengths)

                for key, arrays in nto_arrays.items():
                    if len(arrays) != nframes:
                        continue
                    stacked, lengths = self._stack_result_arrays(arrays)
                    if stacked is None:
                        continue
                    rsp_grp.create_dataset(key, data=stacked)
                    if lengths is not None:
                        rsp_grp.create_dataset(f"{key}_lengths", data=lengths)

    # ── Main compute loop ─────────────────────────────────────────────────────

    def compute(
        self,
        snapshots,
        basis_set,
        scf_options=None,
        property_options=None,
        propagate_guess: bool = True,
        guess_fallback: bool = True,
        qm_charge: int | None = None,
        qm_multiplicity: int | None = None,
        stacked_h5_file: str | Path | None = "trajectory.h5",
        scf_driver=None,
        prop_driver=None,
    ):
        """
        Run SCF (and optionally property) calculations over a trajectory.

        Frames are processed **in the order they appear in** ``snapshots``.
        When ``propagate_guess=True`` (default), the converged MOs from frame
        *N* are reused as the initial density guess for frame *N + 1*,
        exploiting the small geometry change between consecutive AIMD steps to
        reduce the number of SCF iterations.

        On the **first frame** the checkpoint does not yet exist, so the SCF
        driver falls back to its default guess (SAD) automatically.  After
        convergence the checkpoint is written and all subsequent frames use it.

        If a frame fails to converge with the propagated guess and
        ``guess_fallback=True``, the driver retries with a fresh SAD guess.
        The new converged MOs are then used as the checkpoint for the next
        frame.

        :param snapshots:
            Ordered list of snapshot dictionaries (or a single dict) as
            returned by :class:`QMTrajectoryParser`.
        :param basis_set:
            Basis name string (e.g. ``"def2-SVP"``) or a pre-built
            :class:`MolecularBasis` object used for all frames.
        :param scf_options:
            Dictionary of SCF settings.  Recognised key:

            - ``"scf_type"``: ``"scf"`` (restricted, default),
              ``"roscf"`` (restricted open-shell), or ``"uscf"``
              (unrestricted).

            All other keys are forwarded verbatim to the SCF driver's
            ``update_settings`` call (e.g. ``"xcfun"``, ``"eri_thresh"``,
            ``"grid_level"``).
        :param property_options:
            Dictionary of property settings.  The driver to run is selected
            automatically:

            - ``"nstates"`` present → :class:`LinearResponseEigenSolver`.
            - ``"frequencies"`` present → :class:`ComplexResponse`.

            If ``None`` or empty, no property calculation is run.
            All keys are forwarded to the property driver's
            ``update_settings`` call.
        :param propagate_guess:
            If ``True`` (default), propagate the converged MOs from each
            frame as the initial guess for the next frame.  Set to ``False``
            for an independent SAD guess on every frame.
        :param guess_fallback:
            If ``True`` (default) and a frame fails to converge with the
            propagated guess, the driver retries automatically with a SAD
            guess.  The stale checkpoint is removed before the retry so that
            the next frame starts from the recovered geometry.
        :param qm_charge:
            Optional charge override applied to all frames.  If ``None``, the
            per-snapshot ``qm_charge`` value is used.
        :param qm_multiplicity:
            Optional multiplicity override applied to all frames.  If
            ``None``, the per-snapshot ``qm_multiplicity`` value is used.
        :param stacked_h5_file:
            Path for the stacked HDF5 output file.  Defaults to
            ``"trajectory.h5"``.  Set to ``None`` to disable HDF5 writing.
        :param scf_driver:
            Advanced use only — supply a pre-configured SCF driver instance
            directly.  When provided, ``scf_options`` is ignored.
        :param prop_driver:
            Advanced use only — supply a pre-configured property driver
            instance directly.  When provided, ``property_options`` is
            ignored.
        :return:
            Dictionary with keys:

            - ``scf_all``: list of ``(frame, scf_results)`` tuples.
            - ``prop_all``: list of ``(frame, property_results)`` tuples
              (only present when a property driver is used).
        :raises ValueError:
            If ``qm_multiplicity`` is not a positive integer, or if a
            snapshot's charge/multiplicity combination is incompatible.
        """
        if isinstance(snapshots, dict):
            snapshots = [snapshots]

        if qm_charge is not None:
            qm_charge = int(qm_charge)
        if qm_multiplicity is not None:
            qm_multiplicity = int(qm_multiplicity)
            if qm_multiplicity <= 0:
                raise ValueError(
                    "qm_multiplicity must be a positive integer."
                )

        guess_h5 = Path("traj_guess.h5") if propagate_guess else None

        # ── Build SCF driver from scf_options if not provided by caller ───────
        if scf_driver is None:
            scf_opts = dict(scf_options or {})
            scf_type = scf_opts.pop("scf_type", "scf").lower()
            if scf_type == "uscf":
                scf_driver = ScfUnrestrictedDriver(self.comm, self.ostream)
            elif scf_type == "roscf":
                scf_driver = ScfRestrictedOpenDriver(self.comm, self.ostream)
            else:
                scf_driver = ScfRestrictedDriver(self.comm, self.ostream)
            if scf_opts:
                scf_driver.update_settings(scf_opts)

        # ── Build property driver from property_options if not provided ───────
        if prop_driver is None and property_options:
            prop_opts = dict(property_options)
            if "nstates" in prop_opts:
                prop_driver = LinearResponseEigenSolver(
                    self.comm, self.ostream
                )
            elif "frequencies" in prop_opts:
                prop_driver = ComplexResponse(self.comm, self.ostream)
            if prop_driver is not None:
                prop_driver.update_settings(prop_opts)
                if isinstance(prop_driver, ComplexResponse) and "property" in prop_opts:
                    prop_driver.set_cpp_property(prop_opts["property"])

        do_rsp = prop_driver is not None
        rsp_driver = prop_driver
        if do_rsp and hasattr(rsp_driver, "detach_attach"):
            rsp_driver.detach_attach = True

        # ── Build reference molecule and basis from the first snapshot ─────────
        first = snapshots[0]
        labels0 = [str(x) for x in first["qm_elements"]]
        coords0 = np.asarray(first["qm_coords"], dtype=float)
        molecule0 = Molecule(labels0, coords0)
        charge0 = (
            int(qm_charge)
            if qm_charge is not None
            else int(first.get("qm_charge", 0))
        )
        mult0 = (
            int(qm_multiplicity)
            if qm_multiplicity is not None
            else int(first.get("qm_multiplicity", 1))
        )
        molecule0.set_charge(charge0)
        molecule0.set_multiplicity(mult0)

        if isinstance(basis_set, str):
            fixed_basis = MolecularBasis.read(molecule0, basis_set)
        else:
            fixed_basis = basis_set

        scf_all = []
        scf_history_all = []
        rsp_all = [] if do_rsp else None

        # ── Main trajectory loop ───────────────────────────────────────────────
        for snap in snapshots:
            frame = int(snap["frame"])

            labels = [str(x) for x in snap["qm_elements"]]
            coords = np.asarray(snap["qm_coords"], dtype=float)
            snap_charge = (
                qm_charge
                if qm_charge is not None
                else int(snap.get("qm_charge", 0))
            )
            snap_mult = (
                qm_multiplicity
                if qm_multiplicity is not None
                else int(snap.get("qm_multiplicity", 1))
            )

            molecule = Molecule(labels, coords)
            molecule.set_charge(snap_charge)
            molecule.set_multiplicity(snap_mult)

            if not molecule.check_multiplicity():
                raise ValueError(
                    f"Incompatible QM charge ({snap_charge}) and multiplicity "
                    f"({snap_mult}) for frame {frame}."
                )

            # Configure guess propagation.
            # On the first frame guess_h5 does not yet exist; validate_checkpoint
            # inside the SCF driver detects this and falls back to SAD
            # automatically.  On subsequent frames the existing checkpoint file
            # is validated (nuclear charges + basis must match, which they always
            # do for a fixed-topology trajectory) and the converged MOs are used
            # as the starting density.
            if guess_h5 is not None:
                scf_driver.checkpoint_file = str(guess_h5)
                scf_driver.restart = True

            scf_results = scf_driver.compute(molecule, fixed_basis)

            # ── Optional fallback to SAD if propagated guess did not converge ──
            if guess_fallback and not scf_driver.is_converged:
                if self.rank == mpi_master():
                    self.ostream.print_info(
                        f"Frame {frame}: SCF did not converge with the "
                        "propagated guess.  Retrying with SAD initial guess..."
                    )
                    self.ostream.print_blank()

                # Remove the potentially stale/corrupted checkpoint so that the
                # retry and all subsequent frames start from a clean SAD guess
                # rather than a broken one.
                if guess_h5 is not None and guess_h5.is_file():
                    guess_h5.unlink()

                scf_driver.checkpoint_file = (
                    None if guess_h5 is None else str(guess_h5)
                )
                scf_driver.restart = False
                scf_results = scf_driver.compute(molecule, fixed_basis)

                # Restore checkpoint settings for the next frame.
                if guess_h5 is not None:
                    scf_driver.checkpoint_file = str(guess_h5)
                    scf_driver.restart = True

            if scf_results is None:
                if self.rank == mpi_master():
                    self.ostream.print_info(
                        f"Frame {frame}: SCF did not converge — "
                        "skipping this frame in the stacked output."
                    )
                    self.ostream.print_blank()
                continue

            scf_all.append((frame, scf_results))
            scf_history_all.append(
                list(scf_driver._history) if scf_driver._history else []
            )

            if do_rsp:
                rsp_results = rsp_driver.compute(
                    molecule, fixed_basis, scf_results
                )
                rsp_all.append((frame, rsp_results))

        results = {"scf_all": scf_all, "scf_history_all": scf_history_all}
        if do_rsp:
            results["prop_all"] = rsp_all

        if stacked_h5_file is not None:
            _xcfun = getattr(scf_driver, "xcfun", None)
            if hasattr(_xcfun, "get_func_label"):
                h5_dft_func_label = _xcfun.get_func_label()
            elif _xcfun:
                h5_dft_func_label = str(_xcfun)
            else:
                h5_dft_func_label = "HF"
            self._write_stacked_h5(
                snapshots,
                results,
                stacked_h5_file,
                fixed_basis,
                h5_dft_func_label,
                rsp_driver if do_rsp else None,
            )

        return results

