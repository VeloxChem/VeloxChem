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
from .cppsolver import ComplexResponseSolver
from .respchargesdriver import RespChargesDriver
from .veloxchemlib import mpi_master, hartree_in_ev, hartree_in_kjpermol
from .errorhandler import assert_msg_critical


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

        if "charges_resp" in scf_results:
            payload["charges_resp"] = np.asarray(
                scf_results["charges_resp"], dtype=float
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

    # ── HDF5 incremental output ───────────────────────────────────────────────
    #
    # The file is opened before the loop and written frame-by-frame, so it is
    # always valid on disk.  A crash loses at most the in-flight frame.
    #
    # Layout
    # ──────
    # Root-level resizable datasets hold per-frame molecular metadata
    # (nuclear_charges, atom_coordinates, nuclear_repulsion, …) plus a
    # boolean ``converged`` flag.  ``scf/`` and ``rsp/`` groups hold the
    # corresponding quantum-chemical results; non-converged frames store
    # NaN-filled placeholders so frame indices stay aligned.

    def _h5_open_for_writing(
        self,
        h5_path: Path,
        first_snap: dict,
        fixed_basis,
        dft_label: str,
        total_frames: int,
    ):
        """
        Create the trajectory HDF5 file and initialise all resizable datasets.

        Called once before the main loop.  All datasets are created with
        ``maxshape=(None, ...)`` so they can be extended frame-by-frame via
        :meth:`_h5_append_frame`.

        :param h5_path:
            Destination path (created fresh; existing file is removed).
        :param first_snap:
            First snapshot dict, used to infer per-frame dataset shapes.
        :param fixed_basis:
            :class:`MolecularBasis` object (or basis-set label string).
        :param dft_label:
            DFT functional label string.
        :param total_frames:
            Total number of frames in the trajectory (written as metadata).
        """
        if self.rank != mpi_master():
            return

        h5_path.parent.mkdir(parents=True, exist_ok=True)
        if h5_path.exists():
            h5_path.unlink()

        _basis_label = (
            fixed_basis if isinstance(fixed_basis, str)
            else fixed_basis.get_label()
        )
        _dft_label_str = str(dft_label) if dft_label else "HF"

        # Build a reference molecule to infer shapes.
        _lbl0 = [str(x) for x in first_snap["qm_elements"]]
        _crd0 = np.asarray(first_snap["qm_coords"], dtype=float)
        _mol0 = Molecule(_lbl0, _crd0)
        _mol0.set_charge(int(first_snap.get("qm_charge", 0)))
        _mol0.set_multiplicity(int(first_snap.get("qm_multiplicity", 1)))
        MolecularBasis.read(_mol0, _basis_label)  # assign atom basis labels
        _abl_flat = []
        for _a, _e in _mol0.get_atom_basis_labels():
            _abl_flat += [_a, _e]
        _abl_bytes = np.bytes_(_abl_flat)

        natoms = _mol0.number_of_atoms()
        ncoords = natoms * 3

        with h5py.File(h5_path, "w") as hf:
            # ── Scalar metadata written once ──────────────────────────────────
            hf.create_dataset(
                "total_frames", data=np.array([total_frames], dtype=np.int32)
            )
            traj_name = str(first_snap.get("trajectory_name", ""))
            hf.create_dataset("trajectory_name", data=np.bytes_([traj_name]))

            qm_idx = np.asarray(
                first_snap.get("qm_atom_indices", []), dtype=np.int32
            )
            if qm_idx.size:
                hf.create_dataset(
                    "qm_atom_indices",
                    data=qm_idx.reshape(1, -1),
                    maxshape=(None, qm_idx.size),
                )

            # ── Per-frame resizable root datasets ─────────────────────────────
            def _res(name, shape, dtype, fill=np.nan, *, vlen=False):
                """Helper: create a (0, *shape) resizable dataset."""
                if vlen:
                    dt = h5py.special_dtype(vlen=bytes)
                    hf.create_dataset(
                        name,
                        shape=(0,),
                        maxshape=(None,),
                        dtype=dt,
                    )
                else:
                    hf.create_dataset(
                        name,
                        shape=(0,) + tuple(shape),
                        maxshape=(None,) + tuple(shape),
                        dtype=dtype,
                        fillvalue=fill,
                    )

            _res("frame_id",               (),        np.int32,   fill=-1)
            _res("converged",              (),        bool,       fill=False)
            _res("nuclear_repulsion",      (),        np.float64)
            _res("nuclear_charges",        (natoms,), np.int32,   fill=-1)
            _res("atom_coordinates",       (natoms, 3), np.float64)
            _res("number_of_atoms",        (),        np.int32,   fill=0)
            _res("number_of_alpha_electrons", (),     np.int32,   fill=0)
            _res("number_of_beta_electrons",  (),     np.int32,   fill=0)
            _res("molecular_charge",       (),        np.float64)
            _res("spin_multiplicity",      (),        np.int32,   fill=0)
            _res("basis_set",              (),        "S64",      fill=b"")
            _res("dft_func_label",         (),        "S64",      fill=b"")
            _res("potfile_text",           (),        "S1",       fill=b"")

            # atom_basis_labels_flattened has variable per-atom string size;
            # store as a fixed-length byte dataset sized from the first frame.
            hf.create_dataset(
                "atom_basis_labels_flattened",
                data=_abl_bytes.reshape(1, -1),
                maxshape=(None, len(_abl_bytes)),
                dtype=_abl_bytes.dtype,
            )

            # Apply metadata attributes.
            dset_nc = hf["nuclear_charges"]
            self._apply_atomic_metadata(dset_nc, None, "nuclear_charges")
            self._apply_temporal_metadata(hf["nuclear_repulsion"], None, "nuclear_repulsion")

    def _h5_append_frame(
        self,
        h5_path: Path,
        frame: int,
        snap: dict,
        molecule,
        fixed_basis,
        dft_label: str,
        scf_results,
        scf_history,
        rsp_results,
        prop_driver,
    ):
        """
        Append one frame's results to the trajectory HDF5 file.

        Non-converged frames (``scf_results is None``) are written with a
        ``converged=False`` flag and NaN placeholders in all numeric datasets,
        so the frame index is preserved.

        :param h5_path:
            Path to the already-open trajectory file.
        :param frame:
            Integer frame number.
        :param snap:
            Snapshot dictionary for this frame.
        :param molecule:
            :class:`Molecule` object for this frame.
        :param fixed_basis:
            :class:`MolecularBasis` object.
        :param dft_label:
            DFT functional label string.
        :param scf_results:
            SCF results dict, or ``None`` if SCF did not converge.
        :param scf_history:
            SCF history list, or ``[]`` if not available.
        :param rsp_results:
            Property results dict, or ``None`` if not computed.
        :param prop_driver:
            The property driver instance (used for CPP spectrum / LR extras),
            or ``None``.
        """
        if self.rank != mpi_master():
            return

        converged = scf_results is not None
        _basis_label = (
            fixed_basis if isinstance(fixed_basis, str)
            else fixed_basis.get_label()
        )
        _dft_label_str = str(dft_label) if dft_label else "HF"

        with h5py.File(h5_path, "a") as hf:

            def _ext(name):
                """Extend dataset along axis-0 by one and return new index."""
                dset = hf[name]
                n = dset.shape[0]
                dset.resize(n + 1, axis=0)
                return n

            # ── Root-level per-frame metadata ─────────────────────────────────
            i = _ext("frame_id")
            hf["frame_id"][i] = int(snap.get("frame_id", frame))
            _ext("converged");  hf["converged"][i] = converged

            e_nuc = molecule.effective_nuclear_repulsion_energy(fixed_basis)
            _ext("nuclear_repulsion"); hf["nuclear_repulsion"][i] = e_nuc

            nc = np.asarray(molecule.get_element_ids(), dtype=np.int32)
            _ext("nuclear_charges"); hf["nuclear_charges"][i] = nc

            coords = np.asarray(molecule.get_coordinates_in_bohr())
            _ext("atom_coordinates"); hf["atom_coordinates"][i] = coords

            _ext("number_of_atoms");
            hf["number_of_atoms"][i] = molecule.number_of_atoms()
            _ext("number_of_alpha_electrons")
            hf["number_of_alpha_electrons"][i] = molecule.number_of_alpha_electrons()
            _ext("number_of_beta_electrons")
            hf["number_of_beta_electrons"][i] = molecule.number_of_beta_electrons()
            _ext("molecular_charge")
            hf["molecular_charge"][i] = molecule.get_charge()
            _ext("spin_multiplicity")
            hf["spin_multiplicity"][i] = molecule.get_multiplicity()

            bs_bytes = np.bytes_([_basis_label])
            _ext("basis_set"); hf["basis_set"][i] = bs_bytes[0]
            dft_bytes = np.bytes_([_dft_label_str])
            _ext("dft_func_label"); hf["dft_func_label"][i] = dft_bytes[0]
            _ext("potfile_text"); hf["potfile_text"][i] = b""

            abl = hf["atom_basis_labels_flattened"]
            n_abl = abl.shape[0]
            abl.resize(n_abl + 1, axis=0)
            # shape already fixed from init; values are same for all frames
            # (fixed topology), so we just copy the first row or write zeros.
            abl[n_abl] = abl[0] if n_abl > 0 else abl[0]

            # qm_atom_indices (fixed across frames — only extend row count)
            if "qm_atom_indices" in hf:
                qi = hf["qm_atom_indices"]
                qi.resize(qi.shape[0] + 1, axis=0)
                qi[qi.shape[0] - 1] = qi[0]

            # ── SCF group ─────────────────────────────────────────────────────
            if converged:
                payload = self._get_scf_hdf5_payload_local(
                    scf_results, scf_history
                )
                scf_grp = hf.require_group("scf")
                for key, value in payload.items():
                    arr = self._result_value_to_array(value)
                    if arr is None:
                        continue
                    if key not in scf_grp:
                        # All axes unbounded so variable-length history arrays
                        # (scf_history_*) can grow in both dimensions.
                        mshape = (None,) + tuple(None for _ in arr.shape)
                        scf_grp.create_dataset(
                            key,
                            data=arr.reshape((1,) + arr.shape),
                            maxshape=mshape,
                        )
                        self._apply_atomic_metadata(scf_grp[key], "scf", key)
                        self._apply_temporal_metadata(scf_grp[key], "scf", key)
                    else:
                        dset = scf_grp[key]
                        n = dset.shape[0]
                        # Handle variable-length trailing dimension
                        # (e.g. scf_history_* arrays differ per frame).
                        if arr.ndim >= 1 and dset.ndim == arr.ndim + 1:
                            cur_w = dset.shape[1]
                            new_w = arr.shape[0]
                            if new_w > cur_w:
                                # Grow width; back-fill old rows with NaN.
                                dset.resize(cur_w if dset.dtype.kind != 'f'
                                            else cur_w, axis=1)
                                dset.resize(new_w, axis=1)
                                if dset.dtype.kind == 'f':
                                    dset[:n, cur_w:] = np.nan
                            elif new_w < cur_w:
                                # Pad incoming array with NaN to match width.
                                padded = np.full(cur_w, np.nan,
                                                 dtype=dset.dtype)
                                padded[:new_w] = arr
                                arr = padded
                        dset.resize(n + 1, axis=0)
                        dset[n] = arr

            # ── RSP group ─────────────────────────────────────────────────────
            if converged and rsp_results is not None and prop_driver is not None:
                is_cpp = hasattr(prop_driver, "get_spectrum")
                if is_cpp:
                    rsp_payload = self._get_cpp_rsp_hdf5_payload_local(
                        rsp_results, prop_driver
                    )
                else:
                    rsp_payload = self._get_lr_rsp_hdf5_payload_local(
                        rsp_results
                    )

                rsp_grp = hf.require_group("rsp")
                for key, value in rsp_payload.items():
                    arr = self._result_value_to_array(value)
                    if arr is None:
                        continue
                    if key not in rsp_grp:
                        mshape = (None,) + tuple(None for _ in arr.shape)
                        rsp_grp.create_dataset(
                            key,
                            data=arr.reshape((1,) + arr.shape),
                            maxshape=mshape,
                        )
                        self._apply_atomic_metadata(rsp_grp[key], "rsp", key)
                        self._apply_temporal_metadata(rsp_grp[key], "rsp", key)
                    else:
                        dset = rsp_grp[key]
                        n = dset.shape[0]
                        dset.resize(n + 1, axis=0)
                        dset[n] = arr

                # LR-specific: state eigenvectors
                if (
                    not is_cpp
                    and hasattr(prop_driver, "get_full_solution_vector")
                ):
                    basis_for_snap = fixed_basis
                    nocc = molecule.number_of_alpha_occupied_orbitals(
                        basis_for_snap
                    )
                    eigvecs = rsp_results.get("eigenvectors_distributed", {})
                    nstates = int(
                        rsp_results.get("number_of_states", len(eigvecs))
                    )
                    for s in range(nstates):
                        if s not in eigvecs:
                            continue
                        eigvec = prop_driver.get_full_solution_vector(eigvecs[s])
                        if eigvec is None:
                            continue
                        s_key = f"S{s + 1}"
                        ev_arr = np.asarray(eigvec)
                        if s_key not in rsp_grp:
                            rsp_grp.create_dataset(
                                s_key,
                                data=ev_arr.reshape(1, -1),
                                maxshape=(None, ev_arr.size),
                            )
                        else:
                            dset = rsp_grp[s_key]
                            n = dset.shape[0]
                            dset.resize(n + 1, axis=0)
                            dset[n] = ev_arr

                        if getattr(prop_driver, "nto", False) and hasattr(
                            prop_driver, "get_nto"
                        ):
                            if getattr(prop_driver, "core_excitation", False):
                                occ_count = int(
                                    getattr(prop_driver, "num_core_orbitals", 0)
                                )
                                mo_occ = scf_results["C_alpha"][
                                    :, :occ_count
                                ].copy()
                                mo_vir = scf_results["C_alpha"][
                                    :, nocc:
                                ].copy()
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
                                    scf_results["C_alpha"][
                                        :, :num_core
                                    ].copy(),
                                    scf_results["C_alpha"][
                                        :, nocc - num_val : nocc
                                    ].copy(),
                                ))
                                mo_vir = scf_results["C_alpha"][
                                    :, nocc : nocc + num_vir
                                ].copy()
                                z_mat = eigvec[: eigvec.size // 2].reshape(
                                    occ_count, -1
                                )
                                y_mat = eigvec[eigvec.size // 2 :].reshape(
                                    occ_count, -1
                                )
                            else:
                                mo_occ = scf_results["C_alpha"][
                                    :, :nocc
                                ].copy()
                                mo_vir = scf_results["C_alpha"][
                                    :, nocc:
                                ].copy()
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
                            for nto_key, nto_arr in [
                                (
                                    nto_prefix + "alpha_orbitals",
                                    np.asarray(nto_mo.alpha_to_numpy()),
                                ),
                                (
                                    nto_prefix + "alpha_energies",
                                    np.asarray(nto_mo.ea_to_numpy()),
                                ),
                                (
                                    nto_prefix + "alpha_occupations",
                                    np.asarray(nto_mo.occa_to_numpy()),
                                ),
                            ]:
                                if nto_key not in rsp_grp:
                                    rsp_grp.create_dataset(
                                        nto_key,
                                        data=nto_arr.reshape(
                                            (1,) + nto_arr.shape
                                        ),
                                        maxshape=(None,) + nto_arr.shape,
                                    )
                                else:
                                    dset = rsp_grp[nto_key]
                                    n = dset.shape[0]
                                    dset.resize(n + 1, axis=0)
                                    dset[n] = nto_arr

    # ── Main compute loop ─────────────────────────────────────────────────────

    def compute(
        self,
        snapshots,
        basis_set,
        scf_options=None,
        property_options=None,
        resp_options=None,
        propagate_guess: bool = True,
        guess_fallback: bool = True,
        qm_charge: int | None = None,
        qm_multiplicity: int | None = None,
        stacked_h5_file: str | Path | None = "trajectory.h5",
        scf_driver=None,
        prop_driver=None,
        resp_driver=None,
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
        :param resp_options:
            Dictionary of RESP charge options forwarded to
            :class:`RespChargesDriver.update_settings`.  Pass an empty dict
            ``{}`` to enable RESP with default settings.  If ``None``
            (default), RESP charges are not computed.
        :param resp_driver:
            Advanced use only — supply a pre-configured
            :class:`RespChargesDriver` instance directly.  When provided,
            ``resp_options`` is ignored.
        :return:
            Dictionary with keys:

            - ``scf_all``: list of ``(frame, scf_results)`` tuples, where
              ``scf_results`` is ``None`` for frames that did not converge.
            - ``prop_all``: list of ``(frame, property_results)`` tuples
              (only present when a property driver is used; ``property_results``
              is ``None`` for frames that did not converge).
            - ``resp_all``: list of ``(frame, charges)`` tuples where
              ``charges`` is a 1-D NumPy array of length *natoms*, or
              ``None`` for non-converged frames.  Only present when RESP is
              requested.
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
                prop_driver = ComplexResponseSolver(self.comm, self.ostream)
            if prop_driver is not None:
                prop_driver.update_settings(prop_opts)
                if isinstance(prop_driver, ComplexResponseSolver) and "property" in prop_opts:
                    prop_driver.set_cpp_property(prop_opts["property"])

        # ── Build RESP driver from resp_options if not provided ───────────────
        if resp_driver is None and resp_options is not None:
            resp_driver = RespChargesDriver(self.comm, self.ostream)
            if resp_options:
                resp_driver.update_settings(resp_options)

        do_resp = resp_driver is not None
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
        resp_all = [] if do_resp else None

        # ── Determine DFT label once (needed before the loop for H5 init) ─────
        _xcfun = getattr(scf_driver, "xcfun", None)
        if hasattr(_xcfun, "get_func_label"):
            h5_dft_func_label = _xcfun.get_func_label()
        elif _xcfun:
            h5_dft_func_label = str(_xcfun)
        else:
            h5_dft_func_label = "HF"

        # ── Open H5 file before the loop so writes are incremental ────────────
        h5_path = Path(stacked_h5_file) if stacked_h5_file else None
        if h5_path is not None:
            self._h5_open_for_writing(
                h5_path,
                snapshots[0],
                fixed_basis,
                h5_dft_func_label,
                total_frames=len(snapshots),
            )

        # ── Wipe any stale guess from a previous run before starting ──────────
        if guess_h5 is not None and guess_h5.is_file():
            guess_h5.unlink()

        # ── Main trajectory loop ───────────────────────────────────────────────
        n_total = len(snapshots)
        for i_snap, snap in enumerate(snapshots):
            frame = int(snap["frame"])
            if self.rank == mpi_master():
                self.ostream.print_info(
                    f"Frame {i_snap + 1}/{n_total}"
                )
                self.ostream.print_blank()

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
                        f"Frame {i_snap + 1}/{n_total}: SCF did not converge "
                        "with the propagated guess.  Retrying with SAD initial "
                        "guess..."
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
                        f"Frame {i_snap + 1}/{n_total}: SCF did not converge — "
                        "recording as non-converged frame."
                    )
                    self.ostream.print_blank()
                scf_all.append((frame, None))
                scf_history_all.append([])
                if do_rsp:
                    rsp_all.append((frame, None))
                if do_resp:
                    resp_all.append((frame, None))
                if h5_path is not None:
                    self._h5_append_frame(
                        h5_path,
                        frame,
                        snap,
                        molecule,
                        fixed_basis,
                        h5_dft_func_label,
                        scf_results=None,
                        scf_history=[],
                        rsp_results=None,
                        prop_driver=rsp_driver if do_rsp else None,
                    )
                continue

            scf_all.append((frame, scf_results))
            scf_history_all.append(
                list(scf_driver._history) if scf_driver._history else []
            )

            rsp_results = None
            if do_rsp:
                rsp_results = rsp_driver.compute(
                    molecule, fixed_basis, scf_results
                )
                rsp_all.append((frame, rsp_results))

            if do_resp:
                q = resp_driver.compute(molecule, fixed_basis, scf_results)
                # RespChargesDriver also mutates scf_results['charges_resp'],
                # which is then picked up by _get_scf_hdf5_payload_local.
                if self.rank == mpi_master():
                    resp_all.append((frame, np.asarray(q, dtype=float)))
                else:
                    resp_all.append((frame, None))

            if h5_path is not None:
                self._h5_append_frame(
                    h5_path,
                    frame,
                    snap,
                    molecule,
                    fixed_basis,
                    h5_dft_func_label,
                    scf_results=scf_results,
                    scf_history=scf_history_all[-1],
                    rsp_results=rsp_results,
                    prop_driver=rsp_driver if do_rsp else None,
                )

        # ── End-of-run summary ────────────────────────────────────────────────
        if self.rank == mpi_master():
            non_conv = [
                i_snap + 1
                for i_snap, (_, res) in enumerate(scf_all)
                if res is None
            ]
            n_conv = n_total - len(non_conv)
            self.ostream.print_info(
                f"QMTrajectoryDriver: {n_total} frames processed, "
                f"{n_conv} converged, {len(non_conv)} non-converged."
            )
            if non_conv:
                self.ostream.print_info(
                    "Non-converged frames (1-based): "
                    + ", ".join(str(f) for f in non_conv)
                )
            self.ostream.print_blank()

        results = {"scf_all": scf_all, "scf_history_all": scf_history_all}
        if do_rsp:
            results["prop_all"] = rsp_all
        if do_resp:
            results["resp_all"] = resp_all

        return results

    def plot_scf(self, results, y_unit='a.u.', ax=None):
        """
        Plot the SCF energy **relative to the minimum energy frame** as a
        function of trajectory frame.

        Non-converged frames are omitted from the line but highlighted as
        red crosses at the bottom of the plot so that gaps are visible.

        :param results:
            The dictionary returned by :meth:`compute`.
        :param y_unit:
            Unit for the y-axis energy values.  Accepted values (case-
            insensitive): ``'a.u.'`` (Hartree, default), ``'ev'``,
            ``'kj/mol'``.
        :param ax:
            Optional matplotlib ``Axes`` to plot into.  If ``None`` a new
            figure is created.
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError(
                "matplotlib is required for plot_scf. "
                "Install it with: conda install -c conda-forge matplotlib"
            )

        assert_msg_critical(
            y_unit.lower() in ['a.u.', 'ev', 'kj/mol'],
            "plot_scf: Invalid y_unit. Choose 'a.u.', 'ev', or 'kj/mol'."
        )

        _unit = y_unit.lower()
        if _unit == 'ev':
            _factor = hartree_in_ev()
            _ylabel = 'Relative SCF Energy [eV]'
        elif _unit == 'kj/mol':
            _factor = hartree_in_kjpermol()
            _ylabel = 'Relative SCF Energy [kJ mol$^{-1}$]'
        else:
            _factor = 1.0
            _ylabel = 'Relative SCF Energy [Hartree]'

        frames_conv, energies_conv = [], []
        frames_nonconv = []

        for frame, scf_res in results['scf_all']:
            if scf_res is not None:
                frames_conv.append(frame)
                energies_conv.append(float(scf_res['scf_energy']) * _factor)
            else:
                frames_nonconv.append(frame)

        if energies_conv:
            e_min = min(energies_conv)
            energies_conv = [e - e_min for e in energies_conv]

        if ax is None:
            fig, ax = plt.subplots(figsize=(6.5, 4))
        else:
            fig = None

        if frames_conv:
            ax.plot(
                frames_conv,
                energies_conv,
                color='black',
                alpha=0.9,
                linewidth=2.5,
                ls='-',
                zorder=0,
            )
            ax.scatter(
                frames_conv,
                energies_conv,
                facecolors='none',
                edgecolors='darkcyan',
                linewidths=1.5,
                s=40,
                zorder=1,
            )

        if frames_nonconv:
            y_min = 0.0
            ax.scatter(
                frames_nonconv,
                [y_min] * len(frames_nonconv),
                marker='x',
                color='red',
                s=60,
                zorder=2,
                label='non-converged',
            )
            ax.legend(frameon=False)

        ax.set_xlabel('Frame')
        ax.set_ylabel(_ylabel)
        ax.set_title('SCF Energy along Trajectory')
        if fig is not None:
            fig.tight_layout()
            plt.show()

    def show_trajectory(self, frames, mode='animate', stride=1, width=600,
                        height=400, interval=100, loop='forward'):
        """
        Display a list of parsed frames using py3Dmol.

        :param frames:
            List of frame dicts as returned by :class:`QMTrajectoryParser`.
        :param mode:
            Visualisation mode: ``'animate'`` plays frames as a movie;
            ``'superimpose'`` overlays all selected frames simultaneously.
        :param stride:
            Use every ``stride``-th frame (default: 1, i.e. all frames).
        :param width:
            Viewer width in pixels.
        :param height:
            Viewer height in pixels.
        :param interval:
            Playback delay between frames in milliseconds (``'animate'`` only).
        :param loop:
            Animation loop mode: ``'forward'`` or ``'backAndForth'``
            (``'animate'`` only).
        """
        from .qmtrajectoryanalyzer import show_trajectory as _show_trajectory
        _show_trajectory(frames, mode=mode, stride=stride, width=width,
                         height=height, interval=interval, loop=loop)

