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

import sys
from pathlib import Path

import h5py
import numpy as np

from .veloxchemlib import mpi_master, hartree_in_ev, hartree_in_kjpermol
from .veloxchemlib import bohr_in_angstrom
from .molecule import Molecule
from .molecularbasis import MolecularBasis
from .errorhandler import assert_msg_critical
from .spectrumplot import plot_trajectory_spectrum


def show_trajectory(snapshots, mode='animate', stride=1, width=600, height=400,
                    interval=100, loop='forward'):
    """
    Display a list of snapshot dicts as a 3D trajectory using py3Dmol.

    Works with the output of :meth:`QMTrajectoryParser.structures` or
    :meth:`QMTrajectoryAnalyzer.get_snapshots`.

    :param snapshots:
        List of snapshot dicts, each containing ``'qm_elements'`` and
        ``'qm_coords'`` (coordinates in Angstrom).
    :param mode:
        Visualisation mode: ``'animate'`` plays frames as a movie;
        ``'superimpose'`` overlays all selected frames simultaneously.
    :param stride:
        Use every ``stride``-th frame (default: 1, i.e. all frames).
        For example, ``stride=5`` uses frames 0, 5, 10, ...
    :param width:
        Viewer width in pixels.
    :param height:
        Viewer height in pixels.
    :param interval:
        Playback delay between frames in milliseconds (``'animate'`` mode only).
    :param loop:
        Animation loop mode: ``'forward'`` or ``'backAndForth'``
        (``'animate'`` mode only).
    """
    assert_msg_critical(
        len(snapshots) > 0,
        "show_trajectory: snapshots list is empty."
    )
    assert_msg_critical(
        isinstance(stride, int) and stride >= 1,
        "show_trajectory: stride must be a positive integer."
    )
    assert_msg_critical(
        mode in ('animate', 'superimpose'),
        "show_trajectory: mode must be 'animate' or 'superimpose'."
    )
    try:
        import py3Dmol
    except ImportError:
        raise ImportError(
            "py3Dmol is required for show_trajectory. "
            "Install it with: conda install -c conda-forge py3dmol"
        )

    selected = snapshots[::stride]

    # Build a single concatenated multi-frame XYZ string.
    # py3Dmol's addModels() parses this as separate frames that can be
    # animated or superimposed.  Using addModel() in a loop does NOT work
    # for animation — it only superimposes.
    xyz_blocks = []
    for snap in selected:
        elems = snap['qm_elements']
        coords = snap['qm_coords']  # Angstrom
        n = len(elems)
        lines = [str(n), f"frame {snap.get('frame', '')}"]
        for el, coord in zip(elems, coords):
            lines.append(f"{el}  {coord[0]:.6f}  {coord[1]:.6f}  {coord[2]:.6f}")
        xyz_blocks.append('\n'.join(lines))

    viewer = py3Dmol.view(width=width, height=height)

    if mode == 'animate':
        viewer.addModelsAsFrames('\n'.join(xyz_blocks), 'xyz')
        viewer.setStyle({}, {"stick": {"radius": 0.15}, "sphere": {"scale": 0.30}})
        viewer.animate({'loop': loop, 'reps': 0, 'interval': interval})
    else:  # superimpose
        viewer.addModels('\n'.join(xyz_blocks), 'xyz')
        viewer.setStyle({}, {"stick": {"radius": 0.15}, "sphere": {"scale": 0.30}})

    viewer.zoomTo()
    viewer.show()


class QMTrajectoryAnalyzer:
    """
    Post-processing and analysis tool for trajectory HDF5 files produced by
    :class:`QMTrajectoryDriver`.

    Typical workflow::

        analyzer = vlx.QMTrajectoryAnalyzer("enzyme.h5")
        analyzer.summary()
        analyzer.plot_scf(y_unit='kj/mol')

        # Re-run non-converged frames with looser thresholds
        drv = vlx.QMTrajectoryDriver()
        analyzer.recompute_frames(
            frame_ids=analyzer.non_converged_frames,
            scf_options={"conv_thresh": 1e-3, "max_iter": 300},
            driver=drv,
            update_h5=True,
        )

    :param h5_path:
        Path to the trajectory HDF5 file.
    """

    def __init__(self, h5_path):
        self._path = Path(h5_path)
        assert_msg_critical(
            self._path.is_file(),
            f"QMTrajectoryAnalyzer: file not found: {self._path}"
        )
        # Cache lightweight arrays on construction; heavy data read on demand.
        with h5py.File(self._path, 'r') as hf:
            self._frame_ids = np.asarray(hf['frame_id'], dtype=int)
            self._converged = np.asarray(hf['converged'], dtype=bool)
            self._n_frames = len(self._frame_ids)
            self._basis_label = hf['basis_set'][0].decode().strip() if 'basis_set' in hf else ''
            self._dft_label = hf['dft_func_label'][0].decode().strip() if 'dft_func_label' in hf else 'HF'
            self._trajectory_name = hf['trajectory_name'][0].decode().strip() if 'trajectory_name' in hf else ''
            self._has_rsp = 'rsp' in hf
            self._has_resp = 'scf' in hf and 'charges_resp' in hf.get('scf', {})

    # ── Properties ─────────────────────────────────────────────────────────────

    @property
    def frame_ids(self):
        """numpy.ndarray of trajectory frame ids (1 per row written)."""
        return self._frame_ids.copy()

    @property
    def converged_mask(self):
        """Boolean numpy.ndarray, True where SCF converged."""
        return self._converged.copy()

    @property
    def non_converged_frames(self):
        """List of frame ids where SCF did not converge."""
        return list(self._frame_ids[~self._converged])

    @property
    def converged_frames(self):
        """List of frame ids where SCF converged."""
        return list(self._frame_ids[self._converged])

    # ── Data accessors ──────────────────────────────────────────────────────────

    def scf_energies(self, y_unit='a.u.'):
        """
        Return per-frame SCF energies as a numpy array.

        Non-converged frames are filled with ``numpy.nan``.

        :param y_unit:
            Energy unit: ``'a.u.'`` (Hartree), ``'ev'``, or ``'kj/mol'``.

        :return:
            numpy.ndarray of shape ``(n_frames,)``.
        """
        assert_msg_critical(
            y_unit.lower() in ['a.u.', 'ev', 'kj/mol'],
            "scf_energies: Invalid y_unit. Choose 'a.u.', 'ev', or 'kj/mol'."
        )
        if y_unit.lower() == 'ev':
            factor = hartree_in_ev()
        elif y_unit.lower() == 'kj/mol':
            factor = hartree_in_kjpermol()
        else:
            factor = 1.0

        energies = np.full(self._n_frames, np.nan)
        with h5py.File(self._path, 'r') as hf:
            if 'scf' in hf and 'scf_energy' in hf['scf']:
                raw = np.asarray(hf['scf']['scf_energy'], dtype=float).ravel()
                n = min(len(raw), self._n_frames)
                # scf group only has rows for converged frames — align via mask
                conv_indices = np.where(self._converged)[0]
                for j, row_idx in enumerate(conv_indices):
                    if j < len(raw):
                        energies[row_idx] = raw[j]
        return energies * factor

    def get_resp_charges(self, frame_ids=None):
        """
        Return RESP charges for the requested frames.

        Non-converged frames or frames without RESP data are filled with
        ``numpy.nan``.

        :param frame_ids:
            List of frame ids to extract.  Defaults to all
            :attr:`converged_frames`.

        :return:
            numpy.ndarray of shape ``(n_frames, n_atoms)``.
        :raises AssertionError:
            If no RESP charges are stored in the HDF5 file.
        """
        assert_msg_critical(
            self._has_resp,
            "get_resp_charges: no RESP charges found in the H5 file. "
            "Re-run with resp_options={} to compute them."
        )

        if frame_ids is None:
            frame_ids = list(self._frame_ids)

        with h5py.File(self._path, 'r') as hf:
            raw = np.asarray(hf['scf']['charges_resp'], dtype=float)

        # raw rows correspond only to converged frames
        conv_indices = np.where(self._converged)[0]
        natoms = raw.shape[1] if raw.ndim == 2 else raw.shape[0]

        all_charges = np.full((self._n_frames, natoms), np.nan)
        for j, row_idx in enumerate(conv_indices):
            if j < len(raw):
                all_charges[row_idx] = raw[j]

        # Select requested frames
        rows = [self._row_for_frame(fid) for fid in frame_ids]
        return all_charges[rows]

    def get_molecule(self, frame_id):
        """
        Reconstruct a :class:`Molecule` from the stored nuclear charges and
        coordinates for the given frame id.

        :param frame_id:
            The trajectory frame id (as in :attr:`frame_ids`).

        :return:
            :class:`Molecule` object with coordinates in Bohr.
        """
        row = self._row_for_frame(frame_id)
        with h5py.File(self._path, 'r') as hf:
            charges = np.asarray(hf['nuclear_charges'][row], dtype=int)
            coords = np.asarray(hf['atom_coordinates'][row], dtype=float)
            charge = int(np.asarray(hf['molecular_charge'][row]))
            mult = int(np.asarray(hf['spin_multiplicity'][row]))

        mol = Molecule(charges, coords, units='bohr')
        return mol

    def get_snapshots(self, frame_ids=None):
        """
        Reconstruct snapshot dictionaries for the given frame ids.

        The returned list has the same format as :meth:`QMTrajectoryParser.structures`
        output and can be passed directly to :meth:`QMTrajectoryDriver.compute`.

        :param frame_ids:
            List of frame ids to extract.  Defaults to all frames.

        :return:
            List of snapshot dicts.
        """
        if frame_ids is None:
            frame_ids = list(self._frame_ids)

        snapshots = []
        ang_per_bohr = bohr_in_angstrom()
        total = len(frame_ids)

        with h5py.File(self._path, 'r') as hf:
            for seq_idx, fid in enumerate(frame_ids):
                row = self._row_for_frame(fid)
                charges = np.asarray(hf['nuclear_charges'][row], dtype=int)
                coords_bohr = np.asarray(hf['atom_coordinates'][row], dtype=float)
                charge = int(np.asarray(hf['molecular_charge'][row]))
                mult = int(np.asarray(hf['spin_multiplicity'][row]))

                # Molecule() accepts nuclear charges directly (same as valetanalyzer.py)
                mol_tmp = Molecule(charges, coords_bohr, units='bohr')
                labels = np.asarray(
                    [mol_tmp.get_label(i) for i in range(mol_tmp.number_of_atoms())],
                    dtype=object
                )
                # Coordinates in Angstrom (parser convention)
                coords_ang = coords_bohr * ang_per_bohr

                snapshots.append({
                    'frame': int(fid),
                    'frame_index': seq_idx,
                    'trajectory_name': self._trajectory_name,
                    'total_frames': total,
                    'qm_coords': coords_ang,
                    'qm_elements': labels,
                    'qm_charge': charge,
                    'qm_multiplicity': mult,
                    'qm_atom_indices': np.arange(len(labels), dtype=int),
                })
        return snapshots

    def _get_rsp_frame_results(self, frame_ids=None):
        """Collect per-frame response dictionaries from the stacked HDF5 file."""
        assert_msg_critical(
            self._has_rsp,
            'plot_spectrum: no response data found in the H5 file.'
        )

        if frame_ids is None:
            frame_ids = list(self._frame_ids)

        conv_rows = {
            int(fid): i for i, fid in enumerate(self._frame_ids[self._converged])
        }
        frame_results = []

        with h5py.File(self._path, 'r') as hf:
            rsp_grp = hf['rsp']
            keys = list(rsp_grp.keys())
            for fid in frame_ids:
                fid = int(fid)
                if fid not in conv_rows:
                    continue
                row = conv_rows[fid]
                rsp_res = {}
                for key in keys:
                    dset = rsp_grp[key]
                    if dset.ndim == 0:
                        rsp_res[key] = dset[()]
                    elif dset.shape[0] > row:
                        rsp_res[key] = np.asarray(dset[row])
                frame_results.append((fid, rsp_res))

        return frame_results

    # ── Diagnostics ─────────────────────────────────────────────────────────────

    def summary(self):
        """
        Print a human-readable summary of the trajectory H5 file.
        """
        n_conv = int(self._converged.sum())
        n_nonconv = self._n_frames - n_conv
        pct = 100.0 * n_conv / self._n_frames if self._n_frames else 0.0

        energies = self.scf_energies(y_unit='a.u.')
        e_conv = energies[~np.isnan(energies)]

        lines = [
            f"Trajectory H5 file : {self._path}",
            f"Trajectory         : {self._trajectory_name or '(unknown)'}",
            f"Basis set          : {self._basis_label or '(unknown)'}",
            f"DFT functional     : {self._dft_label or 'HF'}",
            f"Total frames       : {self._n_frames}",
            f"Converged          : {n_conv}  ({pct:.1f} %)",
            f"Non-converged      : {n_nonconv}",
        ]
        if n_nonconv:
            nc_list = ', '.join(str(f) for f in self.non_converged_frames)
            lines.append(f"  → frame ids      : {nc_list}")
        if e_conv.size:
            lines += [
                f"Energy min         : {e_conv.min():.6f} Hartree",
                f"Energy max         : {e_conv.max():.6f} Hartree",
                f"Energy range       : {e_conv.max() - e_conv.min():.6f} Hartree",
            ]
        if self._has_rsp:
            lines.append("Response data      : present")
        if self._has_resp:
            lines.append("RESP charges       : present")

        print('\n'.join(lines))

    # ── Plotting ────────────────────────────────────────────────────────────────

    def plot_scf(self, y_unit='a.u.', ax=None):
        """
        Plot the SCF energy **relative to the minimum energy frame** as a
        function of trajectory frame, reading directly from the HDF5 file.

        Non-converged frames are shown as red crosses at the bottom of the
        plot so that gaps are visible.

        :param y_unit:
            Energy unit: ``'a.u.'`` (Hartree, default), ``'ev'``, or
            ``'kj/mol'``.
        :param ax:
            Optional matplotlib ``Axes``.  If ``None`` a new figure is
            created.
        """
        assert_msg_critical(
            y_unit.lower() in ['a.u.', 'ev', 'kj/mol'],
            "plot_scf: Invalid y_unit. Choose 'a.u.', 'ev', or 'kj/mol'."
        )
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError(
                "matplotlib is required for plot_scf. "
                "Install it with: conda install -c conda-forge matplotlib"
            )

        if y_unit.lower() == 'ev':
            _ylabel = 'Relative SCF Energy [eV]'
        elif y_unit.lower() == 'kj/mol':
            _ylabel = 'Relative SCF Energy [kJ mol$^{-1}$]'
        else:
            _ylabel = 'Relative SCF Energy [Hartree]'

        energies = self.scf_energies(y_unit=y_unit)
        frames_conv = self._frame_ids[~np.isnan(energies)].tolist()
        energies_conv = energies[~np.isnan(energies)].tolist()
        frames_nonconv = list(self._frame_ids[np.isnan(energies)])

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

    def plot_spectrum(self,
                      property=None,
                      x_unit='nm',
                      broadening_type='lorentzian',
                      broadening_value=None,
                      x_range=None,
                      frame_ids=None,
                      show_frames=True,
                      show_average=True,
                      show_sticks=True,
                      ax=None):
        """Plot one stored spectrum per frame together with the average."""
        if broadening_value is None:
            from .veloxchemlib import hartree_in_wavenumber
            broadening_value = (1000.0 / hartree_in_wavenumber() *
                                hartree_in_ev())

        frame_results = self._get_rsp_frame_results(frame_ids)
        return plot_trajectory_spectrum(
            frame_results,
            property=property,
            x_unit=x_unit,
            broadening_type=broadening_type,
            broadening_value=broadening_value,
            x_range=x_range,
            show_frames=show_frames,
            show_average=show_average,
            show_sticks=show_sticks,
            ax=ax,
        )

    def show_trajectory(self, frame_range=None, mode='animate', stride=1,
                        width=600, height=400, interval=100, loop='forward'):
        """
        Display the trajectory stored in the HDF5 file using py3Dmol.

        :param frame_range:
            ``(start, end)`` tuple of frame ids (both inclusive) to restrict
            the view.  ``None`` (default) uses all frames.
        :param mode:
            Visualisation mode: ``'animate'`` plays frames as a movie;
            ``'superimpose'`` overlays all selected frames simultaneously.
        :param stride:
            Use every ``stride``-th frame in the selected range (default: 1).
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
        if frame_range is not None:
            start, end = frame_range
            frame_ids = [int(fid) for fid in self._frame_ids
                         if start <= int(fid) <= end]
        else:
            frame_ids = None
        snapshots = self.get_snapshots(frame_ids)
        show_trajectory(snapshots, mode=mode, stride=stride, width=width,
                        height=height, interval=interval, loop=loop)

    # ── Recompute ────────────────────────────────────────────────────────────────

    def recompute_frames(
        self,
        frame_ids=None,
        scf_options=None,
        basis_set=None,
        driver=None,
        update_h5=True,
        resp_options=None,
        resp_driver=None,
    ):
        """
        Re-run SCF (and optionally RESP charges) for the specified frames and
        optionally patch the results back into the HDF5 file.

        :param frame_ids:
            List of frame ids to recompute.  Defaults to
            :attr:`non_converged_frames`.
        :param scf_options:
            Dict of SCF options to override (e.g. looser thresholds).
        :param basis_set:
            Basis set label to use.  Defaults to the value stored in the H5.
        :param driver:
            A :class:`QMTrajectoryDriver` instance to use.  If ``None``,
            a fresh one is created.
        :param update_h5:
            If ``True``, patch the converged results back into the HDF5 file.
        :param resp_options:
            Dict of RESP charge options forwarded to
            :class:`RespChargesDriver`.  Pass ``{}`` to use defaults.
            If ``None`` (default), RESP charges are not recomputed.
        :param resp_driver:
            Pre-configured :class:`RespChargesDriver` instance.  When
            provided, ``resp_options`` is ignored.

        :return:
            The ``results`` dict from :meth:`QMTrajectoryDriver.compute`.
        """
        from .qmtrajectorydriver import QMTrajectoryDriver

        if frame_ids is None:
            frame_ids = self.non_converged_frames
        if not frame_ids:
            print("QMTrajectoryAnalyzer.recompute_frames: no frames to recompute.")
            return {}

        if driver is None:
            driver = QMTrajectoryDriver()

        _basis = basis_set or self._basis_label
        assert_msg_critical(
            _basis,
            "QMTrajectoryAnalyzer.recompute_frames: basis_set not specified "
            "and could not be read from the H5 file."
        )

        snapshots = self.get_snapshots(frame_ids)
        results = driver.compute(
            snapshots,
            basis_set=_basis,
            scf_options=scf_options or {},
            resp_options=resp_options,
            resp_driver=resp_driver,
            stacked_h5_file=None,       # don't create a new H5
            propagate_guess=False,      # fresh SAD for each frame
        )

        if update_h5:
            self._patch_h5(frame_ids, snapshots, results, driver, _basis)
            # Refresh cached arrays
            with h5py.File(self._path, 'r') as hf:
                self._converged = np.asarray(hf['converged'], dtype=bool)
                self._has_resp = 'scf' in hf and 'charges_resp' in hf.get('scf', {})

        return results

    # ── Internal helpers ────────────────────────────────────────────────────────

    def _row_for_frame(self, frame_id):
        """Return the H5 row index for a given frame id."""
        matches = np.where(self._frame_ids == int(frame_id))[0]
        assert_msg_critical(
            len(matches) > 0,
            f"QMTrajectoryAnalyzer: frame id {frame_id} not found in {self._path}"
        )
        return int(matches[0])

    def _patch_h5(self, frame_ids, snapshots, results, driver, basis_label):
        """
        Write newly converged results back into the existing H5 file,
        overwriting the NaN placeholders for those rows.
        """
        from .qmtrajectorydriver import QMTrajectoryDriver

        scf_all = {frame: res for frame, res in results.get('scf_all', [])}
        resp_all = {frame: q for frame, q in results.get('resp_all', [])}
        scf_history_map = {}
        for idx, (frame, _) in enumerate(results.get('scf_all', [])):
            hist = results.get('scf_history_all', [])
            if idx < len(hist):
                scf_history_map[frame] = hist[idx]

        fixed_basis = MolecularBasis.read(
            driver._build_reference_molecule(snapshots[0], basis_label)
            if hasattr(driver, '_build_reference_molecule')
            else self.get_molecule(frame_ids[0]),
            basis_label,
        )

        # Determine DFT label
        _xcfun = getattr(getattr(driver, 'scf_driver', None) or object(), 'xcfun', None)
        if hasattr(_xcfun, 'get_func_label'):
            dft_label = _xcfun.get_func_label()
        elif _xcfun:
            dft_label = str(_xcfun)
        else:
            dft_label = self._dft_label or 'HF'

        with h5py.File(self._path, 'a') as hf:
            # Root-level datasets are frame-aligned, but scf/ datasets are aligned
            # only to the subset of frames that were already converged in the file.
            original_converged = np.asarray(hf['converged'][:], dtype=bool)
            converged_rows = np.flatnonzero(original_converged)
            converged_row_map = {
                int(full_row): idx for idx, full_row in enumerate(converged_rows)
            }

            for snap, fid in zip(snapshots, frame_ids):
                full_row = self._row_for_frame(fid)
                scf_res = scf_all.get(int(fid))
                converged = scf_res is not None

                # Update root-level converged flag using the full trajectory row.
                hf['converged'][full_row] = converged

                if not converged:
                    continue

                # Reuse `row` below as the index into converged-only datasets
                # (for example scf/ and downstream response payloads).
                row = converged_row_map.get(int(full_row))

                # Newly converged frames do not have a pre-existing scf/ row in
                # files written without placeholders, so avoid misaligned writes.
                if row is None:
                    continue
                # Patch scf group — only scf_energy and a few key scalars
                if 'scf' not in hf:
                    continue
                scf_grp = hf['scf']
                payload = QMTrajectoryDriver._get_scf_hdf5_payload_local(
                    scf_res,
                    scf_history_map.get(int(fid), [])
                )
                for key, value in payload.items():
                    if key not in scf_grp:
                        continue
                    dset = scf_grp[key]
                    if row < dset.shape[0]:
                        arr = QMTrajectoryDriver._result_value_to_array(value)
                        if arr is None:
                            continue
                        # Handle variable-width history arrays
                        if arr.ndim >= 1 and dset.ndim == arr.ndim + 1:
                            cur_w = dset.shape[1]
                            new_w = arr.shape[0]
                            if new_w < cur_w:
                                padded = np.full(cur_w, np.nan, dtype=dset.dtype)
                                padded[:new_w] = arr
                                arr = padded
                            elif new_w > cur_w:
                                arr = arr[:cur_w]
                        try:
                            dset[row] = arr
                        except Exception:
                            pass  # shape mismatch — skip silently

                # Patch RESP charges if present in recompute results
                q = resp_all.get(int(fid))
                if q is not None:
                    q_arr = np.asarray(q, dtype=float)
                    scf_grp = hf.require_group('scf')
                    if 'charges_resp' not in scf_grp:
                        natoms = q_arr.shape[0]
                        n_total = len(self._frame_ids)
                        placeholder = np.full((n_total, natoms), np.nan)
                        scf_grp.create_dataset(
                            'charges_resp',
                            data=placeholder,
                            maxshape=(None, natoms),
                        )
                        dset_q = scf_grp['charges_resp']
                        from .qmtrajectorydriver import QMTrajectoryDriver as _D
                        _D._apply_atomic_metadata(
                            _D, dset_q, 'scf', 'charges_resp'
                        )
                    dset_q = scf_grp['charges_resp']
                    if row < dset_q.shape[0]:
                        try:
                            dset_q[row] = q_arr
                        except Exception:
                            pass
