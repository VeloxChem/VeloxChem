#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2026 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice,
#     this list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software
#     without specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
#  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
#  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
#  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#  POSSIBILITY OF SUCH DAMAGE.

"""Post-processing and plotting for relaxed constrained scans.

A relaxed excited-state scan is only trustworthy if three separate things are
checked, and the energy profile alone shows none of them:

* the constraint was actually enforced at every converged point;
* the tracked state stayed on the same physical surface, i.e. the gradient was
  never taken on S0 (see :mod:`openqpstatetracker`);
* the relaxed coordinates stayed chemically sensible instead of escaping the
  constraint -- a torsion whose central angle opens to 180 degrees is the
  classic failure and is what geomeTRIC reports as ``LinearTorsionError``.

:func:`collect_scan_data` extracts all three from an
``OptimizationDriver.compute`` result plus, optionally, the tracking history of
an OpenQP or Serenity excited-state gradient driver.  :func:`plot_scan` renders
them.  Both work from plain arrays, so a scan can also be re-plotted later from
a saved JSON payload without repeating the calculation.
"""

from pathlib import Path
import json

import numpy as np

from .veloxchemlib import bohr_in_angstrom

HARTREE_IN_EV = 27.211386245988
HARTREE_IN_KCALMOL = 627.509474063

# Chart palette: the first three categorical slots, which are the ones
# validated for all pairs in both light and dark mode, plus recessive chrome.
_SERIES = ('#2a78d6', '#eb6834', '#1baf7a')
_INK = '#0b0b0b'
_MUTED = '#898781'
_GRID = '#e1e0d9'
_WARN = '#d03b3b'


def _xyz_to_array(xyz_string):
    """Parses one XYZ block into labels and coordinates in Angstrom."""

    lines = str(xyz_string).strip().splitlines()
    natoms = int(lines[0].split()[0])
    labels = []
    coordinates = []
    for line in lines[2:2 + natoms]:
        fields = line.split()
        labels.append(fields[0])
        coordinates.append([float(value) for value in fields[1:4]])

    return labels, np.array(coordinates, dtype=float)


def _unwrap_degrees(values):
    """Removes the +-180 degree branch cuts from an ordered angle series."""

    return [float(value) for value in
            np.degrees(np.unwrap(np.radians(np.asarray(values, dtype=float))))]


def _internal(coordinate, xyz, atoms):
    """Distance in Angstrom, or angle/dihedral in degrees, from Angstrom xyz."""

    points = [xyz[index] for index in atoms]

    if coordinate == 'distance':
        return float(np.linalg.norm(points[1] - points[0]))

    if coordinate == 'angle':
        u = points[0] - points[1]
        v = points[2] - points[1]
        cosine = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
        return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))

    b0 = points[0] - points[1]
    b1 = points[2] - points[1]
    b2 = points[3] - points[2]
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w),
                                       np.dot(v, w))))


def pyramidalization(xyz, center, substituents):
    """Improper angle at ``center``, in degrees.

    Zero for a planar sp2 center; it grows as the center pyramidalizes, which
    is the coordinate that carries a twisted iminium through its S1/S0 conical
    intersection.  Defined as the angle between the bond to the first
    substituent and the plane of the remaining two.

    :param xyz:
        Coordinates in Angstrom.
    :param center:
        Zero-based index of the central atom.
    :param substituents:
        Three zero-based indices of the bonded atoms.
    """

    a, b, c = (xyz[index] - xyz[center] for index in substituents)
    normal = np.cross(b, c)
    norm = np.linalg.norm(normal) * np.linalg.norm(a)
    if norm < 1.0e-12:
        return 0.0
    sine = np.dot(normal, a) / norm

    return float(abs(np.degrees(np.arcsin(np.clip(sine, -1.0, 1.0)))))


def collect_scan_data(opt_results,
                      scan_atoms,
                      scan_coordinate='dihedral',
                      tracking_history=None,
                      geometry_probes=None):
    """Reduces a scan result to the arrays the plots and checks need.

    :param opt_results:
        The dictionary returned by ``OptimizationDriver.compute`` for a job
        whose constraints contained a ``scan`` entry.
    :param scan_atoms:
        One-based atom indices of the scanned coordinate, as written in the
        constraint string.
    :param scan_coordinate:
        ``'dihedral'``, ``'angle'`` or ``'distance'``.
    :param tracking_history:
        Optional ``gradient_driver.tracking_history``.  Only the entries at the
        converged geometries are used, so the state energies line up with the
        scan points.  Entries must carry the ``scan_point`` stamp that
        ``OptimizationEngine`` adds during a constrained scan; without it the
        state data are omitted rather than misaligned.
    :param geometry_probes:
        Optional ``{label: (coordinate, one-based atom indices)}`` of extra
        internal coordinates to follow along the scan.

    :return:
        A dictionary of plain lists and arrays, JSON-serializable via
        :func:`save_scan_data`.
    """

    geometries = opt_results.get('scan_geometries')
    energies = opt_results.get('scan_energies')
    if not geometries or energies is None:
        raise ValueError(
            'collect_scan_data: the optimization result carries no scan data. '
            'Pass the result of a job whose constraints contained a "scan" '
            'entry.')

    atoms = [index - 1 for index in scan_atoms]
    parsed = [_xyz_to_array(xyz) for xyz in geometries]
    labels = parsed[0][0]
    coordinates = [xyz for _, xyz in parsed]

    raw_values = [_internal(scan_coordinate, xyz, atoms)
                  for xyz in coordinates]

    data = {
        'labels': labels,
        'scan_coordinate': scan_coordinate,
        'scan_atoms': list(scan_atoms),
        # A dihedral comes back in (-180, 180], so a scan that runs to 180
        # degrees ends at -180 and the profile jumps back across the whole
        # axis.  The grid points are in scan order, so unwrapping restores the
        # coordinate the scan actually drove.
        'scan_values': (_unwrap_degrees(raw_values)
                        if scan_coordinate == 'dihedral' else
                        [float(value) for value in raw_values]),
        'energies': [float(value) for value in energies],
        'geometries': list(geometries),
        'probes': {},
        'probe_kinds': {},
    }

    # A torsion is only well defined while both of its central angles stay away
    # from 180 degrees; carrying them makes the LinearTorsionError visible as a
    # trend rather than as a traceback.
    if scan_coordinate == 'dihedral':
        for triple in (atoms[0:3], atoms[1:4]):
            label = 'central angle %d-%d-%d' % tuple(
                index + 1 for index in triple)
            data['probes'][label] = [
                _internal('angle', xyz, triple) for xyz in coordinates
            ]
            data['probe_kinds'][label] = 'angle'

    for label, (coordinate, probe_atoms) in (geometry_probes or {}).items():
        probe_indices = [index - 1 for index in probe_atoms]
        series = [_internal(coordinate, xyz, probe_indices)
                  for xyz in coordinates]
        data['probes'][label] = (_unwrap_degrees(series)
                                 if coordinate == 'dihedral' else series)
        data['probe_kinds'][label] = coordinate

    if tracking_history:
        data.update(_collect_tracking(tracking_history, len(coordinates)))

    return data


def _tracking_field(entry, name):
    """Reads a field from a result object or a restored dictionary."""

    if isinstance(entry, dict):
        return entry.get(name)

    return getattr(entry, name, None)


def _collect_tracking(tracking_history, n_points):
    """Aligns the tracking history with the converged scan points.

    The history has one entry per accepted electronic evaluation, which is many
    more than the number of scan points -- a five-point scan easily produces
    sixty.  The alignment therefore relies on the ``scan_point`` stamp that
    :class:`~veloxchem.optimizationengine.OptimizationEngine` puts on every
    evaluation, taking the *last* entry of each grid point: that is the one at
    the converged geometry.

    Without the stamp the entries cannot be attributed to grid points at all,
    and guessing produces a plot whose state energies belong to different
    geometries than its profile.  In that case nothing is returned.
    """

    grouped = {}
    for entry in tracking_history:
        stamp = _tracking_field(entry, 'scan_point')
        if stamp is None:
            return {}
        grouped[int(stamp)] = entry

    entries = [grouped[key] for key in sorted(grouped)]
    if len(entries) != n_points:
        return {}

    spectra = []
    selected = []
    collisions = []
    for entry in entries:
        energies = np.asarray(_tracking_field(entry, 'state_energies'),
                              dtype=float)
        root = int(_tracking_field(entry, 'selected_raw_root'))
        spectra.append(energies.tolist())
        selected.append(root)
        collisions.append(bool(
            _tracking_field(entry, 'ground_state_collision') or False))

    spectra_array = np.array(spectra, dtype=float)
    ground = spectra_array.min(axis=1)
    tracked = np.array([
        spectra_array[index, root - 1]
        for index, root in enumerate(selected)
    ])

    return {
        'state_energies': spectra_array.tolist(),
        'selected_raw_root': selected,
        'ground_state_collision': collisions,
        'gap_to_s0_ev': ((tracked - ground) * HARTREE_IN_EV).tolist(),
    }


def check_scan(data, linear_torsion_threshold=165.0):
    """Reports whether the scan is usable, as a list of message strings.

    :param data:
        The dictionary from :func:`collect_scan_data`.
    :param linear_torsion_threshold:
        Central-angle value in degrees beyond which a constrained torsion is
        treated as ill-defined.

    :return:
        ``(ok, messages)``.
    """

    messages = []
    ok = True

    # scan_values is already unwrapped by collect_scan_data, so the spacing is
    # the real grid spacing rather than a branch-cut artefact.
    values = np.asarray(data['scan_values'], dtype=float)
    spacing = np.abs(np.diff(values))
    if spacing.size and np.max(spacing) > 1.5 * np.median(spacing) + 1.0:
        messages.append(
            'Scan points are unevenly spaced along the constrained '
            'coordinate; a point may not have converged onto its target.')
        ok = False

    for label, series in data.get('probes', {}).items():
        if not label.startswith('central angle'):
            continue
        worst = float(np.max(series))
        if worst >= linear_torsion_threshold:
            messages.append(
                f'{label} reaches {worst:.1f} degrees: the constrained '
                'torsion is becoming undefined and the molecule is escaping '
                'the constraint.')
            ok = False

    collisions = data.get('ground_state_collision')
    if collisions and any(collisions):
        points = [str(i + 1) for i, flag in enumerate(collisions) if flag]
        messages.append(
            'The tracked state reached the S1/S0 seam at scan point(s) ' +
            ', '.join(points) + '. Energies there are the lowest excited '
            'adiabatic state, not a diabatically continued one.')

    gaps = data.get('gap_to_s0_ev')
    if gaps is not None:
        smallest = float(np.min(gaps))
        if smallest < 0.2:
            messages.append(
                f'Smallest tracked-state/S0 gap is {smallest:.3f} eV; the '
                'Born-Oppenheimer separation is not meaningful there.')
        if smallest <= 1.0e-8:
            messages.append(
                'The tracked state is degenerate with S0: the optimization '
                'was running on the ground state.')
            ok = False

    if ok and not messages:
        messages.append('All checks passed.')

    return ok, messages


def save_scan_data(data, filename):
    """Writes the collected scan data as JSON for later re-plotting."""

    with open(filename, 'w', encoding='utf-8') as handle:
        json.dump(data, handle, indent=2)


def load_scan_data(filename):
    """Reads a payload written by :func:`save_scan_data`."""

    with open(filename, encoding='utf-8') as handle:
        return json.load(handle)


def _style_axis(axis):
    """Applies the recessive grid/axis chrome the whole figure shares."""

    axis.grid(True, color=_GRID, linewidth=0.8, zorder=0)
    axis.set_axisbelow(True)
    for side in ('top', 'right'):
        axis.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        axis.spines[side].set_color(_GRID)
    axis.tick_params(colors=_MUTED, labelsize=9)
    axis.xaxis.label.set_color(_MUTED)
    axis.yaxis.label.set_color(_MUTED)
    axis.title.set_color(_INK)


def plot_scan(data,
              filename=None,
              energy_unit='kcal/mol',
              state_label='tracked state',
              show=False):
    """Plots the relaxed scan: profile, state spectrum, and geometry checks.

    :param data:
        The dictionary from :func:`collect_scan_data`.
    :param filename:
        Where to write the figure; ``None`` skips saving.
    :param energy_unit:
        ``'kcal/mol'``, ``'eV'`` or ``'hartree'`` for the relative profile.
    :param state_label:
        Name used for the followed state in the legend.
    :param show:
        Calls ``matplotlib.pyplot.show()`` when true.

    :return:
        The matplotlib figure.
    """

    import matplotlib
    if not show:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    scale, unit_text = {
        'kcal/mol': (HARTREE_IN_KCALMOL, 'kcal/mol'),
        'ev': (HARTREE_IN_EV, 'eV'),
        'hartree': (1.0, 'Hartree'),
    }[str(energy_unit).lower()]

    values = np.asarray(data['scan_values'], dtype=float)
    energies = np.asarray(data['energies'], dtype=float)
    relative = (energies - energies.min()) * scale

    axis_label = {
        'dihedral': 'dihedral %s (degrees)',
        'angle': 'angle %s (degrees)',
        'distance': 'distance %s (Angstrom)',
    }[data['scan_coordinate']] % '-'.join(
        str(index) for index in data['scan_atoms'])

    has_states = 'state_energies' in data
    n_panels = 3 if has_states else 2
    figure, axes = plt.subplots(n_panels, 1, figsize=(8.0, 3.2 * n_panels),
                                sharex=True)
    figure.patch.set_facecolor('#fcfcfb')

    # ---- Panel 1: the relaxed profile itself -----------------------------
    profile = axes[0]
    _style_axis(profile)
    profile.set_facecolor('#fcfcfb')
    profile.plot(values, relative, color=_SERIES[0], linewidth=2.0,
                 marker='o', markersize=5, zorder=3)
    profile.set_ylabel(f'relative energy ({unit_text})')
    profile.set_title(f'Relaxed scan on the {state_label}', loc='left',
                      fontsize=11, fontweight='bold')

    # Direct-label the extrema rather than every point.  The maximum is
    # labelled below its marker so the text cannot run into the panel title.
    for index, offset in ((int(np.argmax(relative)), -16),
                          (int(np.argmin(relative)), 9)):
        profile.annotate(f'{relative[index]:.1f}',
                         (values[index], relative[index]),
                         textcoords='offset points', xytext=(0, offset),
                         ha='center', fontsize=9, color=_INK)
    profile.margins(y=0.15)

    # ---- Panel 2: the state spectrum, if tracking data was supplied ------
    if has_states:
        spectrum = axes[1]
        _style_axis(spectrum)
        spectrum.set_facecolor('#fcfcfb')
        state_energies = np.asarray(data['state_energies'], dtype=float)
        reference = state_energies.min(axis=1)
        relative_states = (state_energies -
                           reference[:, None]) * HARTREE_IN_EV

        # State index is an ordered quantity, so the untracked roots use one
        # ordinal blue ramp instead of eight competing categorical hues.
        ramp = ['#86b6ef', '#5598e7', '#2a78d6', '#1c5cab', '#104281']
        for root in range(relative_states.shape[1]):
            spectrum.plot(values, relative_states[:, root],
                          color=ramp[min(root, len(ramp) - 1)],
                          linewidth=1.2, alpha=0.55, zorder=2)

        gaps = np.asarray(data['gap_to_s0_ev'], dtype=float)
        spectrum.plot(values, gaps, color=_SERIES[1], linewidth=2.0,
                      marker='o', markersize=5, zorder=4,
                      label=f'{state_label} above S0')
        spectrum.plot(values, np.zeros_like(values), color=_SERIES[2],
                      linewidth=2.0, zorder=3, label='S0')

        collisions = np.asarray(data.get('ground_state_collision',
                                         [False] * values.size), dtype=bool)
        if collisions.any():
            spectrum.scatter(values[collisions], gaps[collisions],
                             s=110, facecolors='none', edgecolors=_WARN,
                             linewidths=2.0, zorder=5,
                             label='S1/S0 seam reached')

        spectrum.set_ylabel('energy above S0 (eV)')
        spectrum.set_title('Adiabatic spectrum along the scan', loc='left',
                           fontsize=11, fontweight='bold')
        spectrum.legend(frameon=False, fontsize=9, labelcolor=_INK)

    # ---- Panel 3: geometry sanity ----------------------------------------
    geometry = axes[-1]
    _style_axis(geometry)
    geometry.set_facecolor('#fcfcfb')

    # Filter by the coordinate type, not by the label text: only bend angles
    # answer "is the constrained torsion still defined", and mixing a wrapped
    # dihedral onto the same axis makes the panel unreadable.
    probes = data.get('probes', {})
    kinds = data.get('probe_kinds', {})
    angle_probes = {label: series for label, series in probes.items()
                    if kinds.get(label, 'angle') == 'angle'}
    # Hues are assigned in fixed order and never cycled into new colors.  Past
    # the third probe the hue order repeats with a different dash pattern, so
    # identity is carried by hue *and* line style rather than by hue alone.
    dashes = ((), (5, 2), (1, 1.6), (6, 2, 1, 2))
    for index, (label, series) in enumerate(sorted(angle_probes.items())):
        line, = geometry.plot(
            values, series, color=_SERIES[index % len(_SERIES)],
            linewidth=2.0, marker='o', markersize=4, zorder=3, label=label)
        pattern = dashes[min(index // len(_SERIES), len(dashes) - 1)]
        if pattern:
            line.set_dashes(pattern)
    if angle_probes:
        geometry.axhline(165.0, color=_WARN, linewidth=1.5, linestyle='--',
                         zorder=2)
        # Sits just under its own line, right-aligned, so it can collide with
        # neither the panel title above nor the data below.
        geometry.annotate('torsion undefined above this line', (0.99, 164.0),
                          xycoords=('axes fraction', 'data'), fontsize=8,
                          color=_WARN, va='top', ha='right')
        # Upper left: the band under the threshold line is empty there, while
        # the data themselves sit low in the panel.
        geometry.legend(frameon=False, fontsize=9, labelcolor=_INK,
                        ncol=2, loc='upper left')
        geometry.margins(y=0.12)

    geometry.set_ylabel('angle (degrees)')
    geometry.set_xlabel(axis_label)
    geometry.set_title('Geometry checks: is the constraint still meaningful?',
                       loc='left', fontsize=11, fontweight='bold')

    figure.tight_layout()

    if filename is not None:
        figure.savefig(filename, dpi=160, facecolor=figure.get_facecolor())

    if show:
        plt.show()

    return figure


def write_scan_trajectory(data, filename):
    """Writes the converged scan geometries as a multi-frame XYZ trajectory.

    Opening this in any viewer is the fastest way to confirm that the relaxed
    structures are chemically sensible along the whole scan.

    :param data:
        The dictionary from :func:`collect_scan_data`.
    :param filename:
        Destination path.
    """

    lines = []
    for value, energy, xyz_string in zip(data['scan_values'],
                                         data['energies'],
                                         data['geometries']):
        labels, xyz = _xyz_to_array(xyz_string)
        lines.append(str(len(labels)))
        lines.append(f'scan coordinate {value:.3f}  energy {energy:.10f} Eh')
        for label, position in zip(labels, xyz):
            lines.append('{:<3s}{:20.12f}{:20.12f}{:20.12f}'.format(
                label, *position))

    Path(filename).write_text('\n'.join(lines) + '\n', encoding='utf-8')
