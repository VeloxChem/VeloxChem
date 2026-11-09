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
"""Rigid He/Ne QM scans and metal-only LJ fitting of a real metal complex.

The first section holds the inputs and saved data: ``MetalComplexLJStudy``
validates one study (complex, target atom, probe, scan, QM and fit settings),
``LJScanRecord`` is one rigid probe position with its QM energies, and
``QMScanDataset`` is the portable JSON file of all records, independent of fit
weights.

``MetalComplexLJDriver`` then runs the scans and the fit. QM observations are
immutable references to a fixed geometry/protocol. The force-field generator
supplies the baseline; VeloxChem's MMDriver evaluates copied host and
disconnected host+probe models without rebuilding topology or writing files.
Public distances/sigma are angstrom, epsilon and MM fitting energies are
kJ/mol, QM energies are hartree.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict, dataclass, field
import json
from pathlib import Path
from typing import Any
import tempfile

import sys

import numpy as np
from mpi4py import MPI

from .._version import __version__
from ..errorhandler import assert_msg_critical
from ..inputparser import parse_input, print_keywords, print_attributes
from ..mmdriver import MMDriver
from ..mmforcefieldgenerator import MMForceFieldGenerator
from ..molecule import Molecule
from ..outputstream import OutputStream
from ..veloxchemlib import (bohr_in_angstrom, chemical_element_identifier,
                            hartree_in_kjpermol, mpi_master)
from .calculator import DEFAULT_BASIS, SINGLE_POINT_SETTINGS, Calculator
from .evaluator import Evaluator
from .framework import check_force_field_atoms, copy_force_field
from .models import (CalculationResult, FitResult, FitTarget, ParameterSet,
                     ParameterSpec, h5_filename, print_settings, read_h5_group,
                     write_h5_groups,
                     ANGSTROM_TO_NM)

# ==========================================================================
# metal complex data
# ==========================================================================

PROBES = frozenset({"He", "Ne"})
FIT_SETTING_NAMES = frozenset({
    "sigma0_angstrom",
    "epsilon0_kjmol",
    "sigma_bounds_angstrom",
    "epsilon_bounds_kjmol",
    "energy_tolerance_kjmol",
    "point_weights",
    "strategy",
    "de_maxiter",
    "de_popsize",
    "seed",
    "max_nfev",
    "ftol",
    "xtol",
    "gtol",
})
DATASET_GROUP = "metal_complex_lj/qm_scan"
# A supplied baseline matches the study molecule when coordinates agree this closely.
COORDINATE_ATOL_ANGSTROM = 1e-6

# ----------------------------------------------------------------- helpers


def to_json(value) -> str:
    """Canonical JSON text (sorted keys, NumPy values converted, NaN refused)."""
    return json.dumps(value,
                      sort_keys=True,
                      allow_nan=False,
                      default=lambda x: x.item()
                      if isinstance(x, np.generic) else x.tolist())


def write_json(path: Path, payload) -> None:
    """Write ``payload`` as indented canonical JSON."""
    path.write_text(json.dumps(json.loads(to_json(payload)), indent=2) + "\n")


def make_molecule(labels, coordinates, charge, multiplicity, basis_labels=None):
    """A VeloxChem Molecule from labels and angstrom coordinates."""
    args = (list(labels), np.asarray(coordinates, dtype=float), "angstrom")
    mol = Molecule(
        *args, list(basis_labels)) if basis_labels is not None else Molecule(
            *args)
    mol.set_charge(charge)
    mol.set_multiplicity(multiplicity)
    return mol


def molecule_signature(molecule) -> dict[str, Any]:
    """Everything that identifies a molecule exactly (labels, basis labels, geometry, state)."""
    return {
        "labels": list(molecule.get_labels()),
        "basis_labels": list(molecule.get_atom_basis_labels()),
        "coordinates_angstrom": molecule.get_coordinates_in_angstrom().tolist(),
        "charge": float(molecule.get_charge()),
        "multiplicity": int(molecule.get_multiplicity())
    }


def same_molecule(left: dict,
                  right: dict,
                  atol_angstrom: float = COORDINATE_ATOL_ANGSTROM) -> bool:
    """Two ``molecule_signature`` dicts describe the same molecule.

    Labels, basis labels, charge and multiplicity must match exactly;
    coordinates within ``atol_angstrom``, because VeloxChem stores bohr and an
    angstrom round trip can change the last digit.
    """
    if any(left[key] != right[key]
           for key in ("labels", "basis_labels", "charge", "multiplicity")):
        return False
    a = np.asarray(left["coordinates_angstrom"], dtype=float)
    b = np.asarray(right["coordinates_angstrom"], dtype=float)
    return a.shape == b.shape and bool(
        np.allclose(a, b, rtol=0.0, atol=atol_angstrom))


# ------------------------------------------------------------------- study


@dataclass
class MetalComplexLJStudy:
    """Inputs for one real-atom LJ study; charges are generated once by default.

    ``probe_parameters`` is an optional explicit override with
    sigma_angstrom, epsilon_kjmol, and source. Otherwise the generator's probe
    defaults are used. ``baseline_force_field`` accepts a previously generated
    force field for the exact input molecule; it is copied, never modified.
    """

    molecule: Any
    target_atom_index: int
    probe: str = "He"
    scan_distances_angstrom: Any = field(
        default_factory=lambda: np.linspace(2.0, 5.0, 13).tolist())
    scan_directions: Any = "auto"
    clash_cutoff_scale: float = 0.85
    qm_settings: dict[str, Any] = field(default_factory=dict)
    fit_settings: dict[str, Any] = field(default_factory=dict)
    probe_parameters: dict[str, Any] | None = None
    baseline_force_field: Any = None

    def __post_init__(self):
        """Copy the molecule and settings, normalize the probe name, and validate every input.

        An invalid target index, probe, scan, clash scale, or unsupported
        QM/fit setting is a critical error.
        """
        sig = molecule_signature(self.molecule)
        self.molecule = make_molecule(sig["labels"],
                                      sig["coordinates_angstrom"],
                                      sig["charge"], sig["multiplicity"],
                                      sig["basis_labels"])
        self._check_target(sig)
        self.probe = self.probe.capitalize()
        assert_msg_critical(self.probe in PROBES, "probe must be He or Ne")
        self.scan_distances_angstrom = self._checked_distances()
        assert_msg_critical(
            np.isfinite(self.clash_cutoff_scale) and
            self.clash_cutoff_scale >= 0,
            "clash_cutoff_scale must be nonnegative and finite")
        if not isinstance(self.scan_directions, str):
            self.scan_directions = self._checked_directions()
        else:
            assert_msg_critical(
                self.scan_directions == "auto",
                "scan_directions must be 'auto' or a sequence of vectors")
        self.qm_settings = deepcopy(self.qm_settings)
        self.qm_settings.setdefault("basis", DEFAULT_BASIS)
        assert_msg_critical(
            not set(self.qm_settings) - SINGLE_POINT_SETTINGS,
            "unsupported QM settings; see Calculator.single_point")
        self.fit_settings = deepcopy(self.fit_settings)
        unknown = sorted(set(self.fit_settings) - FIT_SETTING_NAMES)
        assert_msg_critical(not unknown,
                            f"unsupported fit settings: {unknown}")
        self.probe_parameters = deepcopy(self.probe_parameters)

    def _check_target(self, sig):
        """Validate the target index and require real atoms with finite coordinates."""
        assert_msg_critical(
            not isinstance(self.target_atom_index, bool) and
            int(self.target_atom_index) == self.target_atom_index,
            "target_atom_index must be an integer")
        self.target_atom_index = int(self.target_atom_index)
        assert_msg_critical(
            0 <= self.target_atom_index < len(sig["labels"]),
            "target_atom_index is outside the input molecule")
        assert_msg_critical(
            all(chemical_element_identifier(label.upper()) > 0
                for label in sig["labels"]),
            "metal-complex fitting requires real atoms, without dummy sites")
        assert_msg_critical(
            bool(np.isfinite(sig["coordinates_angstrom"]).all()),
            "complex coordinates must be finite")

    def _checked_distances(self) -> list[float]:
        """Scan distances (angstrom) as a list; must be unique, positive and finite."""
        distances = np.asarray(self.scan_distances_angstrom, dtype=float)
        assert_msg_critical(
            distances.ndim == 1 and bool(distances.size) and
            bool(np.isfinite(distances).all()) and
            not np.any(distances <= 0),
            "scan distances must be a nonempty sequence of positive finite values")
        assert_msg_critical(
            len(set(distances)) == len(distances),
            "scan distances must be unique")
        return distances.tolist()

    def _checked_directions(self) -> list[list[float]]:
        """User scan directions as unit vectors; must be finite, nonzero, shape (n, 3)."""
        directions = np.asarray(self.scan_directions, dtype=float)
        assert_msg_critical(
            directions.ndim == 2 and directions.shape[1] == 3 and
            bool(len(directions)), "scan directions must have shape (n, 3)")
        norms = np.linalg.norm(directions, axis=1)
        assert_msg_critical(
            bool(np.isfinite(directions).all()) and not np.any(norms == 0),
            "scan directions must be finite nonzero vectors")
        return (directions / norms[:, None]).tolist()

    def signature(self):
        """Molecule, target, probe, scan and QM settings that fix the QM reference data."""
        return {
            "molecule": molecule_signature(self.molecule),
            "target_atom_index": self.target_atom_index,
            "probe": self.probe,
            "scan_distances_angstrom": self.scan_distances_angstrom,
            "scan_directions": self.scan_directions,
            "clash_cutoff_scale": self.clash_cutoff_scale,
            "qm_settings": self.qm_settings
        }


# ----------------------------------------------------------------- records


@dataclass
class LJScanRecord:
    """
    One rigid probe position: complex + probe coordinates and QM energies.

    Energies are in hartree; a rejected point keeps its ``reject_reason``.

    :param record_id:
        The record identifier.
    :param target_atom_index:
        The 0-based index of the target metal atom.
    :param probe:
        The probe element (He or Ne).
    :param distance_angstrom:
        The target-probe distance in angstrom.
    :param direction:
        The unit scan direction.
    :param coordinates_angstrom:
        The complex + probe coordinates in angstrom.
    :param accepted:
        Whether the point passed the clash check (and QM).
    :param reject_reason:
        Why the point was rejected.
    :param qm_converged:
        Whether the combined QM calculation converged.
    :param qm_host_energy_hartree:
        The isolated complex energy.
    :param qm_probe_energy_hartree:
        The isolated probe energy.
    :param qm_combined_energy_hartree:
        The complex + probe energy.
    :param qm_interaction_energy_hartree:
        The interaction energy (combined - complex - probe).
    """

    record_id: str
    target_atom_index: int
    probe: str
    distance_angstrom: float
    direction: list[float]
    coordinates_angstrom: list[list[float]]
    accepted: bool = True
    reject_reason: str | None = None
    qm_converged: bool | None = None
    qm_host_energy_hartree: float | None = None
    qm_probe_energy_hartree: float | None = None
    qm_combined_energy_hartree: float | None = None
    qm_interaction_energy_hartree: float | None = None

    @property
    def usable(self):
        """
        Checks whether the point can be fitted.

        :return:
            True when the point was accepted, QM converged, and all energies are
            finite.
        """

        energies = (self.qm_host_energy_hartree, self.qm_probe_energy_hartree,
                    self.qm_combined_energy_hartree,
                    self.qm_interaction_energy_hartree)
        return (self.accepted and self.qm_converged is True and all(
            value is not None and np.isfinite(value) for value in energies))


@dataclass
class QMScanDataset:
    """
    Portable embedded geometries and raw energies, independent of fit weights.

    :param metadata:
        The study signature, units, energy reference and version.
    :param records:
        The scan records.
    """

    metadata: dict[str, Any]
    records: list[LJScanRecord]

    def save(self, filename, group=DATASET_GROUP):
        """
        Writes the dataset to a group of a run file.

        :param filename:
            The run file (``.h5`` added when missing).
        :param group:
            The group of the dataset.

        :return:
            The path of the HDF5 file.
        """

        payload = {
            "version": 1,
            "metadata": self.metadata,
            "records": [asdict(record) for record in self.records]
        }
        return write_h5_groups(filename, {group: json.loads(to_json(payload))})

    @classmethod
    def load(cls, filename, group=DATASET_GROUP):
        """
        Reads a dataset from a run file, trusting its contents.

        :param filename:
            The run file (``.h5`` added when missing).
        :param group:
            The group of the dataset.

        :return:
            The QM scan dataset.
        """

        payload = read_h5_group(filename, group)
        if payload is None:
            raise FileNotFoundError(
                f"no QM scan dataset in {h5_filename(filename)}")
        return cls(payload.get("metadata", {}),
                   [LJScanRecord(**row) for row in payload["records"]])


# ==========================================================================
# metal complex lj
# ==========================================================================

FITTED_PARAMETERS = frozenset({"metal_lj_sigma", "metal_lj_epsilon"})
DEFAULT_SIGMA_BOUNDS_ANGSTROM = (0.5, 6.0)
DEFAULT_EPSILON_BOUNDS_KJMOL = (1e-6, 100.0)
DEFAULT_ENERGY_TOLERANCE_KJMOL = 1.0
OPTIMIZER_DEFAULTS = {
    "strategy": "least_squares",
    "de_maxiter": 40,
    "de_popsize": 8,
    "seed": 42,
    "max_nfev": 2000,
    "ftol": 1e-8,
    "xtol": 1e-8,
    "gtol": 1e-8
}
STRATEGIES = ("least_squares", "differential_evolution+least_squares")
AUTO_DIRECTION_SAMPLES = 128  # Fibonacci-sphere directions tried for the "auto" scan
ILL_CONDITIONED = 1e8  # scaled Jacobian condition number that triggers a warning
MM_ENERGY_NAME = "MM_INTERACTION_ENERGY"
FIT_RESULT_GROUP = "metal_complex_lj/result"


class MetalComplexLJDriver:
    """
    Fits the LJ parameters of a metal in a real-atom complex to QM probe scans.

    Two stages: QM reference energies of a He/Ne probe scanned around the
    target atom (``compute_qm_dataset``), then a direct fit of the MM
    interaction curves (``fit``). Settings are attributes (see
    ``_input_keywords`` and ``update_settings``); the molecule and target atom
    are given to ``setup`` or directly to ``compute_qm_dataset``/``fit``.

    ``calculator`` may implement ``single_point`` for a different QM backend.
    ``force_field_factory`` is an adapter/testing hook; the default calls the
    native generator with its unchanged host defaults.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes the metal-complex LJ driver with default settings.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        # MPI information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream

        # probe scan
        self.probe = 'He'
        self.scan_distances_angstrom = np.linspace(2.0, 5.0, 13).tolist()
        self.scan_directions = 'auto'
        self.clash_cutoff_scale = 0.85
        self.probe_parameters = None
        self.qm_settings = {}

        # fit (None: the built-in default)
        for name in sorted(FIT_SETTING_NAMES):
            setattr(self, name, None)

        # files
        self.filename = None

        # collaborators
        self.calculator = None
        self.evaluator = None
        self.force_field_factory = None

        self.study = None
        self.dataset = None
        self.baseline_force_field = None
        self.probe_atom = None
        self._host_mm = None
        self._combined_mm = None
        self._interaction_mm = None
        self.result = None
        self.last_calculation = None

        self._input_keywords = {
            'metal_complex_lj': {
                'probe': ('str', 'probe atom: He or Ne'),
                'scan_distances_angstrom':
                    ('raw', 'probe distances from the target atom'),
                'scan_directions': ('raw', "'auto' or scan direction vectors"),
                'clash_cutoff_scale':
                    ('float', 'scale of the clash distance check'),
                'probe_parameters':
                    ('raw', 'explicit probe sigma/epsilon/source'),
                'qm_settings': ('raw', 'QM single-point settings'),
                'sigma0_angstrom': ('float', 'initial sigma'),
                'epsilon0_kjmol': ('float', 'initial epsilon'),
                'sigma_bounds_angstrom': ('seq_fixed', 'sigma bounds'),
                'epsilon_bounds_kjmol': ('seq_fixed', 'epsilon bounds'),
                'energy_tolerance_kjmol':
                    ('float', 'tolerance of the interaction energies'),
                'point_weights': ('raw', 'weights of scan points'),
                'strategy': ('str', 'least_squares or '
                             'differential_evolution+least_squares'),
                'de_maxiter': ('int', 'differential evolution iterations'),
                'de_popsize': ('int', 'differential evolution population'),
                'seed': ('int', 'random seed of the global search'),
                'max_nfev': ('int', 'least-squares function evaluations'),
                'ftol': ('float', 'least-squares ftol'),
                'xtol': ('float', 'least-squares xtol'),
                'gtol': ('float', 'least-squares gtol'),
                'filename':
                    ('str', 'run file <filename>.h5: QM dataset and result'),
            },
        }

    def print_keywords(self):
        """
        Prints the input keywords of the metal-complex LJ driver.
        """

        print_keywords(self._input_keywords, self.ostream)

    def print_attributes(self):
        """
        Prints the attributes of the metal-complex LJ driver.
        """

        print_attributes(self._input_keywords, self.ostream)

    def update_settings(self, settings):
        """
        Updates the settings of the metal-complex LJ driver.

        :param settings:
            The dictionary of settings; unknown keys are an error.
        """

        keywords = {
            key: val[0] for key, val in self._input_keywords['metal_complex_lj'].items()
        }
        unknown = sorted(set(settings) - set(keywords))
        assert_msg_critical(
            not unknown, f'MetalComplexLJDriver: unknown settings {unknown}')
        parse_input(self, keywords, settings)

    def get_references(self):
        """
        Gets the references of the methods used by the fit.

        :return:
            The list of reference strings.
        """

        # least-squares and differential-evolution optimizers (SciPy)
        return ['P. Virtanen et al., Nat. Methods 2020, 17, 261']

    def print_references(self):
        """
        Prints the references of the methods used by the fit.
        """

        self.ostream.print_blank()
        for reference in self.get_references():
            self.ostream.print_reference('Reference: ' + reference)
        self.ostream.print_blank()
        self.ostream.flush()

    def _fit_settings(self):
        """
        Gets the fit settings that differ from the built-in defaults.

        :return:
            The fit settings that are set (not None).
        """

        return {
            name: deepcopy(getattr(self, name))
            for name in sorted(FIT_SETTING_NAMES)
            if getattr(self, name) is not None
        }

    def setup(self, molecule, target_atom_index, baseline_force_field=None):
        """
        Sets up the study of one molecule and target atom.

        The molecule and settings are copied and validated; later changes to
        the scan, QM or probe settings are rejected (see ``_check_study``).

        :param molecule:
            The metal complex (real atoms only).
        :param target_atom_index:
            The 0-based index of the metal atom whose LJ parameters are fitted.
        :param baseline_force_field:
            A previously generated force field of the exact molecule (copied,
            never modified); generated when omitted.
        """

        self.study = MetalComplexLJStudy(
            molecule=molecule,
            target_atom_index=target_atom_index,
            probe=self.probe,
            scan_distances_angstrom=self.scan_distances_angstrom,
            scan_directions=self.scan_directions,
            clash_cutoff_scale=self.clash_cutoff_scale,
            qm_settings=self.qm_settings,
            fit_settings=self._fit_settings(),
            probe_parameters=self.probe_parameters,
            baseline_force_field=baseline_force_field)
        if self.calculator is None:
            self.calculator = Calculator(comm=self.comm, ostream=self.ostream)
        if self.evaluator is None:
            self.evaluator = Evaluator(ostream=self.ostream)
        if self.force_field_factory is None:
            self.force_field_factory = (
                lambda: MMForceFieldGenerator(ostream=self.ostream))
        self._study_signature = deepcopy(self.study.signature())
        self._probe_override = deepcopy(self.study.probe_parameters)
        self._settings_snapshot = self._settings()
        self.dataset = None
        self.baseline_force_field = None
        self.probe_atom = None
        self._host_mm = None
        self._combined_mm = None
        self._interaction_mm = None
        self.result = None
        self.last_calculation = None

    def _settings(self):
        """
        Gets the current settings, to detect changes after ``setup``.

        :return:
            The JSON text of the settings.
        """

        # Fit settings may change between fits; they are read at fit time.
        return to_json({
            key: getattr(self, key)
            for key in self._input_keywords['metal_complex_lj']
            if key not in FIT_SETTING_NAMES
        })

    def _setup_if_given(self, molecule, target_atom_index,
                        baseline_force_field):
        """
        Sets up the study when a molecule is given; requires a study otherwise.
        """

        if molecule is not None:
            self.setup(molecule, target_atom_index, baseline_force_field)
        assert_msg_critical(
            self.study is not None,
            'MetalComplexLJDriver: give the molecule and target atom index')

    def _check_study(self):
        """Critical error if the study or scan/QM settings changed since ``setup``."""
        assert_msg_critical(
            to_json(self.study.signature()) == to_json(self._study_signature)
            and self._settings() == self._settings_snapshot,
            "molecule or scan/QM settings changed; create a new study driver")
        assert_msg_critical(
            self.study.probe_parameters == self._probe_override,
            "probe parameters changed; create a new study driver")

    def _probe_molecule(self):
        """A neutral singlet He/Ne atom at the origin."""
        return make_molecule([self.study.probe], [[0, 0, 0]], 0, 1)

    def _combined_molecule(self, record):
        """Host + probe molecule at one record's coordinates, with the study's basis labels."""
        sig = self._study_signature["molecule"]
        return make_molecule(
            sig["labels"] + [self.study.probe], record.coordinates_angstrom,
            sig["charge"], sig["multiplicity"],
            sig["basis_labels"] + [("", self.study.probe.upper())])

    # ================================================================ baseline
    def prepare_baseline(self):
        """
        Generates or copies the default host force field once.

        The charges and topology of the baseline are preserved; the probe atom
        and the MM engines are set up with it.

        :return:
            A copy of the baseline force field.
        """

        self._check_study()
        if self.baseline_force_field is not None:
            return copy_force_field(self.baseline_force_field)
        supplied = self.study.baseline_force_field
        ff = self._host_force_field(supplied)
        self._check_host_force_field(ff)
        probe_atom = self._probe_atom()
        self.baseline_force_field = copy_force_field(ff)
        self.probe_atom = probe_atom
        self._build_mm_engines(probe_atom)
        self.baseline_metadata = {
            "source": ("supplied generator baseline" if supplied is not None
                       else "MMForceFieldGenerator defaults"),
            "probe_parameters": deepcopy(probe_atom),
            "mixing_rule": "Lorentz–Berthelot",
            "backend": "VeloxChem MMDriver",
            "sigma_unit": "nm",
            "epsilon_unit": "kJ/mol"
        }
        return copy_force_field(self.baseline_force_field)

    def _host_force_field(self, supplied):
        """A copy of the supplied baseline, or the generator's defaults for the complex."""
        if supplied is None:
            ff = self.force_field_factory()
            # RESP/SCF scratch files of the generator go to a temporary directory.
            with tempfile.TemporaryDirectory() as directory:
                ff.molecule_name = str(Path(directory) / "baseline")
                ff.create_topology(self.study.molecule)
            ff.molecule_name = "baseline"
            return ff
        assert_msg_critical(
            same_molecule(molecule_signature(supplied.molecule),
                          self._study_signature["molecule"]),
            "supplied baseline does not match the exact input molecule")
        return copy_force_field(supplied)

    def _check_host_force_field(self, ff):
        """Validate the baseline: atom count/order, finite nonnegative LJ, LB mixing, total charge.

        Stops at the first problem found (atom count/order and total charge
        are critical errors; ``check_force_field_atoms`` raises ValueError).
        """
        natoms = self.study.molecule.number_of_atoms()
        assert_msg_critical(
            set(ff.atoms) == set(range(natoms)),
            "baseline atom ordering/count does not match the complex")
        check_force_field_atoms(ff)
        assert_msg_critical(
            bool(np.isclose(sum(atom["charge"] for atom in ff.atoms.values()),
                            self.study.molecule.get_charge(),
                            atol=1e-5)),
            "baseline partial charges do not sum to the molecular charge")

    def _probe_atom(self) -> dict:
        """The neutral probe's force-field atom (nm, kJ/mol): generator default or override."""
        override = self.study.probe_parameters
        if override is None:
            probe_ff = self.force_field_factory()
            try:
                with tempfile.TemporaryDirectory() as directory:
                    probe_ff.molecule_name = str(Path(directory) / "probe")
                    probe_ff.create_topology(self._probe_molecule(), resp=False)
            except Exception as exc:
                raise ValueError(
                    f"no generator defaults for {self.study.probe}; "
                    "supply probe_parameters") from exc
            probe_atom = deepcopy(probe_ff.atoms[0])
            assert_msg_critical(bool(np.isclose(probe_atom["charge"], 0)),
                                "He/Ne probe must have zero charge")
        else:
            assert_msg_critical(bool(override.get("source")),
                                "explicit probe_parameters require a source")
            probe_atom = {
                "type": self.study.probe,
                "name": self.study.probe,
                "mass": float(self._probe_molecule().get_masses()[0]),
                "charge": 0.0,
                "sigma": float(override["sigma_angstrom"]) * ANGSTROM_TO_NM,
                "epsilon": float(override["epsilon_kjmol"]),
                "comment": override["source"]
            }
        assert_msg_critical(
            all(
                np.isfinite(probe_atom[k]) and probe_atom[k] > 0
                for k in ("sigma", "epsilon")),
            "probe sigma/epsilon must be finite and positive; "
            "provide supported parameters")
        return probe_atom

    def _build_mm_engines(self, probe_atom):
        """Host, host+probe, and cross-terms-only MMDriver engines from the baseline."""
        natoms = self.study.molecule.number_of_atoms()
        self._host_mm = MMDriver(ostream=self.ostream)
        self._host_mm.load_force_field(self.baseline_force_field)
        # MMDriver.load_force_field shallow-copies atoms; isolate the baseline.
        self._host_mm.force_field = deepcopy(self._host_mm.force_field)
        self._combined_mm = MMDriver(ostream=self.ostream)
        self._combined_mm.force_field = deepcopy(self._host_mm.force_field)
        self._combined_mm.force_field["atoms"][natoms] = deepcopy(probe_atom)
        # Internal host terms cancel exactly for the rigid disconnected probe.
        # Evaluate only surviving cross terms with the existing MMDriver kernel.
        # This avoids cancellation noise and singular internal torsion gradients.
        self._interaction_mm = MMDriver(ostream=self.ostream)
        self._interaction_mm.force_field = deepcopy(
            self._combined_mm.force_field)
        self._interaction_mm.force_field["pairs"] = {}
        self._interaction_mm.exclusions = [
            (i, j) for i in range(natoms) for j in range(i + 1, natoms)
        ]

    # =================================================================== scans
    def _clash(self, position):
        """Reason string if a probe at ``position`` (angstrom) clashes with the host, else None.

        A clash is a host-probe distance below ``clash_cutoff_scale`` times the sum
        of covalent radii.
        """
        mol = self.study.molecule
        distances = np.linalg.norm(mol.get_coordinates_in_angstrom() - position,
                                   axis=1)
        if np.any(distances < 1e-8):
            return "probe overlaps a host atom"
        radii = mol.covalent_radii_to_numpy() * bohr_in_angstrom()
        probe_radius = self._probe_molecule().covalent_radii_to_numpy(
        )[0] * bohr_in_angstrom()
        bad = np.flatnonzero(distances < self.study.clash_cutoff_scale *
                             (radii + probe_radius))
        return f"clash with host atom {int(bad[0])}" if len(bad) else None

    def _directions(self):
        """User directions, or the single most open clash-free Fibonacci-sphere direction."""
        if not isinstance(self.study.scan_directions, str):
            return self.study.scan_directions
        directions = _fibonacci_sphere(AUTO_DIRECTION_SAMPLES)
        coords = self.study.molecule.get_coordinates_in_angstrom()
        target = coords[self.study.target_atom_index]
        positions = target + min(
            self.study.scan_distances_angstrom) * directions
        neighbors = np.delete(coords, self.study.target_atom_index, axis=0)
        if len(neighbors):
            gaps = np.linalg.norm(positions[:, None] - neighbors, axis=2)
            clearance = np.min(gaps, axis=1)
        else:
            clearance = np.zeros(AUTO_DIRECTION_SAMPLES)
        valid = np.array(
            [self._clash(position) is None for position in positions])
        indices = np.flatnonzero(valid) if valid.any() else np.arange(
            AUTO_DIRECTION_SAMPLES)
        return [directions[indices[np.argmax(clearance[indices])]].tolist()]

    def generate_scans(self):
        """
        Generates the probe positions along each scan direction.

        Clashing points are kept but rejected. An existing dataset is reused.

        :return:
            The scan records.
        """

        self._check_study()
        if self.dataset is not None:
            return self.dataset.records
        coords = self.study.molecule.get_coordinates_in_angstrom()
        target = coords[self.study.target_atom_index]
        records = []
        for d, direction in enumerate(self._directions()):
            for i, distance in enumerate(self.study.scan_distances_angstrom):
                position = target + distance * np.array(direction)
                reason = self._clash(position)
                records.append(
                    LJScanRecord(f"direction-{d:03d}-point-{i:04d}",
                                 self.study.target_atom_index, self.study.probe,
                                 distance, direction,
                                 np.vstack((coords, position)).tolist(), reason
                                 is None, reason))
        self.dataset = QMScanDataset(
            {
                "study": deepcopy(self._study_signature),
                "units": {
                    "coordinates": "angstrom",
                    "qm_energy": "hartree"
                },
                "energy_reference": "combined - complex - probe; no counterpoise",
                "veloxchem_version": __version__
            }, records)
        return records

    # ============================================================ QM reference
    def compute_qm_dataset(self,
                           molecule=None,
                           target_atom_index=None,
                           baseline_force_field=None):
        """
        Calculates missing references and accepted scan points.

        Failure reasons of rejected points are kept.

        :param molecule:
            The metal complex; sets up a new study (see ``setup``).
        :param target_atom_index:
            The 0-based index of the target metal atom.
        :param baseline_force_field:
            An optional force field of the exact molecule.

        :return:
            The QM scan dataset.
        """

        self._setup_if_given(molecule, target_atom_index, baseline_force_field)
        self.generate_scans()
        self._print_qm_setup()
        references = self._compute_references()
        for record in self.dataset.records:
            if record.accepted and not record.usable:
                self._compute_record(record, references)
        self.save_dataset()
        return self.dataset

    def _compute_references(self) -> dict[str, float]:
        """Isolated complex and probe energies (hartree), computed once and saved."""
        settings = self.study.qm_settings
        references = self.dataset.metadata.setdefault("references_hartree", {})
        try:
            for name, molecule in (("host", self.study.molecule),
                                   ("probe", self._probe_molecule())):
                if name not in references:
                    energy = float(
                        self.calculator.single_point(molecule,
                                                     settings=settings))
                    if not np.isfinite(energy):
                        raise RuntimeError(f"nonfinite {name} reference energy")
                    references[name] = energy
        except Exception as exc:
            self.dataset.metadata["reference_failure"] = str(exc)
            self.save_dataset()
            raise RuntimeError(
                "QM monomer reference calculation failed; dataset saved"
            ) from exc
        self.dataset.metadata.pop("reference_failure", None)
        return references

    def _compute_record(self, record, references):
        """Combined energy and interaction energy for one point; a failure rejects it."""
        record.qm_host_energy_hartree = references["host"]
        record.qm_probe_energy_hartree = references["probe"]
        try:
            energy = float(
                self.calculator.single_point(self._combined_molecule(record),
                                             settings=self.study.qm_settings))
            if not np.isfinite(energy):
                raise RuntimeError("nonfinite combined energy")
            record.qm_combined_energy_hartree = energy
            record.qm_interaction_energy_hartree = energy - references[
                "host"] - references["probe"]
            if not np.isfinite(record.qm_interaction_energy_hartree):
                record.qm_interaction_energy_hartree = None
                raise RuntimeError("nonfinite interaction energy")
            record.qm_converged = True
        except Exception as exc:
            record.qm_converged = False
            record.accepted = False
            record.reject_reason = f"qm_failed: {exc}"

    def _print_qm_setup(self):
        """Print the scan, probe and QM settings before the QM calculations."""
        study, records = self.study, self.dataset.records
        distances = study.scan_distances_angstrom
        charge = study.molecule.get_charge()
        multiplicity = study.molecule.get_multiplicity()
        span = f"{min(distances):g}-{max(distances):g}"
        if isinstance(study.scan_directions, str):
            direction_source = "auto, most open"
        else:
            direction_source = "user supplied"
        print_settings(self.ostream, "Metal-Complex QM Scan Setup", [
            ("Target Atom",
             f"{study.molecule.get_labels()[study.target_atom_index]} "
             f"(index {study.target_atom_index})"),
            ("Complex Charge / Multiplicity", f"{charge:g} / {multiplicity}"),
            ("Probe", study.probe),
            ("Scan Distances", f"{len(distances)} points, {span} angstrom"),
            ("Scan Directions",
             f"{len(self._directions())} ({direction_source})"),
            ("Clash Cutoff Scale", study.clash_cutoff_scale),
            ("Accepted / Clashing Points",
             f"{sum(r.accepted for r in records)} / "
             f"{sum(not r.accepted for r in records)}"),
            ("Points Already Computed", sum(r.usable for r in records)),
            *((f"QM {name.replace('_', ' ').title()}", value)
              for name, value in sorted(study.qm_settings.items())),
            ("Energy Reference", self.dataset.metadata["energy_reference"]),
            ("Run File", h5_filename(self.filename) if self.filename else None),
        ])

    def save_dataset(self, filename=None):
        """
        Writes the scan dataset to a run file.

        Baseline metadata is attached when the baseline has been prepared.
        Raises RuntimeError if no scans have been generated or loaded.

        :param filename:
            The run file (default: the ``filename`` setting; nothing is
            written without either).

        :return:
            The path of the HDF5 file, or None.
        """
        self._check_study()
        if self.dataset is None:
            raise RuntimeError("generate or load scans before saving a dataset")
        if self.baseline_force_field is not None:
            self.dataset.metadata["baseline"] = deepcopy(self.baseline_metadata)
        filename = filename or self.filename
        if filename is None or self.rank != mpi_master():
            return h5_filename(filename) if filename is not None else None
        return self.dataset.save(filename)

    def load_dataset(self, filename):
        """
        Uses the dataset of a run file as-is; the caller ensures correctness.

        :param filename:
            The run file.

        :return:
            The QM scan dataset.
        """

        self._check_study()
        dataset = QMScanDataset.load(filename)
        self.dataset = dataset
        self.result = None
        return dataset

    def _usable_records(self):
        """Accepted, converged records with finite QM energies; raises RuntimeError if none."""
        self._check_study()
        rows = [] if self.dataset is None else [
            row for row in self.dataset.records if row.usable
        ]
        if not rows:
            raise RuntimeError(
                "no accepted finite converged QM scan points available")
        return rows

    # ============================================================ MM and score
    def _parameter_specs(self):
        """Sigma (angstrom) and epsilon (kJ/mol) specs; initial values default to the baseline atom.

        A lower search bound that is not positive is a critical error.
        """
        atom = self.baseline_force_field.atoms[self.study.target_atom_index]
        settings = self._fit_settings()
        specs = [
            ParameterSpec(
                "metal_lj_sigma",
                settings.get("sigma0_angstrom", atom["sigma"] / ANGSTROM_TO_NM),
                tuple(
                    settings.get("sigma_bounds_angstrom",
                                 DEFAULT_SIGMA_BOUNDS_ANGSTROM)), "angstrom"),
            ParameterSpec(
                "metal_lj_epsilon",
                settings.get("epsilon0_kjmol", atom["epsilon"]),
                tuple(
                    settings.get("epsilon_bounds_kjmol",
                                 DEFAULT_EPSILON_BOUNDS_KJMOL)), "kJ/mol")
        ]
        assert_msg_critical(
            all(spec.bounds[0] > 0 for spec in specs),
            "sigma and epsilon search bounds must be positive")
        return specs

    def evaluate_mm(self, parameters=None):
        """
        Calculates the MM host-probe interaction energies of every usable record.

        :param parameters:
            The metal_lj_sigma (angstrom) and metal_lj_epsilon (kJ/mol) values
            (default: the initial values of the fit).

        :return:
            The calculation result with the energies in kJ/mol.
        """

        rows = self._usable_records()
        self.prepare_baseline()
        if parameters is None:
            parameters = ParameterSet(
                {spec.name: spec.initial for spec in self._parameter_specs()})
        elif not isinstance(parameters, ParameterSet):
            parameters = ParameterSet(parameters)
        if set(parameters) != FITTED_PARAMETERS:
            raise ValueError(
                "only metal_lj_sigma and metal_lj_epsilon may be fitted")
        sigma = parameters["metal_lj_sigma"]
        epsilon = parameters["metal_lj_epsilon"]
        if not np.isfinite([sigma, epsilon]).all() or min(sigma, epsilon) <= 0:
            raise ValueError(
                "metal sigma/epsilon must be positive finite values")
        for engine in (self._host_mm, self._combined_mm, self._interaction_mm):
            atom = engine.force_field["atoms"][self.study.target_atom_index]
            atom["sigma"], atom["epsilon"] = sigma * ANGSTROM_TO_NM, epsilon
        # E(combined)-E(host)-E(single probe) leaves only host–probe cross
        # nonbonded terms. Host bonded and 1–4 terms remain in the saved models.
        values = [self._cross_energy(record) for record in rows]
        if not np.isfinite(values).all():
            raise RuntimeError(
                "MM interaction calculation produced nonfinite energies")
        self.last_calculation = CalculationResult(
            values={MM_ENERGY_NAME: values},
            metadata={
                MM_ENERGY_NAME: {
                    "unit": "kJ/mol",
                    "quality_status": "usable",
                    "backend": "VeloxChem MMDriver.compute_nonbonded",
                    "energy_reference": "host–probe cross terms after exact internal cancellation"
                }
            },
            raw_data={
                "record_ids": [row.record_id for row in rows],
                "parameters": parameters.to_dict()
            })
        return self.last_calculation

    def _cross_energy(self, record) -> float:
        """Host-probe Coulomb + van der Waals cross energy (kJ/mol) for one record's geometry."""
        engine = self._interaction_mm
        coords_nm = np.array(record.coordinates_angstrom) * ANGSTROM_TO_NM
        engine.gradient = np.zeros_like(coords_nm)
        engine.coulomb_potential = engine.vdw_potential = 0.0
        engine.compute_nonbonded(coords_nm)
        return engine.coulomb_potential + engine.vdw_potential

    def _targets(self):
        """One FitTarget per usable record (QM interaction, kJ/mol); sets point weights."""
        settings = self._fit_settings()
        tolerance = float(
            settings.get("energy_tolerance_kjmol",
                         DEFAULT_ENERGY_TOLERANCE_KJMOL))
        weights = settings.get("point_weights", {})
        assert_msg_critical(
            isinstance(weights, dict),
            "point_weights must map record IDs to nonnegative weights")
        ids = {row.record_id for row in self.dataset.records}
        assert_msg_critical(not set(weights) - ids,
                            "point_weights contain unknown record IDs")
        assert_msg_critical(
            all(
                np.isfinite(weight) and weight >= 0
                for weight in weights.values()),
            "point weights must be finite and nonnegative")
        rows = self._usable_records()
        self._point_weights = np.array(
            [weights.get(row.record_id, 1.0) for row in rows])
        assert_msg_critical(
            bool(np.any(self._point_weights > 0)),
            "at least one usable point must have positive weight")
        return [
            FitTarget(row.record_id,
                      row.qm_interaction_energy_hartree * hartree_in_kjpermol(),
                      tolerance,
                      float(weight),
                      unit="kJ/mol")
            for row, weight in zip(rows, self._point_weights)
        ]

    def evaluate(self, calculation=None):
        """
        Scores MM predictions against the QM interaction energies with the Evaluator.

        Raises ValueError if the predictions do not match the usable records.

        :param calculation:
            The MM calculation (default: the last MM calculation, else a fresh one
            at the initial parameters).

        :return:
            The evaluation result.
        """

        calculation = calculation if calculation is not None else self.last_calculation
        if calculation is None:
            calculation = self.evaluate_mm()
        targets = self._targets()
        ids = calculation.raw_data["record_ids"]
        if ids != [target.name for target in targets]:
            raise ValueError(
                "MM predictions do not match the reference record ordering")
        if len(calculation.values[MM_ENERGY_NAME]) != len(ids):
            raise ValueError("MM predictions have the wrong number of points")
        values = dict(zip(ids, calculation.values[MM_ENERGY_NAME]))
        mapped = CalculationResult(values=values,
                                   metadata=calculation.metadata,
                                   raw_data=calculation.raw_data)
        return self.evaluator.evaluate(mapped, targets)

    # ===================================================================== fit
    def fit(self, molecule=None, target_atom_index=None,
            baseline_force_field=None):
        """
        Fits sigma and epsilon of the target atom to the QM references.

        Optimizer termination and target satisfaction are reported separately.

        :param molecule:
            The metal complex; sets up a new study (see ``setup``).
        :param target_atom_index:
            The 0-based index of the target metal atom.
        :param baseline_force_field:
            An optional force field of the exact molecule.

        :return:
            The fit result.
        """

        self._setup_if_given(molecule, target_atom_index, baseline_force_field)
        if self.dataset is None:
            self.compute_qm_dataset()
        rows = self._usable_records()
        self.prepare_baseline()
        specs = self._parameter_specs()
        targets = self._targets()
        settings = {**OPTIMIZER_DEFAULTS, **self._fit_settings()}
        strategy = settings["strategy"]
        self.print_references()
        assert_msg_critical(strategy in STRATEGIES, "unknown fitting strategy")
        reference = np.array([target.reference for target in targets])
        tolerance = targets[0].tolerance
        self._print_fit_setup(strategy, specs, rows, tolerance, settings)
        optimized, global_result, history = self._optimize(
            strategy, specs, reference, tolerance, settings)
        parameters = ParameterSet(
            dict(zip((spec.name for spec in specs), optimized.x)))
        calculation = self.evaluate_mm(parameters)
        evaluation = self.evaluate(calculation)
        self.result = self._fit_result(parameters, calculation, evaluation,
                                       optimized, global_result, history,
                                       strategy, specs, rows, reference,
                                       tolerance)
        self.save_dataset()
        self._save_fit_outputs()
        return self.result

    def _optimize(self, strategy, specs, reference, tolerance, settings):
        """Optional differential evolution, then bounded least squares on weighted residuals."""
        from scipy.optimize import differential_evolution, least_squares

        # Residuals sqrt(w) * (MM - QM) / tolerance, so the sum of squares is the
        # weighted, tolerance-normalized objective.
        scale = np.sqrt(self._point_weights) / tolerance
        history = []

        def residual(values):
            """Weighted, tolerance-scaled MM - QM residuals; appends each trial to ``history``."""
            parameters = ParameterSet(
                dict(zip((spec.name for spec in specs), values)))
            predicted = np.array(
                self.evaluate_mm(parameters).values[MM_ENERGY_NAME])
            errors = scale * (predicted - reference)
            history.append({
                "iteration": len(history) + 1,
                "parameters": parameters.to_dict(),
                "objective": float(errors @ errors)
            })
            return errors

        bounds = np.array([spec.bounds for spec in specs])
        initial = np.array([spec.initial for spec in specs])
        global_result = None
        if strategy == "differential_evolution+least_squares":
            global_result = differential_evolution(
                lambda values: float(np.sum(residual(values)**2)),
                bounds,
                maxiter=settings["de_maxiter"],
                popsize=settings["de_popsize"],
                seed=settings["seed"],
                polish=False,
                workers=1)
            initial = global_result.x
        optimized = least_squares(residual,
                                  initial,
                                  bounds=(bounds[:, 0], bounds[:, 1]),
                                  max_nfev=settings["max_nfev"],
                                  ftol=settings["ftol"],
                                  xtol=settings["xtol"],
                                  gtol=settings["gtol"],
                                  x_scale="jac")
        return optimized, global_result, history

    def _fit_result(self, parameters, calculation, evaluation, optimized,
                    global_result, history, strategy, specs, rows, reference,
                    tolerance) -> FitResult:
        """Convergence, conditioning and bound warnings, and the full diagnostics."""
        bounds = np.array([spec.bounds for spec in specs])
        mm_values = calculation.values[MM_ENERGY_NAME]
        differences = np.array(mm_values) - reference
        meets_tolerance = bool(
            np.all(np.abs(differences[self._point_weights > 0]) <= tolerance))
        converged = bool(optimized.success and meets_tolerance and
                         not evaluation.failed)
        singular = np.linalg.svd(optimized.jac[:, :] * optimized.x[None, :],
                                 compute_uv=False)
        condition = (float(singular[0] / singular[-1])
                     if len(singular) == 2 and singular[-1] > 0 else None)
        warnings = list(evaluation.warnings)
        if condition is None or condition > ILL_CONDITIONED:
            warnings.append(
                "scan Jacobian poorly constrains sigma/epsilon; inspect parameter sensitivity"
            )
        if np.any(np.isclose(optimized.x, bounds[:, 0])) or np.any(
                np.isclose(optimized.x, bounds[:, 1])):
            warnings.append("fitted parameter lies near a search bound")
        baseline = self.baseline_force_field.atoms[self.study.target_atom_index]
        if converged:
            stop_reason = "target_tolerance"
        elif optimized.success:
            stop_reason = "optimizer_terminated"
        else:
            stop_reason = "optimizer_limit_or_failure"
        diagnostics = {
            "stop_reason": stop_reason,
            "optimizer_success": bool(optimized.success),
            "optimizer_message": str(optimized.message),
            "de_success": None if global_result is None else bool(
                global_result.success),
            "de_message": None if global_result is None else str(
                global_result.message),
            "strategy": strategy,
            "n_evaluations": len(history),
            "n_usable": len(rows),
            "n_rejected": len(self.dataset.records) - len(rows),
            "bounds": bounds.tolist(),
            "parameter_units": {
                spec.name: spec.unit for spec in specs
            },
            "baseline_sigma_angstrom": baseline["sigma"] / ANGSTROM_TO_NM,
            "baseline_epsilon_kjmol": baseline["epsilon"],
            "baseline": self.baseline_metadata,
            "point_weights": self._point_weights.tolist(),
            "energy_tolerance_kjmol": tolerance,
            "jacobian_singular_values_scaled": singular.tolist(),
            "jacobian_condition_scaled": condition,
            "record_ids": [row.record_id for row in rows],
            "scan_distances_angstrom": [row.distance_angstrom for row in rows],
            "scan_directions": [row.direction for row in rows],
            "probe": self.study.probe,
            "qm_interaction_kjmol": reference.tolist(),
            "mm_interaction_kjmol": mm_values,
            "fit_settings": deepcopy(self._fit_settings()),
            "scientific_validation": "not performed",
        }
        return FitResult(
            parameters,
            metrics={
                "rmse_kjmol": float(np.sqrt(np.mean(differences**2))),
                "max_abs_error_kjmol": float(np.max(np.abs(differences)))
            },
            objective=evaluation.objective,
            status="failed" if evaluation.failed else
            ("success" if converged else "warning"),
            converged=converged,
            iterations=int(optimized.nfev),
            history=history,
            warnings=warnings,
            diagnostics=diagnostics)

    def _print_fit_setup(self, strategy, specs, rows, tolerance, settings):
        """Print the optimizer, bounds, probe and tolerance settings before fitting."""
        probe = self.probe_atom
        uses_de = strategy != "least_squares"
        print_settings(self.ostream, "Metal-Complex LJ Fitting Setup", [
            ("Strategy", strategy),
            *((f"Parameter {spec.name}",
               f"[{spec.bounds[0]:g}, {spec.bounds[1]:g}] {spec.unit}, "
               f"initial {spec.initial:.6g}") for spec in specs),
            ("Probe sigma / epsilon",
             f"{probe['sigma'] / ANGSTROM_TO_NM:.4g} angstrom / "
             f"{probe['epsilon']:.4g} kJ/mol"),
            ("Mixing Rule", "Lorentz-Berthelot"),
            ("Usable QM Points", f"{len(rows)} of {len(self.dataset.records)}"),
            ("Energy Tolerance", f"{tolerance:g} kJ/mol"),
            ("Weighted Points", int(np.sum(self._point_weights != 1.0)) or
             None),
            ("DE Max. Iterations / Pop. Size",
             f"{settings['de_maxiter']} / {settings['de_popsize']}"
             if uses_de else None),
            ("Least-Squares Max. Evaluations", settings["max_nfev"]),
            ("ftol / xtol / gtol",
             f"{settings['ftol']:g} / {settings['xtol']:g} / "
             f"{settings['gtol']:g}"),
            ("Random Seed", settings["seed"] if uses_de else None),
        ])

    # ================================================================= outputs
    def fitted_force_field(self, parameters=None):
        """
        Creates a copy of the baseline force field with only the target atom's LJ
        changed.

        :param parameters:
            The metal_lj_sigma (angstrom) and metal_lj_epsilon (kJ/mol) values
            (default: the fitted parameters).

        :return:
            The copied force-field generator.
        """

        self.prepare_baseline()
        if parameters is None:
            if self.result is None:
                raise RuntimeError("fit first or supply parameters")
            parameters = self.result.parameters
        if not isinstance(parameters, ParameterSet):
            parameters = ParameterSet(parameters)
        assert_msg_critical(
            set(parameters) == FITTED_PARAMETERS and all(
                np.isfinite(v) and v > 0 for v in parameters.values.values()),
            "supply positive finite metal_lj_sigma/metal_lj_epsilon only")
        ff = copy_force_field(self.baseline_force_field)
        atom = ff.atoms[self.study.target_atom_index]
        atom["sigma"] = parameters["metal_lj_sigma"] * ANGSTROM_TO_NM
        atom["epsilon"] = parameters["metal_lj_epsilon"]
        return ff

    def export_fitted_model(self, path=None):
        """
        Writes a per-atom force-field XML using the generator's existing writer.

        The fitted generator (see ``fitted_force_field``) can also be loaded
        directly into MMDriver. The writer assigns per-atom classes, preventing
        other same-type atoms from inheriting the target metal's fitted LJ
        parameters. A JSON file with the study, parameters and diagnostics is
        written next to it.

        :param path:
            The output XML file.

        :return:
            The path of the XML file.
        """

        assert_msg_critical(path is not None,
                            "give an output file for the fitted model")
        ff = self.fitted_force_field()
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        if self.rank != mpi_master():
            return path
        ff.generate_residue_xml(path, mol_name="COMPLEX")
        write_json(
            path.with_suffix(".fit.json"), {
                "study": self._study_signature,
                "parameters": self.result.parameters.to_dict(),
                "diagnostics": self.result.diagnostics
            })
        return path

    def _save_fit_outputs(self):
        """
        Writes the fit result (with the QM/MM curves in its diagnostics) to
        the run file; no-op without a ``filename`` setting.
        """

        if self.filename is None or self.rank != mpi_master():
            return
        write_h5_groups(self.filename,
                        {FIT_RESULT_GROUP: json.loads(to_json(asdict(self.result)))})


def _fibonacci_sphere(n: int) -> np.ndarray:
    """``n`` reproducible, nearly uniform unit vectors."""
    i = np.arange(n)
    z = 1 - (2 * i + 1) / n
    phi = i * np.pi * (3 - np.sqrt(5))
    return np.column_stack(
        (np.sqrt(1 - z * z) * np.cos(phi), np.sqrt(1 - z * z) * np.sin(phi), z))
