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
import numpy as np
import sys

from .veloxchemlib import (mpi_master, hartree_in_ev, avogadro_constant,
                           bohr_in_angstrom)
from .outputstream import OutputStream
from .molecule import Molecule
from .molecularbasis import MolecularBasis
from .scfrestdriver import ScfRestrictedDriver
from .scfrestopendriver import ScfRestrictedOpenDriver
from .scfunrestdriver import ScfUnrestrictedDriver
from .lreigensolver import LinearResponseEigenSolver
from .cppsolver import ComplexResponseSolver
from .spectrumaverager import SpectrumAverager
from .ensembleparser import EnsembleParser
from .embeddingpotfilewriter import EmbeddingPotfileWriter
from .errorhandler import assert_msg_critical
from .inputparser import parse_bool, parse_input, parse_seq_range


class EnsembleDriver:
    """
    Implements ensemble calculations over a list of snapshots produced by
    :class EnsembleParser.

    Each snapshot contains teh QM subsysem and (optionally) environment atoms.
    The driver can run, for each snapshot:

    - A SCF calculation (restricted or unrestricted)
    - Polarizable embedding (fika) via per-snapshot potfiles.
    - Non-polarizable embedding (NPE) through the fika potfile.
    - Optional linear response properties using
      :class LinearResponseEigenSolver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables:
        - pe_models: Dictionary of available PE models (e.g., SEP, CP3).
        - npe_models: Dictionary of available NPE models (e.g., tip3p, ff19sb).
        - pe_model: Selected PE model dictctionary.
        - npe_model: Selected NPE model dictionary.
        - grid_level: The accuracy level of DFT grid.
        - xcfun: The XC functional.
        - eri_thresh: The electron repulsion integrals screening threshold.
        - excited_states: If True, run linear response after SCF.
          response is run as well when any response option is set (e.g., :attr: `nstates`).
        - nstates: Number of excited states to compute in linear response.
        - nto: Whether to run NOT analysis in linear response.
        - core_excitation: Whether to compute core excitations in linear response.
        - num_core_orbitals: Number of core orbitals to consider for core excitations in linear response.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes the ensemble driver.

        PE (SEP/CP3) and NPE (TIP3P/ff19sb) parameter tables are read from the
        ``database/environment_parameters`` directory.

        See this reference, Figure 4, for a summary of
        an overview of SEP/CP3 parametrizations:
        https://doi.org/10.1021/acs.jctc.5c01719
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

        self._potfile_writer = EmbeddingPotfileWriter(comm, ostream)

        # Input-driven state (see update_settings/compute_from_input)
        self.description = 'N/A'
        self.trajectory_file = None
        self.topology_file = None
        self.pdb_file = None
        self.qm_region = ''
        self.qm_charge = 0
        self.qm_multiplicity = 1
        self.pe_cutoff = None
        self.npe_cutoff = None
        self.start = None
        self.end = None
        self.num_snapshots = None
        self.sampling_times = None
        self.last_snapshot_only = False
        self.output_dir = Path('pot_frames')

        self.method_dict = {}
        self.scf_dict = {}
        self.rsp_dict = {}
        self.spectrum_dict = {}

    @property
    def pe_models(self):
        return self._potfile_writer.pe_models

    @pe_models.setter
    def pe_models(self, models):
        self._potfile_writer.pe_models = models

    @property
    def npe_models(self):
        return self._potfile_writer.npe_models

    @npe_models.setter
    def npe_models(self, models):
        self._potfile_writer.npe_models = models

    @property
    def pe_model(self):
        return self._potfile_writer.pe_model

    @pe_model.setter
    def pe_model(self, model):
        self._potfile_writer.pe_model = model

    @property
    def npe_model(self):
        return self._potfile_writer.npe_model

    @npe_model.setter
    def npe_model(self, model):
        self._potfile_writer.npe_model = model

    def set_env_models(self,
                       pe_model: str | list[str] | tuple[str, ...] | None = None,
                       npe_model: str | list[str] | tuple[str, ...] | None = None):
        """Set the PE and/or NPE models used for potential files."""
        return self._potfile_writer.set_env_models(
            pe_model=pe_model, npe_model=npe_model)

    def write_pot_file(self, snapshot, filename):
        """Write one embedding snapshot to an explicitly named potfile."""
        return self._potfile_writer.write_pot_file(snapshot, filename)

    def write_pot_files(self, snapshots, outdir: str | Path = "pot_frames"):
        """Write frame-numbered potfiles for embedding snapshots."""
        return self._potfile_writer.write_pot_files(snapshots, outdir=outdir)

    @staticmethod
    def _parse_model_names(model_names):
        """
        Parse a model selection from an input-file value.

        Accepts a single string, a comma/space-separated string, or a list/tuple
        of such strings. Returns a list of model names or None when empty.
        """
        if model_names is None:
            return None

        if isinstance(model_names, str):
            names = model_names.replace(',', ' ').split()
        else:
            names = []
            for item in model_names:
                names.extend(str(item).replace(',', ' ').split())

        return names or None

    def update_settings(self,
                        ensemble_dict,
                        spectrum_dict=None,
                        scf_dict=None,
                        rsp_dict=None,
                        method_dict=None):
        """
        Updates settings for an input-file driven ensemble calculation.

        The ``@ensemble`` input section is mapped onto the arguments of
        :meth:`EnsembleParser.structures`. The ``@scf`` and ``@response``
        sections are merged with the common ``@method settings`` section to
        build the SCF and response option dictionaries.

        :param ensemble_dict:
            The dictionary of ensemble/trajectory input.
        :param spectrum_dict:
            The dictionary of spectrum settings.
        :param scf_dict:
            The dictionary of SCF settings.
        :param rsp_dict:
            The dictionary of response settings.
        :param method_dict:
            The dictionary of method settings.
        """
        ensemble_dict = dict(ensemble_dict or {})

        self.spectrum_dict = dict(spectrum_dict or {})
        self.scf_dict = dict(scf_dict or {})
        self.rsp_dict = dict(rsp_dict or {})
        self.method_dict = dict(method_dict or {})

        self.description = str(ensemble_dict.get('description', 'N/A'))

        self.trajectory_file = ensemble_dict.get('trajectory_file')
        self.topology_file = ensemble_dict.get('topology_file')
        self.pdb_file = ensemble_dict.get('pdb_file')

        self.qm_region = str(ensemble_dict.get('qm_region', '')).strip()
        self.qm_charge = int(ensemble_dict.get('qm_charge', 0))
        self.qm_multiplicity = int(ensemble_dict.get('qm_multiplicity', 1))

        self.pe_cutoff = None
        self.npe_cutoff = None
        self.start = None
        self.end = None
        self.num_snapshots = None
        self.sampling_times = None
        self.last_snapshot_only = False

        if 'pe_cutoff' in ensemble_dict:
            self.pe_cutoff = float(ensemble_dict['pe_cutoff'])
        if 'npe_cutoff' in ensemble_dict:
            self.npe_cutoff = float(ensemble_dict['npe_cutoff'])

        if 'start' in ensemble_dict:
            self.start = float(ensemble_dict['start'])
        if 'end' in ensemble_dict:
            self.end = float(ensemble_dict['end'])
        if 'num_snapshots' in ensemble_dict:
            self.num_snapshots = int(ensemble_dict['num_snapshots'])
        if 'sampling_times' in ensemble_dict:
            self.sampling_times = np.asarray(
                parse_seq_range(ensemble_dict['sampling_times']), dtype=float)
        if 'last_snapshot_only' in ensemble_dict:
            self.last_snapshot_only = parse_bool(
                ensemble_dict['last_snapshot_only'])

        self.output_dir = Path('pot_frames')
        if 'output_dir' in ensemble_dict:
            self.output_dir = Path(str(ensemble_dict['output_dir']))
        elif 'potdir' in ensemble_dict:
            self.output_dir = Path(str(ensemble_dict['potdir']))
        elif 'filename' in ensemble_dict:
            self.output_dir = Path(str(ensemble_dict['filename']) + '_files')

        pe_model = self._parse_model_names(ensemble_dict.get('pe_model'))
        npe_model = self._parse_model_names(ensemble_dict.get('npe_model'))
        if pe_model is not None or npe_model is not None:
            self.set_env_models(pe_model=pe_model, npe_model=npe_model)

    def compute_from_input(self):
        """
        Runs an input-file driven ensemble calculation.

        The trajectory is parsed with :class:`EnsembleParser` using the settings
        stored by :meth:`update_settings`, followed by :meth:`compute`.

        :return:
            The results dictionary returned by :meth:`compute`.
        """
        basis_set = self.method_dict.get('basis')
        basis_path = self.method_dict.get('basis_path', '.')
        if basis_set is None:
            raise ValueError(
                "EnsembleDriver: missing 'basis' in @method settings.")

        method_options = {
            key: val
            for key, val in self.method_dict.items()
            if key not in {'basis', 'basis_path', 'pe_options'}
        }

        # Common method settings take precedence over explicit @scf/@response
        # settings.
        scf_options = {**self.scf_dict, **method_options}
        property_options = {**self.rsp_dict, **method_options}

        snapshots = EnsembleParser(self.comm, self.ostream).structures(
            trajectory_file=self.trajectory_file,
            topology_file=self.topology_file,
            pdb_file=self.pdb_file,
            qm_region=self.qm_region,
            qm_charge=self.qm_charge,
            qm_multiplicity=self.qm_multiplicity,
            pe_cutoff=self.pe_cutoff,
            npe_cutoff=self.npe_cutoff,
            start=self.start,
            end=self.end,
            num_snapshots=self.num_snapshots,
            sampling_times=self.sampling_times,
            last_snapshot_only=self.last_snapshot_only,
        )

        results = self.compute(
            snapshots,
            basis_set=basis_set,
            scf_options=scf_options,
            property_options=property_options,
            potdir=self.output_dir,
            basis_path=basis_path,
        )

        if (self.spectrum_dict and results.get("rsp_all") and
                self.rank == mpi_master()):
            self.output_dir.mkdir(parents=True, exist_ok=True)

            broadening_type = None
            if 'line_profile' in self.spectrum_dict:
                broadening_type = str(
                    self.spectrum_dict['line_profile']).strip().lower()

            broadening_value_ev = None
            if 'broadening_parameter' in self.spectrum_dict:
                broadening_value = float(
                    self.spectrum_dict['broadening_parameter'])
                units = str(
                    self.spectrum_dict.get('units') or 'ev').strip().lower()
                if units in ('ev', 'electronvolt', 'electronvolts'):
                    broadening_value_ev = broadening_value
                elif units in ('au', 'hartree', 'a.u.'):
                    broadening_value_ev = broadening_value * hartree_in_ev()
                else:
                    raise ValueError(
                        "EnsembleDriver: unsupported broadening units "
                        f"'{units}'.")

            source_file = self.trajectory_file or self.pdb_file
            source_stem = (Path(str(source_file)).stem
                           if source_file else 'ensemble')
            spectrum_file = self.output_dir / (
                f'{source_stem}_averaged_spectra.csv')

            self.save_uv_vis_spectra(
                results,
                broadening_type=broadening_type,
                broadening_value_ev=broadening_value_ev,
                averaged_spectra_filename=spectrum_file,
            )

        return results

    def save_uv_vis_spectra(
        self,
        results: dict,
        energy_min_ev=None,
        energy_max_ev=None,
        broadening_type: str | None = None,
        broadening_value_ev: float | None = None,
        averaged_spectra_filename: str | Path = "averaged_spectra.csv",
    ):
        """
        Computes and saves averaged absorption or ECD spectra data to a CSV
        file.

        Unlike :meth:`plot_uv_vis_spectra`, this method does not require
        matplotlib.

        :param results:
            Results dictionary returned by :meth:`compute`.
        :param energy_min_ev:
            Minimum photon energy in eV for the common grid. If None, it is set
            to (min excitation energy - padding).
        :param energy_max_ev:
            Maximum photon energy in eV for the common grid. If None, it is set
            to (max excitation energy + padding).
        :param broadening_type:
            Broadening line profile, 'gaussian' or 'lorentzian'. If None, the
            SpectrumAverager default is used.
        :param broadening_value_ev:
            Broadening parameter in eV. If None, the SpectrumAverager default is
            used.
        :param averaged_spectra_filename:
            Filename for saving the averaged spectra CSV.
        :return:
            The dictionary of averaged-spectrum data from SpectrumAverager.
        :raises KeyError:
            If `rsp_all` is not present in the results dictionary.
        """
        if not isinstance(results, dict):
            raise TypeError(
                "results must be a dictionary as returned by EnsembleDriver.compute()."
            )

        if "rsp_all" not in results or results.get("rsp_all", None) is None:
            raise KeyError(
                "No 'rsp_all' found in results. Make sure you ran a response "
                "calculation (e.g., pass property_options with 'property' and "
                "either 'nstates' or 'frequencies' to EnsembleDriver.compute())."
            )

        if 'property' not in self.rsp_dict:
            raise ValueError(
                "A 'property' is required for response calculations.")

        rsp_all = results["rsp_all"]
        property_name, is_cpp = self._parse_property_type(
            self.rsp_dict['property'])
        if is_cpp:
            return self._save_cpp_spectra(
                rsp_all, averaged_spectra_filename, property_name)

        spec_avg = SpectrumAverager(comm=self.comm, ostream=self.ostream)
        if broadening_type is not None:
            spec_avg.broadening_type = str(broadening_type).lower()
        if broadening_value_ev is not None:
            spec_avg.broadening_value_ev = float(broadening_value_ev)

        if property_name == 'absorption' and not is_cpp:
            spectrum = spec_avg.compute(
                rsp_all,
                energy_min_ev=energy_min_ev,
                energy_max_ev=energy_max_ev,
            )
            data = np.column_stack((spectrum["wavelength_nm"],
                                    spectrum["eps_mean_wl"]))
            header = "Wavelength [nm],epsilon [L mol^-1 cm^-1]"
        elif property_name == 'ecd' and not is_cpp:
            spectrum = spec_avg.compute_ecd(
                rsp_all,
                energy_min_ev=energy_min_ev,
                energy_max_ev=energy_max_ev,
            )
            data = np.column_stack((spectrum["wavelength_nm"],
                                    spectrum["delta_eps_mean_wl"]))
            header = "Wavelength [nm],Delta epsilon [L mol^-1 cm^-1]"
        else:
            raise ValueError(
                f"Unsupported non-CPP property '{property_name}'.")

        np.savetxt(
            averaged_spectra_filename,
            data,
            delimiter=",",
            header=header,
            comments="",
        )

        return spectrum

    def _save_cpp_spectra(self, rsp_all, averaged_spectra_filename,
                          property_name):
        """
        Averages per-snapshot CPP spectra and saves them to CSV.
        """
        cpp_driver = ComplexResponseSolver(self.comm, self.ostream)
        cpp_driver.set_cpp_property(property_name)

        x_nm = None
        y_all = []
        for item in rsp_all:
            rsp_results = item[1] if isinstance(item, tuple) else item
            spectrum = cpp_driver.get_spectrum(rsp_results, 'nm')

            x = np.asarray(spectrum['x_data'], dtype=float)[::-1].copy()
            y = np.asarray(spectrum['y_data'], dtype=float)[::-1].copy()

            if property_name == 'absorption':
                expected_unit = '[a.u.]'
                a0 = bohr_in_angstrom() * 1.0e-10
                NA = avogadro_constant()
                y = y * a0**2 * 1.0e4 * NA / (np.log(10.0) * 1.0e3)
            elif property_name == 'ecd':
                expected_unit = r'[L mol$^{-1}$ cm$^{-1}$]'
            else:
                raise ValueError(
                    f"Unsupported CPP property '{property_name}'.")

            assert_msg_critical(
                expected_unit in spectrum['y_label'],
                'EnsembleDriver: unexpected CPP spectrum unit.')

            if x_nm is None:
                x_nm = x
            elif not np.allclose(x, x_nm):
                raise ValueError(
                    'CPP spectra have inconsistent frequency grids.')
            y_all.append(y)

        y_mean = np.mean(y_all, axis=0)
        if property_name == 'absorption':
            y_label = 'epsilon [L mol^-1 cm^-1]'
        elif property_name == 'ecd':
            y_label = 'delta epsilon [L mol^-1 cm^-1]'
        else:
            raise ValueError(
                f"Unsupported CPP property '{property_name}'.")

        data = np.column_stack((x_nm, y_mean))
        np.savetxt(
            averaged_spectra_filename,
            data,
            delimiter=',',
            header=f'Wavelength [nm],{y_label}',
            comments='',
        )

        return {
            'wavelength_nm': x_nm,
            'y_mean': y_mean,
            'y_all': y_all,
        }

    @staticmethod
    def _apply_options_to_driver(driver, options: dict, skip_keys: set[str] | None = None):
        """
        Apply dictionary options to a driver using its input keyword types when
        available, and setattr with key validation otherwise.
        """
        if skip_keys is None:
            skip_keys = set()

        keyword_types = {}
        for group, keywords in getattr(driver, '_input_keywords', {}).items():
            if group in ('scf', 'response', 'method_settings'):
                keyword_types.update({
                    key: val[0] for key, val in keywords.items()
                })

        parsed_options = {
            key: val
            for key, val in options.items()
            if key not in skip_keys and key in keyword_types
        }
        parse_input(driver, keyword_types, parsed_options)

        for key, val in options.items():
            if key in skip_keys or key in keyword_types:
                continue
            if not hasattr(driver, key):
                raise ValueError(
                    f"Unknown option '{key}' for {type(driver).__name__}."
                )
            setattr(driver, key, val)

    def _build_scf_driver(self, scf_options: dict):
        """
        Build SCF driver from scf_options['scf_type'].
        Allowed values:
            - 'restricted'           : ScfRestrictedDriver
            - 'unrestricted'         : ScfUnrestrictedDriver
            - 'restricted_openshell' : ScfRestrictedOpenDriver
        """
        scf_type = str(scf_options.get("scf_type", "restricted")).lower()
        if scf_type == "restricted":
            return ScfRestrictedDriver(self.comm, self.ostream)
        if scf_type == "unrestricted":
            return ScfUnrestrictedDriver(self.comm, self.ostream)
        if scf_type == "restricted_openshell":
            return ScfRestrictedOpenDriver(self.comm, self.ostream)
        raise ValueError(
            "Invalid scf_type in scf_options. "
            "Expected one of: 'restricted', 'unrestricted', 'restricted_openshell'."
        )

    @staticmethod
    def _parse_property_type(property_name):
        """
        Parse a response property into its final calculation decision.

        :return:
            A tuple ``(property_name, is_cpp)``, where ``property_name`` is
            either ``'absorption'`` or ``'ecd'``.
        """
        prop = str(property_name).lower().strip()
        is_cpp = '(cpp)' in prop

        if 'absorption' in prop or 'uv-vis' in prop:
            return 'absorption', is_cpp
        if 'ecd' in prop or 'circular dichroism' in prop:
            return 'ecd', is_cpp

        raise ValueError(
            "Invalid property in property_options. Expected one of: "
            "'absorption', 'absorption(cpp)', 'ecd', or 'ecd(cpp)'."
        )

    def _build_property_driver(self, property_options: dict):
        """
        Build a response driver for one of the four supported property types.
        """
        if not property_options:
            return None

        if "property" not in property_options:
            raise ValueError("property_options must contain key 'property'.")

        property_name, is_cpp = self._parse_property_type(
            property_options["property"])
        has_nstates = "nstates" in property_options
        has_freqs = "frequencies" in property_options

        if is_cpp:
            if not has_freqs:
                raise ValueError(
                    "A CPP property requires 'frequencies'.")
            if has_nstates:
                raise ValueError(
                    "A CPP property does not accept 'nstates'.")

            drv = ComplexResponseSolver(self.comm, self.ostream)
            drv.set_cpp_property(property_name)
            return drv
        else:
            if not has_nstates:
                raise ValueError(
                    "A non-CPP property requires 'nstates'.")
            if has_freqs:
                raise ValueError(
                    "A non-CPP property does not accept 'frequencies'.")

            return LinearResponseEigenSolver(self.comm, self.ostream)

    def compute(
        self,
        snapshots,
        basis_set: str,
        scf_options: dict | None = None,
        property_options: dict | None = None,
        potdir: str | Path = "pot_frames",
        write_pe_potfiles: bool = True,
        qm_charge: int | None = None,
        qm_multiplicity: int | None = None,
        embedding_method: str = "fika",
        fika_r_on: float = 1.0,
        fika_r_off: float = 1.9,
        basis_path: str | Path = ".",
    ):
        """
        Drives the computation over the ensemble of snapshots.

        For each snapshot, an SCF calculation is performed for the QM subsystem, with
        optional PE/NPE environment terms. Optional response calculations are
        controlled via property_options.

        :param snapshots:
            A list of snapshot dictionaries (or a single snapshot dict).
        :param basis_set: (str)
            Basis set.
        :param scf_options: (dict)
            SCF/method options merged in one dictionary.
            Includes optional key 'scf_type' in {'restricted', 'unrestricted', 'restricted_openshell'}
            to select the SCF driver.
        :param property_options: (dict)
            Property options dictionary. Supported properties are 'absorption',
            'absorption (cpp)', 'ecd', and 'ecd (cpp)'. TD-DFT properties use
            'nstates'; CPP properties use 'frequencies'.
        :param potdir : (str or Path)
            Directory to store/read the per-snapshot potfiles.
        :param write_pe_potfiles: (bool)
            If True, potfiles are (re)generated from the snapshots before the
            loop. If False, existing potfiles are read from potdir.
        :param qm_charge:
            Optional override for the QM-region charge. If None, the value stored
            in each snapshot is used, defaulting to 0 when absent.
        :param qm_multiplicity:
            Optional override for the QM-region multiplicity. If None, the value stored
            in each snapshot is used, defaulting to 1 when absent.
        :param embedding_method: (str)
            Must be 'fika' (the default). The environment is read from the
            per-snapshot potfile, generated from the snapshots or supplied in
            potdir. The fika embedding uses the inner switching function
            (default r_on/r_off of 1.0/1.9 angstrom, see fika_r_on/fika_r_off)
            on the permanent-charge and induced-dipole interactions of pairs of
            different residues, and reduces anisotropic PE polarizabilities to
            their isotropic part (trace/3).
        :param fika_r_on: (float)
            Inner switching radius (angstrom) of the fika embedding: S(r) = 0
            below it.
        :param fika_r_off: (float)
            Outer switching radius (angstrom) of the fika embedding: S(r) = 1
            above it.
        :param basis_path: (str or Path)
            Directory containing custom basis-set files.

        :return:
            Dictionary with keys:
            - scf_all: list of (frame, scf_results)
            - rsp_all: list of (frame, rsp_results), only present if response is run.

        :raises RuntimeError:
            If required PE/NPE models have not been selected.
        :raises ValueError:
            If embedding_method is not 'fika' or a required potfile is missing.
        """

        if isinstance(snapshots, dict):
            snapshots = [snapshots]

        if scf_options is None:
            scf_options = {}
        else:
            scf_options = dict(scf_options)

        if property_options is None:
            property_options = {}

        do_rsp = bool(property_options)
        if do_rsp and "property" not in property_options:
            raise ValueError(
                "property_options must contain key 'property'.")

        if qm_charge is not None:
            qm_charge = int(qm_charge)
        if qm_multiplicity is not None:
            qm_multiplicity = int(qm_multiplicity)
            if qm_multiplicity <= 0:
                raise ValueError("QM multiplicity must be a positive integer.")

        potdir = Path(potdir)

        embedding_method = str(embedding_method).lower()
        if embedding_method != "fika":
            raise ValueError("embedding_method must be 'fika'.")

        # Detect whether PE / NPE models are needed
        has_any_pe = any(np.asarray(s.get("pe_coords", [])).size > 0 for s in snapshots)
        has_any_npe = any(np.asarray(s.get("npe_coords", [])).size > 0 for s in snapshots)

        # Only require environment models when corresponding environment atoms exist
        if has_any_pe and self.pe_model is None:
            raise RuntimeError(
                "Snapshots contain PE atoms but pe_model is not set. "
                "Call set_env_models(pe_model=...)."
            )

        if has_any_npe and self.npe_model is None:
            raise RuntimeError(
                "Snapshots contain NPE atoms but npe_model is not set. "
                "Call set_env_models(npe_model=...)."
            )

        for snap in snapshots:
            self._potfile_writer.validate_partition(snap)

        if write_pe_potfiles and (has_any_pe or has_any_npe):
            self.write_pot_files(snapshots, outdir=potdir)

        # Guard against open-shell response requests for now:
        if do_rsp:
            scf_type = str(scf_options.get("scf_type", "restricted")).lower()
            if scf_type != "restricted":
                raise NotImplementedError(
                    "Open-shell response is not yet supported in EnsembleDriver.compute. "
                    "Use scf_type='restricted' for response calculations."
                )

        scf_all = []
        rsp_all = [] if do_rsp else None

        for snap in snapshots:
            scf_driver = self._build_scf_driver(scf_options)
            self._apply_options_to_driver(scf_driver, scf_options, skip_keys={"scf_type"})

            frame = int(snap["frame"])

            labels = [str(x) for x in snap["qm_elements"]]
            coords = np.asarray(snap["qm_coords"], dtype=float)

            snap_qm_charge = (
                qm_charge
                if qm_charge is not None
                else int(snap.get("qm_charge", 0))
            )
            snap_qm_multiplicity = (
                qm_multiplicity
                if qm_multiplicity is not None
                else int(snap.get("qm_multiplicity", 1))
            )

            molecule = Molecule(labels, coords)
            molecule.set_charge(snap_qm_charge)
            molecule.set_multiplicity(snap_qm_multiplicity)

            if not molecule.check_multiplicity():
                raise ValueError(
                    f"Incompatible QM charge ({snap_qm_charge}) and multiplicity"
                    f"({snap_qm_multiplicity}) for frame {frame}."
                )
            basis = MolecularBasis.read(molecule, basis_set, str(basis_path))

            pe_coords = np.asarray(snap.get("pe_coords", []), dtype=float)
            npe_coords = np.asarray(snap.get("npe_coords", []), dtype=float)

            has_pe = pe_coords.size > 0
            has_npe = npe_coords.size > 0

            # Reset embedding inputs every frame
            scf_driver.potfile = None
            scf_driver.point_charges = None
            scf_driver.embedding = None
            fika_embedding = None

            if has_pe or has_npe:
                potfile = potdir / f"pe_frame_{frame:06d}.pot"
                if not potfile.is_file():
                    raise ValueError(
                        f"Missing potfile for frame {frame}: {potfile}. " +
                        "Set write_pe_potfiles=True or provide the file.")
                fika_embedding = {
                    "settings": {
                        "embedding_method": "fika",
                        "damping": None,
                        "switching": {
                            "r_on": float(fika_r_on),
                            "r_off": float(fika_r_off),
                        },
                    },
                    "inputs": {
                        "potfile": str(potfile),
                    },
                }
                scf_driver.embedding = fika_embedding

            scf_results = scf_driver.compute(molecule, basis)
            scf_all.append((frame, scf_results))

            if do_rsp:
                rsp_driver = self._build_property_driver(property_options)
                # property is routing metadata; do not setattr it on the driver
                self._apply_options_to_driver(
                    rsp_driver,
                    property_options,
                    skip_keys={"property"},
                )
                rsp_driver.embedding = fika_embedding
                rsp_results = rsp_driver.compute(molecule, basis, scf_results)
                rsp_all.append((frame, rsp_results))

        results = {"scf_all": scf_all}
        if do_rsp:
            results["rsp_all"] = rsp_all
        return results

    def plot_uv_vis_spectra(
        self,
        results: dict,
        energy_min_ev=None,
        energy_max_ev=None,
        broadening_type: str | None = None,
        broadening_value_ev: float | None = None,
        show_individual: bool = False,
        show_sticks: bool = True,
        show_std: bool = False,
        title: str = "Absorption Spectrum (Averaged)",
        ax=None,
        xlim_nm=None,
        save_averaged_spectra: bool = True,
        averaged_spectra_filename: str | Path = "averaged_spectra.csv",
    ):
        """
        Convenience wrapper to plot averaged UV/Vis spectra for an ensemble.

        This method expects the `results` dictionary returned by :meth:`compute`.
        It extracts `results["rsp_all"]` and passes it to
        :class:`SpectrumAverager`.

        :param results:
            Results dictionary returned by :meth:`compute`.
        :param energy_min_ev:
            Minimum photon energy in eV for the common grid. If None, it is set to
            (min excitation energy - padding).
        :param energy_max_ev:
            Maximum photon energy in eV for the common grid. If None, it is set to
            (max excitation energy + padding).
        :param broadening_type:
            Broadening line profile, 'gaussian' or 'lorentzian'. If None, the
            SpectrumAverager default is used.
        :param broadening_value_ev:
            Broadening parameter in eV. If None, the SpectrumAverager default is
            used.
        :param show_individual:
            If True, plot individual broadened spectra for each snapshot.
        :param show_sticks:
            If True, plots oscillator strength sticks at transition wavelengths.
        :param show_std:
            If True, shows a shaded area corresponding to +/- one standard deviation.
        :param title:
            Title of the plot.
        :param ax:
            Matplotlib Axes object to plot on. If None, a new figure and axes are
            created.
        :param xlim_nm:
            Tuple (xmin, xmax) to set x-axis limits in nm. If None, automatic limits are used.
        :param save_averaged_spectra:
            If True, saves the averaged spectra data to a CSV.
        :param averaged_spectra_filename:
            Filename for saving the averaged spectra CSV if `save_averaged_spectra` is True.
            Default is "averaged_spectra.csv" in current working directory.
        :return:
            The Matplotlib Axes object containing the plot.
        :raises KeyError:
            If `rsp_all` is not present in the results dictionary.
        """
        if not isinstance(results, dict):
            raise TypeError(
                "results must be a dictionary as returned by EnsembleDriver.compute()."
            )

        if "rsp_all" not in results or results.get("rsp_all", None) is None:
            raise KeyError(
                "No 'rsp_all' found in results. Make sure you ran a response calculation "
                "(e.g., pass property_optios with 'property' and either "
                "'nstates' or 'frequencies' to EnsembleDriver.compute())"
            )

        rsp_all = results["rsp_all"]

        spec_avg = SpectrumAverager(comm=self.comm, ostream=self.ostream)
        if broadening_type is not None:
            spec_avg.broadening_type = str(broadening_type).lower()
        if broadening_value_ev is not None:
            spec_avg.broadening_value_ev = float(broadening_value_ev)
        ax = spec_avg.plot_uv_vis_spectra(
            rsp_all,
            energy_min_ev=energy_min_ev,
            energy_max_ev=energy_max_ev,
            show_individual=show_individual,
            show_sticks=show_sticks,
            show_std=show_std,
            title=title,
            ax=ax,
            xlim_nm=xlim_nm,
            save_averaged_spectra=save_averaged_spectra,
            averaged_spectra_filename=averaged_spectra_filename,
        )
        return ax
