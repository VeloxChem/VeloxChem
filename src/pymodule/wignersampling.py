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
#  modification, are permitted provided that the following conditions
#  are met:
#
#  1. Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in
#     the documentation and/or other materials provided with the
#     distribution.
#  3. Neither the name of the copyright holder nor the names of its
#     contributors may be used to endorse or promote products derived
#     from this software without specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES ARE DISCLAIMED.
#

import numpy as np
import h5py


class WignerSamplingDriver:
    """
    Generates nuclear initial conditions from the harmonic Wigner
    distribution using the results of a VeloxChem VibrationalAnalysis.

    The sampling follows the independent harmonic-oscillator Wigner
    distribution used for initial-condition generation in programs such
    as Newton-X.

    The VeloxChem vibrational-analysis results are expected to contain

        vib_frequencies
            Harmonic frequencies in cm^-1.
            Shape: (n_modes,)

        normal_modes
            Normalized Cartesian normal-mode displacement vectors.
            Shape: (n_modes, n_atoms, 3)

        reduced_masses
            Vibrational reduced masses in atomic mass units.
            Shape: (n_modes,)

        molecule_xyz_string
            Equilibrium geometry used for the vibrational analysis.

    The generated ensemble contains both Cartesian coordinates and
    Cartesian velocities.

    Notes
    -----
    The vibrational normal modes returned by VeloxChem already contain
    only the vibrational degrees of freedom. Therefore translational
    and rotational modes are not removed again in this driver.

    Internal Wigner sampling is carried out in atomic units.
    """

    # ============================================================
    # Physical constants
    # ============================================================

    # Boltzmann constant in Hartree / K
    _kb_au = 3.166811563e-6

    # 1 atomic mass unit in electron masses
    _amu_to_me = 1822.888486209

    # Conversion between Hartree and cm^-1
    #
    # E / Eh = wavenumber / 219474.63136320
    #
    # Since hbar = 1 in atomic units, this is also the angular
    # frequency omega in atomic units.
    _hartree_to_wavenumber = 219474.63136320

    # Bohr -> Angstrom
    _bohr_to_angstrom = 0.529177210903

    # Atomic unit of time -> fs
    _au_time_to_fs = 0.024188843265857

    # Atomic unit velocity -> Angstrom / fs
    _au_velocity_to_angstrom_fs = (
        _bohr_to_angstrom /
        _au_time_to_fs
    )

    def __init__(self, ostream=None):
        """
        Initializes the Wigner sampling driver.
        """

        self.ostream = ostream

        # ========================================================
        # User settings
        # ========================================================
        self.temperature = 0.0
        self.n_samples = 1
        self.seed = None
        # Modes below this frequency can optionally be excluded.
        #
        # Default 0.0 means that all real vibrational modes are used.
        self.min_frequency = 0.0
        # What to do with imaginary frequencies:
        #
        #     "error"
        #     "exclude"
        self.imaginary_mode_policy = "error"
        # Small numerical frequencies with |nu| below this threshold
        # are considered zero.
        self.zero_frequency_threshold = 1.0
        # Store sampled normal coordinates/momenta in the results.
        self.store_normal_coordinates = True
        # ========================================================
        # Internal data
        # ========================================================
        self._rng = None

        self.write_initial_conditions = True

        self.output_directory = "wigner_initial_conditions"
        self.xyz_filename = "initial_conditions.xyz"
        self.hdf5_filename = "initial_conditions.h5"

        self.hdf5_compression = "gzip"
        self.hdf5_compression_level = 4

    # ============================================================
    # Public API
    # ============================================================

    def compute(self, molecule, vib_results):
        """
        Generates a Wigner ensemble.

        Parameters
        ----------
        molecule
            VeloxChem Molecule corresponding to the equilibrium
            geometry.

        vib_results : dict
            Dictionary returned by VibrationalAnalysis.compute().

        Returns
        -------
        dict
            Wigner ensemble.
        """

        self._validate_settings()

        self._rng = np.random.default_rng(
            self.seed
        )

        data = self._extract_vibrational_data(
            molecule,
            vib_results,
        )

        frequencies = data["frequencies_cm-1"]

        normal_modes = data["normal_modes"]

        reduced_masses = data["reduced_masses"]

        equilibrium_coordinates = data[
            "equilibrium_coordinates_bohr"
        ]

        # --------------------------------------------------------
        # Select physically usable modes
        # --------------------------------------------------------

        (
            frequencies,
            normal_modes,
            reduced_masses,
            mode_indices,
        ) = self._select_modes(
            frequencies,
            normal_modes,
            reduced_masses,
        )

        # --------------------------------------------------------
        # Convert frequencies to atomic units
        # --------------------------------------------------------

        omega = self._frequencies_to_au(
            frequencies
        )

        # --------------------------------------------------------
        # Calculate Wigner widths
        # --------------------------------------------------------

        (
            sigma_q,
            sigma_p,
        ) = self._compute_wigner_widths(
            omega
        )

        # --------------------------------------------------------
        # Sample dimensionless/mass-weighted oscillator variables
        # --------------------------------------------------------

        q, p = self._sample_phase_space(
            sigma_q,
            sigma_p,
        )

        # --------------------------------------------------------
        # Transform to Cartesian coordinates and velocities
        # --------------------------------------------------------

        (
            coordinates,
            velocities,
        ) = self._transform_to_cartesian(
            equilibrium_coordinates,
            normal_modes,
            reduced_masses,
            q,
            p,
        )

        # --------------------------------------------------------
        # Build result dictionary
        # --------------------------------------------------------

        results = {
            "coordinates_bohr": coordinates,
            "velocities_au": velocities,

            "coordinates_angstrom":
                coordinates * self._bohr_to_angstrom,

            "velocities_angstrom_fs":
                velocities *
                self._au_velocity_to_angstrom_fs,

            "vib_frequencies":
                frequencies.copy(),

            "reduced_masses":
                reduced_masses.copy(),

            "mode_indices":
                mode_indices.copy(),

            "temperature":
                float(self.temperature),

            "n_samples":
                int(self.n_samples),

            "seed":
                self.seed,

            "sampling_method":
                "harmonic_wigner",
        }

        if self.store_normal_coordinates:

            results["normal_coordinates"] = q

            results["normal_momenta"] = p

        # --------------------------------------------------------
        # Statistical diagnostics
        # --------------------------------------------------------

        results["statistics"] = (
            self._compute_statistics(
                q,
                p,
                omega,
            )
        )

        self._print_summary(
            frequencies,
            reduced_masses,
        )

        self.validate_cartesian_transformation(vib_results=vib_results)

        if self.write_initial_conditions:
            self._write_initial_conditions(molecule, results)

        return results

    # ============================================================
    # Initial-condition output
    # ============================================================

    def _write_initial_conditions(
        self,
        molecule,
        results,
    ):
        """
        Writes the generated Wigner ensemble.

        Three files are generated:

            initial_conditions.xyz

                Multi-frame XYZ containing all sampled geometries.

            initial_conditions.dat

                Human-readable initial-condition file containing
                coordinates and velocities.

            initial_conditions.h5

                HDF5 container containing the complete Wigner ensemble
                in both atomic units and convenient external units.
                This is the preferred file for VeloxChem trajectory
                initialization.
        """

        from pathlib import Path

        output_directory = Path(
            self.output_directory
        )

        output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        labels = self._get_atomic_labels(
            molecule
        )

        # ========================================================
        # XYZ
        # ========================================================

        xyz_path = (
            output_directory /
            self.xyz_filename
        )

        self._write_xyz_ensemble(
            xyz_path,
            labels,
            results,
        )

        # ========================================================
        # Human-readable coordinates + velocities
        # ========================================================

        # ic_path = (
        #     output_directory /
        #     self.initial_conditions_filename
        # )

        # self._write_initial_condition_file(
        #     ic_path,
        #     labels,
        #     results,
        # )

        # ========================================================
        # HDF5
        # ========================================================

        h5_path = (
            output_directory /
            self.hdf5_filename
        )

        self._write_hdf5_initial_conditions(
            h5_path,
            labels,
            results,
        )

        # ========================================================
        # Store paths
        # ========================================================

        results["xyz_file"] = str(
            xyz_path
        )

        # results["initial_conditions_file"] = str(
        #     ic_path
        # )

        results["hdf5_file"] = str(
            h5_path
        )

        self._print_output_files(
            xyz_path,
            # ic_path,
            h5_path,
        )


    # ============================================================
    # Atomic labels
    # ============================================================

    def _get_atomic_labels(
        self,
        molecule,
    ):
        """
        Returns atomic element labels from a VeloxChem Molecule.
        """

        # Different VeloxChem versions may expose the labels through
        # slightly different interfaces.

        if hasattr(
            molecule,
            "get_labels",
        ):

            labels = molecule.get_labels()

        elif hasattr(
            molecule,
            "get_element_labels",
        ):

            labels = molecule.get_element_labels()

        else:

            raise AttributeError(
                "Unable to obtain atomic labels from Molecule."
            )

        return [
            str(label)
            for label in labels
        ]


    # ============================================================
    # XYZ ensemble
    # ============================================================

    def _write_xyz_ensemble(
        self,
        filename,
        labels,
        results,
    ):
        """
        Writes all Wigner geometries as a standard multi-frame XYZ file.

        Coordinates are written in Angstrom.
        """

        coordinates = results[
            "coordinates_angstrom"
        ]

        frequencies = results[
            "vib_frequencies"
        ]

        natoms = len(labels)

        with open(
            filename,
            "w",
            encoding="utf-8",
        ) as handle:

            for sample_index in range(
                self.n_samples
            ):

                handle.write(
                    f"{natoms}\n"
                )

                handle.write(
                    "Wigner initial condition "
                    f"{sample_index + 1} / "
                    f"{self.n_samples}; "
                    f"T = {self.temperature:.2f} K; "
                    f"seed = {self.seed}\n"
                )

                geometry = coordinates[
                    sample_index
                ]

                for atom_index in range(
                    natoms
                ):

                    x, y, z = geometry[
                        atom_index
                    ]

                    handle.write(
                        f"{labels[atom_index]:<3s}"
                        f" {x:20.12f}"
                        f" {y:20.12f}"
                        f" {z:20.12f}\n"
                    )


    # ============================================================
    # Geometry + velocity initial-condition file
    # ============================================================

    def _write_initial_condition_file(
        self,
        filename,
        labels,
        results,
    ):
        """
        Writes coordinates and velocities for every Wigner initial
        condition.

        Coordinates are written in Angstrom.

        Velocities are written in Angstrom/fs.

        The format is intentionally simple and self-contained so that
        it can be parsed directly by a VeloxChem trajectory driver.
        """

        coordinates = results[
            "coordinates_angstrom"
        ]

        velocities = results[
            "velocities_angstrom_fs"
        ]

        natoms = len(labels)

        with open(
            filename,
            "w",
            encoding="utf-8",
        ) as handle:

            # ----------------------------------------------------
            # File header
            # ----------------------------------------------------

            handle.write(
                "# VELOXCHEM WIGNER INITIAL CONDITIONS\n"
            )

            handle.write("#\n")

            handle.write(
                "# Harmonic Wigner sampling\n"
            )

            handle.write(
                f"# Temperature_K = "
                f"{self.temperature:.8f}\n"
            )

            handle.write(
                f"# Number_of_samples = "
                f"{self.n_samples}\n"
            )

            handle.write(
                f"# Number_of_atoms = "
                f"{natoms}\n"
            )

            handle.write(
                f"# Random_seed = "
                f"{self.seed}\n"
            )

            handle.write("#\n")

            handle.write(
                "# Coordinate units: Angstrom\n"
            )

            handle.write(
                "# Velocity units: Angstrom/fs\n"
            )

            handle.write("#\n")

            handle.write(
                "# Columns:\n"
            )

            handle.write(
                "# element "
                "x y z "
                "vx vy vz\n"
            )

            handle.write("#\n")

            # ----------------------------------------------------
            # Initial conditions
            # ----------------------------------------------------

            for sample_index in range(
                self.n_samples
            ):

                handle.write(
                    "BEGIN_INITIAL_CONDITION "
                    f"{sample_index}\n"
                )

                handle.write(
                    f"INDEX {sample_index}\n"
                )

                handle.write(
                    f"NATOMS {natoms}\n"
                )

                handle.write(
                    f"TEMPERATURE_K "
                    f"{self.temperature:.8f}\n"
                )

                handle.write(
                    "ATOMS\n"
                )

                geometry = coordinates[
                    sample_index
                ]

                velocity = velocities[
                    sample_index
                ]

                for atom_index in range(
                    natoms
                ):

                    x, y, z = geometry[
                        atom_index
                    ]

                    vx, vy, vz = velocity[
                        atom_index
                    ]

                    handle.write(
                        f"{labels[atom_index]:<3s}"
                        f" {x:20.12f}"
                        f" {y:20.12f}"
                        f" {z:20.12f}"
                        f" {vx:20.12f}"
                        f" {vy:20.12f}"
                        f" {vz:20.12f}\n"
                    )

                handle.write(
                    "END_INITIAL_CONDITION\n\n"
                )


    def _write_hdf5_initial_conditions(
        self,
        filename,
        labels,
        results,
    ):
        """
        Writes the complete Wigner ensemble to an HDF5 file.

        Atomic units are regarded as the native representation for
        subsequent VeloxChem dynamics calculations.

        HDF5 layout
        -----------

        /metadata
            temperature_K
            n_samples
            n_atoms
            n_modes
            seed
            sampling_method

        /molecule
            labels

        /vibrational_analysis
            frequencies_cm1
            reduced_masses_amu
            mode_indices

        /initial_conditions
            coordinates_bohr
            velocities_au
            coordinates_angstrom
            velocities_angstrom_fs

        /normal_modes
            Q
            P
        """

        try:
            import h5py

        except ImportError as exc:

            raise ImportError(
                "Writing Wigner initial conditions requires h5py."
            ) from exc

        labels = np.asarray(
            labels,
            dtype="S",
        )

        coordinates_bohr = np.asarray(
            results["coordinates_bohr"],
            dtype=np.float64,
        )

        velocities_au = np.asarray(
            results["velocities_au"],
            dtype=np.float64,
        )

        coordinates_angstrom = np.asarray(
            results["coordinates_angstrom"],
            dtype=np.float64,
        )

        velocities_angstrom_fs = np.asarray(
            results["velocities_angstrom_fs"],
            dtype=np.float64,
        )

        frequencies = np.asarray(
            results["vib_frequencies"],
            dtype=np.float64,
        )

        reduced_masses = np.asarray(
            results["reduced_masses"],
            dtype=np.float64,
        )

        mode_indices = np.asarray(
            results["mode_indices"],
            dtype=np.int64,
        )

        nsamples = coordinates_bohr.shape[0]
        natoms = coordinates_bohr.shape[1]
        nmodes = frequencies.size

        # ========================================================
        # Create HDF5 file
        # ========================================================

        with h5py.File(
            filename,
            "w",
        ) as h5file:

            # ====================================================
            # Root attributes
            # ====================================================

            h5file.attrs["file_type"] = (
                "VeloxChem Wigner initial conditions"
            )

            h5file.attrs["format_version"] = "1.0"

            h5file.attrs["sampling_method"] = (
                "harmonic_wigner"
            )

            # ====================================================
            # Metadata
            # ====================================================

            metadata = h5file.create_group(
                "metadata"
            )

            metadata.create_dataset(
                "temperature_K",
                data=np.float64(
                    self.temperature
                ),
            )

            metadata.create_dataset(
                "n_samples",
                data=np.int64(
                    nsamples
                ),
            )

            metadata.create_dataset(
                "n_atoms",
                data=np.int64(
                    natoms
                ),
            )

            metadata.create_dataset(
                "n_modes",
                data=np.int64(
                    nmodes
                ),
            )

            metadata.create_dataset(
                "seed",
                data=np.int64(
                    -1
                    if self.seed is None
                    else self.seed
                ),
            )

            string_dtype = h5py.string_dtype(
                encoding="utf-8"
            )

            metadata.create_dataset(
                "sampling_method",
                data="harmonic_wigner",
                dtype=string_dtype,
            )

            # ====================================================
            # Molecule
            # ====================================================

            molecule_group = h5file.create_group(
                "molecule"
            )

            molecule_group.create_dataset(
                "labels",
                data=labels,
            )

            # ====================================================
            # Vibrational analysis
            # ====================================================

            vib_group = h5file.create_group(
                "vibrational_analysis"
            )

            frequency_dataset = (
                vib_group.create_dataset(
                    "frequencies_cm1",
                    data=frequencies,
                )
            )

            frequency_dataset.attrs["units"] = (
                "cm^-1"
            )

            reduced_mass_dataset = (
                vib_group.create_dataset(
                    "reduced_masses_amu",
                    data=reduced_masses,
                )
            )

            reduced_mass_dataset.attrs["units"] = (
                "amu"
            )

            vib_group.create_dataset(
                "mode_indices",
                data=mode_indices,
            )

            # ====================================================
            # Initial conditions
            # ====================================================

            ic_group = h5file.create_group(
                "initial_conditions"
            )

            # ----------------------------------------------------
            # Native atomic-unit data
            # ----------------------------------------------------

            coordinate_dataset = (
                ic_group.create_dataset(
                    "coordinates_bohr",
                    data=coordinates_bohr,
                    compression=self.hdf5_compression,
                    compression_opts=(
                        self.hdf5_compression_level
                        if self.hdf5_compression == "gzip"
                        else None
                    ),
                    shuffle=True,
                    chunks=(
                        1,
                        natoms,
                        3,
                    ),
                )
            )

            coordinate_dataset.attrs["units"] = (
                "bohr"
            )

            coordinate_dataset.attrs[
                "description"
            ] = (
                "Cartesian nuclear coordinates"
            )

            velocity_dataset = (
                ic_group.create_dataset(
                    "velocities_au",
                    data=velocities_au,
                    compression=self.hdf5_compression,
                    compression_opts=(
                        self.hdf5_compression_level
                        if self.hdf5_compression == "gzip"
                        else None
                    ),
                    shuffle=True,
                    chunks=(
                        1,
                        natoms,
                        3,
                    ),
                )
            )

            velocity_dataset.attrs["units"] = (
                "bohr / atomic_time"
            )

            velocity_dataset.attrs[
                "description"
            ] = (
                "Cartesian nuclear velocities"
            )

            # ----------------------------------------------------
            # Human-friendly units
            # ----------------------------------------------------

            coordinate_angstrom_dataset = (
                ic_group.create_dataset(
                    "coordinates_angstrom",
                    data=coordinates_angstrom,
                    compression=self.hdf5_compression,
                    compression_opts=(
                        self.hdf5_compression_level
                        if self.hdf5_compression == "gzip"
                        else None
                    ),
                    shuffle=True,
                    chunks=(
                        1,
                        natoms,
                        3,
                    ),
                )
            )

            coordinate_angstrom_dataset.attrs[
                "units"
            ] = "angstrom"

            velocity_external_dataset = (
                ic_group.create_dataset(
                    "velocities_angstrom_fs",
                    data=velocities_angstrom_fs,
                    compression=self.hdf5_compression,
                    compression_opts=(
                        self.hdf5_compression_level
                        if self.hdf5_compression == "gzip"
                        else None
                    ),
                    shuffle=True,
                    chunks=(
                        1,
                        natoms,
                        3,
                    ),
                )
            )

            velocity_external_dataset.attrs[
                "units"
            ] = "angstrom / fs"

            # ====================================================
            # Sampled normal coordinates
            # ====================================================

            if (
                "normal_coordinates" in results
                and
                "normal_momenta" in results
            ):

                normal_group = h5file.create_group(
                    "normal_modes"
                )

                q = np.asarray(
                    results[
                        "normal_coordinates"
                    ],
                    dtype=np.float64,
                )

                p = np.asarray(
                    results[
                        "normal_momenta"
                    ],
                    dtype=np.float64,
                )

                q_dataset = (
                    normal_group.create_dataset(
                        "Q",
                        data=q,
                        compression=self.hdf5_compression,
                        compression_opts=(
                            self.hdf5_compression_level
                            if self.hdf5_compression == "gzip"
                            else None
                        ),
                        shuffle=True,
                        chunks=(
                            min(
                                256,
                                nsamples,
                            ),
                            nmodes,
                        ),
                    )
                )

                q_dataset.attrs[
                    "description"
                ] = (
                    "Sampled mass-weighted "
                    "normal coordinates"
                )

                p_dataset = (
                    normal_group.create_dataset(
                        "P",
                        data=p,
                        compression=self.hdf5_compression,
                        compression_opts=(
                            self.hdf5_compression_level
                            if self.hdf5_compression == "gzip"
                            else None
                        ),
                        shuffle=True,
                        chunks=(
                            min(
                                256,
                                nsamples,
                            ),
                            nmodes,
                        ),
                    )
                )

                p_dataset.attrs[
                    "description"
                ] = (
                    "Sampled mass-weighted "
                    "normal momenta"
                )

    # ============================================================
    # Output summary
    # ============================================================

    def _print_output_files(
        self,
        xyz_path,
        # ic_path,
        h5_path,
    ):

        lines = [
            "",
            "Wigner initial-condition files:",
            "",
            f"XYZ ensemble             : {xyz_path}",
            # f"Coordinates + velocities : {ic_path}",
            f"HDF5 ensemble            : {h5_path}",
            "",
        ]

        if self.ostream is None:

            print(
                "\n".join(lines)
            )

        else:

            for line in lines:

                try:

                    self.ostream.print_info(
                        line
                    )

                except AttributeError:

                    print(line)

    @staticmethod
    def read_initial_condition(
        filename,
        index,
    ):
        """
        Reads one Wigner initial condition from a VeloxChem HDF5
        ensemble.

        Returns
        -------
        dict
            Cartesian coordinates in bohr and velocities in atomic units.
        """

        with h5py.File(
            filename,
            "r",
        ) as h5file:

            n_samples = int(
                h5file[
                    "metadata/n_samples"
                ][()]
            )

            if index < 0 or index >= n_samples:

                raise IndexError(
                    f"Initial-condition index {index} "
                    f"is outside the available range "
                    f"0 <= index < {n_samples}."
                )

            coordinates = h5file[
                "initial_conditions/coordinates_bohr"
            ][index].copy()

            velocities = h5file[
                "initial_conditions/velocities_au"
            ][index].copy()

            labels = h5file[
                "molecule/labels"
            ][...]

            labels = [
                label.decode("utf-8")
                if isinstance(label, bytes)
                else str(label)
                for label in labels
            ]

        return {
            "index": int(index),
            "labels": labels,
            "coordinates_bohr": coordinates,
            "velocities_au": velocities,
        }
    # ============================================================
    # Validation
    # ============================================================

    def _validate_settings(self):
        """
        Validates user-defined settings.
        """

        if self.temperature < 0.0:

            raise ValueError(
                "Wigner sampling temperature cannot be negative."
            )

        if self.n_samples < 1:

            raise ValueError(
                "Number of Wigner samples must be >= 1."
            )

        if self.min_frequency < 0.0:

            raise ValueError(
                "Minimum vibrational frequency cannot be negative."
            )

        allowed = (
            "error",
            "exclude",
        )

        if self.imaginary_mode_policy not in allowed:

            raise ValueError(
                "imaginary_mode_policy must be "
                "'error' or 'exclude'."
            )

    # ============================================================
    # VeloxChem vibrational data
    # ============================================================

    def _extract_vibrational_data(
        self,
        molecule,
        vib_results,
    ):
        """
        Extracts the quantities required for Wigner sampling from the
        actual VeloxChem VibrationalAnalysis result dictionary.
        """

        required_keys = (
            "vib_frequencies",
            "normal_modes",
            "reduced_masses",
            "number_of_modes",
        )

        for key in required_keys:

            if key not in vib_results:

                raise KeyError(
                    "Vibrational-analysis results do not contain "
                    f"the required key '{key}'."
                )

        frequencies = np.asarray(
            vib_results["vib_frequencies"],
            dtype=float,
        )

        normal_modes = np.asarray(
            vib_results["normal_modes"],
            dtype=float,
        )

        reduced_masses = np.asarray(
            vib_results["reduced_masses"],
            dtype=float,
        )

        number_of_modes = int(
            vib_results["number_of_modes"]
        )

        # --------------------------------------------------------
        # Equilibrium coordinates
        #
        # Prefer the Molecule object passed to compute().
        # VeloxChem's coordinate API is assumed to return bohr here.
        # --------------------------------------------------------

        equilibrium_coordinates = np.asarray(
            molecule.get_coordinates_in_bohr(),
            dtype=float,
        )

        natoms = equilibrium_coordinates.shape[0]

        # --------------------------------------------------------
        # Validate shapes
        # --------------------------------------------------------

        if frequencies.ndim != 1:

            raise ValueError(
                "vib_frequencies must be one-dimensional."
            )

        if reduced_masses.ndim != 1:

            raise ValueError(
                "reduced_masses must be one-dimensional."
            )

        expected_mode_shape = (
            number_of_modes,
            natoms,
            3,
        )

        if normal_modes.shape != expected_mode_shape:

            raise ValueError(
                "Unexpected normal_modes shape.\n"
                f"Expected: {expected_mode_shape}\n"
                f"Received: {normal_modes.shape}"
            )

        if frequencies.size != number_of_modes:

            raise ValueError(
                "number_of_modes does not agree with "
                "vib_frequencies."
            )

        if reduced_masses.size != number_of_modes:

            raise ValueError(
                "number_of_modes does not agree with "
                "reduced_masses."
            )

        # --------------------------------------------------------
        # Basic reduced-mass validation
        # --------------------------------------------------------

        if np.any(reduced_masses <= 0.0):

            raise ValueError(
                "All vibrational reduced masses must be positive."
            )

        return {
            "frequencies_cm-1":
                frequencies,

            "normal_modes":
                normal_modes,

            "reduced_masses":
                reduced_masses,

            "equilibrium_coordinates_bohr":
                equilibrium_coordinates,
        }

    # ============================================================
    # Mode selection
    # ============================================================

    def _select_modes(
        self,
        frequencies,
        normal_modes,
        reduced_masses,
    ):
        """
        Selects modes used for Wigner sampling.

        VeloxChem has already removed translations and rotations, so
        this routine acts only on the actual vibrational modes.
        """

        keep = []

        for index, frequency in enumerate(
            frequencies
        ):

            # ----------------------------------------------------
            # Imaginary frequency
            # ----------------------------------------------------

            if frequency < -self.zero_frequency_threshold:

                if self.imaginary_mode_policy == "error":

                    raise ValueError(
                        "Imaginary vibrational frequency detected "
                        f"for mode {index + 1}: "
                        f"{frequency:.6f} cm^-1.\n"
                        "A harmonic Wigner distribution is not "
                        "defined for an unstable normal mode."
                    )

                elif self.imaginary_mode_policy == "exclude":

                    continue

            # ----------------------------------------------------
            # Numerical zero
            # ----------------------------------------------------

            if abs(frequency) <= self.zero_frequency_threshold:

                continue

            # ----------------------------------------------------
            # Low-frequency cutoff
            # ----------------------------------------------------

            if frequency < self.min_frequency:

                continue

            keep.append(index)

        if not keep:

            raise ValueError(
                "No vibrational modes remain for Wigner sampling."
            )

        keep = np.asarray(
            keep,
            dtype=int,
        )

        return (
            frequencies[keep],
            normal_modes[keep],
            reduced_masses[keep],
            keep,
        )

    # ============================================================
    # Frequency conversion
    # ============================================================

    def _frequencies_to_au(
        self,
        frequencies_cm,
    ):
        """
        Converts harmonic frequencies from cm^-1 to atomic units.
        """

        return (
            frequencies_cm /
            self._hartree_to_wavenumber
        )

    # ============================================================
    # Thermal Wigner widths
    # ============================================================

    def _compute_wigner_widths(
        self,
        omega,
    ):
        """
        Computes harmonic Wigner standard deviations.

        In atomic units:

            sigma_Q^2 =
                coth(omega / (2 k_B T))
                -----------------------
                       2 omega

            sigma_P^2 =
                omega
                ----- coth(omega / (2 k_B T))
                  2

        At T = 0:

            sigma_Q^2 = 1 / (2 omega)

            sigma_P^2 = omega / 2
        """

        if self.temperature == 0.0:

            thermal_factor = np.ones_like(
                omega
            )

        else:

            x = (
                omega /
                (
                    2.0 *
                    self._kb_au *
                    self.temperature
                )
            )

            thermal_factor = (
                1.0 /
                np.tanh(x)
            )

        sigma_q_squared = (
            thermal_factor /
            (2.0 * omega)
        )

        sigma_p_squared = (
            0.5 *
            omega *
            thermal_factor
        )

        sigma_q = np.sqrt(
            sigma_q_squared
        )

        sigma_p = np.sqrt(
            sigma_p_squared
        )

        return sigma_q, sigma_p

    # ============================================================
    # Phase-space sampling
    # ============================================================

    def _sample_phase_space(
        self,
        sigma_q,
        sigma_p,
    ):
        """
        Draws independent normal coordinates and momenta from the
        harmonic Wigner distribution.
        """

        nmodes = sigma_q.size

        q = self._rng.normal(
            loc=0.0,
            scale=sigma_q,
            size=(
                self.n_samples,
                nmodes,
            ),
        )

        p = self._rng.normal(
            loc=0.0,
            scale=sigma_p,
            size=(
                self.n_samples,
                nmodes,
            ),
        )

        return q, p

    # ============================================================
    # Cartesian transformation
    # ============================================================

    def _transform_to_cartesian(
        self,
        equilibrium_coordinates,
        normal_modes,
        reduced_masses,
        q,
        p,
    ):
        """
        Transforms sampled normal coordinates and momenta to
        Cartesian coordinates and velocities.

        VeloxChem provides normalized Cartesian normal-mode vectors
        e_{i,a} and reduced masses mu_i.

        For mode i:

            delta r_i =
                e_i * Q_i / sqrt(mu_i)

        and

            v_i =
                e_i * P_i / sqrt(mu_i)

        when Q and P are expressed in mass-weighted atomic units.

        The reduced masses are therefore converted from amu to
        electron masses before the transformation.
        """

        nsamples = q.shape[0]

        natoms = equilibrium_coordinates.shape[0]

        nmodes = normal_modes.shape[0]

        # --------------------------------------------------------
        # Reduced masses -> atomic units
        # --------------------------------------------------------

        reduced_masses_au = (
            reduced_masses *
            self._amu_to_me
        )

        inverse_sqrt_mu = (
            1.0 /
            np.sqrt(reduced_masses_au)
        )

        # --------------------------------------------------------
        # Scale each Cartesian normal mode by 1/sqrt(mu_i)
        #
        # normal_modes:
        #     (nmodes, natoms, 3)
        #
        # inverse_sqrt_mu:
        #     (nmodes,)
        # --------------------------------------------------------

        cartesian_modes = (
            normal_modes *
            inverse_sqrt_mu[:, None, None]
        )

        # --------------------------------------------------------
        # Transform Q -> Cartesian displacement
        #
        # q:
        #     (nsamples, nmodes)
        #
        # cartesian_modes:
        #     (nmodes, natoms, 3)
        #
        # result:
        #     (nsamples, natoms, 3)
        # --------------------------------------------------------

        displacement = np.einsum(
            "si,iaj->saj",
            q,
            cartesian_modes,
            optimize=True,
        )

        # --------------------------------------------------------
        # Transform P -> Cartesian velocity
        # --------------------------------------------------------

        velocities = np.einsum(
            "si,iaj->saj",
            p,
            cartesian_modes,
            optimize=True,
        )

        # --------------------------------------------------------
        # Add equilibrium geometry
        # --------------------------------------------------------

        coordinates = (
            equilibrium_coordinates[None, :, :] +
            displacement
        )

        assert coordinates.shape == (
            nsamples,
            natoms,
            3,
        )

        assert velocities.shape == (
            nsamples,
            natoms,
            3,
        )

        return (
            coordinates,
            velocities,
        )

    # ============================================================
    # Statistics
    # ============================================================

    def _compute_statistics(
        self,
        q,
        p,
        omega,
    ):
        """
        Computes ensemble statistics useful for validating the
        generated Wigner distribution.
        """

        mean_q = np.mean(
            q,
            axis=0,
        )

        mean_p = np.mean(
            p,
            axis=0,
        )

        mean_q2 = np.mean(
            q * q,
            axis=0,
        )

        mean_p2 = np.mean(
            p * p,
            axis=0,
        )

        # Harmonic oscillator energy for each sample and mode:
        #
        # E = P^2/2 + omega^2 Q^2/2

        energies = (
            0.5 * p * p +
            0.5 *
            omega[None, :] ** 2 *
            q * q
        )

        mean_energies = np.mean(
            energies,
            axis=0,
        )

        return {
            "mean_Q":
                mean_q,

            "mean_P":
                mean_p,

            "mean_Q2":
                mean_q2,

            "mean_P2":
                mean_p2,

            "mean_mode_energy_hartree":
                mean_energies,

            "total_mean_vibrational_energy_hartree":
                float(
                    np.sum(mean_energies)
                ),
        }

    def validate_cartesian_transformation(
        self,
        vib_results,
    ):

        frequencies = np.asarray(
            vib_results["vib_frequencies"],
            dtype=float,
        )

        normal_modes = np.asarray(
            vib_results["normal_modes"],
            dtype=float,
        )

        reduced_masses = np.asarray(
            vib_results["reduced_masses"],
            dtype=float,
        )

        hessian = np.asarray(
            vib_results["hessian"],
            dtype=float,
        )

        omega = (
            frequencies /
            self._hartree_to_wavenumber
        )

        reduced_masses_au = (
            reduced_masses *
            self._amu_to_me
        )

        print("")
        print(
            "Normal-mode Cartesian transformation validation"
        )
        print(
            "================================================"
        )

        for i in range(len(frequencies)):

            # Use Q = 1 in mass-weighted atomic units

            q = 1.0

            displacement = (
                normal_modes[i].reshape(-1)
                *
                q
                /
                np.sqrt(
                    reduced_masses_au[i]
                )
            )

            energy_cartesian = (
                0.5
                *
                displacement
                @ hessian
                @ displacement
            )

            energy_normal = (
                0.5
                *
                omega[i]**2
                *
                q**2
            )

            ratio = (
                energy_cartesian /
                energy_normal
            )

            print(
                f"{i + 1:4d} "
                f"{frequencies[i]:12.4f} "
                f"{energy_cartesian:16.8e} "
                f"{energy_normal:16.8e} "
                f"{ratio:12.8f}"
            )

    # ============================================================
    # Output
    # ============================================================

    def _print_summary(
        self,
        frequencies,
        reduced_masses,
    ):
        """
        Prints a summary of the Wigner sampling calculation.
        """

        lines = []

        lines.append("")
        lines.append(
            "               Wigner Sampling Driver"
        )
        lines.append(
            "              ========================"
        )
        lines.append("")

        lines.append(
            f"Temperature             : "
            f"{self.temperature:.2f} K"
        )

        lines.append(
            f"Number of samples       : "
            f"{self.n_samples}"
        )

        lines.append(
            f"Number of sampled modes : "
            f"{frequencies.size}"
        )

        lines.append(
            f"Random seed             : "
            f"{self.seed}"
        )

        lines.append(
            f"Minimum frequency       : "
            f"{np.min(frequencies):.2f} cm^-1"
        )

        lines.append(
            f"Maximum frequency       : "
            f"{np.max(frequencies):.2f} cm^-1"
        )

        lines.append("")

        lines.append(
            " Mode       Frequency / cm^-1      "
            "Reduced mass / amu"
        )

        lines.append(
            " -----------------------------------------------"
        )

        for index, (
            frequency,
            reduced_mass,
        ) in enumerate(
            zip(
                frequencies,
                reduced_masses,
            ),
            start=1,
        ):

            lines.append(
                f" {index:4d}"
                f" {frequency:20.4f}"
                f" {reduced_mass:20.6f}"
            )

        lines.append("")

        if self.ostream is None:

            print(
                "\n".join(lines)
            )

        else:

            for line in lines:

                try:
                    self.ostream.print_info(
                        line
                    )

                except AttributeError:

                    print(line)