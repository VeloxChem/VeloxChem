import os
import shutil
import tempfile
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from veloxchem.outputstream import OutputStream


# ---------------------------------------------------------------------------
# EVB: shared helpers.
#
# The committed JSON force fields (ethanol -> ethene + water) let everything
# below the force-field builder run without any QM, which is what makes the
# EVB smoke and cheap-sampling tiers affordable.
# ---------------------------------------------------------------------------
class EvbTestHelper:
    """Namespace for the shared EVB test helpers.

    A class rather than bare module-level functions so that importing this
    module adds exactly one name to a test module's namespace, and pytest
    cannot mistake a helper for a test.
    """

    @staticmethod
    def evb_compare_dict(dict1, dict2, float_tol=1e-2):
        """Recursively compare two force-field-style dicts within a tolerance."""

        assert sorted(list(dict1.keys())) == sorted(list(dict2.keys()))

        for key in dict1:
            if key == 'comment':
                continue

            val1 = dict1[key]
            val2 = dict2[key]

            type1 = type(val1)
            type2 = type(val2)

            if type1 != type2:
                try:
                    val2 = type1(val2)
                except (ValueError, TypeError):
                    raise AssertionError(
                        f"Type mismatch: {type1} != {type2} for key {key}")

            if isinstance(val1, dict):
                EvbTestHelper.evb_compare_dict(val1, val2, float_tol=float_tol)
            elif isinstance(val1, (float, np.float64)):
                assert abs(val1 - val2) < float_tol
            elif isinstance(val1, (list, np.ndarray)):
                assert np.allclose(val1, val2, atol=float_tol)
            else:
                assert val1 == val2

    @staticmethod
    def _evb_round_numbers_in_line(line, decimals=6):
        line = line.split('"')
        for i, part in enumerate(line):
            try:
                num = float(part)
                if num != 0:
                    line[i] = f"{num:.{decimals}f}"
            except ValueError:
                pass
        return '"'.join(line)

    @staticmethod
    def _evb_read_reference_xml(path) -> str:
        """Read a serialized-system reference, transparently gunzipping .gz files."""
        import gzip

        path = str(path)
        if path.endswith(".gz"):
            with gzip.open(path, "rt", encoding="utf-8") as handle:
                return handle.read()
        with open(path, "r") as handle:
            return handle.read()

    @staticmethod
    def evb_compare_serialized_systems(system, path):
        """Compare an OpenMM system against a serialized reference file.

        Compares serialized strings (not the swig proxy objects). Tolerant to a
        differing number of solvent molecules (skips the tail of longer sections)
        and to openmm version / rounding differences. ``path`` may point at a plain
        ``.xml`` or a gzipped ``.xml.gz`` reference.
        """
        import openmm as mm

        sys_string = mm.XmlSerializer.serialize(system)
        ref_string = EvbTestHelper._evb_read_reference_xml(path)

        sys_lines = sys_string.splitlines()
        sys_i_offset = 0
        ref_lines = ref_string.splitlines()
        ref_i_offset = 0

        i = -1
        while (i + sys_i_offset + 1 < len(sys_lines)
               and i + ref_i_offset + 1 < len(ref_lines)):
            i += 1
            sys_line = sys_lines[i + sys_i_offset].strip()
            ref_line = ref_lines[i + ref_i_offset].strip()
            sys_close_sec_line = sys_line.startswith('</')
            ref_close_sec_line = ref_line.startswith('</')
            while sys_close_sec_line and not ref_close_sec_line:
                ref_i_offset += 1
                ref_line = ref_lines[i + ref_i_offset].strip()
                ref_close_sec_line = ref_line.startswith('</')
            while ref_close_sec_line and not sys_close_sec_line:
                sys_i_offset += 1
                sys_line = sys_lines[i + sys_i_offset].strip()
                sys_close_sec_line = sys_line.startswith('</')
            if 'openmmVersion' in sys_line:
                continue

            if (sys_line.startswith('<A x=') and ref_line.startswith('<A x=')
                    or sys_line.startswith('<B x=') and ref_line.startswith('<B x=')
                    or sys_line.startswith('<C x=')
                    and ref_line.startswith('<C x=')):
                cond = (EvbTestHelper._evb_round_numbers_in_line(sys_line,
                                                   0) == EvbTestHelper._evb_round_numbers_in_line(
                                                       ref_line, 0))
            else:
                cond = (EvbTestHelper._evb_round_numbers_in_line(sys_line) ==
                        EvbTestHelper._evb_round_numbers_in_line(ref_line))

            assert cond, f"Line mismatch on line {i}: {sys_line} != {ref_line}"

    @staticmethod
    def evb_assert_energy_finite(system, positions_angstrom, label=""):
        """Assert that the system yields a finite potential energy at the given
        positions.

        Builds a throwaway Context on the Reference platform (deterministic, no GPU),
        sets the positions (interpreted in angstrom), and checks the potential energy
        is neither NaN nor infinite. This is the cheap EVB smoke-tier correctness
        oracle: it catches malformed systems (bad force groups, mismatched particle
        counts, exploding expressions) without running any dynamics.
        """
        import openmm as mm
        import openmm.unit as unit

        integrator = mm.VerletIntegrator(0.001 * unit.picoseconds)
        platform = mm.Platform.getPlatformByName("Reference")
        context = mm.Context(system, integrator, platform)
        try:
            context.setPositions(np.asarray(positions_angstrom) * unit.angstrom)
            state = context.getState(getEnergy=True)
            energy = state.getPotentialEnergy().value_in_unit(
                unit.kilojoule_per_mole)
        finally:
            del context
            del integrator

        assert np.isfinite(energy), (
            f"Non-finite potential energy ({energy}) for system {label}")
        return energy

    @staticmethod
    def evb_data_dir():
        """The tests/data folder holding the committed EVB inputs and references."""

        return Path(__file__).parent / "data"

    @staticmethod
    def evb_ff_pair():
        """Load the committed ethanol reactant and ethene+water product force fields
        with their molecules attached (no QM).

        Freshly loaded on every call: the builders annotate and mutate what they are
        handed, so tests must not share these objects.
        """

        from veloxchem.molecule import Molecule
        from veloxchem.mmforcefieldgenerator import MMForceFieldGenerator

        data_dir = EvbTestHelper.evb_data_dir()

        reactant = MMForceFieldGenerator.load_forcefield_from_json_file(
            str(data_dir / "evb_ethanol_ff_data.json"), ostream=OutputStream(None))
        reactant.molecule = Molecule.read_xyz_file(
            str(data_dir / "evb_ethanol.xyz"))

        product = MMForceFieldGenerator.load_forcefield_from_json_file(
            str(data_dir / "evb_ethene_H2O_ff_data.json"),
            ostream=OutputStream(None))
        product.molecule = Molecule.read_xyz_file(
            str(data_dir / "evb_ethene_H2O.xyz"))

        return SimpleNamespace(reactant=reactant, product=product)

    @staticmethod
    def evb_petase_ms_states():
        """Load the committed three-state PETase force fields (no QM).

        int1 -> int2 -> int3 of the PETase mechanism, already matched into a single
        atom ordering by build_many_force_fields. Between them the five reacting
        pairs cover every case a multi-step path can produce: bonded on both sides
        of a step, bonded on neither, and bonded in exactly one state.

        Freshly loaded on every call: the builders annotate and mutate what they are
        handed, so tests must not share these objects.
        """

        from veloxchem.molecule import Molecule
        from veloxchem.mmforcefieldgenerator import MMForceFieldGenerator

        data_dir = EvbTestHelper.evb_data_dir()

        states = []
        for index in range(3):
            state = MMForceFieldGenerator.load_forcefield_from_json_file(
                str(data_dir / f"evb_petase_ms_state_{index}_ff_data.json"),
                ostream=OutputStream(None))
            state.molecule = Molecule.read_xyz_file(
                str(data_dir / f"evb_petase_ms_state_{index}.xyz"))
            states.append(state)
        return states

    @staticmethod
    @contextmanager
    def evb_chdir_tmp():
        """Run a block of a test inside a fresh temporary directory.

        build_systems / load_initialisation and the FEP driver write their folders
        and output files into the current working directory, so tests that touch
        them must be isolated here.
        """

        previous_cwd = Path.cwd()
        tmp_dir = Path(tempfile.mkdtemp(prefix="evb_test_"))
        os.chdir(tmp_dir)
        try:
            yield tmp_dir
        finally:
            os.chdir(previous_cwd)
            shutil.rmtree(tmp_dir, ignore_errors=True)
