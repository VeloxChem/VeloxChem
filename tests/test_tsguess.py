import pytest
import sys

from veloxchem.tsguesser import TransitionStateGuesser
from veloxchem.optimizationdriver import OptimizationDriver
from veloxchem.molecule import Molecule


class TestTransitionStateGuesser:

    @staticmethod
    def _make_scan_dict():
        """Minimal synthetic scan dict: 3 lambda points, 2 conformers each."""
        return {
            0.0: [{
                'v': 10.0,
                'e1': 10.0,
                'e2': 30.0,
                'e_int': 10.0,
                'xyz': ''
            }, {
                'v': 15.0,
                'e1': 15.0,
                'e2': 35.0,
                'e_int': 15.0,
                'xyz': ''
            }],
            0.5: [{
                'v': 20.0,
                'e1': 12.0,
                'e2': 28.0,
                'e_int': 20.0,
                'xyz': ''
            }, {
                'v': 18.0,
                'e1': 11.0,
                'e2': 25.0,
                'e_int': 18.0,
                'xyz': ''
            }],
            1.0: [{
                'v': 12.0,
                'e1': 20.0,
                'e2': 12.0,
                'e_int': 12.0,
                'xyz': ''
            }, {
                'v': 14.0,
                'e1': 22.0,
                'e2': 14.0,
                'e_int': 14.0,
                'xyz': ''
            }]
        }

    @staticmethod
    def _make_qm_scan_dict():
        """Scan dict with qm_energy populated, including NaN entries."""
        return {
            0.0: [{
                'v': 10.0,
                'e1': 10.0,
                'e2': 30.0,
                'e_int': 0.0,
                'xyz': '',
                'qm_energy': -100.0
            }, {
                'v': 15.0,
                'e1': 15.0,
                'e2': 35.0,
                'e_int': 0.0,
                'xyz': '',
                'qm_energy': float('nan')
            }],
            0.5: [{
                'v': 20.0,
                'e1': 12.0,
                'e2': 28.0,
                'e_int': 0.0,
                'xyz': '',
                'qm_energy': -95.0
            }, {
                'v': 18.0,
                'e1': 11.0,
                'e2': 25.0,
                'e_int': 0.0,
                'xyz': '',
                'qm_energy': -98.0
            }],
            1.0: [{
                'v': 12.0,
                'e1': 20.0,
                'e2': 12.0,
                'e_int': 0.0,
                'xyz': '',
                'qm_energy': float('nan')
            }, {
                'v': 14.0,
                'e1': 22.0,
                'e2': 14.0,
                'e_int': 0.0,
                'xyz': '',
                'qm_energy': float('nan')
            }]
        }

    @pytest.mark.skipif(('openmm' not in sys.modules),
                        reason='openmm not available')
    def test_ts_guesser(self):

        rea = Molecule.read_str("""
        C              1.340494446986         1.458194010980         2.002452328500
        Cl            -0.248712060729         0.603810082277         1.884782585397
        H              1.176565798320         2.498243131144         2.310601520444
        H              1.973676790569         0.954059040738         2.742914700148
        H              1.839289207801         1.441849856114         1.025499053148
        Br             4.316943499142         3.074995880494         2.227046197863"""
                                )

        pro = Molecule.read_str("""
        C              1.644879612388        -0.508351390798        -0.616977056715
        Cl             3.923435352318        -2.563707406826        -1.580365047563
        H              2.132647869077        -0.956075482316         0.257650347399
        H              2.384228922941         0.036800182006        -1.216446463266
        H              1.178719053826        -1.293388510743        -1.224930701103
        Br             0.244603286036         0.757235977330        -0.003426182246"""
                                )

        ref_mm_ts = Molecule.read_str("""
        C              0.287401892289         0.189412663028        -0.181389105630
        Cl            -0.714922847473         1.441402502475        -1.633913908687
        H             -0.506578391806        -0.570890866041        -0.274740014262
        H              1.179217362421         0.137862529013        -0.828734078861
        H              0.201866111086         0.985209360229         0.577838564927
        Br             1.320755564351        -1.098139598962         1.310400014989"""
                                      )
        ref_qm_ts = Molecule.read_str("""
        C             -0.300250466161         0.158721109824         0.132425421981
        Cl            -1.874349934839        -1.112845926481         0.335515890867
        H              0.374806873996        -0.387013493518         0.813286305729
        H             -0.184403283782         0.023551529485        -0.956320205450
        H             -0.827026276510         1.053156785986         0.506230455730
        Br             1.523332262274         1.629693555505        -0.102017809159"""
                                      )
        # check icrmsd between reference and result
        rea.set_charge(-1)
        pro.set_charge(-1)
        ts_guesser = TransitionStateGuesser()
        ts_guesser.ostream.mute()

        ts_guesser.qm_xcfun = "HF"
        ts_guesser.qm_basis = "STO-3G"
        ts_guesser.do_qm_scan = True

        ts_guesser.save_results_file = False
        ts_guesser.mm_steps = 200
        ts_guesser.lambda_vector = [
            0, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 1.0
        ]
        results = ts_guesser.find_transition_state(rea, pro)
        mm_icrmsd = OptimizationDriver.get_ic_rmsd(
            Molecule.read_xyz_string(results['max_mm_xyz']), ref_mm_ts)
        qm_icrmsd = OptimizationDriver.get_ic_rmsd(
            Molecule.read_xyz_string(results['max_qm_xyz']), ref_qm_ts)

        # Generous values based on averages from 100 runs
        max_bond_rms = 0.1
        max_bond_max = 0.2
        max_angle_rms = 10
        max_angle_max = 20

        assert results['max_mm_lambda'] == 0.55
        assert results['max_qm_lambda'] == 0.45

        for icrmsd in [mm_icrmsd, qm_icrmsd]:
            assert icrmsd['bonds']['rms'] < max_bond_rms
            assert icrmsd['bonds']['max'] < max_bond_max
            assert icrmsd['angles']['rms'] < max_angle_rms
            assert icrmsd['angles']['max'] < max_angle_max

    def test_check_discontinuities_smooth(self):
        guesser = TransitionStateGuesser()
        E1 = [10.0, 12.0, 14.0, 16.0]  # monotonically increasing — no drop
        E2 = [16.0, 14.0, 12.0, 10.0]  # monotonically decreasing — no rise
        assert guesser._check_discontinuities(E1, E2) == set()

    def test_check_discontinuities_e1_drop(self):
        guesser = TransitionStateGuesser()
        E1 = [10.0, 12.0, 8.0, 16.0]  # drops between index 1→2
        E2 = [16.0, 14.0, 12.0, 10.0]  # smooth
        result = guesser._check_discontinuities(E1, E2)
        assert {1, 2}.issubset(result)

    def test_check_discontinuities_e2_rise(self):
        guesser = TransitionStateGuesser()
        E1 = [10.0, 12.0, 14.0, 16.0]  # smooth
        E2 = [16.0, 14.0, 18.0, 10.0]  # rises between index 1→2
        result = guesser._check_discontinuities(E1, E2)
        assert {1, 2}.issubset(result)

    def test_check_discontinuities_both(self):
        guesser = TransitionStateGuesser()
        E1 = [10.0, 14.0, 8.0, 16.0]  # drops at 1→2
        E2 = [16.0, 10.0, 15.0, 8.0]  # rises at 1→2
        assert guesser._check_discontinuities(E1, E2) == {1, 2}

    def test_check_discontinuities_returns_set_no_duplicates(self):
        guesser = TransitionStateGuesser()
        # E1 drops at 0→1 and 1→2 — index 1 is a violation from both pairs
        E1 = [15.0, 10.0, 8.0, 12.0]
        E2 = [10.0, 12.0, 11.0, 9.0]  # smooth
        result = guesser._check_discontinuities(E1, E2)
        assert isinstance(result, set)
        assert result == {0, 1, 2}

    def test_get_best_mm_E_selects_lowest_v(self):
        # lambda=0.0: lowest v=10.0 at conf 0 → e1=10.0, e2=30.0
        # lambda=0.5: lowest v=18.0 at conf 1 → e1=11.0, e2=25.0
        # lambda=1.0: lowest v=12.0 at conf 0 → e1=20.0, e2=12.0
        V, E1, E2, conf_indices = TransitionStateGuesser._get_best_mm_E_from_scan_dict(
            self._make_scan_dict())
        assert V == [10.0, 18.0, 12.0]
        assert E1 == [10.0, 11.0, 20.0]
        assert E2 == [30.0, 25.0, 12.0]
        assert conf_indices == [0, 1, 0]

    def test_get_best_mm_E_single_conformer(self):
        scan = {
            0.0: [{
                'v': 5.0,
                'e1': 1.0,
                'e2': 9.0,
                'e_int': 5.0,
                'xyz': ''
            }],
            0.5: [{
                'v': 8.0,
                'e1': 3.0,
                'e2': 7.0,
                'e_int': 8.0,
                'xyz': ''
            }],
            1.0: [{
                'v': 3.0,
                'e1': 7.0,
                'e2': 3.0,
                'e_int': 3.0,
                'xyz': ''
            }],
        }
        V, _, _, conf_indices = TransitionStateGuesser._get_best_mm_E_from_scan_dict(
            scan)
        assert V == [5.0, 8.0, 3.0]
        assert conf_indices == [0, 0, 0]

    def test_get_best_mm_E_length_matches_lambda_count(self):
        scan = {
            l: [{
                'v': float(i),
                'e1': 0.0,
                'e2': 0.0,
                'e_int': 0.0,
                'xyz': ''
            }]
            for i, l in enumerate([0.0, 0.25, 0.5, 0.75, 1.0])
        }
        V, E1, E2, conf_indices = TransitionStateGuesser._get_best_mm_E_from_scan_dict(
            scan)
        assert len(V) == len(E1) == len(E2) == len(conf_indices) == 5

    # --- _get_best_qm_E_from_scan_dict ---

    def test_get_best_qm_E_skips_nan(self):
        qm_energies, conf_indices = TransitionStateGuesser._get_best_qm_E_from_scan_dict(
            self._make_qm_scan_dict())
        # lambda=0.0: only conf 0 is valid (-100.0); conf 1 is nan
        assert qm_energies[0] == -100.0
        assert conf_indices[0] == 0

    def test_get_best_qm_E_selects_lowest_valid(self):
        qm_energies, conf_indices = TransitionStateGuesser._get_best_qm_E_from_scan_dict(
            self._make_qm_scan_dict())
        # lambda=0.5: conf 0 has -95.0, conf 1 has -98.0 → pick -98.0
        assert qm_energies[1] == -98.0
        assert conf_indices[1] == 1

    def test_get_best_qm_E_all_nan_returns_none(self):
        qm_energies, _ = TransitionStateGuesser._get_best_qm_E_from_scan_dict(
            self._make_qm_scan_dict())
        # lambda=1.0: both conformers are nan → sentinel None
        assert qm_energies[2] is None

    def test_get_best_qm_E_no_qm_energy_key(self):
        # Plain MM-only scan dict — no 'qm_energy' key present
        qm_energies, _ = TransitionStateGuesser._get_best_qm_E_from_scan_dict(
            self._make_scan_dict())
        assert all(e is None for e in qm_energies)

    # -----------------------------------------------------------------------
    # Helpers for conformational-TS tests
    # -----------------------------------------------------------------------

    @staticmethod
    def _make_conformer_pair():
        """Return (anti, gauche) Molecule objects for 1,2-difluoroethane.

        Anti:   F-C-C-F dihedral ≈  180°
        Gauche: F-C-C-F dihedral ≈   60°

        Atom order (0-indexed):
          0 C, 1 C, 2 F (on C0), 3 H, 4 H, 5 F (on C1), 6 H, 7 H
        F-C-C-F torsion: indices (2, 0, 1, 5) / 1-based (3, 1, 2, 6)
        """
        anti = Molecule.read_str("""
        C    0.000   0.000   0.762
        C    0.000   0.000  -0.762
        F    1.101   0.000   1.415
        H   -0.520   0.901   1.153
        H   -0.520  -0.901   1.153
        F   -1.101   0.000  -1.415
        H    0.520   0.901  -1.153
        H    0.520  -0.901  -1.153""")

        gauche = Molecule.read_str("""
        C    0.000   0.000   0.762
        C    0.000   0.000  -0.762
        F    1.101   0.000   1.415
        H   -0.520   0.901   1.153
        H   -0.520  -0.901   1.153
        F    0.551  -0.954  -1.415
        H   -1.040  -0.001  -1.153
        H    0.520   0.901  -1.153""")

        return anti, gauche

    # -----------------------------------------------------------------------
    # _detect_active_dihedral
    # -----------------------------------------------------------------------

    @pytest.mark.skipif(('openmm' not in sys.modules),
                        reason='openmm not available')
    def test_detect_active_dihedral_auto(self):
        """Auto-detection picks the C-C bond with |Δφ| ≈ 120°."""
        anti, gauche = self._make_conformer_pair()
        tsguesser = TransitionStateGuesser()
        tsguesser.ostream.mute()
        tsguesser.save_results_file = False
        tsguesser.build_forcefields(anti, gauche)

        active_torsion, phi_rea, phi_pro = tsguesser._detect_active_dihedral()

        # Central bond atoms (1-indexed j, k in the torsion) must be C–C (0,1)
        central = tuple(sorted([active_torsion[1], active_torsion[2]]))
        assert central == (0, 1), (
            f"Expected C-C central bond (0,1), got {central}")

        # Shortest-path angle difference should be ≈ 120°
        delta = phi_pro - phi_rea
        if delta > 180.0:
            delta -= 360.0
        elif delta <= -180.0:
            delta += 360.0
        assert abs(abs(delta) - 120.0) < 5.0, (
            f"|Δφ| = {abs(delta):.1f}°, expected ≈ 120°")

    @pytest.mark.skipif(('openmm' not in sys.modules),
                        reason='openmm not available')
    def test_detect_active_dihedral_threshold(self):
        """Raises ValueError when no dihedral changes by more than 20°."""
        anti, _ = self._make_conformer_pair()
        # Use the anti geometry for both reactant and product (Δφ = 0°)
        anti2 = Molecule.read_str("""
        C    0.000   0.000   0.762
        C    0.000   0.000  -0.762
        F    1.101   0.000   1.415
        H   -0.520   0.901   1.153
        H   -0.520  -0.901   1.153
        F   -1.101   0.000  -1.415
        H    0.520   0.901  -1.153
        H    0.520  -0.901  -1.153""")
        tsguesser = TransitionStateGuesser()
        tsguesser.ostream.mute()
        tsguesser.save_results_file = False
        tsguesser.build_forcefields(anti, anti2)

        with pytest.raises(ValueError, match='20'):
            tsguesser._detect_active_dihedral()

    # -----------------------------------------------------------------------
    # active_torsion explicit override
    # -----------------------------------------------------------------------

    @pytest.mark.skipif(('openmm' not in sys.modules),
                        reason='openmm not available')
    def test_active_torsion_explicit_override(self):
        """Explicit active_torsion bypasses auto-detection and is stored correctly."""
        anti, gauche = self._make_conformer_pair()
        tsguesser = TransitionStateGuesser()
        tsguesser.ostream.mute()
        tsguesser.save_results_file = False

        # F-C-C-F torsion: 1-indexed (3, 1, 2, 6) → 0-indexed (2, 0, 1, 5)
        tsguesser.active_torsion = (3, 1, 2, 6)

        tsguesser.build_forcefields(anti, gauche)
        tsguesser.build_systems()

        # Stored torsion should be converted to 0-indexed
        assert tsguesser._conformer_active_torsion == (2, 0, 1, 5)

        # Stored angles should match a direct measurement on the molecules
        expected_phi_rea = tsguesser.reactant.molecule.get_dihedral_in_degrees(
            [3, 1, 2, 6])
        expected_phi_pro = tsguesser.product.molecule.get_dihedral_in_degrees(
            [3, 1, 2, 6])
        assert abs(tsguesser._conformer_phi_reactant - expected_phi_rea) < 0.1
        assert abs(tsguesser._conformer_phi_product - expected_phi_pro) < 0.1

        # Reactant should be close to 180° (anti) and product to ±60° (gauche)
        assert abs(abs(tsguesser._conformer_phi_reactant) - 180.0) < 5.0
        assert abs(abs(tsguesser._conformer_phi_product) - 60.0) < 5.0

    # -----------------------------------------------------------------------
    # Implicit solvation
    # -----------------------------------------------------------------------

    @pytest.mark.skipif(('openmm' not in sys.modules),
                        reason='openmm not available')
    def test_implicit_solvent_system_has_gb_force(self):
        """Integration systems contain a CustomGBForce when implicit_solvent_model is set."""
        rea = Molecule.read_str("""
        C              1.340494446986         1.458194010980         2.002452328500
        Cl            -0.248712060729         0.603810082277         1.884782585397
        H              1.176565798320         2.498243131144         2.310601520444
        H              1.973676790569         0.954059040738         2.742914700148
        H              1.839289207801         1.441849856114         1.025499053148
        Br             4.316943499142         3.074995880494         2.227046197863""")

        pro = Molecule.read_str("""
        C              1.644879612388        -0.508351390798        -0.616977056715
        Cl             3.923435352318        -2.563707406826        -1.580365047563
        H              2.132647869077        -0.956075482316         0.257650347399
        H              2.384228922941         0.036800182006        -1.216446463266
        H              1.178719053826        -1.293388510743        -1.224930701103
        Br             0.244603286036         0.757235977330        -0.003426182246""")

        rea.set_charge(-1)
        pro.set_charge(-1)

        tsguesser = TransitionStateGuesser()
        tsguesser.ostream.mute()
        tsguesser.save_results_file = False
        tsguesser.implicit_solvent_model = 'gbn2'

        tsguesser.build_forcefields(rea, pro)
        tsguesser.build_systems()

        # Every integration system (keyed by lambda float) should contain a GB force
        sample_system = tsguesser.systems[0.5]
        force_type_names = [type(f).__name__ for f in sample_system.getForces()]
        assert 'CustomGBForce' in force_type_names, (
            f"Expected CustomGBForce in system forces, got: {force_type_names}")
