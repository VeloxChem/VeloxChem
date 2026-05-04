import importlib.util
from pathlib import Path
from types import SimpleNamespace
import sys

import numpy as np
import pytest

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver

ANALYZER_PATH = (Path(__file__).resolve().parents[1] / 'src' / 'pymodule' /
                 'orbitalanalyzerdriver.py')
ANALYZER_SPEC = importlib.util.spec_from_file_location(
    'veloxchem.orbitalanalyzerdriver', ANALYZER_PATH)
ANALYZER_MODULE = importlib.util.module_from_spec(ANALYZER_SPEC)
sys.modules['veloxchem.orbitalanalyzerdriver'] = ANALYZER_MODULE
ANALYZER_SPEC.loader.exec_module(ANALYZER_MODULE)
_build_orbital_candidates = ANALYZER_MODULE._build_orbital_candidates

NBODRIVER_PATH = (Path(__file__).resolve().parents[1] / 'src' / 'pymodule' /
                  'nbodriver.py')
NBODRIVER_SPEC = importlib.util.spec_from_file_location(
    'veloxchem.nbodriver', NBODRIVER_PATH)
NBODRIVER_MODULE = importlib.util.module_from_spec(NBODRIVER_SPEC)
sys.modules['veloxchem.nbodriver'] = NBODRIVER_MODULE
NBODRIVER_SPEC.loader.exec_module(NBODRIVER_MODULE)

VBDRIVER_PATH = (Path(__file__).resolve().parents[1] / 'src' / 'pymodule' /
                 'vbdriver.py')
VBDRIVER_SPEC = importlib.util.spec_from_file_location(
    'veloxchem.vbdriver', VBDRIVER_PATH)
VBDRIVER_MODULE = importlib.util.module_from_spec(VBDRIVER_SPEC)
sys.modules['veloxchem.vbdriver'] = VBDRIVER_MODULE
VBDRIVER_SPEC.loader.exec_module(VBDRIVER_MODULE)


def _pd_ligand_molecule(ligand):
    if ligand == 'N':
        xyz = """
Pd  0.0000  0.0000  0.0000
N   2.0500  0.0000  0.0000
H   2.4500  0.9300  0.0000
H   2.4500 -0.4650  0.8050
H   2.4500 -0.4650 -0.8050
"""
    elif ligand == 'P':
        xyz = """
Pd  0.0000  0.0000  0.0000
P   2.2500  0.0000  0.0000
H   2.8500  1.1500  0.0000
H   2.8500 -0.5750  0.9960
H   2.8500 -0.5750 -0.9960
"""
    else:
        raise ValueError(ligand)
    molecule = Molecule.read_str(xyz)
    molecule.set_charge(0)
    molecule.set_multiplicity(1)
    return molecule


def _pd_ligand_molecule_at_distance(ligand, distance):
    if ligand == 'N':
        hx = distance + 0.4000
        xyz = f"""
Pd  0.0000  0.0000  0.0000
N   {distance:.4f}  0.0000  0.0000
H   {hx:.4f}  0.9300  0.0000
H   {hx:.4f} -0.4650  0.8050
H   {hx:.4f} -0.4650 -0.8050
"""
    elif ligand == 'P':
        hx = distance + 0.6000
        xyz = f"""
Pd  0.0000  0.0000  0.0000
P   {distance:.4f}  0.0000  0.0000
H   {hx:.4f}  1.1500  0.0000
H   {hx:.4f} -0.5750  0.9960
H   {hx:.4f} -0.5750 -0.9960
"""
    else:
        raise ValueError(ligand)
    molecule = Molecule.read_str(xyz)
    molecule.set_charge(0)
    molecule.set_multiplicity(1)
    return molecule


def _mock_nao_data(ligand, pi_coupling):
    # Atom 0 is Pd, atom 1 is the ligand donor. Hydrogens are omitted from
    # this minimal NAO payload because the metal-ligand candidate builder only
    # requires atom maps for the centers that carry diagnostic orbitals.
    atom_map = np.array([0, 0, 0, 0, 1, 1, 1, 1], dtype=int)
    angular_map = np.array([2, 2, 2, 0, 0, 1, 1, 2], dtype=int)
    populations = np.array([1.75, 1.70, 1.65, 0.25, 2.00, 1.90, 0.20, 0.20])
    density = np.diag(populations)

    # Ligand sigma donation into a low-occupation Pd acceptor.
    density[5, 3] = density[3, 5] = 0.12
    # Metal-to-ligand pi back-donation. This is deliberately larger for the
    # phosphine mock payload, reflecting a stronger ligand pi-acceptor channel.
    density[0, 7] = density[7, 0] = pi_coupling
    density[1, 7] = density[7, 1] = 0.5 * pi_coupling

    return SimpleNamespace(
        populations=populations,
        density=density,
        atom_map=atom_map,
        angular_momentum_map=angular_map,
    )


def _metal_ligand_candidates(ligand, pi_coupling):
    molecule = _pd_ligand_molecule(ligand)
    nao_data = _mock_nao_data(ligand, pi_coupling)
    candidates = _build_orbital_candidates(
        molecule,
        nao_data,
        include_metal_ligand_candidates=True,
    )
    return [candidate for candidate in candidates if candidate.get('type') == 'ML']


def _mock_nonbonding_pi_payload():
    atom_map = np.array([0] * 20 + [1] * 4, dtype=int)
    angular_map = np.array([2] * 20 + [0, 1, 1, 2], dtype=int)
    populations = np.array(
        [1.80] * 18 + [0.05, 0.03] + [1.95, 1.45, 1.20, 0.03])
    density = np.diag(populations)
    density[0, 22] = density[22, 0] = 0.18
    density[1, 21] = density[21, 1] = 0.09
    density[21, 19] = density[19, 21] = 0.04
    return SimpleNamespace(
        populations=populations,
        density=density,
        atom_map=atom_map,
        angular_momentum_map=angular_map,
    )


def test_metal_ligand_sigma_and_pi_channels_are_separated():
    candidates = _metal_ligand_candidates('P', pi_coupling=0.20)
    subtypes = {candidate['subtype'] for candidate in candidates}

    assert 'sigma-acceptor' in subtypes
    assert 'pi-donor' in subtypes
    assert all(candidate['metal_atom'] == 0 for candidate in candidates)
    assert all(candidate['ligand_atom'] == 1 for candidate in candidates)

    sigma = next(candidate for candidate in candidates
                 if candidate['subtype'] == 'sigma-acceptor')
    pi = next(candidate for candidate in candidates
              if candidate['subtype'] == 'pi-donor')

    assert sigma['channel'] == 'ligand-to-metal-sigma-donation'
    assert sigma['donor_atom'] == 1
    assert sigma['acceptor_atom'] == 0
    assert sigma['donation_strength'] > 0.0

    assert pi['channel'] == 'metal-to-ligand-pi-back-donation'
    assert pi['donor_atom'] == 0
    assert pi['acceptor_atom'] == 1
    assert pi['back_donation_strength'] > 0.0


def test_metal_ligand_channels_remain_available_for_dissociation_scan():
    molecule = _pd_ligand_molecule_at_distance('N', 5.0)
    candidates = _build_orbital_candidates(
        molecule,
        _mock_nao_data('N', pi_coupling=0.04),
        include_metal_ligand_candidates=True,
    )
    subtypes = {
        candidate['subtype']
        for candidate in candidates
        if candidate.get('type') == 'ML'
    }

    assert 'sigma-acceptor' in subtypes
    assert 'pi-donor' in subtypes


def test_phosphine_mock_has_larger_pi_back_donation_than_ammine_mock():
    nh3_candidates = _metal_ligand_candidates('N', pi_coupling=0.04)
    ph3_candidates = _metal_ligand_candidates('P', pi_coupling=0.20)

    nh3_pi = next(candidate for candidate in nh3_candidates
                  if candidate['subtype'] == 'pi-donor')
    ph3_pi = next(candidate for candidate in ph3_candidates
                  if candidate['subtype'] == 'pi-donor')

    assert ph3_pi['back_donation_strength'] > nh3_pi['back_donation_strength']


def test_back_donation_uses_metal_d_and_ligand_nonbonding_space():
    candidates = _build_orbital_candidates(
        _pd_ligand_molecule('P'),
        _mock_nonbonding_pi_payload(),
        include_metal_ligand_candidates=True,
    )
    pi = next(candidate for candidate in candidates
              if candidate.get('type') == 'ML' and
              candidate.get('subtype') == 'pi-donor')

    assert pi['donor_atom'] == 0
    assert pi['acceptor_atom'] == 1
    assert pi['back_donation_strength'] > 0.0
    assert 'nonbonding' in pi['donor_acceptor_role']


def test_vb_ao_hamiltonian_reconstructs_ecp_hf_energy():
    molecule = _pd_ligand_molecule_at_distance('N', 1.90)
    basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)

    scf = ScfRestrictedDriver()
    scf.ostream.mute()
    scf.xcfun = 'hf'
    scf_results = scf.compute(molecule, basis)

    analysis = ANALYZER_MODULE.OrbitalAnalyzer(
        molecule,
        basis,
        mol_orbs=scf.mol_orbs,
        options=ANALYZER_MODULE.OrbitalAnalyzerOptions(
            include_mo_analysis=False),
    ).run()
    h_ao, s_ao, eri_ao, e_nuc = VBDRIVER_MODULE.VbDriver()._ao_integrals(
        molecule,
        basis,
    )
    density = analysis.density
    j_matrix = np.einsum('lm,uvlm->uv', density, eri_ao, optimize=True)
    k_matrix = np.einsum('lm,ulvm->uv', density, eri_ao, optimize=True)
    g_matrix = j_matrix - 0.5 * k_matrix
    reconstructed_energy = (
        float(np.einsum('uv,uv->', density, h_ao, optimize=True)) +
        0.5 * float(np.einsum('uv,uv->', density, g_matrix, optimize=True)) +
        float(e_nuc)
    )

    assert reconstructed_energy == pytest.approx(
        float(scf_results['scf_energy']), abs=1.0e-8)


def test_nbo_keeps_metal_ligand_records_candidate_only():
    molecule = _pd_ligand_molecule('P')
    candidates = _build_orbital_candidates(
        molecule,
        _mock_nao_data('P', pi_coupling=0.20),
        include_metal_ligand_candidates=True,
    )
    primary = NBODRIVER_MODULE._build_primary_assignment(
        molecule,
        candidates,
        NBODRIVER_MODULE.NboConstraints(),
    )

    assert any(candidate.get('type') == 'ML' for candidate in candidates)
    assert primary['counts'].get('ML', 0) == 0
    assert all(candidate.get('type') != 'ML'
               for candidate in primary['nbo_list'])


def test_nbo_result_payload_exposes_metal_ligand_diagnostics(monkeypatch):
    molecule = _pd_ligand_molecule('P')
    nao_data = _mock_nao_data('P', pi_coupling=0.20)
    candidates = _build_orbital_candidates(
        molecule,
        nao_data,
        include_metal_ligand_candidates=True,
    )
    nao_data.transform = np.eye(len(nao_data.populations))
    nao_data.overlap = np.eye(len(nao_data.populations))
    spin_data = SimpleNamespace(
        alpha_density=0.5 * nao_data.density,
        beta_density=0.5 * nao_data.density,
        spin_density=np.zeros_like(nao_data.density),
        spin_populations=np.zeros_like(nao_data.populations),
        unpaired_electrons=0.0,
    )
    analysis = SimpleNamespace(
        overlap=np.eye(len(nao_data.populations)),
        density=nao_data.density,
        nao_data=nao_data,
        spin_data=spin_data,
        mo_analysis={},
        orbital_candidates=candidates,
    )

    class FakeOrbitalAnalyzer:
        def __init__(self, *args, **kwargs):
            pass

        def run(self):
            return analysis

    monkeypatch.setattr(NBODRIVER_MODULE, 'OrbitalAnalyzer',
                        FakeOrbitalAnalyzer)
    results = NBODRIVER_MODULE.NboDriver().compute(
        molecule,
        basis=SimpleNamespace(),
        mol_orbs=SimpleNamespace(),
        options=NBODRIVER_MODULE.NboComputeOptions(
            include_diagnostics=False,
            include_mo_analysis=False,
            max_alternatives=0,
        ),
    )

    diagnostics = results['metal_ligand_diagnostics']
    subtypes = {record['subtype'] for record in diagnostics}

    assert subtypes == {'sigma-acceptor', 'pi-donor'}
    assert all(record['metal_atom'] == 1 for record in diagnostics)
    assert all(record['ligand_atom'] == 2 for record in diagnostics)
    assert all(candidate.get('type') != 'ML'
               for candidate in results['nbo_list'])


def test_vb_partitions_metal_ligand_records_with_metadata():
    candidates = _build_orbital_candidates(
        _pd_ligand_molecule('P'),
        _mock_nao_data('P', pi_coupling=0.20),
        include_metal_ligand_candidates=True,
    )
    partitions = VBDRIVER_MODULE.VbDriver()._candidate_partition_diagnostics(
        candidates)
    metal_ligand = partitions['metal_ligand']
    subtypes = {record['subtype'] for record in metal_ligand}

    assert subtypes == {'sigma-acceptor', 'pi-donor'}
    assert all(record['metal_atom'] == 0 for record in metal_ligand)
    assert all(record['ligand_atom'] == 1 for record in metal_ligand)
    assert any(record['channel'] == 'ligand-to-metal-sigma-donation'
               for record in metal_ligand)
    assert any(record['channel'] == 'metal-to-ligand-pi-back-donation'
               for record in metal_ligand)


def test_vb_can_select_metal_ligand_backdonation_candidate():
    candidates = _build_orbital_candidates(
        _pd_ligand_molecule('P'),
        _mock_nao_data('P', pi_coupling=0.20),
        include_metal_ligand_candidates=True,
    )
    driver = VBDRIVER_MODULE.VbDriver()
    selected = driver._select_active_bond_candidate(
        candidates,
        VBDRIVER_MODULE.VbComputeOptions(
            active_candidate_subtype='pi-donor'),
    )

    assert selected is not None
    assert selected['type'] == 'ML'
    assert selected['subtype'] == 'pi-donor'


def test_vb_builds_combined_metal_ligand_active_space(monkeypatch):
    molecule = _pd_ligand_molecule('P')
    nao_data = _mock_nao_data('P', pi_coupling=0.20)
    candidates = _build_orbital_candidates(
        molecule,
        nao_data,
        include_metal_ligand_candidates=True,
    )
    nao_data.transform = np.eye(len(nao_data.populations))
    analysis = SimpleNamespace(
        nao_data=nao_data,
        orbital_candidates=candidates,
    )

    monkeypatch.setattr(
        VBDRIVER_MODULE.OrbitalAnalyzer,
        'ao_shell_map',
        staticmethod(lambda _molecule, _basis: (
            np.arange(len(nao_data.populations)),
            np.zeros(len(nao_data.populations), dtype=int),
        )),
    )
    monkeypatch.setattr(
        VBDRIVER_MODULE.VbDriver,
        '_s_normalize',
        lambda self, vector, overlap: vector / np.linalg.norm(vector),
    )
    monkeypatch.setitem(sys.modules, 'veloxchem', SimpleNamespace(
        OverlapDriver=lambda: SimpleNamespace(
            compute=lambda _molecule, _basis: SimpleNamespace(
                to_numpy=lambda: np.eye(len(nao_data.populations))))))

    active_space = VBDRIVER_MODULE.VbDriver()._build_metal_ligand_active_space(
        molecule,
        basis=SimpleNamespace(get_dimensions_of_basis=lambda: len(
            nao_data.populations)),
        candidates=candidates,
        nao_data=nao_data,
        options=VBDRIVER_MODULE.VbComputeOptions(
            active_metal_ligand_channels=('sigma-acceptor', 'pi-donor')),
    )

    assert active_space.metadata['metal_ligand_channels'] == (
        'sigma-acceptor', 'pi-donor')
    assert active_space.metadata['electron_count'] == 4
    assert len(active_space.active_orbitals) == 4
    assert len(active_space.structures) == 36


def test_vb_metal_ligand_active_space_marks_backdonation_switch(monkeypatch):
    molecule = _pd_ligand_molecule('P')
    nao_data = _mock_nao_data('P', pi_coupling=0.20)
    candidates = _build_orbital_candidates(
        molecule,
        nao_data,
        include_metal_ligand_candidates=True,
    )
    nao_data.transform = np.eye(len(nao_data.populations))

    monkeypatch.setattr(
        VBDRIVER_MODULE.OrbitalAnalyzer,
        'ao_shell_map',
        staticmethod(lambda _molecule, _basis: (
            np.arange(len(nao_data.populations)),
            np.zeros(len(nao_data.populations), dtype=int),
        )),
    )
    monkeypatch.setattr(
        VBDRIVER_MODULE.VbDriver,
        '_s_normalize',
        lambda self, vector, overlap: vector / np.linalg.norm(vector),
    )
    monkeypatch.setitem(sys.modules, 'veloxchem', SimpleNamespace(
        OverlapDriver=lambda: SimpleNamespace(
            compute=lambda _molecule, _basis: SimpleNamespace(
                to_numpy=lambda: np.eye(len(nao_data.populations))))))

    driver = VBDRIVER_MODULE.VbDriver()
    basis = SimpleNamespace(get_dimensions_of_basis=lambda: len(
        nao_data.populations))
    sigma_only = driver._build_metal_ligand_active_space(
        molecule,
        basis=basis,
        candidates=candidates,
        nao_data=nao_data,
        options=VBDRIVER_MODULE.VbComputeOptions(
            active_metal_ligand_channels=('sigma-acceptor',)),
    )
    sigma_plus_pi = driver._build_metal_ligand_active_space(
        molecule,
        basis=basis,
        candidates=candidates,
        nao_data=nao_data,
        options=VBDRIVER_MODULE.VbComputeOptions(
            active_metal_ligand_channels=('sigma-acceptor', 'pi-donor')),
    )

    assert sigma_only.metadata['metal_ligand_model'] == 'sigma-only'
    assert sigma_only.metadata['backdonation_blocked'] is True
    assert sigma_only.metadata['backdonation_enabled'] is False
    assert sigma_only.metadata['metal_ligand_channels'] == ('sigma-acceptor',)
    assert sigma_only.metadata['electron_count'] == 2
    assert len(sigma_only.active_orbitals) == 2
    assert len(sigma_only.structures) == 4
    assert [record['subtype'] for record in
            sigma_only.metadata['selected_metal_ligand_records']] == [
                'sigma-acceptor'
            ]

    assert sigma_plus_pi.metadata['metal_ligand_model'] == (
        'sigma-plus-backdonation')
    assert sigma_plus_pi.metadata['backdonation_blocked'] is False
    assert sigma_plus_pi.metadata['backdonation_enabled'] is True
    assert sigma_plus_pi.metadata['metal_ligand_channels'] == (
        'sigma-acceptor', 'pi-donor')
    assert sigma_plus_pi.metadata['electron_count'] == 4
    assert len(sigma_plus_pi.active_orbitals) == 4
    assert len(sigma_plus_pi.structures) == 36
    assert {record['subtype'] for record in
            sigma_plus_pi.metadata['selected_metal_ligand_records']} == {
                'sigma-acceptor', 'pi-donor'
            }

    diagnostics = driver._active_space_diagnostics(sigma_plus_pi)
    assert diagnostics['metal_ligand_model'] == 'sigma-plus-backdonation'
    assert diagnostics['metal_ligand_backdonation_enabled'] is True
    assert diagnostics['active_orbital_count'] == 4
    assert diagnostics['determinant_count'] == 36


def test_vb_result_payload_exposes_metal_ligand_partition(monkeypatch):
    molecule = _pd_ligand_molecule('P')
    nao_data = _mock_nao_data('P', pi_coupling=0.20)
    candidates = _build_orbital_candidates(
        molecule,
        nao_data,
        include_metal_ligand_candidates=True,
    )
    analysis = SimpleNamespace(
        nao_data=nao_data,
        orbital_candidates=candidates,
    )

    class FakeOrbitalAnalyzer:
        def __init__(self, *args, **kwargs):
            pass

        def classify_nao_data(self, _nao_data):
            return analysis

    monkeypatch.setattr(VBDRIVER_MODULE, 'OrbitalAnalyzer',
                        FakeOrbitalAnalyzer)
    results = VBDRIVER_MODULE.VbDriver().compute(
        molecule,
        basis=SimpleNamespace(get_dimensions_of_basis=lambda: len(
            nao_data.populations)),
        options=VBDRIVER_MODULE.VbComputeOptions(use_active_space=False),
        nao_data=nao_data,
    )
    metal_ligand = results['diagnostics']['candidate_partitions'][
        'metal_ligand']

    assert {record['subtype'] for record in metal_ligand} == {
        'sigma-acceptor',
        'pi-donor',
    }
    assert any(record['channel'] == 'ligand-to-metal-sigma-donation'
               for record in metal_ligand)
    assert any(record['channel'] == 'metal-to-ligand-pi-back-donation'
               for record in metal_ligand)


def _mock_metal_ligand_active_space(channel_count):
    driver = VBDRIVER_MODULE.VbDriver()
    n_orbitals = 2 * channel_count
    n_ao = 2 * n_orbitals
    orbitals = []
    for index in range(n_orbitals):
        coefficients = np.zeros(n_ao)
        coefficients[index] = 1.0
        orbitals.append(
            VBDRIVER_MODULE.VbOrbital(
                label=f'active_ml_{index + 1}',
                coefficients=coefficients,
                center=index,
                kind='active',
            )
        )
    structures = driver._determinant_active_structures(
        n_orbitals,
        channel_count,
        channel_count,
    )
    channels = ('sigma-acceptor',) if channel_count == 1 else (
        'sigma-acceptor', 'pi-donor')
    model = 'sigma-only' if channel_count == 1 else 'sigma-plus-backdonation'
    records = [
        {
            'label': f'ML_{channel}',
            'subtype': channel,
            'channel': channel,
            'metal_atom': 0,
            'ligand_atom': 1,
            'donation_strength': 0.1,
            'back_donation_strength': 0.2 if channel == 'pi-donor' else 0.0,
            'interaction_strength': 0.2,
            'occupation': 2.0,
        }
        for channel in channels
    ]
    active_space = VBDRIVER_MODULE.VbActiveSpace(
        active_bond=(0, 1),
        active_candidate_label='metal_ligand_' + '+'.join(channels),
        active_candidate={
            'type': 'ML_SYSTEM',
            'subtype': '+'.join(channels),
            'atoms': (0, 1),
        },
        active_orbitals=tuple(orbitals),
        structures=tuple(structures),
        metadata={
            'source': 'mock',
            'model': 'metal-ligand-channel-determinant-ci',
            'metal_ligand_model': model,
            'metal_ligand_channels': channels,
            'backdonation_enabled': channel_count > 1,
            'backdonation_blocked': channel_count == 1,
            'selected_metal_ligand_records': records,
            'electron_count': 2 * channel_count,
            'spin': 'singlet',
            'n_alpha': channel_count,
            'n_beta': channel_count,
            'determinant_ci': True,
            'determinant_count': len(structures),
            'active_pi_atoms': tuple(range(n_orbitals)),
        },
    )
    return active_space, orbitals, structures, n_ao


def test_metal_ligand_bovb_lowers_sigma_only_determinant_limit(monkeypatch):
    active_space, orbitals, structures, n_ao = _mock_metal_ligand_active_space(1)
    driver = VBDRIVER_MODULE.VbDriver()

    def fake_ao_integrals(self, molecule, basis):
        h_ao = np.zeros((n_ao, n_ao))
        h_ao[2, 2] = -1.0
        h_ao[3, 3] = -1.0
        return h_ao, np.eye(n_ao), np.zeros((n_ao, n_ao, n_ao, n_ao)), 0.0

    def fake_breathing(self, molecule, basis, center, active_vector,
                       occupied_vectors, overlap):
        vector = np.zeros(n_ao)
        vector[int(np.argmax(np.abs(active_vector))) + len(orbitals)] = 1.0
        return vector

    monkeypatch.setattr(VBDRIVER_MODULE.VbDriver, '_ao_integrals', fake_ao_integrals)
    monkeypatch.setattr(VBDRIVER_MODULE.VbDriver, '_center_breathing_vector', fake_breathing)

    result = driver._compute_metal_ligand_bovb(
        molecule=SimpleNamespace(),
        basis=SimpleNamespace(),
        structures=list(structures),
        orbitals=list(orbitals),
        options=VBDRIVER_MODULE.VbComputeOptions(mode='bovb'),
        active_space=active_space,
    )
    diagnostics = result['diagnostics']

    assert np.isfinite(result['energy'])
    assert diagnostics['metal_ligand_bovb_model'] == (
        'channel-local-determinant-ci-breathing-orbital-relaxation')
    assert diagnostics['metal_ligand_bovb_has_external_breathing_space'] is True
    assert result['energy'] < diagnostics['metal_ligand_bovb_initial_energy'] - 1.0e-4
    assert diagnostics['metal_ligand_bovb_energy_lowering'] > 1.0e-4
    assert diagnostics['metal_ligand_bovb_used_fixed_orbital_limit'] is False
    assert len(diagnostics['metal_ligand_bovb_breathing_records']) == 2


def test_metal_ligand_sigma_vbscf_uses_spin_adapted_active_space(monkeypatch):
    active_space, orbitals, _, n_ao = _mock_metal_ligand_active_space(1)
    driver = VBDRIVER_MODULE.VbDriver()
    options = VBDRIVER_MODULE.VbComputeOptions(
        mode='vbscf',
        optimize_orbitals=True,
    )
    sigma_space = driver._spin_adapted_metal_ligand_sigma_space(
        active_space,
        options,
    )

    def fake_ao_integrals(self, molecule, basis):
        h_ao = np.zeros((n_ao, n_ao))
        h_ao[0, 1] = h_ao[1, 0] = -0.2
        return h_ao, np.eye(n_ao), np.zeros((n_ao, n_ao, n_ao, n_ao)), 0.0

    monkeypatch.setattr(VBDRIVER_MODULE.VbDriver, '_ao_integrals', fake_ao_integrals)

    result = driver.compute_vbscf_h2(
        molecule=SimpleNamespace(),
        basis=SimpleNamespace(),
        structures=list(sigma_space.structures),
        orbitals=list(orbitals),
        reference_orbitals=None,
        options=options,
        active_space=sigma_space,
    )
    diagnostics = result['diagnostics']

    assert sigma_space.metadata['model'] == 'metal-ligand-sigma-spin-adapted-vbscf'
    assert sigma_space.metadata['spin_adapted_sigma_vbscf'] is True
    assert len(sigma_space.structures) == 3
    assert np.isfinite(result['energy'])
    assert diagnostics['retained_overlap_rank'] == 3
    assert np.isclose(np.sum(result['lowdin_weights']), 1.0, atol=1.0e-8)


def test_metal_ligand_sigma_plus_backdonation_vbscf_keeps_switch_metadata(monkeypatch):
    active_space, orbitals, structures, n_ao = _mock_metal_ligand_active_space(2)
    driver = VBDRIVER_MODULE.VbDriver()

    def fake_ao_integrals(self, molecule, basis):
        h_ao = np.zeros((n_ao, n_ao))
        h_ao[4, 4] = -1.0
        h_ao[5, 5] = -1.0
        h_ao[6, 6] = -0.5
        h_ao[7, 7] = -0.5
        return h_ao, np.eye(n_ao), np.zeros((n_ao, n_ao, n_ao, n_ao)), 0.0

    monkeypatch.setattr(VBDRIVER_MODULE.VbDriver, '_ao_integrals', fake_ao_integrals)

    result = driver._compute_metal_ligand_vbscf_determinant(
        molecule=SimpleNamespace(),
        basis=SimpleNamespace(),
        structures=list(structures),
        orbitals=list(orbitals),
        active_space=active_space,
    )
    diagnostics = result['diagnostics']

    assert diagnostics['metal_ligand_model'] == 'sigma-plus-backdonation'
    assert diagnostics['metal_ligand_backdonation_enabled'] is True
    assert diagnostics['metal_ligand_orbital_relaxation_method'] == 'vbscf'
    assert diagnostics['metal_ligand_channels'] == ('sigma-acceptor', 'pi-donor')
    assert diagnostics['metal_ligand_vbscf_energy_lowering'] == 0.0
    assert diagnostics['metal_ligand_vbscf_has_external_relaxation_space'] is False
    assert len(diagnostics['metal_ligand_vbscf_orbital_amplitudes']) == 4


def test_metal_ligand_bovb_preserves_backdonation_switch_metadata(monkeypatch):
    active_space, orbitals, structures, n_ao = _mock_metal_ligand_active_space(2)
    driver = VBDRIVER_MODULE.VbDriver()

    def fake_ao_integrals(self, molecule, basis):
        h_ao = np.zeros((n_ao, n_ao))
        h_ao[4, 4] = -1.0
        h_ao[5, 5] = -1.0
        h_ao[6, 6] = -0.5
        h_ao[7, 7] = -0.5
        return h_ao, np.eye(n_ao), np.zeros((n_ao, n_ao, n_ao, n_ao)), 0.0

    def fake_breathing(self, molecule, basis, center, active_vector,
                       occupied_vectors, overlap):
        vector = np.zeros(n_ao)
        vector[int(np.argmax(np.abs(active_vector))) + len(orbitals)] = 1.0
        return vector

    monkeypatch.setattr(VBDRIVER_MODULE.VbDriver, '_ao_integrals', fake_ao_integrals)
    monkeypatch.setattr(VBDRIVER_MODULE.VbDriver, '_center_breathing_vector', fake_breathing)

    result = driver._compute_metal_ligand_bovb(
        molecule=SimpleNamespace(),
        basis=SimpleNamespace(),
        structures=list(structures),
        orbitals=list(orbitals),
        options=VBDRIVER_MODULE.VbComputeOptions(mode='bovb'),
        active_space=active_space,
    )
    diagnostics = result['diagnostics']
    diagnostics.update(driver._active_space_diagnostics(active_space))

    assert diagnostics['metal_ligand_model'] == 'sigma-plus-backdonation'
    assert diagnostics['metal_ligand_backdonation_enabled'] is True
    assert diagnostics['metal_ligand_backdonation_blocked'] is False
    assert diagnostics['metal_ligand_channels'] == ('sigma-acceptor', 'pi-donor')
    assert diagnostics['metal_ligand_bovb_has_external_breathing_space'] is True
    assert result['energy'] <= diagnostics['metal_ligand_bovb_initial_energy'] + 1.0e-8
    assert len(diagnostics['metal_ligand_bovb_breathing_amplitudes']) == 4
