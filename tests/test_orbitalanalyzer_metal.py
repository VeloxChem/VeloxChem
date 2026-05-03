import importlib.util
from pathlib import Path
from types import SimpleNamespace
import sys

import numpy as np

from veloxchem.molecule import Molecule

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
