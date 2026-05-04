import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest
from mpi4py import MPI

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver

NBODRIVER_PATH = (Path(__file__).resolve().parents[1] / 'src' / 'pymodule' /
                  'nbodriver.py')
ANALYZER_PATH = (Path(__file__).resolve().parents[1] / 'src' / 'pymodule' /
                 'orbitalanalyzerdriver.py')
ANALYZER_SPEC = importlib.util.spec_from_file_location(
    'veloxchem.orbitalanalyzerdriver', ANALYZER_PATH)
ANALYZER_MODULE = importlib.util.module_from_spec(ANALYZER_SPEC)
sys.modules['veloxchem.orbitalanalyzerdriver'] = ANALYZER_MODULE
ANALYZER_SPEC.loader.exec_module(ANALYZER_MODULE)
NBODRIVER_SPEC = importlib.util.spec_from_file_location('veloxchem.nbodriver',
                                                        NBODRIVER_PATH)
NBODRIVER_MODULE = importlib.util.module_from_spec(NBODRIVER_SPEC)
sys.modules['veloxchem.nbodriver'] = NBODRIVER_MODULE
NBODRIVER_SPEC.loader.exec_module(NBODRIVER_MODULE)
NboDriver = NBODRIVER_MODULE.NboDriver

MOLECULE_XYZ = {
    'water': """
O   0.000000   0.000000   0.000000
H   0.758602   0.000000   0.504284
H  -0.758602   0.000000   0.504284
""",
    'methane': """
C   0.000000   0.000000   0.000000
H   0.629118   0.629118   0.629118
H  -0.629118  -0.629118   0.629118
H  -0.629118   0.629118  -0.629118
H   0.629118  -0.629118  -0.629118
""",
    'ethylene': """
C  -0.669500   0.000000   0.000000
C   0.669500   0.000000   0.000000
H  -1.232100   0.923000   0.000000
H  -1.232100  -0.923000   0.000000
H   1.232100   0.923000   0.000000
H   1.232100  -0.923000   0.000000
""",
    'allyl_cation': """
C  -1.280000   0.000000   0.000000
C   0.000000   0.000000   0.000000
C   1.280000   0.000000   0.000000
H  -1.880000   0.920000   0.000000
H  -1.880000  -0.920000   0.000000
H   0.000000   1.080000   0.000000
H   1.880000   0.920000   0.000000
H   1.880000  -0.920000   0.000000
""",
    'allyl_anion': """
C  -1.280000   0.000000   0.000000
C   0.000000   0.000000   0.000000
C   1.280000   0.000000   0.000000
H  -1.880000   0.920000   0.000000
H  -1.880000  -0.920000   0.000000
H   0.000000   1.080000   0.000000
H   1.880000   0.920000   0.000000
H   1.880000  -0.920000   0.000000
""",
    'allyl_radical': """
C  -1.280000   0.000000   0.000000
C   0.000000   0.000000   0.000000
C   1.280000   0.000000   0.000000
H  -1.880000   0.920000   0.000000
H  -1.880000  -0.920000   0.000000
H   0.000000   1.080000   0.000000
H   1.880000   0.920000   0.000000
H   1.880000  -0.920000   0.000000
""",
    'benzene': """
C   1.396792   0.000000   0.000000
C   0.698396   1.209951   0.000000
C  -0.698396   1.209951   0.000000
C  -1.396792   0.000000   0.000000
C  -0.698396  -1.209951   0.000000
C   0.698396  -1.209951   0.000000
H   2.490290   0.000000   0.000000
H   1.245145   2.156660   0.000000
H  -1.245145   2.156660   0.000000
H  -2.490290   0.000000   0.000000
H  -1.245145  -2.156660   0.000000
H   1.245145  -2.156660   0.000000
""",
    'formate': """
C   0.000000   0.000000   0.000000
O   1.250000   0.000000   0.000000
O  -1.250000   0.000000   0.000000
H   0.000000   1.100000   0.000000
""",
    'ozone': """
O  -1.080000   0.000000   0.000000
O   0.000000   0.000000   0.000000
O   1.080000   0.000000   0.000000
""",
    'nitrate': """
N   0.000000   0.000000   0.000000
O   1.240000   0.000000   0.000000
O  -0.620000   1.073872   0.000000
O  -0.620000  -1.073872   0.000000
""",
}

MOLECULE_CHARGES = {
    'allyl_cation': 1,
    'allyl_anion': -1,
    'formate': -1,
    'nitrate': -1,
}

MOLECULE_MULTIPLICITIES = {
    'allyl_radical': 2,
}


def _run_nbo_with_driver(name,
                         constraints=None,
                         options=None,
                         basis_label='sto-3g'):

    molecule = Molecule.read_str(MOLECULE_XYZ[name])
    if name in MOLECULE_CHARGES:
        molecule.set_charge(MOLECULE_CHARGES[name])
    if name in MOLECULE_MULTIPLICITIES:
        molecule.set_multiplicity(MOLECULE_MULTIPLICITIES[name])
    basis = MolecularBasis.read(molecule, basis_label, ostream=None)

    if molecule.get_multiplicity() == 1:
        scf_drv = ScfRestrictedDriver()
    else:
        scf_drv = ScfUnrestrictedDriver()
    scf_drv.ostream.mute()
    scf_drv.xcfun = 'hf'
    scf_drv.compute(molecule, basis)

    nbo_drv = NboDriver()
    nbo_drv.verbose = False
    nbo_drv.ostream.mute()
    results = nbo_drv.compute(molecule,
                              basis,
                              scf_drv.mol_orbs,
                              constraints=constraints,
                              options=options)

    return molecule, results, nbo_drv


def _run_nbo(name, constraints=None, options=None, basis_label='sto-3g'):

    molecule, results, _ = _run_nbo_with_driver(name,
                                               constraints=constraints,
                                               options=options,
                                               basis_label=basis_label)
    return molecule, results


def _mock_pd_ligand_molecule(ligand='NH3'):

    if ligand == 'NH3':
        xyz = """
Pd  0.0000  0.0000  0.0000
N   2.0500  0.0000  0.0000
H   2.4500  0.9300  0.0000
H   2.4500 -0.4650  0.8050
H   2.4500 -0.4650 -0.8050
"""
    elif ligand == 'PH3':
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


def _mock_metal_ligand_analysis(ligand='NH3', sigma_coupling=0.12,
                                pi_coupling=0.04):

    molecule = _mock_pd_ligand_molecule(ligand)
    atom_map = np.array([0, 0, 0, 0, 1, 1, 1, 1], dtype=int)
    angular_map = np.array([2, 2, 2, 0, 0, 1, 1, 2], dtype=int)
    populations = np.array([1.75, 1.70, 1.65, 0.25,
                            2.00, 1.90, 0.20, 0.20])
    density = np.diag(populations)
    density[5, 3] = density[3, 5] = sigma_coupling
    density[0, 7] = density[7, 0] = pi_coupling
    density[1, 7] = density[7, 1] = 0.5 * pi_coupling
    density[2, 6] = density[6, 2] = 0.25 * pi_coupling

    nao_data = SimpleNamespace(
        populations=populations,
        density=density,
        atom_map=atom_map,
        angular_momentum_map=angular_map,
        transform=np.eye(len(populations)),
        overlap=np.eye(len(populations)),
        equivalent_atom_groups=[],
        local_frames=[],
    )
    spin_data = SimpleNamespace(
        alpha_density=0.5 * density,
        beta_density=0.5 * density,
        spin_density=np.zeros_like(density),
        spin_populations=np.zeros_like(populations),
        unpaired_electrons=0.0,
    )
    candidates = ANALYZER_MODULE._build_orbital_candidates(
        molecule,
        nao_data,
        spin_data=None,
        include_metal_ligand_candidates=True,
    )
    analysis = SimpleNamespace(
        overlap=np.eye(len(populations)),
        density=density,
        nao_data=nao_data,
        spin_data=spin_data,
        mo_analysis={'max_normalization_error': 0.0},
        orbital_candidates=candidates,
    )
    return molecule, analysis


def _run_mock_metal_ligand_nbo(monkeypatch, ligand='NH3',
                               sigma_coupling=0.12, pi_coupling=0.04):

    molecule, analysis = _mock_metal_ligand_analysis(
        ligand,
        sigma_coupling=sigma_coupling,
        pi_coupling=pi_coupling,
    )

    class MockOrbitalAnalyzer:

        def __init__(self, *args, **kwargs):
            pass

        def run(self):
            return analysis

    monkeypatch.setattr(NBODRIVER_MODULE, 'OrbitalAnalyzer',
                        MockOrbitalAnalyzer)
    driver = NboDriver()
    driver.verbose = False
    driver.ostream.mute()
    results = driver.compute(
        molecule,
        basis=SimpleNamespace(),
        mol_orbs=SimpleNamespace(),
        options=NBODRIVER_MODULE.NboComputeOptions(
            include_diagnostics=True,
            include_mo_analysis=False,
            max_alternatives=0,
        ),
    )
    return molecule, results, driver


def _nbo_type_counts(results):

    counts = {'CR': 0, 'BD': 0, 'LP': 0, 'BD*': 0, 'RY': 0}
    for item in results['nbo_list']:
        item_type = item['type']
        counts[item_type] = counts.get(item_type, 0) + 1
    return counts


def _selected_bond_pairs(results, subtype=None):

    pairs = set()
    for item in results['nbo_list']:
        if item.get('type') != 'BD' or len(item.get('atoms', ())) != 2:
            continue
        if subtype is not None and item.get('subtype') != subtype:
            continue
        pairs.add(tuple(atom + 1 for atom in sorted(item['atoms'])))
    return pairs


class _RecordingOutputStream:

    def __init__(self):
        self.lines = []
        self.flush_count = 0

    def print_header(self, line):
        self.lines.append(line)

    def print_line(self, line):
        self.lines.append(line)

    def print_blank(self):
        self.lines.append('')

    def flush(self):
        self.flush_count += 1


def _alternative_pi_sets(results):

    return {
        tuple(sorted(tuple(sorted(pair))
                     for pair in alt.get('pi_bonds', [])))
        for alt in results.get('alternatives', [])
    }


def _candidate_indices(candidates):

    return {candidate['index'] for candidate in candidates}


def _candidate_electron_count(candidate):

    return float(candidate.get('electron_count', 2.0))


def _candidate_vector(candidate, size):

    vector = np.zeros(size)
    for item in candidate.get('coefficients', []):
        vector[int(item['nao_index']) - 1] = float(item['coefficient'])
    return vector


def _assert_lewis_accounting(molecule, assignment):

    atom_electrons = np.array(assignment['atom_electron_count'])
    atom_valence = np.array(assignment['atom_valence_electron_count'])
    formal_charges = np.array(assignment['formal_charges'])
    octet_excess = np.array(assignment['octet_excess'])
    selected_electrons = assignment['selected_lewis_electron_count']

    assert len(atom_electrons) == molecule.number_of_atoms()
    assert len(atom_valence) == molecule.number_of_atoms()
    assert len(formal_charges) == molecule.number_of_atoms()
    assert len(octet_excess) == molecule.number_of_atoms()
    assert abs(float(np.sum(atom_electrons)) - selected_electrons) < 1.0e-12
    assert abs(selected_electrons - sum(
        _candidate_electron_count(candidate)
        for candidate in assignment['nbo_list'])) < 1.0e-12
    assert assignment['electron_ownership_error'] < 1.0e-12
    assert abs(float(np.sum(formal_charges)) - molecule.get_charge()) < 1.0e-12
    assert np.all(octet_excess < 1.0e-12)
    assert assignment['valence_warnings'] == []


def _assert_score_terms(assignment):

    terms = assignment['score_terms']
    expected_keys = {
        'occupation_reward',
        'allowed_pi_bonus',
        'electron_pair_penalty',
        'formal_charge_penalty',
        'formal_charge_error',
        'target_formal_charges',
        'octet_penalty',
        'octet_error',
        'valence_penalty',
        'nonclassical_penalty',
        'nonclassical_count',
        'total',
    }

    assert expected_keys.issubset(terms)
    assert terms['electron_pair_penalty'] >= 0.0
    assert terms['formal_charge_penalty'] >= 0.0
    assert terms['octet_penalty'] >= 0.0
    assert terms['valence_penalty'] >= 0.0
    assert terms['nonclassical_penalty'] >= 0.0
    assert abs(assignment['score'] - terms['total']) < 1.0e-12


def _assert_spin_data(molecule, results):

    alpha_density = results['nao_alpha_density_matrix']
    beta_density = results['nao_beta_density_matrix']
    total_density = results['nao_density_matrix']
    spin_density = results['nao_spin_density_matrix']
    spin_populations = np.array(results['nao_spin_populations'])

    assert alpha_density.shape == total_density.shape
    assert beta_density.shape == total_density.shape
    assert spin_density.shape == total_density.shape
    np.testing.assert_allclose(alpha_density + beta_density,
                               total_density,
                               atol=1.0e-10,
                               rtol=0.0)
    np.testing.assert_allclose(alpha_density - beta_density,
                               spin_density,
                               atol=1.0e-12,
                               rtol=0.0)
    assert len(spin_populations) == total_density.shape[0]
    assert abs(float(np.sum(spin_populations)) -
               results['unpaired_electrons']) < 1.0e-10
    assert abs(results['unpaired_electrons'] -
               (molecule.get_multiplicity() - 1)) < 1.0e-8


def _assert_active_partition(molecule, assignment):

    fixed_ids = _candidate_indices(assignment['fixed_nbo_list'])
    active_ids = _candidate_indices(assignment['active_nbo_list'])
    nbo_ids = _candidate_indices(assignment['nbo_list'])
    active_pi_ids = _candidate_indices(assignment['active_pi_nbo_list'])
    active_lp_ids = _candidate_indices(assignment['active_lone_pair_nbo_list'])

    assert fixed_ids.isdisjoint(active_ids)
    assert fixed_ids | active_ids == nbo_ids
    assert active_pi_ids | active_lp_ids == active_ids
    assert active_pi_ids.isdisjoint(active_lp_ids)
    assert assignment['active_electron_pairs'] == len(assignment['active_nbo_list'])
    assert assignment['active_pi_electron_pairs'] == len(assignment['active_pi_nbo_list'])
    assert assignment['active_lone_pair_electron_pairs'] == len(assignment['active_lone_pair_nbo_list'])
    assert abs(assignment['selected_lewis_electron_count'] -
               molecule.number_of_electrons()) < 1.0e-12


def _assert_ranked_by_descending_weight(alternatives):

    weights = [float(alternative.get('weight', 0.0))
               for alternative in alternatives]
    ranks = [int(alternative.get('rank', 0))
             for alternative in alternatives]

    assert ranks == list(range(1, len(alternatives) + 1))
    assert all(left >= right - 1.0e-14
               for left, right in zip(weights, weights[1:]))


def _pi_bond_set(alternative):

    return {
        tuple(sorted(pair))
        for pair in alternative.get('pi_bonds', [])
    }


@pytest.mark.solvers
class TestNboDriver:

    def test_compute_flushes_verbose_report_output(self):

        molecule = Molecule.read_str(MOLECULE_XYZ['water'])
        basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = 'hf'
        scf_drv.compute(molecule, basis)

        ostream = _RecordingOutputStream()
        nbo_drv = NboDriver(ostream=ostream)
        nbo_drv.compute(molecule, basis, scf_drv.mol_orbs)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert ostream.flush_count == 1
            assert any('Natural Population Analysis' in line
                       for line in ostream.lines)
            assert any('Natural Bond Orbital (NBO) Primary Summary' in line
                       for line in ostream.lines)

    @pytest.mark.parametrize('name', ['water', 'methane', 'ethylene', 'benzene'])
    def test_nao_foundation_invariants(self, name):

        molecule, results = _run_nbo(name)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            diagnostics = results['diagnostics']
            invariants = diagnostics['foundation_invariants']

            assert abs(diagnostics['electron_count'] -
                       molecule.number_of_electrons()) < 1.0e-8
            assert diagnostics['orthonormality_error'] < 1.0e-8
            assert diagnostics['density_formula_error'] < 1.0e-8
            assert diagnostics['density_symmetry_error'] < 1.0e-10
            assert diagnostics['electron_count_error'] < 1.0e-8
            assert diagnostics['charge_conservation_error'] < 1.0e-8
            assert diagnostics['population_trace_error'] < 1.0e-10
            assert diagnostics['mo_nao_max_normalization_error'] < 1.0e-8
            assert abs(
                np.sum(results['natural_charges']) -
                molecule.get_charge()) < 1.0e-8
            assert diagnostics['foundation_invariants_passed']
            assert set(invariants) == {
                'orthonormality',
                'density_transform',
                'density_symmetry',
                'electron_count',
                'charge_conservation',
                'population_trace',
                'mo_normalization',
            }
            assert all(item['passed'] for item in invariants.values())

    def test_nao_foundation_is_deterministic_for_water(self):

        _, first = _run_nbo('water')
        _, second = _run_nbo('water')

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            np.testing.assert_allclose(first['nao_transform'],
                                       second['nao_transform'],
                                       atol=1.0e-12,
                                       rtol=0.0)
            np.testing.assert_allclose(first['nao_density_matrix'],
                                       second['nao_density_matrix'],
                                       atol=1.0e-12,
                                       rtol=0.0)
            np.testing.assert_allclose(first['nao_populations'],
                                       second['nao_populations'],
                                       atol=1.0e-12,
                                       rtol=0.0)
            np.testing.assert_allclose(first['natural_charges'],
                                       second['natural_charges'],
                                       atol=1.0e-12,
                                       rtol=0.0)
            assert [candidate['index'] for candidate in first['nbo_candidates']] == [
                candidate['index'] for candidate in second['nbo_candidates']
            ]
            assert [candidate['type'] for candidate in first['nbo_candidates']] == [
                candidate['type'] for candidate in second['nbo_candidates']
            ]

    @pytest.mark.parametrize(
        'name, expected',
        [
            ('water', {
                'BD': 2,
                'LP': 2,
                'CR': 1
            }),
            ('methane', {
                'BD': 4,
                'LP': 0,
                'CR': 1
            }),
            ('ethylene', {
                'BD': 6,
                'LP': 0,
                'CR': 2
            }),
            ('benzene', {
                'BD': 15,
                'LP': 0,
                'CR': 6
            }),
        ],
    )
    def test_nbo_counts_for_closed_shell_examples(self, name, expected):

        _, results = _run_nbo(name)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            observed = _nbo_type_counts(results)

            assert observed['BD'] == expected['BD']
            assert observed['LP'] == expected['LP']
            assert observed['CR'] == expected['CR']

    def test_antibonding_complements_are_candidate_only(self):

        _, results = _run_nbo('ethylene')

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            norb = len(results['nao_populations'])
            candidates = results['nbo_candidates']
            primary_ids = _candidate_indices(results['primary']['nbo_list'])
            candidates_by_index = {
                int(candidate['index']): candidate for candidate in candidates
            }
            bd_candidates = [
                candidate for candidate in candidates
                if candidate['type'] == 'BD'
            ]
            bdstars = [
                candidate for candidate in candidates
                if candidate['type'] == 'BD*'
            ]

            assert len(bdstars) == len(bd_candidates)
            assert not any(candidate['index'] in primary_ids
                           for candidate in bdstars)
            assert not any(candidate['type'] == 'BD*'
                           for candidate in results['nbo_list'])

            for alternative in results['alternatives']:
                assert not any(candidate['type'] == 'BD*'
                               for candidate in alternative['nbo_list'])

            for complement in bdstars:
                parent = candidates_by_index[complement['parent_index']]
                parent_vector = _candidate_vector(parent, norb)
                complement_vector = _candidate_vector(complement, norb)

                assert parent['type'] == 'BD'
                assert complement['parent_type'] == 'BD'
                assert complement['subtype'] == parent['subtype']
                assert complement['atoms'] == parent['atoms']
                assert complement['electron_count'] == 0.0
                assert abs(float(np.dot(parent_vector,
                                        complement_vector))) < 1.0e-10
                assert abs(complement['parent_overlap']) < 1.0e-10

            assert results['primary']['counts'] == {'CR': 2, 'BD': 6}
            assert abs(results['primary']['selected_lewis_electron_count'] -
                       16.0) < 1.0e-12

    def test_rydberg_complements_and_donor_acceptor_diagnostics_are_candidate_only(self):

        _, results, drv = _run_nbo_with_driver('water', basis_label='def2-svp')

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            candidates = results['nbo_candidates']
            primary = results['primary']
            primary_ids = _candidate_indices(primary['nbo_list'])
            acceptors = [
                candidate for candidate in candidates
                if candidate['type'] in {'BD*', 'RY'}
            ]
            rydbergs = [
                candidate for candidate in candidates
                if candidate['type'] == 'RY'
            ]
            acceptor_ids = _candidate_indices(acceptors)
            diagnostics = results['donor_acceptor_diagnostics']
            summary = results['diagnostics']['donor_acceptor_summary']

            assert rydbergs
            assert primary['counts'] == {'CR': 1, 'LP': 2, 'BD': 2}
            assert not any(candidate['index'] in primary_ids
                           for candidate in acceptors)
            assert not any(candidate['type'] in {'BD*', 'RY'}
                           for candidate in results['nbo_list'])

            for alternative in results['alternatives']:
                assert not any(candidate['type'] in {'BD*', 'RY'}
                               for candidate in alternative['nbo_list'])

            for candidate in rydbergs:
                assert candidate['subtype'] == 'one-center'
                assert len(candidate['atoms']) == 1
                assert candidate['electron_count'] == 0.0
                assert candidate['occupation'] <= 0.5 + 1.0e-12
                assert candidate['source'] == 'one-center low-occupation NAO complement'

            assert diagnostics['model'] == 'nao_density_coupling'
            assert diagnostics['donor_count'] == len(primary['nbo_list'])
            assert diagnostics['acceptor_count'] == len(acceptors)
            assert diagnostics['interaction_count'] > 0
            assert summary['model'] == diagnostics['model']
            assert summary['donor_count'] == diagnostics['donor_count']
            assert summary['acceptor_count'] == diagnostics['acceptor_count']
            assert summary['interaction_count'] == diagnostics['interaction_count']

            interaction_acceptor_types = {
                item['acceptor_type']
                for item in diagnostics['interactions']
            }
            assert {'BD*', 'RY'}.issubset(interaction_acceptor_types)

            report = drv.nbo_report(level='full', return_text=True)
            assert 'Candidate-only acceptors (BD*/RY)' in report
            assert 'diagnostics only' in report
            assert 'BD*' in report
            assert 'RY(' in report
            assert '|Dens. coup.|' in report

            for interaction in diagnostics['interactions']:
                assert interaction['donor_index'] in primary_ids
                assert interaction['acceptor_index'] in acceptor_ids
                assert interaction['acceptor_type'] in {'BD*', 'RY'}
                assert interaction['density_coupling_squared'] >= 0.0
                assert interaction['abs_density_coupling'] == abs(
                    interaction['density_coupling'])
                assert interaction['abs_overlap'] == abs(interaction['overlap'])

    def test_organic_system_has_no_metal_ligand_diagnostics_or_report_section(self):

        _, results, drv = _run_nbo_with_driver('water', basis_label='def2-svp')

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert results['metal_ligand_diagnostics'] == []
            assert results['diagnostics']['metal_ligand_summary'] == {
                'candidate_count': 0,
                'channels': [],
            }
            assert not any(candidate.get('type') == 'ML'
                           for candidate in results['nbo_candidates'])
            assert not any(candidate.get('type') == 'ML'
                           for candidate in results['nbo_list'])
            for alternative in results['alternatives']:
                assert not any(candidate.get('type') == 'ML'
                               for candidate in alternative['nbo_list'])

            report = drv.nbo_report(level='full', return_text=True)
            assert 'Metal-ligand diagnostics (ML)' not in report

    def test_metal_ligand_diagnostics_are_reported_but_not_selected(self,
                                                                    monkeypatch):

        _, results, drv = _run_mock_metal_ligand_nbo(monkeypatch, ligand='NH3')

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            records = results['metal_ligand_diagnostics']
            subtypes = {record['subtype'] for record in records}
            channels = {record['channel'] for record in records}

            assert subtypes == {'sigma-acceptor', 'pi-donor'}
            assert channels == {
                'ligand-to-metal-sigma-donation',
                'metal-to-ligand-pi-back-donation',
            }
            assert not any(candidate.get('type') == 'ML'
                           for candidate in results['nbo_list'])
            assert not any(candidate.get('type') == 'ML'
                           for alternative in results['alternatives']
                           for candidate in alternative['nbo_list'])
            assert results['diagnostics']['metal_ligand_summary'] == {
                'candidate_count': len(records),
                'channels': sorted(channels),
            }

            sigma = next(record for record in records
                         if record['subtype'] == 'sigma-acceptor')
            pi = next(record for record in records
                      if record['subtype'] == 'pi-donor')
            assert sigma['metal_atom'] == 1
            assert sigma['ligand_atom'] == 2
            assert sigma['donor_atom'] == 2
            assert sigma['acceptor_atom'] == 1
            assert sigma['donation_strength'] > 0.0
            assert sigma['back_donation_strength'] == 0.0
            assert pi['donor_atom'] == 1
            assert pi['acceptor_atom'] == 2
            assert pi['back_donation_strength'] > 0.0
            assert pi['donation_strength'] == 0.0

            report = drv.nbo_report(level='full', return_text=True)
            assert 'Metal-ligand diagnostics (ML)' in report
            assert 'diagnostics only' in report
            assert 'Pd1' in report
            assert 'N2' in report
            assert 'sigma-acceptor' in report
            assert 'pi-donor' in report
            assert 'ligand-to-metal-sigma-donation' in report
            assert 'metal-to-ligand-pi-back-donation' in report
            assert 'Role: ligand sigma donor / metal acceptor' in report

    @pytest.mark.parametrize('name', ['water', 'methane', 'ethylene', 'benzene'])
    def test_lewis_accounting_for_closed_shell_examples(self, name):

        molecule, results = _run_nbo(name)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            primary = results['primary']
            _assert_lewis_accounting(molecule, primary)
            _assert_score_terms(primary)

            assert abs(primary['selected_lewis_electron_count'] -
                       molecule.number_of_electrons()) < 1.0e-12

            for alternative in results['alternatives']:
                _assert_lewis_accounting(molecule, alternative)
                _assert_score_terms(alternative)
                assert abs(alternative['selected_lewis_electron_count'] -
                           molecule.number_of_electrons()) < 1.0e-12

    def test_lewis_score_penalizes_formal_charge_separation(self):

        ring_pi_bonds = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1)]
        cross_ring_pi_bonds = [(1, 4), (2, 5), (3, 6)]
        _, results = _run_nbo(
            'benzene',
            constraints={
                'allowed_pi_bonds': ring_pi_bonds + cross_ring_pi_bonds
            },
            options={
                'max_alternatives': 24,
                'pi_min_occupation': 0.05,
                'conjugated_pi_max_path': 2,
            },
        )

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            for assignment in [results['primary']] + results['alternatives']:
                formal_charges = np.array(assignment['formal_charges'])
                formal_charge_error = float(np.sum(formal_charges**2))
                terms = assignment['score_terms']

                assert terms['target_formal_charges'] == [0.0] * len(formal_charges)
                assert abs(terms['formal_charge_error'] -
                           formal_charge_error) < 1.0e-12
                assert abs(terms['formal_charge_penalty'] -
                           formal_charge_error) < 1.0e-12

            assert results['primary']['score_terms']['formal_charge_penalty'] > 0.0
            assert any(alt['score_terms']['formal_charge_penalty'] == 0.0
                       for alt in results['alternatives'])

    def test_required_and_forbidden_bond_constraints(self):

        constraints = {
            'required_bonds': [(1, 2), (3, 4)],
            'forbidden_bonds': [(1, 6)],
        }
        _, results = _run_nbo('benzene', constraints=constraints)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            selected_pairs = _selected_bond_pairs(results)
            summary = results['diagnostics']['constraint_summary']

            assert (1, 2) in selected_pairs
            assert (3, 4) in selected_pairs
            assert (1, 6) not in selected_pairs
            assert summary['required_bonds'] == 2
            assert summary['forbidden_bonds'] == 1
            assert any('Forbidden bonds were excluded' in warning
                       for warning in results['primary']['warnings'])

    def test_allowed_and_required_pi_bond_constraints(self):

        ring_pi_bonds = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1)]
        cross_ring_pi_bonds = [(1, 4), (2, 5), (3, 6)]
        options = {
            'include_diagnostics': True,
            'max_alternatives': 24,
            'pi_min_occupation': 0.05,
            'conjugated_pi_max_path': 2,
        }

        _, allowed_results = _run_nbo(
            'benzene',
            constraints={
                'allowed_pi_bonds': ring_pi_bonds + cross_ring_pi_bonds
            },
            options=options,
        )
        _, required_results = _run_nbo(
            'benzene',
            constraints={'required_pi_bonds': [(1, 2), (3, 4), (5, 6)]},
            options=options,
        )

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            allowed_pi_sets = _alternative_pi_sets(allowed_results)
            required_pi_pairs = _selected_bond_pairs(required_results,
                                                     subtype='pi')

            assert any((1, 4) in pi_set or (2, 5) in pi_set or (3, 6) in pi_set
                       for pi_set in allowed_pi_sets)
            assert {(1, 2), (3, 4), (5, 6)}.issubset(required_pi_pairs)
            assert (required_results['diagnostics']['constraint_summary']
                    ['required_pi_bonds'] == 3)

    def test_lewis_alternatives_carry_sigma_pi_split(self):

        ring_pi_bonds = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1)]
        cross_ring_pi_bonds = [(1, 4), (2, 5), (3, 6)]
        _, results = _run_nbo(
            'benzene',
            constraints={
                'allowed_pi_bonds': ring_pi_bonds + cross_ring_pi_bonds
            },
            options={
                'max_alternatives': 24,
                'pi_min_occupation': 0.05,
                'conjugated_pi_max_path': 2,
            },
        )

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            primary = results['primary']
            assert primary['sigma_electron_pairs'] + primary['pi_electron_pairs'] == primary['electron_pairs']
            assert all(candidate.get('subtype') != 'pi'
                       for candidate in primary['sigma_nbo_list'])
            assert all(candidate.get('type') == 'BD' and
                       candidate.get('subtype') == 'pi'
                       for candidate in primary['pi_nbo_list'])

            for alternative in results['alternatives']:
                sigma_ids = _candidate_indices(alternative['sigma_nbo_list'])
                pi_ids = _candidate_indices(alternative['pi_nbo_list'])
                nbo_ids = _candidate_indices(alternative['nbo_list'])

                assert sigma_ids.isdisjoint(pi_ids)
                assert sigma_ids | pi_ids == nbo_ids
                assert alternative['sigma_electron_pairs'] == len(alternative['sigma_nbo_list'])
                assert alternative['pi_electron_pairs'] == len(alternative['pi_nbo_list'])
                assert alternative['pi_electron_pairs'] == len(alternative['pi_bonds'])
                assert alternative['sigma_electron_pairs'] + alternative['pi_electron_pairs'] == alternative['electron_pairs']

    def test_benzene_resonance_classes_group_kekule_alternatives(self):

        ring_pi_bonds = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1)]
        _, results, drv = _run_nbo_with_driver(
            'benzene',
            constraints={'allowed_pi_bonds': ring_pi_bonds},
            options={
                'max_alternatives': 12,
                'pi_min_occupation': 0.05,
                'conjugated_pi_max_path': 2,
            },
        )

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            alternatives = results['alternatives']
            classes = results['resonance_classes']

            assert len(alternatives) == 2
            assert len(classes) == 1
            assert classes[0]['class_degeneracy'] == 2
            assert classes[0]['member_ranks'] == [1, 2]
            assert abs(classes[0]['class_weight_sum'] - 1.0) < 1.0e-12
            assert classes[0]['class_label'].startswith('pi:')

            class_ids = {alternative['resonance_class_id']
                         for alternative in alternatives}
            assert class_ids == {classes[0]['class_id']}
            for alternative in alternatives:
                assert 'resonance_signature' in alternative
                assert 'resonance_label' in alternative
                assert 'class_signature' in alternative
                assert alternative['class_degeneracy'] == 2

            report = drv.nbo_report(level='full', return_text=True)
            assert 'Resonance classes' in report
            assert 'Weight sum' in report
            assert classes[0]['class_id'] in report
            assert 'Lewis/resonance alternatives' in report

    def test_lewis_alternatives_are_ranked_by_score_weight(self):

        ring_pi_bonds = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1)]
        cross_ring_pi_bonds = [(1, 4), (2, 5), (3, 6)]
        _, results = _run_nbo(
            'benzene',
            constraints={
                'allowed_pi_bonds': ring_pi_bonds + cross_ring_pi_bonds
            },
            options={
                'max_alternatives': 24,
                'pi_min_occupation': 0.05,
                'conjugated_pi_max_path': 2,
            },
        )

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            alternatives = results['alternatives']

            assert alternatives
            _assert_ranked_by_descending_weight(alternatives)
            class_member_ranks = sorted(
                rank for record in results['resonance_classes']
                for rank in record['member_ranks']
            )
            assert class_member_ranks == [alternative['rank']
                                          for alternative in alternatives]

    def test_ozone_polar_resonance_metadata_supports_visualization(self):

        molecule, results, drv = _run_nbo_with_driver(
            'ozone',
            options={
                'max_alternatives': 24,
                'lone_pair_min_occupation': 1.0,
                'pi_min_occupation': 0.05,
                'conjugated_pi_max_path': 2,
            },
        )

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            alternatives = results['alternatives']
            polar_alternatives = [
                alternative for alternative in alternatives
                if alternative.get('active_positive_atoms') and
                alternative.get('active_lone_pair_atoms')
            ]

            assert polar_alternatives
            _assert_ranked_by_descending_weight(alternatives)
            assert {(1, 2), (2, 3)}.issubset(set().union(*[
                _pi_bond_set(alternative) for alternative in polar_alternatives
            ]))

            for alternative in polar_alternatives:
                assert alternative['active_positive_atoms'] == [2]
                assert set(alternative['active_lone_pair_atoms']).issubset({1, 3})
                assert alternative['resonance_class_id']
                assert alternative['resonance_label'].startswith('pi:')
                assert 'lp:' in alternative['resonance_label']
                assert 'pos:' in alternative['resonance_label']
                _assert_lewis_accounting(molecule, alternative)
                _assert_active_partition(molecule, alternative)

            pytest.importorskip('py3Dmol')
            view = drv.show_structures(results=results,
                                       molecule=molecule,
                                       ranks=[alternative['rank']
                                              for alternative in polar_alternatives[:2]],
                                       display=False)
            assert view is not None

    @pytest.mark.parametrize(
        'name, options',
        [
            ('formate', {
                'max_alternatives': 24,
                'pi_min_occupation': 0.05,
                'conjugated_pi_max_path': 2,
            }),
            ('ozone', {
                'max_alternatives': 24,
                'lone_pair_min_occupation': 1.0,
                'pi_min_occupation': 0.05,
                'conjugated_pi_max_path': 2,
            }),
        ],
    )
    def test_lone_pair_donation_extends_pi_active_space(self, name, options):

        molecule, results = _run_nbo(name, options=options)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            donating_alternatives = [
                alternative for alternative in results['alternatives']
                if alternative['active_pi_electron_pairs'] > 0 and
                alternative['active_lone_pair_electron_pairs'] > 0
            ]
            charge_separated = [
                alternative for alternative in donating_alternatives
                if (max(alternative['formal_charges']) > 0.5 and
                    min(alternative['formal_charges']) < -0.5)
            ]

            assert donating_alternatives
            assert charge_separated

            for alternative in donating_alternatives:
                _assert_active_partition(molecule, alternative)

    @pytest.mark.parametrize(
        'name, expected_pi_bonds, options',
        [
            ('allyl_cation', {(1, 2), (2, 3)}, {}),
            ('allyl_anion', {(1, 2), (2, 3)}, {}),
            ('formate', {(1, 2), (1, 3)}, {}),
            ('nitrate', {(1, 2), (1, 3), (1, 4)}, {}),
            ('ozone', {(1, 2), (2, 3)}, {
                'lone_pair_min_occupation': 1.0,
            }),
        ],
    )
    def test_nra_structure_pool_invariants_for_representative_systems(
            self, name, expected_pi_bonds, options):

        run_options = {
            'include_nra': True,
            'nra_subspace': 'pi',
            'nra_fit_metric': 'frobenius',
            'max_alternatives': 24,
            'pi_min_occupation': 0.05,
            'conjugated_pi_max_path': 2,
        }
        run_options.update(options)
        molecule, results = _run_nbo(name, options=run_options)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            alternatives = results['alternatives']
            nra = results['nra']
            weights = np.array(nra['weights'])
            observed_pi_bonds = set().union(*[
                _pi_bond_set(alternative) for alternative in alternatives
            ])

            assert alternatives
            assert expected_pi_bonds.issubset(observed_pi_bonds)
            assert nra['subspace'] == 'pi'
            assert nra['fit_metric'] == 'frobenius'
            assert len(nra['structures']) == len(alternatives)
            assert len(weights) == len(alternatives)
            assert np.all(weights >= -1.0e-12)
            assert abs(float(np.sum(weights)) - 1.0) < 1.0e-10
            assert np.isfinite(nra['residual_norm'])
            assert np.isfinite(nra['relative_residual'])
            assert all(structure['label'] == 'sigma-only'
                       or structure['label'].startswith('pi:')
                       for structure in nra['structures'])

            for alternative in alternatives:
                _assert_lewis_accounting(molecule, alternative)
                _assert_score_terms(alternative)
                _assert_active_partition(molecule, alternative)

            for structure, alternative, weight in zip(nra['structures'],
                                                      alternatives,
                                                      weights):
                assert structure['rank'] == alternative['rank']
                assert structure['pi_bonds'] == alternative['pi_bonds']
                assert abs(structure['score_weight'] -
                           alternative['weight']) < 1.0e-12
                assert abs(structure['nra_weight'] - weight) < 1.0e-12
                assert np.isfinite(structure['residual_norm'])

    def test_user_guided_allyl_cation_has_two_pi_structures(self):

        constraints = {
            'allowed_pi_bonds': [(1, 2), (2, 3)],
        }
        options = {
            'include_nra': True,
            'nra_subspace': 'pi',
            'nra_prior_weights': {'pi:1-2': 0.5, 'pi:2-3': 0.5},
            'nra_prior_strength': 5.0,
            'max_alternatives': 12,
            'pi_min_occupation': 0.05,
            'conjugated_pi_max_path': 2,
        }
        _, results = _run_nbo('allyl_cation',
                              constraints=constraints,
                              options=options)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            structures = results['nra']['structures']
            labels = [structure['label'] for structure in structures]
            weights = np.array(results['nra']['weights'])

            assert labels == ['pi:1-2', 'pi:2-3']
            np.testing.assert_allclose(weights,
                                       np.array([0.5, 0.5]),
                                       atol=1.0e-8,
                                       rtol=0.0)

    def test_open_shell_allyl_radical_uses_somo_and_spin_resolved_nra(self):

        molecule, results = _run_nbo(
            'allyl_radical',
            options={
                'include_nra': True,
                'nra_subspace': 'pi',
                'max_alternatives': 24,
                'pi_min_occupation': 0.05,
                'conjugated_pi_max_path': 2,
            },
        )

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            primary = results['primary']
            somo_candidates = [
                candidate for candidate in results['nbo_candidates']
                if candidate['type'] == 'SOMO'
            ]
            radical_lone_pair_candidates = [
                candidate for candidate in results['nbo_candidates']
                if candidate['type'] == 'LP' and
                candidate.get('subtype') == 'radical-lone-pair'
            ]
            radical_pi_candidates = [
                candidate for candidate in results['nbo_candidates']
                if candidate['type'] == 'BD' and
                candidate.get('subtype') == 'pi' and
                candidate.get('electron_count') == 1.0
            ]
            selected_somos = [
                candidate for candidate in primary['nbo_list']
                if candidate['type'] == 'SOMO'
            ]

            _assert_spin_data(molecule, results)
            _assert_lewis_accounting(molecule, primary)
            _assert_score_terms(primary)

            assert somo_candidates
            assert radical_lone_pair_candidates
            assert radical_pi_candidates
            assert all(candidate['electron_count'] == 1.0
                       for candidate in radical_lone_pair_candidates)
            assert all(candidate['electron_count'] == 1.0
                       for candidate in radical_pi_candidates)
            assert all(candidate['spin_occupation'] > 0.0
                       for candidate in radical_lone_pair_candidates)
            assert all(candidate['spin_occupation'] > 0.0
                       for candidate in radical_pi_candidates)
            assert len(selected_somos) == 1
            assert selected_somos[0]['electron_count'] == 1.0
            assert primary['counts']['SOMO'] == 1
            assert abs(primary['selected_lewis_electron_count'] -
                       molecule.number_of_electrons()) < 1.0e-12
            assert any(candidate.get('type') == 'BD' and
                       candidate.get('subtype') == 'pi'
                       for candidate in primary['nbo_list'])

            radical_alternatives = [
                alternative for alternative in results['alternatives']
                if alternative['active_one_electron_count'] == 1
            ]
            assert radical_alternatives
            for alternative in radical_alternatives:
                assert abs(alternative['selected_lewis_electron_count'] -
                           molecule.number_of_electrons()) < 1.0e-12
                assert abs(alternative['active_electron_count'] -
                           sum(_candidate_electron_count(candidate)
                               for candidate in alternative['active_nbo_list'])) < 1.0e-12
                assert alternative['active_one_electron_count'] == len(
                    alternative['active_one_electron_nbo_list'])
                assert all(_candidate_electron_count(candidate) == 1.0
                           for candidate in alternative['active_one_electron_nbo_list'])
                assert any(candidate.get('type') in {'SOMO', 'LP'}
                           for candidate in alternative['active_one_electron_nbo_list'])

            nra = results['nra']
            weights = np.array(nra['weights'])
            assert nra['spin_fit'] == 'total_spin'
            assert len(weights) == len(results['alternatives'])
            assert len(nra['structures']) == len(results['alternatives'])
            assert np.all(weights >= -1.0e-12)
            assert abs(float(np.sum(weights)) - 1.0) < 1.0e-10
            assert np.isfinite(nra['residual_norm'])
            assert np.isfinite(nra['relative_residual'])
            assert np.isfinite(nra['total_residual_norm'])
            assert np.isfinite(nra['relative_total_residual'])
            assert np.isfinite(nra['spin_residual_norm'])
            assert np.isfinite(nra['relative_spin_residual'])
            assert any('spin-density residuals' in warning
                       for warning in nra['warnings'])

            for structure, alternative, weight in zip(nra['structures'],
                                                      results['alternatives'],
                                                      weights):
                assert structure['rank'] == alternative['rank']
                assert abs(structure['nra_weight'] - weight) < 1.0e-12
                assert np.isfinite(structure['residual_norm'])
                assert np.isfinite(structure['spin_residual_norm'])

    def test_nra_is_optional_and_fits_benzene_alternatives(self):

        ring_pi_bonds = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1)]
        cross_ring_pi_bonds = [(1, 4), (2, 5), (3, 6)]
        constraints = {
            'allowed_pi_bonds': ring_pi_bonds + cross_ring_pi_bonds,
        }

        _, default_results = _run_nbo('benzene', constraints=constraints)
        _, nra_results = _run_nbo(
            'benzene',
            constraints=constraints,
            options={
                'include_nra': True,
                'nra_subspace': 'pi',
                'nra_fit_metric': 'frobenius',
                'max_alternatives': 24,
                'pi_min_occupation': 0.05,
                'conjugated_pi_max_path': 2,
            },
        )

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            nra = nra_results['nra']
            weights = np.array(nra['weights'])
            alternating_ring_weights = []
            labels = [structure['label'] for structure in nra['structures']]
            for structure in nra['structures']:
                pi_set = tuple(
                    sorted(
                        tuple(sorted(pair)) for pair in structure['pi_bonds']))
                if pi_set in {
                    ((1, 2), (3, 4), (5, 6)),
                    ((1, 6), (2, 3), (4, 5)),
                }:
                    alternating_ring_weights.append(structure['nra_weight'])

            assert 'nra' not in default_results
            assert nra['subspace'] == 'pi'
            assert nra['fit_metric'] == 'frobenius'
            assert len(weights) == len(nra_results['alternatives'])
            assert np.all(weights >= -1.0e-12)
            assert abs(float(np.sum(weights)) - 1.0) < 1.0e-10
            assert np.isfinite(nra['residual_norm'])
            assert np.isfinite(nra['relative_residual'])
            assert len(alternating_ring_weights) == 2
            assert abs(alternating_ring_weights[0] - alternating_ring_weights[1]) < 1.0e-8
            assert len(set(labels)) == len(labels)
            assert all(label.startswith('pi:') for label in labels)
            assert 'pi:1-2,3-4,5-6' in labels
            assert 'pi:1-6,2-3,4-5' in labels

    def test_nra_prior_regularization_is_optional_and_reported(self):

        ring_pi_bonds = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1)]
        cross_ring_pi_bonds = [(1, 4), (2, 5), (3, 6)]
        constraints = {
            'allowed_pi_bonds': ring_pi_bonds + cross_ring_pi_bonds,
        }
        base_options = {
            'include_nra': True,
            'nra_subspace': 'pi',
            'nra_fit_metric': 'frobenius',
            'max_alternatives': 24,
            'pi_min_occupation': 0.05,
            'conjugated_pi_max_path': 2,
        }

        _, baseline_results = _run_nbo('benzene',
                                       constraints=constraints,
                                       options=base_options)
        alternative_count = len(baseline_results['alternatives'])
        zero_prior_options = dict(base_options)
        zero_prior_options.update({
            'nra_prior_weights': [0.0] * alternative_count,
            'nra_prior_strength': 25.0,
        })
        _, zero_prior_results = _run_nbo('benzene',
                                         constraints=constraints,
                                         options=zero_prior_options)

        regularized_prior = [0.0] * alternative_count
        regularized_prior[0] = 1.0
        prior_options = dict(base_options)
        prior_options.update({
            'nra_prior_weights': regularized_prior,
            'nra_prior_strength': 25.0,
        })
        _, prior_results, prior_drv = _run_nbo_with_driver(
            'benzene',
            constraints=constraints,
            options=prior_options,
        )

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            baseline_nra = baseline_results['nra']
            zero_prior_nra = zero_prior_results['nra']
            prior_nra = prior_results['nra']
            baseline_weights = np.array(baseline_nra['weights'])
            zero_prior_weights = np.array(zero_prior_nra['weights'])
            prior_weights = np.array(prior_nra['weights'])

            np.testing.assert_allclose(zero_prior_weights,
                                       baseline_weights,
                                       atol=1.0e-12,
                                       rtol=0.0)
            assert not zero_prior_nra['prior']['active']
            assert zero_prior_nra['prior']['strength'] == 25.0
            assert any('all zero' in warning
                       for warning in zero_prior_nra['warnings'])

            assert prior_nra['prior']['active']
            assert prior_nra['prior']['mode'] == 'regularized'
            assert prior_nra['prior']['strength'] == 25.0
            assert prior_nra['prior']['fixed'] is False
            assert prior_nra['prior']['input_weights'] == regularized_prior
            assert prior_nra['prior']['normalized_weights'] == regularized_prior
            assert abs(float(np.sum(prior_weights)) - 1.0) < 1.0e-10
            assert np.all(prior_weights >= -1.0e-12)
            assert abs(prior_nra['structures'][0]['prior_weight'] - 1.0) < 1.0e-12
            assert all(structure['prior_weight'] is not None
                       for structure in prior_nra['structures'])

            report = prior_drv.nra_report(return_text=True)
            assert 'Prior mode: regularized' in report
            assert 'Prior strength: 25' in report
            assert 'Prior weight' in report
            assert 'Score weight' in report
            assert 'NRA weight' in report

    def test_nra_report_for_benzene_alternatives(self):

        ring_pi_bonds = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1)]
        cross_ring_pi_bonds = [(1, 4), (2, 5), (3, 6)]
        constraints = {
            'allowed_pi_bonds': ring_pi_bonds + cross_ring_pi_bonds,
        }

        molecule, default_results, default_drv = _run_nbo_with_driver(
            'benzene',
            constraints=constraints,
        )
        _, nra_results, nra_drv = _run_nbo_with_driver(
            'benzene',
            constraints=constraints,
            options={
                'include_nra': True,
                'nra_subspace': 'pi',
                'nra_fit_metric': 'frobenius',
                'max_alternatives': 24,
                'pi_min_occupation': 0.05,
                'conjugated_pi_max_path': 2,
            },
        )

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            report = nra_drv.nra_report(return_text=True)
            full_report = nra_drv.nra_report(level='full', return_text=True)
            weights = np.array(nra_results['nra']['weights'])

            assert 'Natural Resonance Analysis' in report
            assert 'NRA weight' in report
            assert 'Score weight' in report
            assert 'Signature' in report
            assert 'pi:' in report
            assert 'Residual' in report
            assert 'Score weights are Lewis ranking weights' in report
            assert f'{nra_results["nra"]["structures"][0]["nra_weight"]:>12.6f}' in report
            assert f'{nra_results["nra"]["structures"][0]["score_weight"]:>13.6f}' in report
            assert 'Structure details' in full_report
            assert abs(float(np.sum(weights)) - 1.0) < 1.0e-10

            with pytest.raises(Exception, match='include_nra=True'):
                default_drv.format_nra_report(molecule, default_results)
