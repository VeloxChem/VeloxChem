import numpy as np
import pytest
from mpi4py import MPI

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.nbodriver import NboDriver
from veloxchem.scfrestdriver import ScfRestrictedDriver

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
}


def _run_nbo(name, constraints=None, options=None):

    molecule = Molecule.read_str(MOLECULE_XYZ[name])
    basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)

    scf_drv = ScfRestrictedDriver()
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

    return molecule, results


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


def _alternative_pi_sets(results):

    return {
        tuple(sorted(tuple(sorted(pair))
                     for pair in alt.get('pi_bonds', [])))
        for alt in results.get('alternatives', [])
    }


@pytest.mark.solvers
class TestNboDriver:

    def test_nao_invariants_for_water(self):

        molecule, results = _run_nbo('water')

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            diagnostics = results['diagnostics']

            assert abs(diagnostics['electron_count'] -
                       molecule.number_of_electrons()) < 1.0e-8
            assert diagnostics['orthonormality_error'] < 1.0e-8
            assert diagnostics['mo_nao_max_normalization_error'] < 1.0e-8
            assert abs(
                np.sum(results['natural_charges']) -
                molecule.get_charge()) < 1.0e-8

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
        dewar_para_pi_bonds = [(1, 4), (2, 5), (3, 6)]
        options = {
            'include_diagnostics': True,
            'max_alternatives': 24,
            'pi_min_occupation': 0.05,
            'conjugated_pi_max_path': 2,
        }

        _, allowed_results = _run_nbo(
            'benzene',
            constraints={
                'allowed_pi_bonds': ring_pi_bonds + dewar_para_pi_bonds
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

    def test_nra_is_optional_and_fits_benzene_alternatives(self):

        ring_pi_bonds = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1)]
        dewar_para_pi_bonds = [(1, 4), (2, 5), (3, 6)]
        constraints = {
            'allowed_pi_bonds': ring_pi_bonds + dewar_para_pi_bonds,
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
            kekule_weights = []
            for structure in nra['structures']:
                pi_set = tuple(
                    sorted(
                        tuple(sorted(pair)) for pair in structure['pi_bonds']))
                if pi_set in {
                    ((1, 2), (3, 4), (5, 6)),
                    ((1, 6), (2, 3), (4, 5)),
                }:
                    kekule_weights.append(structure['nra_weight'])

            assert 'nra' not in default_results
            assert nra['subspace'] == 'pi'
            assert nra['fit_metric'] == 'frobenius'
            assert len(weights) == len(nra_results['alternatives'])
            assert np.all(weights >= -1.0e-12)
            assert abs(float(np.sum(weights)) - 1.0) < 1.0e-10
            assert np.isfinite(nra['residual_norm'])
            assert np.isfinite(nra['relative_residual'])
            assert len(kekule_weights) == 2
            assert abs(kekule_weights[0] - kekule_weights[1]) < 1.0e-8
