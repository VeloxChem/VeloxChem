from collections import namedtuple
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.dispersionmodel import DispersionModel
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfhessiandriver import ScfHessianDriver
import veloxchem.batchsize as batchsize


class TestHessianAtomPair:

    @pytest.mark.solvers
    @pytest.mark.parametrize('ri', [False, True])
    def test_hessian_atom_pair(self, ri):

        mol = Molecule.read_xyz_string("""15

        O       99.814000000   100.835000000   101.232000000
        H       99.329200000    99.976800000   101.063000000
        H       99.151600000   101.561000000   101.414000000
        O       98.804000000    98.512200000    97.758100000
        H       99.782100000    98.646900000    97.916700000
        H       98.421800000    99.326500000    97.321300000
        O       98.070300000    98.516900000   100.438000000
        H       97.172800000    98.878600000   100.690000000
        H       98.194000000    98.592200000    99.448100000
        O      102.360000000   101.551000000    99.964500000
        H      102.675000000   102.370000000   100.444000000
        H      101.556000000   101.180000000   100.430000000
        O      101.665000000    98.316100000    98.319400000
        H      101.904000000    99.233800000    98.002000000
        H      102.224000000    97.640900000    97.837700000""")

        basis = MolecularBasis.read(mol, 'sto-3g')

        atom_pairs = [(9, 6), (0, 3)]

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_drv.acc_type = 'l2_diis'
        scf_results = scf_drv.compute(mol, basis)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        hess_drv = ScfHessianDriver(scf_drv)
        hess_drv.atom_pairs = atom_pairs
        hess_drv.update_settings({}, cphf_dict={'conv_thresh': 1.0e-6})
        hess_drv.compute(mol, basis)

        pair_hess = hess_drv.hessian

        if scf_drv.rank == mpi_master():
            cwd = Path(__file__).parent
            path = str(
                cwd / 'data' / {
                    False: 'water_hessian_pair_1-4_7-10.txt',
                    True: 'water_hessian_pair_1-4_7-10_ri.txt'
                }[ri])
            reference_hess = np.loadtxt(path)
            np.testing.assert_allclose(pair_hess,
                                       reference_hess,
                                       rtol=1e-8,
                                       atol=5e-7)

    @pytest.mark.solvers
    @pytest.mark.skipif(not DispersionModel.is_available(),
                        reason='dftd4-python not available')
    @pytest.mark.parametrize('ri', [False, True])
    def test_hessian_atom_pair_with_d4_matches_full_hessian(self, ri):

        mol = Molecule.read_xyz_string("""15

        O       99.814000000   100.835000000   101.232000000
        H       99.329200000    99.976800000   101.063000000
        H       99.151600000   101.561000000   101.414000000
        O       98.804000000    98.512200000    97.758100000
        H       99.782100000    98.646900000    97.916700000
        H       98.421800000    99.326500000    97.321300000
        O       98.070300000    98.516900000   100.438000000
        H       97.172800000    98.878600000   100.690000000
        H       98.194000000    98.592200000    99.448100000
        O      102.360000000   101.551000000    99.964500000
        H      102.675000000   102.370000000   100.444000000
        H      101.556000000   101.180000000   100.430000000
        O      101.665000000    98.316100000    98.319400000
        H      101.904000000    99.233800000    98.002000000
        H      102.224000000    97.640900000    97.837700000
        """)
        basis = MolecularBasis.read(mol, 'sto-3g')
        atom_pairs = [(9, 6), (0, 3)]

        def run_hessian_with_dispersion(atom_pairs_input=None):

            scf_drv = ScfRestrictedDriver()
            scf_drv.ri = ri
            scf_drv.ostream.mute()
            scf_drv.dispersion = True
            scf_results = scf_drv.compute(mol, basis)

            if scf_drv.rank == mpi_master():
                if ri:
                    assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                    assert not scf_results['ri_fallback']

            hess_drv = ScfHessianDriver(scf_drv)
            hess_drv.ostream.mute()
            hess_drv.atom_pairs = atom_pairs_input
            hess_drv.update_settings({}, cphf_dict={'conv_thresh': 1.0e-6})
            hess_drv.compute(mol, basis)

            return scf_drv, hess_drv.hessian

        scf_drv, full_hessian = run_hessian_with_dispersion()
        _, pair_hessian = run_hessian_with_dispersion(atom_pairs)

        if scf_drv.rank == mpi_master():
            requested_pairs = set(atom_pairs)
            involved_atoms = sorted(
                {idx for pair in atom_pairs for idx in pair})
            requested_pairs.update((i, i) for i in involved_atoms)

            for i in range(mol.number_of_atoms()):
                for j in range(mol.number_of_atoms()):
                    block = pair_hessian[i * 3:(i + 1) * 3, j * 3:(j + 1) * 3]
                    if ((i, j) in requested_pairs or (j, i) in requested_pairs):
                        np.testing.assert_allclose(
                            pair_hessian[i * 3:(i + 1) * 3, j * 3:(j + 1) * 3],
                            full_hessian[i * 3:(i + 1) * 3, j * 3:(j + 1) * 3],
                            rtol=1.0e-7,
                            atol=1.0e-9)
                    else:
                        np.testing.assert_allclose(block,
                                                   0.0,
                                                   rtol=1.0e-7,
                                                   atol=1.0e-9)

    @pytest.mark.solvers
    def test_batches_of_the_atoms_in_pairs(self, monkeypatch):
        """The exchange-correlation part over the atoms in pairs takes every
        one of them, also where they do not fill the last batch: four atoms,
        three of them in pairs, two to a batch."""

        mol = Molecule.read_xyz_string("""4

        N     0.000000000     0.000000000     0.000000000
        H    -0.653401663     0.309213352     0.817609879
        H     0.695693936     0.071283622    -0.702632331
        H    -0.468651181    -0.844282500    -0.433152944""")

        basis = MolecularBasis.read(mol, 'sto-3g', ostream=None)
        nao = basis.get_dimensions_of_basis()

        atom_pairs = [(0, 1), (1, 2)]

        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = 'b3lyp'
        scf_drv.ostream.mute()
        scf_drv.compute(mol, basis)

        def run():
            hess_drv = ScfHessianDriver(scf_drv)
            hess_drv.atom_pairs = atom_pairs
            hess_drv.update_settings({}, cphf_dict={'conv_thresh': 1.0e-8})
            hess_drv.compute(mol, basis)
            return hess_drv.hessian

        whole = run()

        # NOTE: get_batch_size gives 6 Fock matrices for the 9 of the three
        # atoms in pairs (an estimate of nao^2 doubles a thread for each, with
        # four threads, and 6.5 of them free to the one rank on the node), so
        # two atoms to a batch. Only the name psutil in batchsize is replaced,
        # so the RI memory budgets, which read the free memory elsewhere, are
        # not touched.
        memory = namedtuple('memory', ['total', 'available'])
        available = int(6.5 * nao * nao * 8 * 4)
        monkeypatch.setattr(
            batchsize, 'psutil',
            SimpleNamespace(
                virtual_memory=lambda: memory(64 * 1024**3, available)))
        monkeypatch.setenv('OMP_NUM_THREADS', '4')
        monkeypatch.setenv('SLURM_NTASKS_PER_NODE', '1')

        assert batchsize.get_batch_size(None, 9, nao, scf_drv.comm) == 6

        batched = run()

        if scf_drv.rank == mpi_master():
            assert np.max(np.abs(whole - batched)) < 1.0e-8
