from mpi4py import MPI
from pathlib import Path
import h5py
import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.errorhandler import VeloxChemError
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lrsolver import LinearResponseSolver
from veloxchem.fikapdbreader import FikaPdbReader
from veloxchem.veloxchemlib import (FikaClassicalSystem, FikaForceField,
                                    fika_polarizable_sites)

from test_fika_pdbreader import get_pdb, get_waters


def write_droplet(directory):
    """
    Acrolein (MOL) and 40 waters in a PDB file in a directory shared by all
    ranks.
    """

    path = None
    if MPI.COMM_WORLD.Get_rank() == mpi_master():
        path = Path(directory) / 'fika_droplet.pdb'
        path.write_text(get_pdb(get_waters(40)))
        path = str(path)
    return MPI.COMM_WORLD.bcast(path, root=mpi_master())


def fika_embedding(pdb_file, damping=None):

    return {
        'settings': {
            'embedding_method': 'fika',
            'damping': damping,
        },
        'inputs': {
            'pdb_file': pdb_file,
            'solvent': {
                'HOH': 'ahlstrom'
            },
        },
    }


def run_scf(molecule, basis, embedding):

    scf_drv = ScfRestrictedDriver()
    scf_drv.ostream.mute()
    scf_drv.conv_thresh = 1.0e-8
    scf_drv.embedding = embedding
    return scf_drv, scf_drv.compute(molecule, basis)


def run_lrs(molecule, basis, embedding, scf_results):

    lrs_drv = LinearResponseSolver()
    lrs_drv.ostream.mute()
    lrs_drv.conv_thresh = 1.0e-8
    lrs_drv.frequencies = [0.0, 0.04]
    lrs_drv.embedding = embedding
    results = lrs_drv.compute(molecule, basis, scf_results)
    if MPI.COMM_WORLD.Get_rank() != mpi_master():
        return None
    return np.array([
        results['response_functions'][(a, b, w)]
        for w in lrs_drv.frequencies
        for a in 'xyz'
        for b in 'xyz'
    ])


@pytest.mark.solvers
class TestFikaScf:

    def test_scf_and_response(self, tmp_path):

        pdb_file = write_droplet(tmp_path)
        reader = FikaPdbReader()
        residues = reader.read(pdb_file)
        molecule = reader.get_solute(residues)
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)

        ref_file = Path(__file__).parent / 'data' / 'fika_pe_reference.h5'
        with h5py.File(ref_file, 'r') as hf:
            ref_f_emb = np.array(hf.get('droplet_f_emb'))
            ref_response = np.array(hf.get('droplet_response'))

        fika = fika_embedding(pdb_file)
        fika_drv, fika_results = run_scf(molecule, basis, fika)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert fika_results['scf_energy'] == pytest.approx(
                -190.61770854930165, abs=1.0e-8)
            assert fika_results['E_emb'] == pytest.approx(
                -0.0032606747801258718, abs=1.0e-8)
            assert np.max(np.abs(fika_results['F_emb'] -
                                 ref_f_emb)) < 1.0e-7

            # Induced dipoles of the final density, one per oxygen.
            positions, _, _ = fika_polarizable_sites(
                fika_drv._embedding_drv.classical_system)
            assert fika_results['induced_dipoles'].shape == positions.shape

            # Each SCF step starts the dipoles from those of the previous one.
            iterations = fika_drv._embedding_drv.iterations
            assert max(iterations[len(iterations) // 2:]) < iterations[0]

        fika_prop = run_lrs(molecule, basis, fika, fika_results)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert np.max(np.abs(fika_prop - ref_response)) < 1.0e-6

    def test_excitations_and_complex_response(self, tmp_path):

        from veloxchem.lreigensolver import LinearResponseEigenSolver
        from veloxchem.tdaeigensolver import TdaEigenSolver
        from veloxchem import ComplexResponse

        pdb_file = write_droplet(tmp_path)
        reader = FikaPdbReader()
        residues = reader.read(pdb_file)
        molecule = reader.get_solute(residues)
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)

        ref_file = Path(__file__).parent / 'data' / 'fika_pe_reference.h5'
        with h5py.File(ref_file, 'r') as hf:
            ref_lr_eigenvalues = np.array(hf.get('droplet_lr_eigenvalues'))
            ref_tda_eigenvalues = np.array(hf.get('droplet_tda_eigenvalues'))
            ref_complex_response = np.array(
                hf.get('droplet_complex_response'))

        embedding = fika_embedding(pdb_file)
        scf_results = run_scf(molecule, basis, embedding)[1]

        def solve(solver, **settings):
            solver.ostream.mute()
            solver.conv_thresh = 1.0e-6
            solver.embedding = embedding
            for key, value in settings.items():
                setattr(solver, key, value)
            return solver.compute(molecule, basis, scf_results)

        rank_is_master = MPI.COMM_WORLD.Get_rank() == mpi_master()

        for solver_class, ref_eigenvalues in [
            (LinearResponseEigenSolver, ref_lr_eigenvalues),
            (TdaEigenSolver, ref_tda_eigenvalues),
        ]:
            results = solve(solver_class(), nstates=4)
            if rank_is_master:
                assert np.max(np.abs(results['eigenvalues'] -
                                     ref_eigenvalues)) < 1.0e-8

        frequencies = [0.0, 0.1, 0.2]
        results = solve(ComplexResponse(),
                        frequencies=frequencies,
                        damping=0.004)
        if rank_is_master:
            values = np.array([
                results['response_functions'][(a, a, w)]
                for w in frequencies
                for a in 'xyz'
            ])
            assert np.max(np.abs(values - ref_complex_response)) < (
                1.0e-7 * np.max(np.abs(ref_complex_response)))

    def test_response_with_subcommunicators(self, tmp_path):

        # With subcommunicators each one embeds its own trial densities.
        from veloxchem.lreigensolver import LinearResponseEigenSolver

        pdb_file = write_droplet(tmp_path)
        reader = FikaPdbReader()
        molecule = reader.get_solute(reader.read(pdb_file))
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)
        embedding = fika_embedding(pdb_file)
        _, scf_results = run_scf(molecule, basis, embedding)

        def excitation_energies(use_subcomms):
            solver = LinearResponseEigenSolver()
            solver.ostream.mute()
            solver.conv_thresh = 1.0e-6
            solver.nstates = 4
            solver.use_subcomms = use_subcomms
            solver.embedding = embedding
            results = solver.compute(molecule, basis, scf_results)
            return results['eigenvalues'] if MPI.COMM_WORLD.Get_rank(
            ) == mpi_master() else None

        reference = excitation_energies(False)
        distributed = excitation_energies(True)
        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert np.max(np.abs(distributed - reference)) < 1.0e-10

    def test_unrestricted(self, tmp_path):

        from veloxchem.molecule import Molecule
        from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
        from veloxchem.scfrestopendriver import ScfRestrictedOpenDriver
        from veloxchem import (ComplexResponse, ComplexResponseUnrestrictedSolver,
                               LinearResponseUnrestrictedSolver,
                               LinearResponseUnrestrictedEigenSolver,
                               TdaEigenSolver, TdaUnrestrictedEigenSolver)

        pdb_file = write_droplet(tmp_path)
        reader = FikaPdbReader()
        molecule = reader.get_solute(reader.read(pdb_file))
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)
        embedding = fika_embedding(pdb_file)
        master = MPI.COMM_WORLD.Get_rank() == mpi_master()

        def scf(driver_class, mol):
            driver = driver_class()
            driver.ostream.mute()
            driver.conv_thresh = 1.0e-8
            driver.embedding = embedding
            return driver.compute(mol, basis)

        def solve(solver, scf_results, mol, **settings):
            solver.ostream.mute()
            solver.conv_thresh = 1.0e-6
            solver.embedding = embedding
            for key, value in settings.items():
                setattr(solver, key, value)
            return solver.compute(mol, basis, scf_results)

        def functions(results, frequencies):
            return np.array([
                results['response_functions'][(a, b, w)]
                for w in frequencies
                for a in 'xyz'
                for b in 'xyz'
            ])

        # A closed shell through the unrestricted path gives the restricted
        # results (TDA: closed-shell unrestricted RPA has a triplet
        # instability).
        restricted = run_scf(molecule, basis, embedding)[1]
        unrestricted = scf(ScfUnrestrictedDriver, molecule)
        frequencies = [0.0, 0.05]
        alpha = solve(LinearResponseSolver(), restricted, molecule,
                      frequencies=frequencies)
        alpha_u = solve(LinearResponseUnrestrictedSolver(), unrestricted,
                        molecule,
                        frequencies=frequencies)
        cpp = solve(ComplexResponse(), restricted, molecule,
                    frequencies=[0.1], damping=0.004)
        cpp_u = solve(ComplexResponseUnrestrictedSolver(), unrestricted,
                      molecule,
                      frequencies=[0.1],
                      damping=0.004)
        tda = solve(TdaEigenSolver(), restricted, molecule, nstates=3)
        tda_u = solve(TdaUnrestrictedEigenSolver(), unrestricted, molecule,
                      nstates=8)
        if master:
            assert unrestricted['scf_energy'] == pytest.approx(
                restricted['scf_energy'], abs=1.0e-9)
            assert np.max(
                np.abs(functions(alpha_u, frequencies) -
                       functions(alpha, frequencies))) < 1.0e-7
            assert np.max(
                np.abs(functions(cpp_u, [0.1]) - functions(cpp, [0.1]))) < 1.0e-7
            for energy in tda['eigenvalues']:
                assert np.min(np.abs(tda_u['eigenvalues'] - energy)) < 1.0e-8

        # An open shell: the cation doublet with UHF and ROHF, and the
        # unrestricted solvers.
        cation = Molecule(molecule)
        cation.set_charge(1)
        cation.set_multiplicity(2)
        uhf = scf(ScfUnrestrictedDriver, cation)
        rohf = scf(ScfRestrictedOpenDriver, cation)
        polarizability = solve(LinearResponseUnrestrictedSolver(), uhf, cation,
                               frequencies=[0.0])
        rpa = solve(LinearResponseUnrestrictedEigenSolver(), uhf, cation,
                    nstates=3,
                    max_iter=400)
        if master:
            assert uhf['E_emb'] < 0.0 and rohf['E_emb'] < 0.0
            assert uhf['scf_energy'] < rohf['scf_energy']
            assert all(
                polarizability['response_functions'][(a, a, 0.0)] < 0.0
                for a in 'xyz')
            assert np.all(rpa['eigenvalues'] > 0.0)

    def test_damping_and_objects(self, tmp_path):

        pdb_file = write_droplet(tmp_path)
        reader = FikaPdbReader()
        residues = reader.read(pdb_file)
        molecule = reader.get_solute(residues)
        basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)

        # Thole damping is the default.
        thole = fika_embedding(pdb_file)
        del thole['settings']['damping']
        _, thole_results = run_scf(molecule, basis, thole)
        _, undamped_results = run_scf(molecule, basis, fika_embedding(pdb_file))

        # The same classical system given as an object.
        system, _ = reader.get_classical_system(residues, {'HOH': 'ahlstrom'})
        objects = {
            'settings': {
                'embedding_method': 'FIKA'
            },
            'inputs': {
                'objects': {
                    'classical_system': system
                }
            },
        }
        _, object_results = run_scf(molecule, basis, objects)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert abs(thole_results['scf_energy'] -
                       undamped_results['scf_energy']) > 1.0e-6
            assert object_results['scf_energy'] == pytest.approx(
                thole_results['scf_energy'], abs=1.0e-9)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='input errors abort MPI runs')
    def test_errors(self, tmp_path):

        pdb_file = write_droplet(tmp_path)
        reader = FikaPdbReader()
        residues = reader.read(pdb_file)
        molecule = reader.get_solute(residues)
        basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)

        def fails(embedding, message, molecule=molecule, basis=basis):
            with pytest.raises(VeloxChemError, match=message):
                run_scf(molecule, basis, embedding)

        # The molecule must be the solute of the PDB file.
        moved = molecule.get_coordinates_in_bohr()
        moved[0, 0] += 0.1
        from veloxchem.molecule import Molecule
        other = Molecule([int(z) for z in molecule.get_element_ids()], moved,
                         'bohr')
        fails(fika_embedding(pdb_file), 'not the solute of the PDB file',
              molecule=other,
              basis=MolecularBasis.read(other, 'sto-3g', ostream=None))

        embedding = fika_embedding(pdb_file)
        embedding['settings']['damping'] = 'exponential'
        fails(embedding, "damping must be 'thole' or None")

        embedding = fika_embedding(pdb_file)
        embedding['settings']['solver'] = 'jidiis'
        fails(embedding, 'unknown settings')

        embedding = fika_embedding(pdb_file)
        del embedding['inputs']['solvent']
        fails(embedding, "inputs need 'solvent'")

        embedding = fika_embedding(pdb_file)
        embedding['inputs']['objects'] = {}
        fails(embedding, "either 'pdb_file', 'potfile' or 'objects'")

        # Gradients are not available yet.
        from veloxchem.scfgradientdriver import ScfGradientDriver
        scf_drv, scf_results = run_scf(molecule, basis,
                                       fika_embedding(pdb_file))
        grad_drv = ScfGradientDriver(scf_drv)
        grad_drv.ostream.mute()
        with pytest.raises(VeloxChemError,
                           match='gradients with the fika embedding'):
            grad_drv.compute(molecule, basis, scf_results)

    def test_response_on_fika_reference(self, tmp_path):

        pdb_file = write_droplet(tmp_path)
        reader = FikaPdbReader()
        molecule = reader.get_solute(reader.read(pdb_file))
        basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)

        _, fika_results = run_scf(molecule, basis, fika_embedding(pdb_file))
        _, gas_results = run_scf(molecule, basis, None)
        # the SCF results are filled on the master rank only
        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert fika_results['embedding_method'] == 'fika'
            assert 'embedding_method' not in gas_results

        # Linear response on a fika reference needs the fika embedding; a
        # reused driver accepts a gas-phase reference afterwards.
        lrs_drv = LinearResponseSolver()
        lrs_drv.ostream.mute()
        lrs_drv.frequencies = [0.0]
        lrs_drv.embedding = fika_embedding(pdb_file)
        lrs_drv.compute(molecule, basis, fika_results)
        lrs_drv.embedding = None
        lrs_drv.compute(molecule, basis, gas_results)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_response_on_fika_reference_errors(self, tmp_path):

        pdb_file = write_droplet(tmp_path)
        reader = FikaPdbReader()
        molecule = reader.get_solute(reader.read(pdb_file))
        basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)

        _, fika_results = run_scf(molecule, basis, fika_embedding(pdb_file))
        _, gas_results = run_scf(molecule, basis, None)

        # Linear response on a fika reference needs the fika embedding; a
        # reused driver accepts a gas-phase reference afterwards.
        lrs_drv = LinearResponseSolver()
        lrs_drv.ostream.mute()
        lrs_drv.frequencies = [0.0]
        with pytest.raises(VeloxChemError,
                           match='SCF reference used the fika embedding'):
            lrs_drv.compute(molecule, basis, fika_results)
        lrs_drv.compute(molecule, basis, gas_results)

        # Nonlinear response rejects the fika embedding, whether it comes
        # from the SCF reference or is set on the driver.
        from veloxchem.quadraticresponsedriver import QuadraticResponseDriver
        from veloxchem.cubicresponsedriver import CubicResponseDriver

        def nonlinear_driver(driver_type):
            driver = driver_type()
            driver.ostream.mute()
            for label in 'abcd':
                if hasattr(driver, f'{label}_component'):
                    setattr(driver, f'{label}_component', 'z')
            return driver

        for driver_type in (QuadraticResponseDriver, CubicResponseDriver):
            driver = nonlinear_driver(driver_type)
            with pytest.raises(VeloxChemError,
                               match='fika embedding is not supported'):
                driver.compute(molecule, basis, fika_results)

        driver = nonlinear_driver(QuadraticResponseDriver)
        driver.embedding = fika_embedding(pdb_file)
        with pytest.raises(VeloxChemError,
                           match='fika embedding is not supported'):
            driver.compute(molecule, basis, gas_results)
