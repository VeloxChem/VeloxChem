from pathlib import Path
from unittest.mock import patch
from dataclasses import dataclass
from pytest import approx
import numpy as np
import pytest

from veloxchem.veloxchemlib import (ElectronRepulsionIntegralsDriver,
                                    GridDriver, is_single_node, denmat)
from veloxchem.aodensitymatrix import AODensityMatrix
from veloxchem.aofockmatrix import AOFockMatrix
from veloxchem.main import main, select_scf_driver
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestopendriver import ScfRestrictedOpenDriver
from veloxchem.qqscheme import get_qq_scheme
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.molecule import Molecule


@patch('veloxchem.main.ScfRestrictedOpenDriver')
@patch('veloxchem.main.MpiTask')
@patch('veloxchem.cli.argparse')
class TestROSetup:

    @pytest.fixture
    def input_dict(self):
        return {
            'jobs': {
                'task': 'roscf',
            },
            'method_settings': {
                'basis': 'DEF2-SVP',
            },
            'filename': None,
        }

    def test_scftype(self, mock_argparse, mock_mpi, mock_roscf, input_dict,
                     tmpdir):

        task = mock_mpi()
        task.input_dict = input_dict

        with patch('veloxchem.main.select_scf_driver') as mock_select:
            mock_select().is_converged = False
            main()
        mock_select.assert_called_with(task, 'restricted_openshell')

    def test_roscfdriverreturned(self, mock_argparse, mock_mpi, mock_roscf,
                                 input_dict, tmpdir):

        task = mock_mpi()
        task.input_dict = input_dict

        assert select_scf_driver(task, 'restricted_openshell') is mock_roscf()


@dataclass
class ROSCF_Helper:
    scf_drv: ScfRestrictedOpenDriver
    mol: Molecule
    bas: MolecularBasis
    mats: dict

    def __iter__(self):
        i = 0
        while not self.scf_drv.is_converged:
            e, gn = self.energy(), self.gradient()
            self.scf_drv.store_diis_data(0, self.mats['fock_mat'],
                                         self.mats['ao_density'])
            self.scf_drv.add_iter_data({'energy': e, 'gradient': gn})
            yield e, gn
            self.new_orbitals()
            self.new_density()
            i += 1

    def energy(self):
        e1 = self.e1()
        e2 = self.e2()
        e_Z = self.nuclear()
        return e1 + e2 + e_Z

    def e1(self):
        density = sum(self.ao_density())
        h1 = self.h1()
        return np.trace(np.matmul(h1, density))

    def ao_density(self):
        da = self.mats['ao_density'].alpha_to_numpy(0)
        db = self.mats['ao_density'].beta_to_numpy(0)
        return da, db

    def e2(self):
        self.comp_2e_fock()
        return self.comp_2e_energy()

    def comp_2e_fock(self):

        eri_drv = ElectronRepulsionIntegralsDriver(self.scf_drv.comm)
        qq_data = eri_drv.compute(get_qq_scheme(self.scf_drv.qq_type),
                                  self.scf_drv.eri_thresh, self.mol, self.bas)

        e_grad = None

        fock_mat = AOFockMatrix(self.mats['ao_density'])
        if self.scf_drv.dft:
            self.scf_drv.update_fock_type(fock_mat)

        vxc_mat, e_pe, V_pe = self.scf_drv.comp_2e_fock(fock_mat,
                                                        self.mats['ao_density'],
                                                        self.mol, self.bas,
                                                        qq_data, e_grad)

        self.mats['fock_mat'] = fock_mat
        self.mats['vxc_mat'] = vxc_mat
        self.mats['e_pe'] = e_pe
        self.mats['V_pe'] = V_pe

    def comp_2e_energy(self):
        e_ee = self.comp_ee_energy()
        e_xc = self.comp_xc_energy()
        return e_ee + e_xc

    def comp_ee_energy(self):
        e_ee = self.mats['fock_mat'].get_energy(0, self.mats['ao_density'], 0)
        return e_ee

    def comp_xc_energy(self):
        try:
            e_xc = self.mats['vxc_mat'].get_energy()
        except AttributeError:
            e_xc = 0.0

        return e_xc

    def gradient(self):
        self.e2()
        mats = (
            self.mats[k]
            for k in ['fock_mat', 'vxc_mat', 'V_pe', 'kinetic', 'potential'])
        self.scf_drv.comp_full_fock(*mats)

        oao_mat = self.oao_mat()
        self.scf_drv.scf_tensors = {'S': self.overlap()}

        e_grad, max_grad = self.scf_drv.comp_gradient(self.mats['fock_mat'],
                                                      self.mats['overlap'],
                                                      self.mats['ao_density'],
                                                      oao_mat)
        return e_grad

    def overlap(self):
        return self.mats['overlap'].to_numpy()

    def nuclear(self):
        return self.mol.nuclear_repulsion_energy()

    def h1(self):
        t = self.mats['kinetic'].to_numpy()
        v = self.mats['potential'].to_numpy()
        return t - v

    def oao_mat(self):
        return self.mats['overlap'].get_ortho_matrix(self.scf_drv.ovl_thresh)

    def fock_eff(self):
        self.gradient()  # sets fock_mat to full fock

        fock_mat = self.mats['fock_mat']
        fa = fock_mat.alpha_to_numpy(0)
        fb = fock_mat.beta_to_numpy(0)

        da, db = self.ao_density()
        s_mat = self.mats['overlap'].to_numpy()
        f_proj = self.scf_drv.get_projected_fock(fa, fb, da, db, s_mat)
        return f_proj

    def new_orbitals(self):
        fock_eff = self.fock_eff()
        oao = self.oao_mat()
        mol_orbs = self.scf_drv.gen_molecular_orbitals(fock_eff, oao)
        self.scf_drv.mol_orbs = mol_orbs
        return mol_orbs.alpha_to_numpy()

    def new_density(self):
        mol_orbs = self.scf_drv.mol_orbs
        mo = mol_orbs.alpha_to_numpy()
        na = 2
        nb = 1
        da = np.matmul(mo[:, :na], mo[:, :na].T)
        db = np.matmul(mo[:, :nb], mo[:, :nb].T)
        ao_density = AODensityMatrix([da, db], denmat.unrest)
        self.mats['ao_density'] = ao_density


@pytest.mark.skipif(not is_single_node(),
                    reason='This test only runs on single node')
class TestROSCF:

    @pytest.fixture
    def initial_density(self):
        Da = np.array([[1.017257008959035147e+00, -1.324945784447343899e-01],
                       [-1.324945784447343899e-01, 1.017257008959035369e+00]])
        Db = np.array([[1.004974038596785801e+00, -2.073733392173061560e-02],
                       [-2.073733392173061560e-02, 4.279085843668227365e-04]])
        ao_density = AODensityMatrix([Da, Db], denmat.unrest)
        return ao_density

    @pytest.fixture
    def scf_setup(self, initial_density):
        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'heh.inp')
        task = MpiTask([inpfile, None])
        task.input_dict['scf']['checkpoint_file'] = None

        scf_drv = ScfRestrictedOpenDriver(task.mpi_comm, task.ostream)

        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])

        mol = task.molecule
        bas = task.ao_basis

        S, T, V, _ = scf_drv.comp_one_ints(mol, bas)

        scf_drv.comp_guess_density = lambda *args: initial_density
        ao_density = scf_drv.comp_guess_density()
        mats = {
            'ao_density': ao_density,
            'overlap': S,
            'kinetic': T,
            'potential': V,
        }

        return ROSCF_Helper(scf_drv, mol, bas, mats)

    @pytest.fixture
    def roothan_setup(self, scf_setup):
        scf_setup.scf_drv.max_err_vecs = 1
        return scf_setup

    @pytest.fixture
    def diis_setup(self, scf_setup):
        scf_setup.scf_drv.max_err_vecs = 2
        return scf_setup

    def test_initial_density(self, scf_setup):

        alpha_density, beta_density = scf_setup.ao_density()
        overlap = scf_setup.overlap()
        assert np.trace(np.matmul(overlap, alpha_density)) == approx(2.0)
        assert np.trace(np.matmul(overlap, beta_density)) == approx(1.0)

    def test_initial_one_el(self, scf_setup):
        assert scf_setup.e1() == approx(-5.486532987513215)

    def test_nuclear_potential(self, scf_setup):
        e_Z = scf_setup.nuclear()
        assert e_Z == approx(0.5601421515833619)

    def test_initial_two_el(self, scf_setup):
        e_ee = scf_setup.e2()
        assert e_ee == approx(1.580015911515688)

    def test_initial_total_energy(self, scf_setup):
        total_energy = scf_setup.energy()
        assert total_energy == approx(-3.3463749244141647)

    def test_initial_gradient(self, scf_setup):
        assert scf_setup.gradient() == approx(0.05872228054511208)

    def test_initial_fockeff(self, scf_setup):
        expected = [[-9.172419301346693699e-01, -1.335506295099394558e-01],
                    [-1.335506295099394281e-01, -1.936856073872699480e-01]]
        np.testing.assert_allclose(scf_setup.fock_eff(), expected)

    def test_new_molorbs(self, scf_setup):
        expected = [
            [-9.972950283347190581e-01, -1.505311775609534386e-01],
            [-1.935429905478982734e-02, 1.008405880619075656e+00],
        ]
        np.testing.assert_allclose(scf_setup.new_orbitals(), expected)

    def test_roothan(self, roothan_setup):

        scf_drv = roothan_setup.scf_drv
        mol = roothan_setup.mol
        bas = roothan_setup.bas

        scf_drv.restart = False
        scf_drv.acc_type = 'DIIS'
        scf_drv.max_err_vecs = 1
        scf_drv.conv_thresh = 1e-5

        scf_drv.compute(mol, bas)

        assert scf_drv.is_converged
        assert scf_drv.iter_data[-1]['energy'] == approx(-3.347480513475661)

    def test_step_1(self, roothan_setup):
        iscf = iter(roothan_setup)
        energy, gradient = next(iscf)
        assert energy == approx(-3.3463749244141647)
        assert gradient == approx(0.05872228054511208)

    def test_step_2(self, roothan_setup):
        iscf = iter(roothan_setup)
        energy, gradient = next(iscf)
        energy, gradient = next(iscf)
        assert energy == approx(-3.347477184661499)
        assert gradient == approx(0.0032244743264249445)

    def test_step_3(self, roothan_setup):
        iscf = iter(roothan_setup)
        for _ in range(3):
            energy, gradient = next(iscf)
        assert energy == approx(-3.347480505231504)
        assert gradient == approx(0.00016067000893129152)

    def test_step_4(self, roothan_setup):
        iscf = iter(roothan_setup)
        for _ in range(4):
            energy, gradient = next(iscf)
        assert energy == approx(-3.347480513475661)
        assert gradient == approx(8.066484831554014e-06)

    def test_diis_data_after_1(self, diis_setup):
        next(iter(diis_setup))
        assert len(diis_setup.scf_drv.fock_matrices) == 1

    def test_diis_data_after_2(self, diis_setup):
        iscf = iter(diis_setup)
        energy, gradient = next(iscf)
        energy, gradient = next(iscf)
        assert len(diis_setup.scf_drv.fock_matrices) == 2
        assert len(diis_setup.scf_drv.den_matrices) == 2

    def test_diis_2(self, scf_setup):

        scf_drv = scf_setup.scf_drv
        mol = scf_setup.mol
        bas = scf_setup.bas

        scf_drv.acc_type = 'DIIS'
        scf_drv.max_err_vecs = 2

        scf_drv.compute(mol, bas)

        assert scf_drv.is_converged
        assert scf_drv.iter_data[-1]['energy'] == approx(-3.347480513475661)

    def test_c2diis_5(self, scf_setup):

        scf_drv = scf_setup.scf_drv
        mol = scf_setup.mol
        bas = scf_setup.bas

        scf_drv.max_err_vecs = 5

        scf_drv.compute(mol, bas)

        assert scf_drv.is_converged
        assert scf_drv.iter_data[-1]['energy'] == approx(-3.347480513475661)


@pytest.mark.skipif(not is_single_node(),
                    reason='This test only runs on single node')
class TestRODFT:

    @pytest.fixture
    def initial_density(self):

        Da = np.array([[1.007149386323507922e+00, -8.485576025415243751e-02],
                       [-8.485576025415243751e-02, 1.007149386323507478e+00]])
        Db = np.array([[1.000587154488399744e+00, -3.559606189326910367e-03],
                       [-3.559606189326910367e-03, 1.266336087391910148e-05]])
        ao_density = AODensityMatrix([Da, Db], denmat.unrest)
        return ao_density

    @pytest.fixture
    def scf_setup(self, request, initial_density):

        try:
            xcfun_label, ref_energy = request.param
        except AttributeError:
            xcfun_label, ref_energy = None, 0.

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'heh.inp')
        task = MpiTask([inpfile, None])
        task.input_dict['scf']['checkpoint_file'] = None
        task.ao_basis = MolecularBasis.read(task.molecule, 'STO-3G', '.',
                                            task.ostream)

        scf_drv = ScfRestrictedOpenDriver(task.mpi_comm, task.ostream)
        scf_drv.grid_level = 6

        if xcfun_label:
            grid_drv = GridDriver(scf_drv.comm)
            grid_drv.set_level(scf_drv.grid_level)
            task.input_dict['method_settings']['xcfun'] = xcfun_label
            task.input_dict['method_settings']['grid'] = 6

        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])

        mol = task.molecule
        bas = task.ao_basis

        S, T, V, _ = scf_drv.comp_one_ints(mol, bas)

        if xcfun_label:
            scf_drv.molgrid = grid_drv.generate(mol)

        scf_drv.comp_guess_density = lambda *args: initial_density
        ao_density = scf_drv.comp_guess_density()
        mats = {
            'ao_density': ao_density,
            'overlap': S,
            'kinetic': T,
            'potential': V,
        }

        return ROSCF_Helper(scf_drv, mol, bas, mats), ref_energy

    def test_initial_density(self, scf_setup):

        scf, _ = scf_setup

        alpha_density, beta_density = scf.ao_density()
        overlap = scf.overlap()

        assert np.trace(np.matmul(overlap, alpha_density)) == approx(2.0)
        assert np.trace(np.matmul(overlap, beta_density)) == approx(1.0)

    @pytest.mark.parametrize(
        'scf_setup',
        [
            (None, -3.269190923863),
            ('slater', -3.064460757711),
            ('b3lyp', -3.315910830402),
        ],
        indirect=True,
        ids=['rohf', 'slater', 'b3lyp'],
    )
    def test_initial_total_energy(self, scf_setup):
        xcfun = scf_setup[0].scf_drv.xcfun
        if xcfun and xcfun.get_func_label() == "B3LYP":
            pytest.skip('missing accuracy')
        scf, ref_energy = scf_setup
        total_energy = scf.energy()
        assert total_energy == approx(ref_energy)

    @pytest.mark.parametrize(
        'scf_setup',
        [
            (None, -5.442239744622594),
            ('slater', -5.442239744622594),
        ],
        indirect=True,
        ids=['rohf', 'slater'],
    )
    def test_initial_oneel_energy(self, scf_setup):
        scf, ref_energy = scf_setup
        total_energy = scf.e1()
        assert total_energy == approx(ref_energy)

    @pytest.mark.parametrize(
        'scf_setup',
        [
            (None, 1.612906672206136),  # 1.0*K
            ('slater', 3.060388872497171),  # 0.0*K
            ('b3lyp', 2.7708924339156367),  # 0.2*K
        ],
        indirect=True,
        ids=['rohf', 'slater', 'b3lyp'],
    )
    def test_initial_ee_energy(self, scf_setup):
        scf, ref_energy = scf_setup
        scf.e2()
        e_ee = scf.comp_ee_energy()
        assert e_ee == approx(ref_energy)

    @pytest.mark.parametrize(
        'scf_setup',
        [('slater', -1.242752026437), ('slda', -1.43873045822),
         ('b3lyp', -1.20470566061)],
        indirect=True,
        ids=['slater', 'slda', 'b3lyp'],
    )
    def test_initial_xc_energy(self, scf_setup):
        xcfun = scf_setup[0].scf_drv.xcfun
        if xcfun and xcfun.get_func_label() == "B3LYP":
            pytest.skip('missing accuracy')
        scf, ref_energy = scf_setup
        scf.e2()
        e_xc = scf.comp_xc_energy()
        assert e_xc == approx(ref_energy)
