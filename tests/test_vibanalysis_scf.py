from contextlib import redirect_stdout
import io
from pathlib import Path
import sys
import types
import numpy as np
import h5py
import pytest

from mpi4py import MPI
from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.vibrationalanalysis import VibrationalAnalysis


class TestScfVibrationalAnalysisDriver:

    class RecordingOutput:

        def __init__(self):
            self.headers = []
            self.infos = []
            self.warnings = []
            self.blank_count = 0
            self.flushed = False

        def print_header(self, text):
            self.headers.append(text)

        def print_info(self, text):
            self.infos.append(text)

        def print_warning(self, text):
            self.warnings.append(text)

        def print_blank(self):
            self.blank_count += 1

        def flush(self):
            self.flushed = True

    class FakePy3DmolView:

        def __init__(self, width, height):
            self.width = width
            self.height = height
            self.model = None
            self.view_style = None
            self.style = None
            self.animation = None
            self.rotation = None
            self.zoomed = False
            self.shown = False

        def addModel(self, model, fmt, options):
            self.model = (model, fmt, options)

        def setViewStyle(self, style):
            self.view_style = style

        def setStyle(self, style):
            self.style = style

        def animate(self, options):
            self.animation = options

        def rotate(self, angle, axis):
            self.rotation = (angle, axis)

        def zoomTo(self):
            self.zoomed = True

        def show(self):
            self.shown = True

    @staticmethod
    def _get_water_molecule():

        return Molecule.read_xyz_string("""3
        water
        O    0.000000    0.000000    0.000000
        H    0.000000   -0.757160    0.586260
        H    0.000000    0.757160    0.586260
        """)

    def _get_synthetic_vibanalysis(self, base_name=None):

        vib_drv = VibrationalAnalysis(ScfRestrictedDriver())
        vib_drv.ostream = self.RecordingOutput()

        vib_drv.update_settings(
            method_dict={'xcfun': 'hf'},
            vib_dict={
                'do_ir': 'yes',
                'do_raman': 'yes',
                'filename': base_name,
                'print_depolarization_ratio': 'yes',
            },
            hessian_dict={'numerical': 'no'},
            cphf_dict={'conv_thresh': '1.0e-6'},
            rsp_dict={'conv_thresh': '1.0e-5'},
            polgrad_dict={},
        )

        natm = 3
        nmodes = 6
        mode_grid = np.arange(1.0, natm * 3 * nmodes + 1.0).reshape(
            nmodes, natm * 3)

        vib_drv.hessian = np.diag(np.linspace(0.1, 0.9, natm * 3))
        vib_drv.raw_normal_modes = mode_grid
        vib_drv.normal_modes = mode_grid / np.linalg.norm(mode_grid,
                                                          axis=1)[:, np.newaxis]
        vib_drv.vib_frequencies = np.array(
            [120.0, 240.0, 360.0, 480.0, 600.0, 720.0])
        vib_drv.reduced_masses = np.array([1.1, 1.2, 1.3, 1.4, 1.5, 1.6])
        vib_drv.force_constants = np.array([0.5, 0.7, 0.9, 1.1, 1.3, 1.5])
        vib_drv.dipole_gradient = np.arange(27.0).reshape(3, natm * 3) / 10.0
        vib_drv.ir_intensities = np.array([4.0, 11.0, 9.0, 15.0, 13.0, 7.0])
        vib_drv.free_energy_summary = 'Gibbs free energy summary'
        vib_drv.elec_energy = -76.0
        vib_drv.frequencies = (0.0, 0.2)
        vib_drv.print_depolarization_ratio = True

        base_gradient = np.arange(1.0, 82.0).reshape(3, 3, natm * 3) / 50.0
        vib_drv.polarizability_gradient = {
            0.0: base_gradient,
            0.2: base_gradient * 1.5,
        }
        (vib_drv.raman_activities, vib_drv.int_pol, vib_drv.int_depol,
         vib_drv.depol_ratio) = vib_drv.calculate_raman_activity(
             vib_drv.raw_normal_modes)

        return vib_drv

    @pytest.mark.timeconsuming
    def test_scf_vibrational_analysis_driver_numerical(self):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'water_hessian_scf.inp')
        h5file = str(here / 'data' / 'water_vib_scf.h5')

        task = MpiTask([inpfile, None])

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.acc_type = 'l2_c2diis'
        scf_drv.compute(task.molecule, task.ao_basis)

        vib_settings = {
            'do_ir': 'yes',
            'do_raman': 'yes',
            'numerical_hessian': 'yes',
            'numerical_raman': 'yes'
        }
        method_settings = {}
        vibanalysis_drv = VibrationalAnalysis(scf_drv)
        vibanalysis_drv.update_settings(method_settings, vib_settings)
        vibanalysis_drv.ostream.mute()
        vibanalysis_drv.compute(task.molecule, task.ao_basis)

        if task.mpi_rank == mpi_master():

            hf = h5py.File(h5file)
            ref_hessian = np.array(hf.get('hessian'))
            ref_frequencies = np.array(hf.get('frequencies'))
            ref_ir_intensities = np.array(hf.get('ir'))
            ref_raman_activities = np.array(hf.get('raman'))
            hf.close()

            diff_hessian = np.max(np.abs(vibanalysis_drv.hessian - ref_hessian))
            rel_diff_freq = np.max(
                np.abs(vibanalysis_drv.vib_frequencies / ref_frequencies - 1.0))
            rel_diff_ir = np.max(
                np.abs(vibanalysis_drv.ir_intensities / ref_ir_intensities -
                       1.0))
            rel_diff_raman = np.max(
                np.abs(vibanalysis_drv.raman_activities[0] /
                       ref_raman_activities - 1.0))

            assert diff_hessian < 1.0e-5
            assert rel_diff_freq < 1.0e-3
            assert rel_diff_ir < 1.0e-3
            assert rel_diff_raman < 1.0e-3

        task.finish()

    @pytest.mark.timeconsuming
    def test_scf_resonance_raman_numerical(self):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'water_hessian_scf.inp')
        h5file = str(here / 'data' / 'water_vib_scf.h5')

        task = MpiTask([inpfile, None])
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)

        scf_drv.compute(task.molecule, task.ao_basis)

        vib_settings = {
            'do_ir': 'no',
            'do_resonance_raman': 'yes',
            'numerical_hessian': 'yes',
            'numerical_raman': 'yes',
            'frequencies': (0.4,),
            'rr_damping': 0.05
        }
        method_settings = {}
        vibanalysis_drv = VibrationalAnalysis(scf_drv)
        vibanalysis_drv.update_settings(method_settings, vib_settings)
        vibanalysis_drv.ostream.mute()
        vibanalysis_drv.compute(task.molecule, task.ao_basis)

        if task.mpi_rank == mpi_master():

            hf = h5py.File(h5file)
            ref_frequencies = np.array(hf.get('frequencies'))
            ref_raman_activities = np.array(hf.get('resonance_raman'))
            hf.close()

            rel_diff_freq = np.max(
                np.abs(vibanalysis_drv.vib_frequencies / ref_frequencies - 1.0))
            rel_diff_raman_dyn = np.max(
                np.abs(vibanalysis_drv.raman_activities[0] /
                       ref_raman_activities - 1.0))

            assert rel_diff_freq < 1.0e-3
            assert rel_diff_raman_dyn < 1.0e-3

    def test_vibrational_analysis_writes_outputs_and_hdf5(self, tmp_path):

        molecule = self._get_water_molecule()
        basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)
        synced_base_name = MPI.COMM_WORLD.bcast(
            str(tmp_path / 'synthetic-vib')
            if MPI.COMM_WORLD.Get_rank() == mpi_master() else None,
            root=mpi_master())
        direct_h5_name = MPI.COMM_WORLD.bcast(
            str(tmp_path / 'synthetic-direct.h5')
            if MPI.COMM_WORLD.Get_rank() == mpi_master() else None,
            root=mpi_master())
        wrapped_h5_name = MPI.COMM_WORLD.bcast(
            str(tmp_path / 'synthetic-vib.h5')
            if MPI.COMM_WORLD.Get_rank() == mpi_master() else None,
            root=mpi_master())

        vib_drv = self._get_synthetic_vibanalysis(synced_base_name)

        assert vib_drv.vib_results_txt_file == f'{synced_base_name}-vib-results.out'
        assert vib_drv.method_dict == {'xcfun': 'hf'}
        assert vib_drv.hessian_dict == {'numerical': 'no'}
        assert vib_drv.cphf_dict == {'conv_thresh': '1.0e-6'}
        assert vib_drv.rsp_dict == {'conv_thresh': '1.0e-5'}

        dominant_modes = vib_drv.get_dominant_modes(3)
        assert dominant_modes == [1, 3, 4]

        vib_drv.print_header()
        assert any(
            'Resonance Raman' not in line for line in vib_drv.ostream.headers)

        vib_drv.print_vibrational_analysis(molecule)

        txt_file = Path(vib_drv.vib_results_txt_file)
        assert txt_file.is_file()
        txt_text = txt_file.read_text()
        assert 'Vibrational Mode' in txt_text
        assert 'Raman activity:' in txt_text
        assert 'Depol. ratio:' in txt_text
        assert any(
            'dominant normal modes' in info for info in vib_drv.ostream.infos)

        if vib_drv.rank == mpi_master():
            direct_h5_file = Path(direct_h5_name)
            direct_h5_file.touch()
            vib_drv.write_vib_results_to_hdf5(molecule, str(direct_h5_file),
                                              basis)

            wrapped_h5_file = Path(wrapped_h5_name)
            wrapped_h5_file.touch()
            vib_drv._write_final_hdf5(molecule, basis)

            with h5py.File(direct_h5_file, 'r') as hf:
                assert hf['vib/number_of_modes'][0] == len(
                    vib_drv.vib_frequencies)
                assert hf['vib/hessian'].shape == vib_drv.hessian.shape
                assert hf['vib/normal_modes'].shape == (6, 3, 3)
                assert hf['vib/external_frequencies'].shape == (2,)
                assert hf[
                    'vib/raman_activities'].shape == vib_drv.raman_activities.shape
                assert hf['vib/polarizability_gradient'].shape == (2, 3, 3, 9)
                assert hf['vib/raman_type'][0] == b'normal'

            assert wrapped_h5_file.is_file()

        vib_drv.do_resonance_raman = True
        vib_drv.do_raman = False
        vib_drv.print_resonance_raman()
        vib_drv.print_resonance_raman_to_file()

        txt_text = txt_file.read_text()
        assert 'Resonance Raman' in txt_text
        assert 'Frequency' in txt_text

    def test_vibrational_analysis_plot_and_print_info_helpers(
            self, monkeypatch, capsys):

        plt = pytest.importorskip('matplotlib.pyplot')

        vib_drv = self._get_synthetic_vibanalysis()
        molecule = self._get_water_molecule()
        vib_results = {
            'molecule_xyz_string': molecule.get_xyz_string(),
            'normal_modes': vib_drv.normal_modes,
            'vib_frequencies': vib_drv.vib_frequencies,
            'ir_intensities': vib_drv.ir_intensities,
            'raman_activities': vib_drv.raman_activities,
            'external_frequencies': vib_drv.frequencies,
        }

        monkeypatch.setattr(plt, 'show', lambda: None)

        vib_drv.plot_ir(vib_results,
                        broadening_type='gaussian',
                        broadening_value=12,
                        scaling_factor=0.98,
                        invert_axes=True)
        vib_drv.plot_raman(vib_results,
                           broadening_type='lorentzian',
                           broadening_value=18,
                           scaling_factor=1.02,
                           invert_axes=True)
        vib_drv.plot(vib_results,
                     plot_type='vibrational',
                     broadening_type='gaussian',
                     broadening_value=10)

        vib_drv.print_info(vib_results, info_type='free_energy')
        vib_drv.print_info(vib_results, info_type='ir')
        vib_drv.print_info(vib_results, info_type='raman')
        vib_drv.print_info(vib_results, info_type='vibrational')

        output = capsys.readouterr().out
        assert 'Gibbs free energy summary' in output
        assert '3 modes with the highest IR intensity:' in output
        assert '3 modes with the highest Raman activity:' in output
        assert 'Number of normal modes: 6' in output

        xvals = np.array([100.0, 300.0])
        yvals = np.array([2.0, 5.0])
        lorentz_x, lorentz_y = vib_drv.lorentzian_broadening(
            xvals, yvals, 0, 500, 50, 20)
        gauss_x, gauss_y = vib_drv.gaussian_broadening(xvals, yvals, 0, 500, 50,
                                                       20)
        assert lorentz_x.shape == lorentz_y.shape
        assert gauss_x.shape == gauss_y.shape
        assert np.max(lorentz_y) > 0.0
        assert np.max(gauss_y) > 0.0

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_vibrational_analysis_validation_and_animation_paths(
            self, monkeypatch):

        vib_drv = self._get_synthetic_vibanalysis()
        molecule = self._get_water_molecule()

        empty_hessian_drv = VibrationalAnalysis(ScfRestrictedDriver())
        empty_hessian_drv.ostream = self.RecordingOutput()
        empty_hessian_drv.print_hessian(molecule)
        assert empty_hessian_drv.ostream.warnings == [
            'Hessian is not available.'
        ]

        fake_view = self.FakePy3DmolView(width=420, height=320)
        fake_module = types.SimpleNamespace(
            view=lambda width, height: fake_view)
        monkeypatch.setitem(sys.modules, 'py3Dmol', fake_module)

        vib_results = {
            'molecule_xyz_string': molecule.get_xyz_string(),
            'normal_modes': {
                '1': vib_drv.normal_modes[0]
            },
        }
        vib_drv.animate(vib_results,
                        mode=1,
                        frames=9,
                        amplitude=0.3,
                        width=420,
                        height=320)

        assert fake_view.model is not None
        assert fake_view.model[1] == 'xyz'
        assert fake_view.animation == {'loop': 'backAndForth'}
        assert fake_view.rotation == (-90, 'x')
        assert fake_view.zoomed is True
        assert fake_view.shown is True

        with redirect_stdout(io.StringIO()):
            with pytest.raises(AssertionError,
                               match='No IR intensities available'):
                vib_drv.print_info(
                    {
                        'molecule_xyz_string': molecule.get_xyz_string(),
                        'normal_modes': vib_drv.normal_modes,
                        'vib_frequencies': vib_drv.vib_frequencies,
                        'ir_intensities': None,
                        'raman_activities': vib_drv.raman_activities,
                    },
                    info_type='ir')

        with redirect_stdout(io.StringIO()):
            with pytest.raises(AssertionError,
                               match='No Raman activities available'):
                vib_drv.print_info(
                    {
                        'molecule_xyz_string': molecule.get_xyz_string(),
                        'normal_modes': vib_drv.normal_modes,
                        'vib_frequencies': vib_drv.vib_frequencies,
                        'ir_intensities': vib_drv.ir_intensities,
                        'raman_activities': None,
                    },
                    info_type='raman')

        with redirect_stdout(io.StringIO()):
            with pytest.raises(AssertionError, match='Invalid plot type'):
                vib_drv.print_info(
                    {
                        'molecule_xyz_string': molecule.get_xyz_string(),
                        'normal_modes': vib_drv.normal_modes,
                        'vib_frequencies': vib_drv.vib_frequencies,
                        'ir_intensities': vib_drv.ir_intensities,
                        'raman_activities': vib_drv.raman_activities,
                    },
                    info_type='unsupported')

        with pytest.raises(AssertionError,
                           match='molecule only has 1 normal modes'):
            vib_drv.animate(vib_results, mode=2)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_vibrational_analysis_plot_rejects_invalid_type(self):

        pytest.importorskip('matplotlib.pyplot')

        vib_drv = self._get_synthetic_vibanalysis()
        molecule = self._get_water_molecule()

        with pytest.raises(AssertionError, match='Invalid plot type'):
            vib_drv.plot(
                {
                    'molecule_xyz_string': molecule.get_xyz_string(),
                    'normal_modes': vib_drv.normal_modes,
                    'vib_frequencies': vib_drv.vib_frequencies,
                    'ir_intensities': vib_drv.ir_intensities,
                },
                plot_type='unsupported')

    @pytest.mark.solvers
    def test_scf_vibrational_analysis_ir(self):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'water_hessian_scf.inp')
        h5file = str(here / 'data' / 'water_vib_scf.h5')

        task = MpiTask([inpfile, None])
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)

        scf_drv.compute(task.molecule, task.ao_basis)

        vib_settings = {
            'do_ir': 'yes',
            'do_raman': 'no',
            'numerical_hessian': 'no',
            'numerical_raman': 'no'
        }
        method_settings = {}
        vibanalysis_drv = VibrationalAnalysis(scf_drv)
        vibanalysis_drv.update_settings(method_settings, vib_settings)
        vibanalysis_drv.ostream.mute()
        vibanalysis_drv.compute(task.molecule, task.ao_basis)

        if task.mpi_rank == mpi_master():

            hf = h5py.File(h5file)
            ref_hessian = np.array(hf.get('hessian'))
            ref_frequencies = np.array(hf.get('frequencies'))
            ref_ir_intensities = np.array(hf.get('ir'))
            hf.close()

            diff_hessian = np.max(np.abs(vibanalysis_drv.hessian - ref_hessian))
            rel_diff_freq = np.max(
                np.abs(vibanalysis_drv.vib_frequencies / ref_frequencies - 1.0))
            rel_diff_ir = np.max(
                np.abs(vibanalysis_drv.ir_intensities / ref_ir_intensities -
                       1.0))

            assert diff_hessian < 1.0e-5
            assert rel_diff_freq < 1.0e-3
            assert rel_diff_ir < 1.0e-3

        task.finish()

    def test_scf_vibrational_analysis_raman(self):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'water_hessian_scf.inp')
        h5file = str(here / 'data' / 'water_vib_scf.h5')

        task = MpiTask([inpfile, None])
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)

        scf_drv.compute(task.molecule, task.ao_basis)

        vib_settings = {
            'do_ir': 'yes',
            'do_raman': 'yes',
            'numerical_hessian': 'no',
            'numerical_raman': 'no'
        }
        method_settings = {}
        vibanalysis_drv = VibrationalAnalysis(scf_drv)
        vibanalysis_drv.update_settings(method_settings, vib_settings)
        vibanalysis_drv.ostream.mute()
        vibanalysis_drv.compute(task.molecule, task.ao_basis)

        if task.mpi_rank == mpi_master():

            hf = h5py.File(h5file)
            ref_hessian = np.array(hf.get('hessian'))
            ref_frequencies = np.array(hf.get('frequencies'))
            ref_ir_intensities = np.array(hf.get('ir'))
            ref_raman_activities = np.array(hf.get('raman'))
            hf.close()

            diff_hessian = np.max(np.abs(vibanalysis_drv.hessian - ref_hessian))
            rel_diff_freq = np.max(
                np.abs(vibanalysis_drv.vib_frequencies / ref_frequencies - 1.0))
            rel_diff_ir = np.max(
                np.abs(vibanalysis_drv.ir_intensities / ref_ir_intensities -
                       1.0))
            rel_diff_raman = np.max(
                np.abs(vibanalysis_drv.raman_activities[0] /
                       ref_raman_activities - 1.0))

            assert diff_hessian < 1.0e-5
            assert rel_diff_freq < 1.0e-3
            assert rel_diff_ir < 1.0e-3
            assert rel_diff_raman < 1.0e-3

        task.finish()

    def test_scf_resonance_raman_analytical(self):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'water_hessian_scf.inp')
        h5file = str(here / 'data' / 'water_vib_scf.h5')

        task = MpiTask([inpfile, None])
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)

        scf_drv.compute(task.molecule, task.ao_basis)

        vib_settings = {
            'do_ir': 'no',
            'do_resonance_raman': 'yes',
            'numerical_hessian': 'no',
            'numerical_raman': 'no',
            'frequencies': (0.4,),
            'rr_damping': 0.05
        }
        method_settings = {}
        vibanalysis_drv = VibrationalAnalysis(scf_drv)
        vibanalysis_drv.update_settings(method_settings, vib_settings)
        vibanalysis_drv.ostream.mute()
        vibanalysis_drv.compute(task.molecule, task.ao_basis)

        if task.mpi_rank == mpi_master():

            hf = h5py.File(h5file)
            ref_frequencies = np.array(hf.get('frequencies'))
            ref_raman_activities = np.array(hf.get('resonance_raman'))
            hf.close()

            rel_diff_freq = np.max(
                np.abs(vibanalysis_drv.vib_frequencies / ref_frequencies - 1.0))
            rel_diff_raman_dyn = np.max(
                np.abs(vibanalysis_drv.raman_activities[0] /
                       ref_raman_activities - 1.0))

            assert rel_diff_freq < 1.0e-3
            assert rel_diff_raman_dyn < 1.0e-3
