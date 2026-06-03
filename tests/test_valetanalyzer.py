import io
from contextlib import redirect_stdout
import numpy as np
import pytest
from mpi4py import MPI

import veloxchem.valetanalyzer as valet_module
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.valetanalyzer import ValetAnalyzer
from veloxchem.veloxchemlib import mpi_master
from veloxchem.errorhandler import VeloxChemError

skip_multi_rank = pytest.mark.skipif(
    MPI.COMM_WORLD.Get_size() > 1,
    reason='real ValetAnalyzer coverage tests run only in serial')


@pytest.mark.solvers
class TestValetAnalyzer:

    _water_lr_case = None

    def _get_h2_molecule_and_basis(self):

        molecule = Molecule.read_molecule_string("""
            H  0.0  0.0  0.0
            H  0.0  0.0  1.4
        """,
                                                 units='au')
        basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)

        return molecule, basis

    def _initialize_silently(self, analyzer, molecule, basis, subgroups=None):

        with redirect_stdout(io.StringIO()):
            analyzer.initialize(molecule, basis, subgroups)

    @classmethod
    def _get_water_lr_case(cls):

        if cls._water_lr_case is None:
            xyz_string = """3
            xyz
            O   -0.1858140  -1.1749469   0.7662596
            H   -0.1285513  -0.8984365   1.6808606
            H   -0.0582782  -0.3702550   0.2638279
            """
            molecule = Molecule.read_xyz_string(xyz_string)
            basis = MolecularBasis.read(molecule, '6-31g', ostream=None)

            scf_drv = ScfRestrictedDriver()
            scf_drv.ostream.mute()
            scf_results = scf_drv.compute(molecule, basis)

            lr_drv = LinearResponseEigenSolver()
            lr_drv.ostream.mute()
            lr_drv.nstates = 3
            lr_drv.detach_attach = True
            lr_drv.detach_attach_charges = True
            rsp_results = lr_drv.compute(molecule, basis, scf_results)

            cls._water_lr_case = (
                molecule,
                basis,
                scf_results,
                rsp_results,
                lr_drv,
            )

        return cls._water_lr_case

    def test_initialize_set_subgroups(self, capsys):

        analyzer = ValetAnalyzer()
        molecule, basis = self._get_h2_molecule_and_basis()

        analyzer.initialize(molecule, basis, [('pair', [1, 2])])
        init_output = capsys.readouterr().out

        assert analyzer._molecule is molecule
        assert analyzer._basis is basis
        assert analyzer._subgroups == [('pair', [1, 2])]
        assert isinstance(analyzer._density_viewer, valet_module.DensityViewer)
        assert 'Molecule: 2 atoms' in init_output
        assert f'Basis: {basis.get_label()}' in init_output

        analyzer.set_subgroups([('atom1', 1)])
        subgroup_output = capsys.readouterr().out

        assert analyzer._subgroups == [('atom1', 1)]
        assert "Subgroups updated: ['atom1']" in subgroup_output

    @skip_multi_rank
    def test_compute_detach_attach_densities_validates_state_and_shell(self):

        analyzer = ValetAnalyzer()
        molecule, basis = self._get_h2_molecule_and_basis()

        scf_results = {
            'scf_type': 'restricted',
            'C_alpha': np.eye(2),
        }
        rsp_results = {
            'eigenvalues': np.array([1.0]),
            'eigenvectors': np.array([[1.0]]),
        }

        with pytest.raises(VeloxChemError,
                           match='ValetAnalyzer is not yet initialized'):
            analyzer.compute_detach_attach_densities(scf_results, rsp_results)

        self._initialize_silently(analyzer, molecule, basis)

        with pytest.raises(VeloxChemError,
                           match='does not yet support open-shell'):
            analyzer.compute_detach_attach_densities(
                {
                    'scf_type': 'unrestricted',
                    'C_alpha': np.eye(2),
                }, rsp_results)

        with pytest.raises(VeloxChemError, match='Invalid state_index'):
            analyzer.compute_detach_attach_densities(scf_results,
                                                     rsp_results,
                                                     state_index=2)

    @skip_multi_rank
    def test_real_detach_attach_density_and_charge_paths(self):

        molecule, basis, scf_results, rsp_results, _ = self._get_water_lr_case()
        analyzer = ValetAnalyzer()
        self._initialize_silently(analyzer, molecule, basis)

        densities = analyzer.compute_detach_attach_densities(scf_results,
                                                             rsp_results,
                                                             state_index=1)
        detachment = densities['detachment_density_matrix_AO']
        attachment = densities['attachment_density_matrix_AO']

        assert detachment.shape == attachment.shape == scf_results['S'].shape
        assert np.allclose(detachment, detachment.T)
        assert np.allclose(attachment, attachment.T)

        charges = analyzer.compute_atom_charges(scf_results,
                                                rsp_results,
                                                state_index=1)
        ref_detachment = np.array([
            0.9998357388039689,
            0.00008213057844294671,
            0.00008213061758819769,
        ])
        ref_attachment = np.array([
            0.28456175730239014,
            0.35771907934092895,
            0.35771916335668086,
        ])

        assert np.allclose(charges['detachment_charges'], ref_detachment)
        assert np.allclose(charges['attachment_charges'], ref_attachment)
        assert np.isclose(np.sum(charges['detachment_charges']), 1.0)
        assert np.isclose(np.sum(charges['attachment_charges']), 1.0)

    @skip_multi_rank
    def test_real_transition_data_uses_grouped_atomic_charges(self):

        molecule, basis, scf_results, rsp_results, _ = self._get_water_lr_case()
        analyzer = ValetAnalyzer()
        self._initialize_silently(analyzer, molecule, basis)

        charges = analyzer.compute_atom_charges(scf_results,
                                                rsp_results,
                                                state_index=1)
        transition_data = analyzer.compute_transition_data(
            scf_results,
            rsp_results, [('oxygen', 1), ('hydrogens', [2, 3])],
            state_index=1)

        expected_detachment = [
            0.0,
            charges['detachment_charges'][0],
            charges['detachment_charges'][1] + charges['detachment_charges'][2],
        ]
        expected_attachment = [
            0.0,
            charges['attachment_charges'][0],
            charges['attachment_charges'][1] + charges['attachment_charges'][2],
        ]

        assert transition_data['subgroup_names'] == [
            'Unassigned', 'oxygen', 'hydrogens'
        ]
        assert transition_data['subgroup_colors'][0] == '#727272'
        assert transition_data['atom_to_subgroup_map'] == {0: 1, 1: 2, 2: 2}
        assert np.allclose(transition_data['detachment_charges'],
                           expected_detachment)
        assert np.allclose(transition_data['attachment_charges'],
                           expected_attachment)
        assert np.allclose(
            transition_data['transition_matrix'],
            ValetAnalyzer.compute_transition_matrix(expected_detachment,
                                                    expected_attachment))

    @skip_multi_rank
    def test_compute_transition_data_validates_subgroup_atom_indices(self):

        molecule, basis, scf_results, rsp_results, _ = self._get_water_lr_case()
        analyzer = ValetAnalyzer()
        self._initialize_silently(analyzer, molecule, basis)

        with pytest.raises(VeloxChemError,
                           match='subgroup atom indices must be in the range'):
            analyzer.compute_transition_data(scf_results,
                                             rsp_results,
                                             [('bad', 0)],
                                             state_index=1)

        with pytest.raises(VeloxChemError,
                           match='subgroup atom indices must be in the range'):
            analyzer.compute_transition_data(scf_results,
                                             rsp_results,
                                             [('bad', [1, 4])],
                                             state_index=1)

    def test_compute_transition_matrix(self):

        matrix = ValetAnalyzer.compute_transition_matrix([0.7, 0.3], [0.2, 0.8])

        assert np.allclose(matrix, [[0.2, 0.0], [0.5, 0.3]])

    def test_plot_transition_diagram_renders_bars_connectors_and_title(self):

        pytest.importorskip('matplotlib.pyplot')

        analyzer = ValetAnalyzer()
        transition_data = {
            'subgroup_names': ['Unassigned', 'Donor', 'Acceptor'],
            'subgroup_colors': ['#727272', '#123456', '#abcdef'],
            'detachment_charges': [0.0, 0.7, 0.3],
            'attachment_charges': [0.0, 0.2, 0.8],
            'transition_matrix': [[0.0, 0.0, 0.0], [0.0, 0.2, 0.0],
                                  [0.0, 0.5, 0.3]],
            'atom_to_subgroup_map': {
                0: 1,
                1: 2
            },
        }

        fig, ax = analyzer.plot_transition_diagram(transition_data,
                                                   title='State 1',
                                                   width=400,
                                                   height=300,
                                                   dpi=100)

        texts = [text.get_text() for text in ax.texts]

        assert ax.get_xlim() == (0.0, 400.0)
        assert ax.get_ylim() == (0.0, 300.0)
        assert len(ax.patches) >= 5
        assert 'State 1' in texts
        assert '50.00%' in texts
        assert 'Donor' in texts

        valet_module.plt.close(fig)
