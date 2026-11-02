from mpi4py import MPI
from pathlib import Path
import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tdaeigensolver import TdaEigenSolver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.excitedstateanalysisdriver import ExcitedStateAnalysisDriver
from veloxchem.resultsio import read_molecule_and_basis


@pytest.mark.skipif(MPI.COMM_WORLD.Get_size() != 1,
                    reason='runs only on a single MPI rank')
class TestCTNumbers:

    @pytest.mark.parametrize('ri', [False, True])
    def test_tda(self, ri):

        ct_matrix = {
            False: np.array([
                [0.34837228, 0.19662888],
                [0.1367462, 0.31825265],
            ]),
            True: np.array([[0.34835735, 0.19661779], [0.13675409, 0.31827078]])
        }[ri]

        hole_participation_ratio = {
            False: 1.9839293483578253,
            True: 1.9839477783725588
        }[ri]
        particle_participation_ratio = {
            False: 1.9982298906088756,
            True: 1.9982282163446563
        }[ri]
        avg_participation_ratio = {
            False: 1.9910796194833504,
            True: 1.9910879973586075
        }[ri]

        avg_hole_position = {
            False: np.array([9.48439790e-02, 4.63701220e-02, -9.30375228e-05]),
            True: np.array([0.09474389, 0.04637279, -0.00009304])
        }[ri]
        avg_particle_position = {
            False: np.array([-1.05663691e-02, 1.95196319e-02, -2.61361122e-05]),
            True: np.array([-0.01058115, 0.01951824, -0.00002614])
        }[ri]
        avg_difference_vec = {
            False: np.array([-1.05410348e-01, -2.68504901e-02, 6.69014106e-05]),
            True: np.array([-0.10532504, -0.02685454, 0.00006690])
        }[ri]

        ct_length = {False: 0.10877635215676436, True: 0.10869469104546430}[ri]

        thienyl_thiazole_str = """S   1.5860   -1.4206    0.0001
        C   0.6775    0.0310   -0.0004
        C   1.4857    1.1466   -0.0004
        C   2.8690    0.8051   -0.0001
        C   3.0653   -0.5583    0.0002
        H   1.1367    2.1716   -0.0005
        H   3.6790    1.5228    0.0000
        H   4.0068   -1.0893    0.0006
        S  -1.6588    1.4476    0.0009
        N  -1.4203   -1.1264    0.0001
        C  -0.7409   -0.0017    0.0002
        C  -2.7697   -0.8319    0.0000
        C  -3.0939    0.5086   -0.0006
        H  -4.0731    0.9658   -0.0010
        H  -3.4860   -1.6422   -0.0003
        """

        fragment_dict = {
            "Donor": [1, 2, 3, 4, 5, 6, 7, 8],
            "Acceptor": [9, 10, 11, 12, 13, 14, 15],
        }

        molecule = Molecule.read_molecule_string(thienyl_thiazole_str)
        basis = MolecularBasis.read(molecule, "STO-3G", ostream=None)

        # SCF settings and calculation
        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        # Solve the TDDFT linear response equation for the first 2 excited states
        tda_drv = TdaEigenSolver()
        tda_drv.ostream.mute()
        tda_drv.nstates = 2
        tda_drv.initial_guess_multiplier = 1
        tda_results = tda_drv.compute(molecule, basis, scf_results)

        exc_drv = ExcitedStateAnalysisDriver()
        exc_drv.analysis_type = "mulliken"
        exc_drv.ostream.mute()

        # add fragment dictionary to ExcitedStateAnalysisDriver
        exc_drv.fragment_dict = fragment_dict

        descriptor_dict_s1 = exc_drv.compute(molecule, basis, scf_results,
                                             tda_results, 1)

        assert np.max(np.abs(descriptor_dict_s1["ct_matrix"] -
                             ct_matrix)) < 1.0e-06
        assert np.max(
            np.abs(descriptor_dict_s1["hole_participation_ratio"] -
                   hole_participation_ratio)) < 1.0e-06
        assert np.max(
            np.abs(descriptor_dict_s1["particle_participation_ratio"] -
                   particle_participation_ratio)) < 1.0e-06
        assert np.max(
            np.abs(descriptor_dict_s1["avg_participation_ratio"] -
                   avg_participation_ratio)) < 1.0e-06
        assert np.max(
            np.abs(descriptor_dict_s1["avg_hole_position"] -
                   avg_hole_position)) < 1.0e-06
        assert np.max(
            np.abs(descriptor_dict_s1["avg_particle_position"] -
                   avg_particle_position)) < 1.0e-06
        assert np.max(
            np.abs(descriptor_dict_s1["avg_difference_vector"] -
                   avg_difference_vec)) < 1.0e-06
        assert np.max(np.abs(descriptor_dict_s1["ct_length"] -
                             ct_length)) < 1.0e-06

    @pytest.mark.parametrize('ri', [False, True])
    def test_rpa(self, ri):
        ct_matrix = {
            False: np.array([[0.37055813, 0.19016109], [0.13581409,
                                                        0.32858347]]),
            True: np.array([[0.37054215, 0.19015108], [0.13582063, 0.32860462]])
        }[ri]

        hole_participation_ratio = {
            False: 1.982496948226768,
            True: 1.982516338995640
        }[ri]
        particle_participation_ratio = {
            False: 1.9997087109644545,
            True: 1.9997077427444081
        }[ri]
        avg_participation_ratio = {
            False: 1.9911028295956112,
            True: 1.9911120408700242
        }[ri]

        avg_hole_position = {
            False: np.array([0.10244102, 0.04395257, -0.00010709]),
            True: np.array([0.10234146, 0.04395201, -0.00010709])
        }[ri]
        avg_particle_position = {
            False: np.array([1.55512189e-02, 1.76200117e-02, -2.91052404e-05]),
            True: np.array([0.01552509, 0.01761981, -0.00002911])
        }[ri]
        avg_difference_vec = {
            False: np.array([-8.68898023e-02, -2.63325597e-02, 7.79813570e-05]),
            True: np.array([-0.08681637, -0.02633220, 0.00007798])
        }[ri]

        ct_length = {False: 0.09079233183655822, True: 0.09072195686301675}[ri]

        thienyl_thiazole_str = """
        S   1.5860   -1.4206    0.0001
        C   0.6775    0.0310   -0.0004
        C   1.4857    1.1466   -0.0004
        C   2.8690    0.8051   -0.0001
        C   3.0653   -0.5583    0.0002
        H   1.1367    2.1716   -0.0005
        H   3.6790    1.5228    0.0000
        H   4.0068   -1.0893    0.0006
        S  -1.6588    1.4476    0.0009
        N  -1.4203   -1.1264    0.0001
        C  -0.7409   -0.0017    0.0002
        C  -2.7697   -0.8319    0.0000
        C  -3.0939    0.5086   -0.0006
        H  -4.0731    0.9658   -0.0010
        H  -3.4860   -1.6422   -0.0003
        """

        fragment_dict = {
            "Donor": [1, 2, 3, 4, 5, 6, 7, 8],
            "Acceptor": [9, 10, 11, 12, 13, 14, 15],
        }

        molecule = Molecule.read_molecule_string(thienyl_thiazole_str)
        basis = MolecularBasis.read(molecule, "STO-3G", ostream=None)

        # SCF settings and calculation
        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        # Solve the TDDFT linear response equation for the first 2 excited states
        lreig_drv = LinearResponseEigenSolver()
        lreig_drv.ostream.mute()
        lreig_drv.nstates = 2
        lreig_drv.initial_guess_multiplier = 1
        lreig_drv.max_subspace_dim = 1000
        lreig_results = lreig_drv.compute(molecule, basis, scf_results)

        exc_drv = ExcitedStateAnalysisDriver()
        exc_drv.analysis_type = "mulliken"
        exc_drv.ostream.mute()

        # add fragment dictionary to ExcitedStateAnalysisDriver
        exc_drv.fragment_dict = fragment_dict

        descriptor_dict_s1 = exc_drv.compute(molecule, basis, scf_results,
                                             lreig_results, 1)

        assert np.max(np.abs(descriptor_dict_s1["ct_matrix"] -
                             ct_matrix)) < 1.0e-06
        assert np.max(
            np.abs(descriptor_dict_s1["hole_participation_ratio"] -
                   hole_participation_ratio)) < 1.0e-06
        assert np.max(
            np.abs(descriptor_dict_s1["particle_participation_ratio"] -
                   particle_participation_ratio)) < 1.0e-06
        assert np.max(
            np.abs(descriptor_dict_s1["avg_participation_ratio"] -
                   avg_participation_ratio)) < 1.0e-06
        assert np.max(
            np.abs(descriptor_dict_s1["avg_hole_position"] -
                   avg_hole_position)) < 1.0e-06
        assert np.max(
            np.abs(descriptor_dict_s1["avg_particle_position"] -
                   avg_particle_position)) < 1.0e-06
        assert np.max(
            np.abs(descriptor_dict_s1["avg_difference_vector"] -
                   avg_difference_vec)) < 1.0e-06
        assert np.max(np.abs(descriptor_dict_s1["ct_length"] -
                             ct_length)) < 1.0e-06

    @pytest.mark.parametrize('ri', [False, True])
    def test_core_excited_rpa(self, ri):
        molecule_xyz = """8

        C   -0.158859000000   1.034900000000    1.205196000000
        C   -1.156627000000  -0.031075000000    1.507521000000
        O   -2.065612000000   0.179139000000    2.354094000000
        O   -1.110771000000  -1.239690000000    0.820160000000
        H    0.760665000000   0.583534000000    0.777146000000
        H    0.108925000000   1.580853000000    2.134378000000
        H   -0.589966000000   1.747965000000    0.472205000000
        H   -1.791528000000  -1.967324000000    1.002895000000
        """
        # reference hole and particle positions
        avg_hole_position = np.array([-2.06430183, 0.17889688, 2.35290839])
        avg_particle_position = {
            False: np.array([-1.45474337, 0.02338098, 1.75900387]),
            True: np.array([-1.45473892, 0.02338250, 1.75900196])
        }[ri]
        avg_difference_vec = {
            False: np.array([0.60955845, -0.1555159, -0.59390452]),
            True: np.array([0.60956306, -0.15551441, -0.59390657])
        }[ri]

        # Molecule and basis
        molecule = Molecule.from_xyz_string(molecule_xyz)
        fragment_dict = {
            "CH3": [1, 5, 6, 7],
            "COOH": [2, 3, 4, 8],
        }
        basis = MolecularBasis.read(molecule, "def2-svp", ostream=None)

        # SCF settings and calculation
        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        # Solve the TDDFT linear response equation for the first 2 excited states
        lreig_drv = LinearResponseEigenSolver()
        lreig_drv.ostream.mute()
        lreig_drv.core_excitation = True
        # Core-excited states, O 1s
        lreig_drv.num_core_orbitals = 2
        lreig_drv.nstates = 2
        lreig_results = lreig_drv.compute(molecule, basis, scf_results)

        exc_drv = ExcitedStateAnalysisDriver()
        exc_drv.analysis_type = "mulliken"
        exc_drv.ostream.mute()
        exc_drv.fragment_dict = fragment_dict

        descriptor_dict_s1 = exc_drv.compute(
            molecule=molecule,
            basis=basis,
            scf_results=scf_results,
            rsp_results=lreig_results,
            state_index=1,
            num_core_orbitals=lreig_drv.num_core_orbitals)
        assert np.max(
            np.abs(descriptor_dict_s1["avg_hole_position"] -
                   avg_hole_position)) < 1.0e-06
        assert np.max(
            np.abs(descriptor_dict_s1["avg_particle_position"] -
                   avg_particle_position)) < 1.0e-06
        assert np.max(
            np.abs(descriptor_dict_s1["avg_difference_vector"] -
                   avg_difference_vec)) < 1.0e-06

    @pytest.mark.parametrize('ri', [False, True])
    def test_core_excited_tda(self, ri):
        molecule_xyz = """8

        C   -0.158859000000   1.034900000000    1.205196000000
        C   -1.156627000000  -0.031075000000    1.507521000000
        O   -2.065612000000   0.179139000000    2.354094000000
        O   -1.110771000000  -1.239690000000    0.820160000000
        H    0.760665000000   0.583534000000    0.777146000000
        H    0.108925000000   1.580853000000    2.134378000000
        H   -0.589966000000   1.747965000000    0.472205000000
        H   -1.791528000000  -1.967324000000    1.002895000000
        """
        # reference hole and particle positions
        avg_hole_position = np.array([-2.06430185, 0.17889688, 2.35290842])
        avg_particle_position = {
            False: np.array([-1.45471821, 0.02337554, 1.7589794]),
            True: np.array([-1.45471647, 0.02337748, 1.75898005])
        }[ri]
        avg_difference_vec = {
            False: np.array([0.60958364, -0.15552134, -0.59392902]),
            True: np.array([0.60958554, -0.15551943, -0.59392851])
        }[ri]

        # Molecule and basis
        molecule = Molecule.from_xyz_string(molecule_xyz)
        fragment_dict = {
            "CH3": [1, 5, 6, 7],
            "COOH": [2, 3, 4, 8],
        }
        basis = MolecularBasis.read(molecule, "def2-svp", ostream=None)

        # SCF settings and calculation
        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        # Solve the TDDFT linear response equation for the first 2 excited states
        tda_drv = TdaEigenSolver()
        tda_drv.ostream.mute()
        tda_drv.core_excitation = True
        # Core-excited states, O 1s
        tda_drv.num_core_orbitals = 2
        tda_drv.nstates = 2
        tda_results = tda_drv.compute(molecule, basis, scf_results)

        exc_drv = ExcitedStateAnalysisDriver()
        exc_drv.analysis_type = "mulliken"
        exc_drv.ostream.mute()
        exc_drv.fragment_dict = fragment_dict

        descriptor_dict_s1 = exc_drv.compute(
            molecule=molecule,
            basis=basis,
            scf_results=scf_results,
            rsp_results=tda_results,
            state_index=1,
            num_core_orbitals=tda_drv.num_core_orbitals)
        assert np.max(
            np.abs(descriptor_dict_s1["avg_hole_position"] -
                   avg_hole_position)) < 1.0e-06
        assert np.max(
            np.abs(descriptor_dict_s1["avg_particle_position"] -
                   avg_particle_position)) < 1.0e-06
        assert np.max(
            np.abs(descriptor_dict_s1["avg_difference_vector"] -
                   avg_difference_vec)) < 1.0e-06

    def test_analysis_from_file(self):
        here = Path(__file__).parent
        filename = str(here / 'data' / 'acetic_acid.h5')
        fragment_dict = {
            "CH3": [1, 5, 6, 7],
            "COOH": [2, 3, 4, 8],
        }

        # reference values
        ct_matrix = np.array([
            [0.00806936, 0.05671282],
            [0.12028262, 0.8164166],
        ])
        hole_participation_ratio = 1.1376616910096484
        particle_participation_ratio = 1.2877855259522788
        avg_participation_ratio = 1.2127236084809636
        avg_hole_position = np.array([-1.77749224, 0.04863906, 2.05871356])
        avg_particle_position = np.array([-1.28883208, 0.11053515, 1.59960834])
        avg_difference_vec = np.array([0.48866016, 0.06189609, -0.45910522])
        ct_length = 0.6733479630377631

        # perform analysis using scf and linear response from h5
        exc_drv = ExcitedStateAnalysisDriver()
        exc_drv.analysis_type = "mulliken"
        exc_drv.ostream.mute()

        exc_drv.fragment_dict = fragment_dict

        molecule, basis = read_molecule_and_basis(filename)
        scf_res, rsp_res = exc_drv.read_from_h5(filename)

        descriptor_dict_s1 = exc_drv.compute(molecule,
                                             basis,
                                             scf_res,
                                             rsp_res,
                                             state_index=1)

        assert np.max(np.abs(descriptor_dict_s1["ct_matrix"] -
                             ct_matrix)) < 1.0e-06
        assert np.max(
            np.abs(descriptor_dict_s1["hole_participation_ratio"] -
                   hole_participation_ratio)) < 1.0e-06
        assert np.max(
            np.abs(descriptor_dict_s1["particle_participation_ratio"] -
                   particle_participation_ratio)) < 1.0e-06
        assert np.max(
            np.abs(descriptor_dict_s1["avg_participation_ratio"] -
                   avg_participation_ratio)) < 1.0e-06
        assert np.max(
            np.abs(descriptor_dict_s1["avg_hole_position"] -
                   avg_hole_position)) < 1.0e-06
        assert np.max(
            np.abs(descriptor_dict_s1["avg_particle_position"] -
                   avg_particle_position)) < 1.0e-06
        assert np.max(
            np.abs(descriptor_dict_s1["avg_difference_vector"] -
                   avg_difference_vec)) < 1.0e-06
        assert np.max(np.abs(descriptor_dict_s1["ct_length"] -
                             ct_length)) < 1.0e-06

    def test_lowdin_analysis_from_file(self):
        here = Path(__file__).parent
        filename = str(here / 'data' / 'acetic_acid.h5')
        fragment_dict = {
            "CH3": [1, 5, 6, 7],
            "COOH": [2, 3, 4, 8],
        }

        # reference values
        ct_matrix = np.array([
            [0.0077972, 0.05471233],
            [0.11902532, 0.81994655]
        ])
        hole_participation_ratio = 1.1325571
        particle_participation_ratio = 1.2840218
        avg_participation_ratio = 1.2082895
        avg_hole_position = np.array([-1.77883271, 0.03979253, 2.05341404])
        avg_particle_position = np.array([-1.29253462, 0.08706147, 1.60096523])
        avg_difference_vec = np.array([0.48629809, 0.04726894, -0.45244881])
        ct_length = 0.66590549

        # perform analysis using scf and linear response from h5
        exc_drv = ExcitedStateAnalysisDriver()
        exc_drv.analysis_type = "lowdin"
        exc_drv.ostream.mute()

        exc_drv.fragment_dict = fragment_dict

        molecule, basis = read_molecule_and_basis(filename)
        scf_res, rsp_res = exc_drv.read_from_h5(filename)

        descriptor_dict_s1 = exc_drv.compute(molecule,
                                             basis,
                                             scf_res,
                                             rsp_res,
                                             state_index=1)

        assert np.max(np.abs(descriptor_dict_s1["ct_matrix"] -
                             ct_matrix)) < 1.0e-06
        assert np.max(
            np.abs(descriptor_dict_s1["hole_participation_ratio"] -
                   hole_participation_ratio)) < 1.0e-06
        assert np.max(
            np.abs(descriptor_dict_s1["particle_participation_ratio"] -
                   particle_participation_ratio)) < 1.0e-06
        assert np.max(
            np.abs(descriptor_dict_s1["avg_participation_ratio"] -
                   avg_participation_ratio)) < 1.0e-06
        assert np.max(
            np.abs(descriptor_dict_s1["avg_hole_position"] -
                   avg_hole_position)) < 1.0e-06
        assert np.max(
            np.abs(descriptor_dict_s1["avg_particle_position"] -
                   avg_particle_position)) < 1.0e-06
        assert np.max(
            np.abs(descriptor_dict_s1["avg_difference_vector"] -
                   avg_difference_vec)) < 1.0e-06
        assert np.max(np.abs(descriptor_dict_s1["ct_length"] -
                             ct_length)) < 1.0e-06
