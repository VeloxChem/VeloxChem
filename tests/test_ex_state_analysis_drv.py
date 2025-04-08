import numpy as np

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tdaeigensolver import TdaEigenSolver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.excitedstateanalysisdriver import ExcitedStateAnalysisDriver


class TestCTNumbers:
    def test_tda(self):
        ct_matrix = np.array([[0.34837228, 0.19662888], [0.1367462,  0.31825265]])

        hole_participation_ratio = 1.9839293483578253
        particle_participation_ratio = 1.9982298906088756
        avg_participation_ratio = 1.9910796194833504
        
        avg_hole_position = np.array([ 1.39951260e-01,  1.33518801e-01, -7.95624567e-05])
        avg_particle_position = np.array([-1.46077909e-01,  1.08171427e-01, -8.18080571e-05])
        avg_difference_vec = np.array([-2.86029169e-01, -2.53473739e-02, -2.24560031e-06])
        
        ct_length = 0.2871500913516519

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
        scf_drv.ostream.mute()
        scf_drv.grid_level = 3
        scf_results = scf_drv.compute(molecule, basis)

        # Solve the TDDFT linear response equation for the first 2 excited states
        tda_drv = TdaEigenSolver()
        tda_drv.ostream.mute()
        tda_drv.nstates = 2
        tda_results = tda_drv.compute(molecule, basis, scf_results)

        exc_drv = ExcitedStateAnalysisDriver()
        exc_drv.ostream.mute()

        # add fragment dictionary to ExcitedStateAnalysisDriver
        exc_drv.fragment_dict = fragment_dict

        descriptor_dict_s1 = exc_drv.compute(molecule, basis, scf_results, tda_results, 1)

        assert np.max(np.abs(descriptor_dict_s1['ct_matrix'] - ct_matrix)) <1.0e-06
        assert np.max(np.abs(descriptor_dict_s1['hole_participation_ratio'] - hole_participation_ratio)) <1.0e-06
        assert np.max(np.abs(descriptor_dict_s1['particle_participation_ratio'] - particle_participation_ratio)) <1.0e-06
        assert np.max(np.abs(descriptor_dict_s1['avg_participation_ratio'] - avg_participation_ratio)) <1.0e-06
        assert np.max(np.abs(descriptor_dict_s1['avg_hole_position'] - avg_hole_position)) <1.0e-06
        assert np.max(np.abs(descriptor_dict_s1['avg_particle_position'] - avg_particle_position)) <1.0e-06
        assert np.max(np.abs(descriptor_dict_s1['avg_difference_vector'] - avg_difference_vec)) <1.0e-06
        assert np.max(np.abs(descriptor_dict_s1['ct_length'] - ct_length)) <1.0e-06

        
    def test_rpa(self):
        ct_matrix = np.array([[0.37844034, 0.17007633], [0.11831815, 0.33537697]])
        
        hole_participation_ratio = 1.9822558674307456
        particle_participation_ratio = 1.999849478491969
        avg_participation_ratio = 1.9910526729613571
        
        avg_hole_position = np.array([ 1.50961023e-01,  1.34494466e-01, -7.94760197e-05])
        avg_particle_position = np.array([-9.57159520e-02,  1.12634410e-01, -8.14126679e-05])
        avg_difference_vec = np.array([-2.46676975e-01, -2.18600555e-02, -1.93664826e-06])
        
        ct_length =  0.2476436794645711

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
        scf_drv.ostream.mute()
        scf_drv.grid_level = 3
        scf_results = scf_drv.compute(molecule, basis)

        # Solve the TDDFT linear response equation for the first 2 excited states
        lreig_drv = LinearResponseEigenSolver()
        lreig_drv.ostream.mute()
        lreig_drv.nstates = 2
        lreig_results = lreig_drv.compute(molecule, basis, scf_results)

        exc_drv = ExcitedStateAnalysisDriver()
        exc_drv.ostream.mute()

        # add fragment dictionary to ExcitedStateAnalysisDriver
        exc_drv.fragment_dict = fragment_dict

        descriptor_dict_s1 = exc_drv.compute(molecule, basis, scf_results, lreig_results, 1)

        assert np.max(np.abs(descriptor_dict_s1['ct_matrix'] - ct_matrix)) <1.0e-06
        assert np.max(np.abs(descriptor_dict_s1['hole_participation_ratio'] - hole_participation_ratio)) <1.0e-06
        assert np.max(np.abs(descriptor_dict_s1['particle_participation_ratio'] - particle_participation_ratio)) <1.0e-06
        assert np.max(np.abs(descriptor_dict_s1['avg_participation_ratio'] - avg_participation_ratio)) <1.0e-06
        assert np.max(np.abs(descriptor_dict_s1['avg_hole_position'] - avg_hole_position)) <1.0e-06
        assert np.max(np.abs(descriptor_dict_s1['avg_particle_position'] - avg_particle_position)) <1.0e-06
        assert np.max(np.abs(descriptor_dict_s1['avg_difference_vector'] - avg_difference_vec)) <1.0e-06
        assert np.max(np.abs(descriptor_dict_s1['ct_length'] - ct_length)) <1.0e-06
