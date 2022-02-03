import pytest

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.quadraticresponsedriver import QuadraticResponseDriver
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


@pytest.mark.solvers
class TestQrf:

    def run_scf(self):

        molecule_string = """
        O   0.0   0.0   0.0
        H   0.0   1.4   1.1
        H   0.0  -1.4   1.1
        """

        basis_set_label = 'cc-pVDZ'

        scf_settings = {'conv_thresh': 1.0e-6}

        molecule = Molecule.read_str(molecule_string, units='au')
        molecule.set_charge(0)
        molecule.set_multiplicity(1)

        method_settings = {'xcfun': 'b3lyp', 'grid_level': 4}

        basis = MolecularBasis.read(molecule, basis_set_label)

        scf_drv = ScfRestrictedDriver()
        scf_drv.update_settings(scf_settings, method_settings)
        scf_drv.compute(molecule, basis)

        return scf_drv.scf_tensors,molecule,basis

    def run_qrf(self, ref_result):

        scf_tensors,molecule,ao_basis = self.run_scf()

        method_settings = {'xcfun': 'b3lyp', 'grid_level': 4}

        rsp_settings = {
            'conv_thresh': 1.0e-6, 
            'b_frequencies': [-0.1],
            'c_frequencies': [0.3],
            'damping': 0.1,
            'a_components':'x',
            'b_components':'x',
            'c_components':'z'
            }

        qrf_prop = QuadraticResponseDriver()

        qrf_prop.update_settings(rsp_settings, method_settings)

        qrf_result_xxz = qrf_prop.compute(molecule, ao_basis, scf_tensors,method_settings)

        rsp_settings = {'conv_thresh': 1.0e-6, 
        'b_frequencies': [-0.1],
        'c_frequencies': [0.3],
        'damping': 0.1,
        'a_components':'y',
        'b_components':'z',
        'c_components':'y'
        }

        qrf_prop.update_settings(rsp_settings, method_settings)

        qrf_result_yzy = qrf_prop.compute(molecule, ao_basis, scf_tensors,method_settings)

        rsp_settings = {'conv_thresh': 1.0e-6, 
        'b_frequencies': [-0.1],
        'c_frequencies': [0.3],
        'damping': 0.1,
        'a_components':'z',
        'b_components':'z',
        'c_components':'z'
        }

        qrf_prop.update_settings(rsp_settings, method_settings)
        
        qrf_result_zzz = qrf_prop.compute(molecule, ao_basis, scf_tensors,method_settings)
        
        thresh = 1.0e-4

        # x-component 

        if is_mpi_master():
             
            assert abs(qrf_result_xxz[(-0.1,0.3)].real - ref_result['xxz'].real) < thresh

            assert abs(qrf_result_xxz[(-0.1,0.3)].imag - ref_result['xxz'].imag) < thresh
            
            # y-component 
            
            assert abs(qrf_result_yzy[(-0.1,0.3)].real - ref_result['yzy'].real) < thresh

            assert abs(qrf_result_yzy[(-0.1,0.3)].imag - ref_result['yzy'].imag) < thresh

            # z-component 
            
            assert abs(qrf_result_zzz[(-0.1,0.3)].real - ref_result['zzz'].real) < thresh

            assert abs(qrf_result_zzz[(-0.1,0.3)].imag - ref_result['zzz'].imag) < thresh

    def test_qrf(self):

        ref_result = { 
            'xxz': 1.893733806672778+2.4615609240441643j,
            'yzy': 24.52654292764868+11.22800902489504j,
            'zzz': 13.51752879736258+18.681411175542237j,
        }

        self.run_qrf(ref_result)

