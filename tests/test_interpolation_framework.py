from pathlib import Path
import numpy as np
import pytest

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.interpolationdriver import InterpolationDriver
from veloxchem.interpolationdatapoint import InterpolationDatapoint
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfgradientdriver import ScfGradientDriver
from veloxchem.scfhessiandriver import ScfHessianDriver


class TestInterpolationSetup:
    
    def test_b_b2_matrix(self):

        def compute_numerical_b2_matrix(structure, z_matrix, delta=1e-4):

            b0 = compute_numerical_B_matrix(structure, z_matrix, delta)
            n_internal, n_cart = b0.shape
            n_atoms = structure.shape[0]

            b2 = np.zeros((n_internal, n_cart, n_cart))

            for atom in range(n_atoms):
                for coord in range(3):
                    k = atom * 3 + coord

                    structure_forward = structure.copy()
                    structure_backward = structure.copy()
                    structure_forward[atom, coord] += delta
                    structure_backward[atom, coord] -= delta
                    
                    b_forward = compute_numerical_B_matrix(structure_forward, z_matrix, delta)
                    b_backward = compute_numerical_B_matrix(structure_backward, z_matrix, delta)
                    
                    db_dx = (b_forward - b_backward) / (2 * delta)

                    b2[:, :, k] = db_dx
                    
            return b2

        def compute_numerical_B_matrix(structure, z_matrix, delta=1e-4):

            interpolation_datapoint = InterpolationDatapoint(z_matrix)
            interpolation_datapoint.reset_coordinates(structure)

            q0 = interpolation_datapoint.internal_coordinates_values
            n_q = len(q0)
            n_atoms = structure.shape[0]
            
            b = np.zeros((n_q, n_atoms * 3))
            
            for atom in range(n_atoms):
                for coord in range(3):

                    structure_forward = structure.copy()
                    structure_backward = structure.copy()
                    
                    structure_forward[atom, coord] += delta
                    structure_backward[atom, coord] -= delta
                    
                    interpolation_datapoint.reset_coordinates(structure_forward)
                    q_forward = interpolation_datapoint.internal_coordinates_values
                    interpolation_datapoint.reset_coordinates(structure_backward)
                    q_backward = interpolation_datapoint.internal_coordinates_values
                    
                    derivative = (q_forward - q_backward) / (2 * delta)
                    for i, elem in enumerate(z_matrix):
                        if len(elem) == 4:                
                            derivative[i] = (derivative[i] + np.pi) % (2 * np.pi) - np.pi

                    col_index = atom * 3 + coord
                    b[:, col_index] = derivative
                    
            return b
        
        molecule_xyz = '''7

        C             -0.666468000000        -2.187340000000        -0.571285000000
        C              0.777456000000        -2.504713000000        -0.361732000000
        O              1.429604000000        -1.880916000000         0.459616000000
        H             -0.774975000000        -1.132031000000        -0.898359000000
        H             -1.225231000000        -2.339168000000         0.375900000000
        H             -1.088468000000        -2.854107000000        -1.352043000000
        H              1.251710000000        -3.291432000000        -0.939414000000'''

        molecule = Molecule.from_xyz_string(molecule_xyz)
        z_matrix = [(0, 1), (0, 3), (0, 4), (0, 5), (1, 2), (1, 6), (1, 0, 3), (1, 0, 4), (1, 0, 5), (3, 0, 4), (3, 0, 5), (4, 0, 5), (0, 1, 2), (0, 1, 6), (2, 1, 6), (3, 0, 1, 2), (3, 0, 1, 6), (4, 0, 1, 2), (4, 0, 1, 6), (5, 0, 1, 2), (5, 0, 1, 6)]
        
        num_b = compute_numerical_B_matrix(molecule.get_coordinates_in_bohr(), z_matrix)
        num_b2 = compute_numerical_b2_matrix(molecule.get_coordinates_in_bohr(), z_matrix)

        interpolation_datapoint = InterpolationDatapoint(z_matrix)
        interpolation_datapoint.reset_coordinates(molecule.get_coordinates_in_bohr())
        
        b_matrix = interpolation_datapoint.b_matrix
        b2_matrix = interpolation_datapoint.b2_matrix

        # Check if the B- and B2-Matrix are correctly determined by comparing to the numerical versions

        assert np.max(np.abs(num_b - b_matrix)) < 1.0e-5
        assert np.max(np.abs(num_b2 - b2_matrix)) < 1.0e-5

    def test_internal_transformation(self):
        
        molecule_xyz = '''7

        C             -0.666468000000        -2.187340000000        -0.571285000000
        C              0.777456000000        -2.504713000000        -0.361732000000
        O              1.429604000000        -1.880916000000         0.459616000000
        H             -0.774975000000        -1.132031000000        -0.898359000000
        H             -1.225231000000        -2.339168000000         0.375900000000
        H             -1.088468000000        -2.854107000000        -1.352043000000
        H              1.251710000000        -3.291432000000        -0.939414000000'''
        molecule = Molecule.from_xyz_string(molecule_xyz)
        z_matrix = [(0, 1), (0, 3), (0, 4), (0, 5), (1, 2), (1, 6), (1, 0, 3), (1, 0, 4), (1, 0, 5), (3, 0, 4), (3, 0, 5), (4, 0, 5), (0, 1, 2), (0, 1, 6), (2, 1, 6), (3, 0, 1, 2), (3, 0, 1, 6), (4, 0, 1, 2), (4, 0, 1, 6), (5, 0, 1, 2), (5, 0, 1, 6)]
        interpolation_datapoint = InterpolationDatapoint(z_matrix)
        interpolation_datapoint.cartesian_coordinates = molecule.get_coordinates_in_bohr()
        qm_driver = ScfRestrictedDriver()
        basis = MolecularBasis.read(molecule, 'def2-svp')
        qm_driver.ostream.mute()
        scf_results = qm_driver.compute(molecule, basis)
        qm_energy = qm_driver.scf_energy
        
        qm_grad_driver = ScfGradientDriver(qm_driver)
        qm_grad_driver.compute(molecule, basis, scf_results)
        qm_gradient = qm_grad_driver.gradient

        qm_hess_driver = ScfHessianDriver(qm_driver)
        qm_hess_driver.compute(molecule, basis)
        qm_hessian = qm_hess_driver.hessian

        interpolation_datapoint.energy = qm_energy
        interpolation_datapoint.gradient = qm_gradient
        interpolation_datapoint.hessian = qm_hessian

        interpolation_datapoint.transform_gradient_and_hessian()
        
        cartesian_gradient = interpolation_datapoint.backtransform_internal_gradient_to_cartesian_coordinates()
        cartesian_hessian = interpolation_datapoint.backtransform_internal_hessian_to_cartesian_coordinates()
        
        # Check that the if the QM Gradient is being transformed and backtransformed correclty between
        # internal and Cartesian coordinates

        assert np.max(np.abs(qm_gradient - cartesian_gradient)) < 1.0e-5
        assert np.max(np.abs(cartesian_hessian - cartesian_hessian)) < 1.0e-5

    def test_interpolation_scheme(self):

        here = Path(__file__).parent
        interpolationdatafile = str(here / 'data' / 'test_interpolation_database.h5')

        molecule_xyz_1 = '''7

        C             -0.666468000000        -2.187340000000        -0.571285000000
        C              0.777456000000        -2.504713000000        -0.361732000000
        O              1.429604000000        -1.880916000000         0.459616000000
        H             -0.774975000000        -1.132031000000        -0.898359000000
        H             -1.225231000000        -2.339168000000         0.375900000000
        H             -1.088468000000        -2.854107000000        -1.352043000000
        H              1.251710000000        -3.291432000000        -0.93941400000
        '''

        molecule_xyz_2 = '''7

        C   -0.6883    -2.1533     -0.6223
        C   0.7228    -2.5875     -0.2626
        O   1.4118    -2.1096     0.6195
        H   -0.8346    -1.0939     -0.3143
        H   -1.4433    -2.6479     0.0119
        H   -0.7988    -2.3508     -1.6990
        H   1.0347    -3.5709     -0.8075
        '''

        molecule_xyz_strings = [molecule_xyz_1, molecule_xyz_2]

        energies_mol_1 = [-153.7121168801204, -153.7121168801204]
        energies_mol_2 = [-153.70999371824846, -153.7099979081999]

        gradient_mol_1 = [np.array(
                [[ 0.01259014,  0.00962321,  0.01356831],
                [-0.00900223, -0.03474453, -0.03609507],
                [ 0.00290219,  0.01188851,  0.01228797],
                [ 0.0002292 ,  0.00487737, -0.00379572],
                [-0.00273938, -0.00308223,  0.00460716],
                [ 0.00214565, -0.00506522, -0.00404103],
                [-0.00612558,  0.01650289,  0.01346838]]), 
                np.array(
                [[ 0.01259014,  0.00962321,  0.01356831],
                [-0.00900223, -0.03474453, -0.03609507],
                [ 0.00290219,  0.01188851,  0.01228797],
                [ 0.0002292 ,  0.00487737, -0.00379572],
                [-0.00273938, -0.00308223,  0.00460716],
                [ 0.00214565, -0.00506522, -0.00404103],
                [-0.00612558,  0.01650289,  0.01346838]])] 

        gradient_mol_2 = [np.array(
                [[-0.01506613, -0.01803113, -0.00699575],
                 [-0.00249755,  0.02186579, -0.02299591],
                 [ 0.01118285,  0.00521398,  0.02315084],
                 [-0.00404425,  0.00367571,  0.0109869 ],
                 [-0.00264075,  0.01235758,  0.00635172],
                 [ 0.01268892, -0.00372708, -0.00254009],
                 [ 0.00037691, -0.02135485, -0.00795771]]), 
                np.array(
                [[-0.01505409, -0.01802155, -0.00699913],
                 [-0.0024999 ,  0.02185242, -0.02297884],
                 [ 0.01118164,  0.00520655,  0.02315804],
                 [-0.00403725,  0.00368312,  0.01096718],
                 [-0.00264169,  0.01234974,  0.00635032],
                 [ 0.01266616, -0.00373323, -0.00252843],
                 [ 0.00038513, -0.02133706, -0.00796914]])]

        energies = [energies_mol_1, energies_mol_2]

        gradients = [gradient_mol_1, gradient_mol_2]

        for i, mol_xyz in enumerate(molecule_xyz_strings):
            molecule = Molecule.from_xyz_string(mol_xyz)
            
            correct_shep_interpolation_energy = energies[i][0]
            correct_shep_interpolation_gradient = gradients[i][0]

            interpolation_settings = { 'interpolation_type':'shepard', 
                                'exponent_p':'2',
                                'exponent_q':'2', 
                                'confidence_radius':'0.5',
                                'use_inverse_bond_length':True
                            }

            interpolation_driver = InterpolationDriver()
            interpolation_driver.update_settings(interpolation_settings)
            interpolation_driver.imforcefield_file = interpolationdatafile
            labels, z_matrix = interpolation_driver.read_labels()
            interpolation_driver.impes_coordinate.z_matrix = z_matrix
            sorted_labels = sorted(labels, key=lambda x: int(x.split('_')[1]))
            org_z_matrix = [(0, 1), (0, 3), (0, 4), (0, 5), (1, 2), (1, 6), (1, 0, 3), (1, 0, 4), (1, 0, 5), (3, 0, 4), (3, 0, 5), (4, 0, 5), (0, 1, 2), (0, 1, 6), (2, 1, 6), (3, 0, 1, 2), (3, 0, 1, 6), (4, 0, 1, 2), (4, 0, 1, 6), (5, 0, 1, 2), (5, 0, 1, 6)]
            
            # Check if the internal coordinates are being defined correctly

            assert sorted(map(tuple, org_z_matrix)) == sorted(map(tuple, z_matrix))

            interpolation_driver.impes_coordinate.z_matrix = z_matrix
            
            im_datapoints = []
            for ordered_label in sorted_labels:
                im_datapoint = InterpolationDatapoint(z_matrix)
                im_datapoint.read_hdf5(interpolationdatafile, ordered_label)
                im_datapoints.append(im_datapoint)

            interpolation_driver.compute(molecule, im_datapoints)
            im_shepard_energy = interpolation_driver.impes_coordinate.energy
            im_shepard_gradient = interpolation_driver.impes_coordinate.gradient

            # Check if the interpolation energy computed matches the stored interpolation energy for shepard type

            assert np.max(np.abs(im_shepard_energy - correct_shep_interpolation_energy)) < 1.0e-5
            assert np.max(np.abs(im_shepard_gradient - correct_shep_interpolation_gradient)) < 1.0e-5

            interpolation_settings = { 'interpolation_type':'simple', 
                                'exponent_p':'2',
                                'exponent_q':'2', 
                                'confidence_radius':'0.5',
                                'use_inverse_bond_length':True
                            }

            interpolation_driver_simple = InterpolationDriver()
            interpolation_driver_simple.update_settings(interpolation_settings)
            interpolation_driver_simple.imforcefield_file = interpolationdatafile
            interpolation_driver_simple.impes_coordinate.z_matrix = z_matrix

            interpolation_driver_simple.compute(molecule, im_datapoints)
            im_simple_energy = interpolation_driver_simple.impes_coordinate.energy
            im_simple_gradient = interpolation_driver_simple.impes_coordinate.gradient

            correct_simple_interpolation_energy = energies[i][1]
            correct_simple_interpolation_gradient = gradients[i][1]

            # Check if the interpolation energy computed matches the stored interpolation energy for simple type

            assert np.max(np.abs(im_simple_energy - correct_simple_interpolation_energy)) < 1.0e-5
            assert np.max(np.abs(im_simple_gradient - correct_simple_interpolation_gradient)) < 1.0e-5
