import numpy as np

from veloxchem.mmdriver import MMDriver
from veloxchem.mmgradientdriver import MMGradientDriver
from veloxchem.molecule import Molecule
from veloxchem.forcefieldgenerator import ForceFieldGenerator


class TestMMDriver:

    def run_mmdriver(self, molecule, partial_charges, expected_gradient,
                     expected_energy):

        ff_gen = ForceFieldGenerator()
        ff_gen.ostream.mute()
        ff_gen.eq_param = False
        ff_gen.partial_charges = partial_charges
        ff_gen.create_topology(molecule)

        mmdriver = MMDriver()
        mmdriver.load_force_field(ff_gen)
        mmdriver.compute(molecule)

        # Test against the numerical gradient
        grad_drv = MMGradientDriver(mmdriver)
        grad_drv.delta_h = 0.0001
        grad_drv.compute_numerical(molecule)
        assert np.max(np.abs(mmdriver.gradient - grad_drv.gradient)) < 1.0e-4

        # Test against the OpenMM gradient
        assert np.max(np.abs(mmdriver.gradient - expected_gradient)) < 1.0e-6

        # Test total energies against OpenMM
        assert np.abs(mmdriver.energy - expected_energy) < 1.0e-6

    def test_mmdriver_COO(self):

        xyz_string = """9

        C              1.067790000000         0.045770000000         0.065120000000
        C              2.582380000000         0.043650000000         0.060130000000
        O              3.063500000000         0.231070000000        -1.262860000000
        H              0.681060000000        -0.084160000000         1.079730000000
        H              0.677390000000        -0.758940000000        -0.566530000000
        H              0.682590000000         0.986030000000        -0.342860000000
        H              2.974890000000        -0.902100000000         0.445850000000
        H              2.968390000000         0.856610000000         0.681920000000
        H              2.724970000000        -0.502170000000        -1.804450000000
        """

        expected_gradient = [[0.007601, 0.000432, 0.004564],
                             [-0.015, -0.00187, 0.006224],
                             [0.00102, -0.000419, 0.001837],
                             [0.000768, 0.000154, -0.001983],
                             [0.000234, 0.002531, -0.000958],
                             [0.000263, -0.002411, -0.001634],
                             [0.002318, 0.001745, -0.002779],
                             [0.001997, -0.000441, -0.003949],
                             [0.000798, 0.00028, -0.001323]]

        expected_energy = -0.011381738932848186

        molecule = Molecule.read_xyz_string(xyz_string)

        partial_charges = [
            -0.13161353, 0.37116662, -0.7019936, 0.02540256, 0.02540256,
            0.02540256, -0.01666828, -0.01666828, 0.41956939
        ]

        self.run_mmdriver(molecule, partial_charges, expected_gradient,
                          expected_energy)

    def test_water_dimer(self):

        xyz_string = """6

        O             -0.471000000000         0.049100000000         0.409300000000
        H             -0.553100000000        -0.081200000000        -0.580200000000
        H              0.498200000000         0.049600000000         0.662700000000
        O             -1.995300000000         1.134700000000         2.419500000000
        H             -2.816300000000         0.573800000000         2.515900000000
        H             -1.430000000000         0.764800000000         1.669500000000
        """

        expected_gradient = [[-0.012024, 0.000819, 0.009577],
                             [-0.018689, -0.004335, -0.036627],
                             [0.032503, 0.002214, 0.024549],
                             [0.009845, 0.013545, 0.004913],
                             [-0.038857, -0.013445, 0.016343],
                             [0.027223, 0.001203, -0.018755]]

        expected_energy = -0.005108847272598356

        molecule = Molecule.read_xyz_string(xyz_string)

        partial_charges = [
            -0.91582993, 0.45791496, 0.45791496, -0.91582993, 0.45791496,
            0.45791496
        ]

        self.run_mmdriver(molecule, partial_charges, expected_gradient,
                          expected_energy)
