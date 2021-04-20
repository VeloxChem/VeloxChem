from pathlib import Path
import unittest
import tempfile

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.cubicgrid import CubicGrid


class TestCheckpoint(unittest.TestCase):

    def write_cube(self, fname):

        text = """
            test mo=9
            MO coefficients
              -6   -0.100000   -0.100000   -0.100000    1
               2    0.300000    0.000000    0.000000
               2    0.000000    0.300000    0.000000
               3    0.000000    0.000000    0.200000
               6    6.000000    0.000000    0.000000   -1.243600
               6    6.000000    0.000000    0.000000    1.243600
               1    1.000000    0.000000    1.727200   -2.323000
               1    1.000000    0.000000   -1.727200   -2.323000
               1    1.000000    0.000000    1.727200    2.323000
               1    1.000000    0.000000   -1.727200    2.323000
               1    9
             4.77717E-03 -4.77717E-03 -1.52606E-02
             4.62157E-03 -4.62157E-03 -1.47339E-02
            -9.24233E-03  9.24233E-03  2.94651E-02
            -8.94621E-03  8.94621E-03  2.84649E-02
        """

        with open(fname, 'w') as f_cube:
            for line in text.splitlines():
                if line.strip():
                    print(line.lstrip(), file=f_cube)

    def test_cube(self):

        with tempfile.TemporaryDirectory() as temp_dir:
            if not is_mpi_master():
                return

            fname = str(Path(temp_dir, 'test.cube'))
            self.write_cube(fname)

            grid_1 = CubicGrid.read_cube(fname)
            grid_2 = CubicGrid.read_cube(fname)
            diff = grid_1.compare(grid_2)
            self.assertTrue(diff < 1.0e-12)


if __name__ == "__main__":
    unittest.main()
