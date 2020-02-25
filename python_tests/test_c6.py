from mpi4py import MPI
import numpy as np
import unittest
import os

from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.c6solver import C6Solver
from veloxchem.rspc6 import C6


class TestC6(unittest.TestCase):

    def run_c6(self, inpfile, potfile, xcfun_label, data_lines, ref_c6_value):

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)
        task.input_dict['scf']['checkpoint_file'] = None

        if potfile is not None:
            task.input_dict['method_settings']['potfile'] = potfile
            try:
                import cppe
            except ImportError:
                return

        if xcfun_label is not None:
            task.input_dict['method_settings']['xcfun'] = xcfun_label

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        ref_freqs = []
        for line in data_lines:
            if float(line.split()[3]) not in ref_freqs:
                ref_freqs.append(float(line.split()[3]))
        ref_freqs = np.sort(ref_freqs)
        ref_n_points = int(len(ref_freqs) - 1)
        ref_prop_real = [float(line.split()[4]) for line in data_lines]

        c6_solver = C6Solver(task.mpi_comm, task.ostream)
        c6_solver.update_settings({'n_points': ref_n_points},
                                   task.input_dict['method_settings'])
        c6_results = c6_solver.compute(task.molecule, task.ao_basis,
                                         scf_drv.scf_tensors)

        if task.mpi_rank == mpi_master():
            keys = [iw for iw in c6_results['solutions'].keys()]
            freqs = []
            for key in keys:
                if key[1] not in freqs:
                    freqs.append(key[1])
            prop_freqs = np.sort(freqs)
            self.assertTrue(np.max(prop_freqs - ref_freqs) < 1.0e-6)

            prop = np.array([
                -c6_results['response_functions'][(a, b, iw)]
                for iw in prop_freqs
                for (a, b) in ['xx', 'yy', 'zz', 'xy', 'xz', 'yz']
            ])
            self.assertTrue(np.max(np.abs(prop.real - ref_prop_real)) < 1.0e-4)

            points, weights = np.polynomial.legendre.leggauss(ref_n_points)
            c6_value = C6.integrate(self, c6_results, prop_freqs[1:][::-1], points,
                                    weights, 0.3)
            print(c6_value, ref_c6_value, prop_freqs[1:])
            self.assertTrue(np.abs(c6_value - ref_c6_value) < 1.0e-4)

    def test_c6_hf(self):

        inpfile = os.path.join('inputs', 'water.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        potfile = None

        xcfun_label = None

        #   --------------------------------------------------
        #   No    A-oper    B-oper   Frequency       Real part
        #   --------------------------------------------------
        raw_data = """
             1   XDIPLEN   XDIPLEN    0.000000        7.251840
             2   YDIPLEN   YDIPLEN    0.000000        8.724528
             3   ZDIPLEN   ZDIPLEN    0.000000        7.880569
             4   XDIPLEN   YDIPLEN    0.000000        0.000000
             5   XDIPLEN   ZDIPLEN    0.000000        0.000000
             6   YDIPLEN   ZDIPLEN    0.000000       -0.000000
             7   XDIPLEN   XDIPLEN    0.006077        7.250980
             8   YDIPLEN   YDIPLEN    0.006077        8.723851
             9   ZDIPLEN   ZDIPLEN    0.006077        7.879847
            10   XDIPLEN   YDIPLEN    0.006077        0.000000
            11   XDIPLEN   ZDIPLEN    0.006077        0.000000
            12   YDIPLEN   ZDIPLEN    0.006077       -0.000000
            13   XDIPLEN   XDIPLEN    0.033952        7.225189
            14   YDIPLEN   YDIPLEN    0.033952        8.703479
            15   ZDIPLEN   ZDIPLEN    0.033952        7.858125
            16   XDIPLEN   YDIPLEN    0.033952        0.000000
            17   XDIPLEN   ZDIPLEN    0.033952        0.000000
            18   YDIPLEN   ZDIPLEN    0.033952       -0.000000
            19   XDIPLEN   XDIPLEN    0.093305        7.060265
            20   YDIPLEN   YDIPLEN    0.093305        8.568694
            21   ZDIPLEN   ZDIPLEN    0.093305        7.716148
            22   XDIPLEN   YDIPLEN    0.093305        0.000000
            23   XDIPLEN   ZDIPLEN    0.093305        0.000000
            24   YDIPLEN   ZDIPLEN    0.093305       -0.000000
            25   XDIPLEN   XDIPLEN    0.206999        6.468356
            26   YDIPLEN   YDIPLEN    0.206999        8.019561
            27   ZDIPLEN   ZDIPLEN    0.206999        7.164827
            28   XDIPLEN   YDIPLEN    0.206999        0.000000
            29   XDIPLEN   ZDIPLEN    0.206999        0.000000
            30   YDIPLEN   ZDIPLEN    0.206999       -0.000000
            31   XDIPLEN   XDIPLEN    0.434785        5.092646
            32   YDIPLEN   YDIPLEN    0.434785        6.390553
            33   ZDIPLEN   ZDIPLEN    0.434785        5.691196
            34   XDIPLEN   YDIPLEN    0.434785        0.000000
            35   XDIPLEN   ZDIPLEN    0.434785        0.000000
            36   YDIPLEN   ZDIPLEN    0.434785        0.000000
            37   XDIPLEN   XDIPLEN    0.964575        2.916907
            38   YDIPLEN   YDIPLEN    0.964575        3.423636
            39   ZDIPLEN   ZDIPLEN    0.964575        3.175088
            40   XDIPLEN   YDIPLEN    0.964575        0.000000
            41   XDIPLEN   ZDIPLEN    0.964575        0.000000
            42   YDIPLEN   ZDIPLEN    0.964575        0.000000
            43   XDIPLEN   XDIPLEN    2.650817        0.803666
            44   YDIPLEN   YDIPLEN    2.650817        0.842586
            45   ZDIPLEN   ZDIPLEN    2.650817        0.827618
            46   XDIPLEN   YDIPLEN    2.650817        0.000000
            47   XDIPLEN   ZDIPLEN    2.650817        0.000000
            48   YDIPLEN   ZDIPLEN    2.650817       -0.000000
            49   XDIPLEN   XDIPLEN   14.809490        0.035343
            50   YDIPLEN   YDIPLEN   14.809490        0.035496
            51   ZDIPLEN   ZDIPLEN   14.809490        0.035432
            52   XDIPLEN   YDIPLEN   14.809490        0.000000
            53   XDIPLEN   ZDIPLEN   14.809490        0.000000
            54   YDIPLEN   ZDIPLEN   14.809490       -0.000000
        """
        data_lines = raw_data.split(os.linesep)[1:-1]

        ref_c6_value = 36.230454

        self.run_c6(inpfile, potfile, xcfun_label, data_lines, ref_c6_value)

    def test_c6_dft(self):

        inpfile = os.path.join('inputs', 'water.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        potfile = None

        xcfun_label = 'b3lyp'

        #   --------------------------------------------------
        #   No    A-oper    B-oper   Frequency       Real part
        #   --------------------------------------------------
        raw_data = """
             1   XDIPLEN   XDIPLEN    0.000000        8.769176
             2   YDIPLEN   YDIPLEN    0.000000        9.696176
             3   ZDIPLEN   ZDIPLEN    0.000000        9.066348
             4   XDIPLEN   YDIPLEN    0.000000        0.000000
             5   XDIPLEN   ZDIPLEN    0.000000       -0.000000
             6   YDIPLEN   ZDIPLEN    0.000000       -0.000000
             7   XDIPLEN   XDIPLEN    0.006077        8.767399
             8   YDIPLEN   YDIPLEN    0.006077        9.695275
             9   ZDIPLEN   ZDIPLEN    0.006077        9.065171
            10   XDIPLEN   YDIPLEN    0.006077        0.000000
            11   XDIPLEN   ZDIPLEN    0.006077       -0.000000
            12   YDIPLEN   ZDIPLEN    0.006077       -0.000000
            13   XDIPLEN   XDIPLEN    0.033952        8.714428
            14   YDIPLEN   YDIPLEN    0.033952        9.668164
            15   ZDIPLEN   ZDIPLEN    0.033952        9.029850
            16   XDIPLEN   YDIPLEN    0.033952        0.000000
            17   XDIPLEN   ZDIPLEN    0.033952       -0.000000
            18   YDIPLEN   ZDIPLEN    0.033952       -0.000000
            19   XDIPLEN   XDIPLEN    0.093305        8.388972
            20   YDIPLEN   YDIPLEN    0.093305        9.489670
            21   ZDIPLEN   ZDIPLEN    0.093305        8.803011
            22   XDIPLEN   YDIPLEN    0.093305        0.000000
            23   XDIPLEN   ZDIPLEN    0.093305       -0.000000
            24   YDIPLEN   ZDIPLEN    0.093305       -0.000000
            25   XDIPLEN   XDIPLEN    0.206999        7.367332
            26   YDIPLEN   YDIPLEN    0.206999        8.777294
            27   ZDIPLEN   ZDIPLEN    0.206999        7.978451
            28   XDIPLEN   YDIPLEN    0.206999        0.000000
            29   XDIPLEN   ZDIPLEN    0.206999        0.000000
            30   YDIPLEN   ZDIPLEN    0.206999       -0.000000
            31   XDIPLEN   XDIPLEN    0.434785        5.465732
            32   YDIPLEN   YDIPLEN    0.434785        6.781963
            33   ZDIPLEN   ZDIPLEN    0.434785        6.057655
            34   XDIPLEN   YDIPLEN    0.434785        0.000000
            35   XDIPLEN   ZDIPLEN    0.434785        0.000000
            36   YDIPLEN   ZDIPLEN    0.434785       -0.000000
            37   XDIPLEN   XDIPLEN    0.964575        2.970734
            38   YDIPLEN   YDIPLEN    0.964575        3.481899
            39   ZDIPLEN   ZDIPLEN    0.964575        3.229928
            40   XDIPLEN   YDIPLEN    0.964575        0.000000
            41   XDIPLEN   ZDIPLEN    0.964575        0.000000
            42   YDIPLEN   ZDIPLEN    0.964575       -0.000000
            43   XDIPLEN   XDIPLEN    2.650817        0.798664
            44   YDIPLEN   YDIPLEN    2.650817        0.839972
            45   ZDIPLEN   ZDIPLEN    2.650817        0.825605
            46   XDIPLEN   YDIPLEN    2.650817        0.000000
            47   XDIPLEN   ZDIPLEN    2.650817        0.000000
            48   YDIPLEN   ZDIPLEN    2.650817       -0.000000
            49   XDIPLEN   XDIPLEN   14.809490        0.035158
            50   YDIPLEN   YDIPLEN   14.809490        0.035549
            51   ZDIPLEN   ZDIPLEN   14.809490        0.035467
            52   XDIPLEN   YDIPLEN   14.809490        0.000000
            53   XDIPLEN   ZDIPLEN   14.809490        0.000000
            54   YDIPLEN   ZDIPLEN   14.809490        0.000000
        """
        data_lines = raw_data.split(os.linesep)[1:-1]

        ref_c6_value = 42.315098

        self.run_c6(inpfile, potfile, xcfun_label, data_lines, ref_c6_value)

    def test_c6_dft_slda(self):

        inpfile = os.path.join('inputs', 'water.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        potfile = None

        xcfun_label = 'slda'

        #   --------------------------------------------------
        #   No    A-oper    B-oper   Frequency       Real part
        #   --------------------------------------------------
        raw_data = """
             1   XDIPLEN   XDIPLEN    0.000000        9.287658 
             2   YDIPLEN   YDIPLEN    0.000000        9.980757 
             3   ZDIPLEN   ZDIPLEN    0.000000        9.449695 
             4   XDIPLEN   YDIPLEN    0.000000        0.000000 
             5   XDIPLEN   ZDIPLEN    0.000000       -0.000000 
             6   YDIPLEN   ZDIPLEN    0.000000        0.000000 
             7   XDIPLEN   XDIPLEN    0.006077        9.285571 
             8   YDIPLEN   YDIPLEN    0.006077        9.979787 
             9   ZDIPLEN   ZDIPLEN    0.006077        9.448365 
            10   XDIPLEN   YDIPLEN    0.006077        0.000000 
            11   XDIPLEN   ZDIPLEN    0.006077       -0.000000 
            12   YDIPLEN   ZDIPLEN    0.006077        0.000000 
            13   XDIPLEN   XDIPLEN    0.033952        9.223464 
            14   YDIPLEN   YDIPLEN    0.033952        9.950588 
            15   ZDIPLEN   ZDIPLEN    0.033952        9.408517 
            16   XDIPLEN   YDIPLEN    0.033952        0.000000 
            17   XDIPLEN   ZDIPLEN    0.033952       -0.000000 
            18   YDIPLEN   ZDIPLEN    0.033952        0.000000 
            19   XDIPLEN   XDIPLEN    0.093305        8.844825 
            20   YDIPLEN   YDIPLEN    0.093305        9.758587 
            21   ZDIPLEN   ZDIPLEN    0.093305        9.153665 
            22   XDIPLEN   YDIPLEN    0.093305        0.000000 
            23   XDIPLEN   ZDIPLEN    0.093305       -0.000000 
            24   YDIPLEN   ZDIPLEN    0.093305        0.000000 
            25   XDIPLEN   XDIPLEN    0.206999        7.685719 
            26   YDIPLEN   YDIPLEN    0.206999        8.996280 
            27   ZDIPLEN   ZDIPLEN    0.206999        8.241285 
            28   XDIPLEN   YDIPLEN    0.206999        0.000000 
            29   XDIPLEN   ZDIPLEN    0.206999        0.000000 
            30   YDIPLEN   ZDIPLEN    0.206999        0.000000 
            31   XDIPLEN   XDIPLEN    0.434785        5.615926 
            32   YDIPLEN   YDIPLEN    0.434785        6.892611 
            33   ZDIPLEN   ZDIPLEN    0.434785        6.180450 
            34   XDIPLEN   YDIPLEN    0.434785        0.000000 
            35   XDIPLEN   ZDIPLEN    0.434785        0.000000 
            36   YDIPLEN   ZDIPLEN    0.434785        0.000000 
            37   XDIPLEN   XDIPLEN    0.964575        3.007170 
            38   YDIPLEN   YDIPLEN    0.964575        3.500276 
            39   ZDIPLEN   ZDIPLEN    0.964575        3.253097 
            40   XDIPLEN   YDIPLEN    0.964575        0.000000 
            41   XDIPLEN   ZDIPLEN    0.964575        0.000000 
            42   YDIPLEN   ZDIPLEN    0.964575       -0.000000 
            43   XDIPLEN   XDIPLEN    2.650817        0.800788
            44   YDIPLEN   YDIPLEN    2.650817        0.840180 
            45   ZDIPLEN   ZDIPLEN    2.650817        0.826151 
            46   XDIPLEN   YDIPLEN    2.650817        0.000000 
            47   XDIPLEN   ZDIPLEN    2.650817        0.000000 
            48   YDIPLEN   ZDIPLEN    2.650817        0.000000 
            49   XDIPLEN   XDIPLEN   14.809490        0.035212 
            50   YDIPLEN   YDIPLEN   14.809490        0.035603 
            51   ZDIPLEN   ZDIPLEN   14.809490        0.035510 
            52   XDIPLEN   YDIPLEN   14.809490        0.000000 
            53   XDIPLEN   ZDIPLEN   14.809490        0.000000 
            54   YDIPLEN   ZDIPLEN   14.809490        0.000000 
        """
        data_lines = raw_data.split(os.linesep)[1:-1]

        ref_c6_value = 44.494604

        self.run_c6(inpfile, potfile, xcfun_label, data_lines, ref_c6_value)

#    def test_c6_hf_pe(self):
#
#        inpfile = os.path.join('inputs', 'pe_water.inp')
#        if not os.path.isfile(inpfile):
#            inpfile = os.path.join('python_tests', inpfile)
#
#        potfile = os.path.join('inputs', 'pe_water.pot')
#        if not os.path.isfile(potfile):
#            potfile = os.path.join('python_tests', potfile)
#
#        xcfun_label = None
#
#        #   --------------------------------------------------
#        #   No    A-oper    B-oper   Frequency       Real part
#        #   --------------------------------------------------
#        raw_data = """
#             1   XDIPLEN   XDIPLEN    0.000000       10.314350
#             2   YDIPLEN   YDIPLEN    0.000000       11.223308
#             3   ZDIPLEN   ZDIPLEN    0.000000        9.892232
#             4   XDIPLEN   YDIPLEN    0.000000        0.199246
#             5   XDIPLEN   ZDIPLEN    0.000000        0.233161
#             6   YDIPLEN   ZDIPLEN    0.000000       -0.390806
#             7   XDIPLEN   XDIPLEN    0.006077       10.312720
#             8   YDIPLEN   YDIPLEN    0.006077       11.222201
#             9   ZDIPLEN   ZDIPLEN    0.006077        9.891287
#            10   XDIPLEN   YDIPLEN    0.006077        0.199106
#            11   XDIPLEN   ZDIPLEN    0.006077        0.233080
#            12   YDIPLEN   ZDIPLEN    0.006077       -0.390662
#            13   XDIPLEN   XDIPLEN    0.033952       10.263898
#            14   YDIPLEN   YDIPLEN    0.033952       11.188871
#            15   ZDIPLEN   ZDIPLEN    0.033952        9.862879
#            16   XDIPLEN   YDIPLEN    0.033952        0.194927
#            17   XDIPLEN   ZDIPLEN    0.033952        0.230674
#            18   YDIPLEN   ZDIPLEN    0.033952       -0.386346
#            19   XDIPLEN   XDIPLEN    0.093305        9.953362
#            20   YDIPLEN   YDIPLEN    0.093305       10.969883
#            21   ZDIPLEN   ZDIPLEN    0.093305        9.677249
#            22   XDIPLEN   YDIPLEN    0.093305        0.169257
#            23   XDIPLEN   ZDIPLEN    0.093305        0.215332
#            24   YDIPLEN   ZDIPLEN    0.093305       -0.358584
#            25   XDIPLEN   XDIPLEN    0.206999        8.864034
#            26   YDIPLEN   YDIPLEN    0.206999       10.102401
#            27   ZDIPLEN   ZDIPLEN    0.206999        8.956787
#            28   XDIPLEN   YDIPLEN    0.206999        0.092132
#            29   XDIPLEN   ZDIPLEN    0.206999        0.161361
#            30   YDIPLEN   ZDIPLEN    0.206999       -0.258727
#            31   XDIPLEN   XDIPLEN    0.434785        6.509070
#            32   YDIPLEN   YDIPLEN    0.434785        7.710062
#            33   ZDIPLEN   ZDIPLEN    0.434785        7.028854
#            34   XDIPLEN   YDIPLEN    0.434785       -0.005216
#            35   XDIPLEN   ZDIPLEN    0.434785        0.051524
#            36   YDIPLEN   ZDIPLEN    0.434785       -0.062669
#            37   XDIPLEN   XDIPLEN    0.964575        3.336900
#            38   YDIPLEN   YDIPLEN    0.964575        3.839739
#            39   ZDIPLEN   ZDIPLEN    0.964575        3.753804
#            40   XDIPLEN   YDIPLEN    0.964575       -0.023506
#            41   XDIPLEN   ZDIPLEN    0.964575       -0.031440
#            42   YDIPLEN   ZDIPLEN    0.964575        0.044480
#            43   XDIPLEN   XDIPLEN    2.650817        0.825249
#            44   YDIPLEN   YDIPLEN    2.650817        0.878798
#            45   ZDIPLEN   ZDIPLEN    2.650817        0.895013
#            46   XDIPLEN   YDIPLEN    2.650817       -0.005278
#            47   XDIPLEN   ZDIPLEN    2.650817       -0.010019
#            48   YDIPLEN   ZDIPLEN    2.650817        0.010647
#            49   XDIPLEN   XDIPLEN   14.809490        0.034853
#            50   YDIPLEN   YDIPLEN   14.809490        0.035934
#            51   ZDIPLEN   ZDIPLEN   14.809490        0.036380
#            52   XDIPLEN   YDIPLEN   14.809490       -0.000186
#            53   XDIPLEN   ZDIPLEN   14.809490       -0.000140
#            54   YDIPLEN   ZDIPLEN   14.809490        0.000165
#        """
#        data_lines = raw_data.split(os.linesep)[1:-1]
#
#        ref_c6_value =
#
#        self.run_c6(inpfile, potfile, xcfun_label, data_lines, ref_c6_value)

#    def test_c6_dft_pe(self):
#
#        inpfile = os.path.join('inputs', 'pe_water.inp')
#        if not os.path.isfile(inpfile):
#            inpfile = os.path.join('python_tests', inpfile)
#
#        potfile = os.path.join('inputs', 'pe_water.pot')
#        if not os.path.isfile(potfile):
#            potfile = os.path.join('python_tests', potfile)
#
#        xcfun_label = 'b3lyp'
#
#        #   --------------------------------------------------
#        #   No    A-oper    B-oper   Frequency       Real part
#        #   --------------------------------------------------
#        raw_data = """
#             1   XDIPLEN   XDIPLEN    0.000000       13.747730
#             2   YDIPLEN   YDIPLEN    0.000000       13.210995
#             3   ZDIPLEN   ZDIPLEN    0.000000       11.668275
#             4   XDIPLEN   YDIPLEN    0.000000        0.705252
#             5   XDIPLEN   ZDIPLEN    0.000000        0.057873
#             6   YDIPLEN   ZDIPLEN    0.000000       -0.873707
#             7   XDIPLEN   XDIPLEN    0.006077       13.743406
#             8   YDIPLEN   YDIPLEN    0.006077       13.209203
#             9   ZDIPLEN   ZDIPLEN    0.006077       11.666562
#            10   XDIPLEN   YDIPLEN    0.006077        0.704596
#            11   XDIPLEN   ZDIPLEN    0.006077        0.057975
#            12   YDIPLEN   ZDIPLEN    0.006077       -0.873325
#            13   XDIPLEN   XDIPLEN    0.033952       13.615031
#            14   YDIPLEN   YDIPLEN    0.033952       13.155435
#            15   ZDIPLEN   ZDIPLEN    0.033952       11.615240
#            16   XDIPLEN   YDIPLEN    0.033952        0.685204
#            17   XDIPLEN   ZDIPLEN    0.033952        0.060974
#            18   YDIPLEN   ZDIPLEN    0.033952       -0.861917
#            19   XDIPLEN   XDIPLEN    0.093305       12.846446
#            20   YDIPLEN   YDIPLEN    0.093305       12.808454
#            21   ZDIPLEN   ZDIPLEN    0.093305       11.288109
#            22   XDIPLEN   YDIPLEN    0.093305        0.573038
#            23   XDIPLEN   ZDIPLEN    0.093305        0.076666
#            24   YDIPLEN   ZDIPLEN    0.093305       -0.790688
#            25   XDIPLEN   XDIPLEN    0.206999       10.629163
#            26   YDIPLEN   YDIPLEN    0.206999       11.515804
#            27   ZDIPLEN   ZDIPLEN    0.206999       10.128299
#            28   XDIPLEN   YDIPLEN    0.206999        0.294992
#            29   XDIPLEN   ZDIPLEN    0.206999        0.097340
#            30   YDIPLEN   ZDIPLEN    0.206999       -0.559007
#            31   XDIPLEN   XDIPLEN    0.434785        7.110153
#            32   YDIPLEN   YDIPLEN    0.434785        8.352163
#            33   ZDIPLEN   ZDIPLEN    0.434785        7.528445
#            34   XDIPLEN   YDIPLEN    0.434785        0.032203
#            35   XDIPLEN   ZDIPLEN    0.434785        0.050712
#            36   YDIPLEN   ZDIPLEN    0.434785       -0.175284
#            37   XDIPLEN   XDIPLEN    0.964575        3.404707
#            38   YDIPLEN   YDIPLEN    0.964575        3.925597
#            39   ZDIPLEN   ZDIPLEN    0.964575        3.823051
#            40   XDIPLEN   YDIPLEN    0.964575       -0.024632
#            41   XDIPLEN   ZDIPLEN    0.964575       -0.026606
#            42   YDIPLEN   ZDIPLEN    0.964575        0.027592
#            43   XDIPLEN   XDIPLEN    2.650817        0.821089
#            44   YDIPLEN   YDIPLEN    2.650817        0.876836
#            45   ZDIPLEN   ZDIPLEN    2.650817        0.894607
#            46   XDIPLEN   YDIPLEN    2.650817       -0.005882
#            47   XDIPLEN   ZDIPLEN    2.650817       -0.010986
#            48   YDIPLEN   ZDIPLEN    2.650817        0.010017
#            49   XDIPLEN   XDIPLEN   14.809490        0.034741
#            50   YDIPLEN   YDIPLEN   14.809490        0.036012
#            51   ZDIPLEN   ZDIPLEN   14.809490        0.036532
#            52   XDIPLEN   YDIPLEN   14.809490       -0.000192
#            53   XDIPLEN   ZDIPLEN   14.809490       -0.000235
#            54   YDIPLEN   ZDIPLEN   14.809490        0.000156
#        """
#        data_lines = raw_data.split(os.linesep)[1:-1]
#
#        ref_c6_value =
#
#        self.run_c6(inpfile, potfile, xcfun_label, data_lines, ref_c6_value)


if __name__ == "__main__":
    unittest.main()
