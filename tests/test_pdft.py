from mpi4py import MPI
import numpy as np

from veloxchem.veloxchemlib import XCIntegrator, XCPairDensityFunctional
from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.griddriver import GridDriver
from veloxchem.scfrestopendriver import ScfRestrictedOpenDriver
from veloxchem.aodensitymatrix import AODensityMatrix


class TestPDFT:

    def test_pfunc(self):

        pfunc = XCPairDensityFunctional('TSLATER', ['TSLATER'], [1.0])
        assert pfunc.get_func_label() == 'TSLATER'
        assert pfunc.get_func_type() == 'PLDA'

        pfunc = XCPairDensityFunctional('TLDA', ['TSLATER', 'TVWN'], [1.0, 1.0])
        assert pfunc.get_func_label() == 'TLDA'
        assert pfunc.get_func_type() == 'PLDA'

        pfunc = XCPairDensityFunctional('TPBE', ['TPBE_X', 'TPBE_C'],
                                        [1.0, 1.0])
        assert pfunc.get_func_label() == 'TPBE'
        assert pfunc.get_func_type() == 'PGGA'

        pfunc = XCPairDensityFunctional('TBLYP', ['TSLATER', 'TB88', 'TLYP'],
                                        [1.0, 1.0, 1.0])
        assert pfunc.get_func_label() == 'TBLYP'
        assert pfunc.get_func_type() == 'PGGA'

        pfunc2 = XCPairDensityFunctional('TSLATER', ['TSLATER'], [1.0])
        assert pfunc != pfunc2

        pfunc2 = XCPairDensityFunctional('TBLYP', ['TSLATER', 'TB88', 'TLYP'],
                                         [1.0, 1.0, 1.0])
        assert pfunc == pfunc2
