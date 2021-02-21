import numpy as np
import unittest
import pytest
import sys
from pathlib import Path
try:
    import cppe
except ImportError:
    pass

from veloxchem.veloxchemlib import ElectronRepulsionIntegralsDriver
from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import denmat
from veloxchem.mpitask import MpiTask
from veloxchem.aodensitymatrix import AODensityMatrix
from veloxchem.aofockmatrix import AOFockMatrix
from veloxchem.scfdriver import ScfDriver
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scffirstorderprop import ScfFirstOrderProperties
from veloxchem.qqscheme import get_qq_scheme


class TestScfRestricted(unittest.TestCase):

    def run_scf(self, inpfile, potfile, xcfun_label, ref_e_scf, ref_dip):

        task = MpiTask([inpfile, None])
        task.input_dict['scf']['checkpoint_file'] = None

        if potfile is not None:
            task.input_dict['method_settings']['potfile'] = potfile

        if xcfun_label is not None:
            task.input_dict['method_settings']['xcfun'] = xcfun_label

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        scf_prop = ScfFirstOrderProperties(task.mpi_comm, task.ostream)
        scf_prop.compute(task.molecule, task.ao_basis, scf_drv.scf_tensors)

        if task.mpi_rank == mpi_master():
            e_scf = scf_drv.get_scf_energy()
            tol = 1.0e-5 if xcfun_label is not None else 1.0e-6
            self.assertTrue(np.max(np.abs(e_scf - ref_e_scf)) < tol)

            dip = scf_prop.get_property('dipole moment')
            self.assertTrue(np.max(np.abs(dip - ref_dip)) < 1.0e-5)

    def test_scf_hf(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water.inp')

        potfile = None

        xcfun_label = None

        #    Final HF energy:             -76.041697549811
        ref_e_scf = -76.041697549811

        ref_dip = np.array([0.000000, 0.000000, 0.786770])

        self.run_scf(inpfile, potfile, xcfun_label, ref_e_scf, ref_dip)

    def test_scf_dft(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water.inp')

        potfile = None

        xcfun_label = 'b3lyp'

        #    Final DFT energy:            -76.443545741524
        ref_e_scf = -76.443545741524

        ref_dip = np.array([0.000000, 0.000000, 0.731257])

        self.run_scf(inpfile, potfile, xcfun_label, ref_e_scf, ref_dip)

    def test_scf_dft_slda(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water.inp')

        potfile = None

        xcfun_label = 'slda'

        #    Final DFT energy:            -76.074208234637
        ref_e_scf = -76.074208234637

        ref_dip = np.array([0.000000, 0.000000, 0.731291])

        self.run_scf(inpfile, potfile, xcfun_label, ref_e_scf, ref_dip)

    @pytest.mark.skipif('cppe' not in sys.modules, reason='cppe not available')
    def test_scf_hf_pe(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'pe_water.inp')
        potfile = str(here / 'inputs' / 'pe_water.pot')

        xcfun_label = None

        #    Final HF energy:             -76.067159426565
        ref_e_scf = -76.067159426565

        ref_dip = np.array([-0.039715, 0.098785, 0.945128])

        self.run_scf(inpfile, potfile, xcfun_label, ref_e_scf, ref_dip)

    @pytest.mark.skipif('cppe' not in sys.modules, reason='cppe not available')
    def test_scf_dft_pe(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'pe_water.inp')
        potfile = str(here / 'inputs' / 'pe_water.pot')

        xcfun_label = 'b3lyp'

        #    Final DFT energy:            -76.468733754150
        ref_e_scf = -76.468733754150

        ref_dip = np.array([-0.048195, 0.098715, 0.902822])

        self.run_scf(inpfile, potfile, xcfun_label, ref_e_scf, ref_dip)

    def test_comp_fock_split_comm(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water.inp')
        task = MpiTask([inpfile, None])

        mol = task.molecule
        bas = task.ao_basis
        nao = bas.get_dimensions_of_basis(mol)

        dmat = np.diag(np.ones(nao))
        dens = AODensityMatrix([dmat], denmat.rest)
        fock = AOFockMatrix(dens)

        eri_drv = ElectronRepulsionIntegralsDriver(task.mpi_comm)
        screening = eri_drv.compute(get_qq_scheme('QQ_DEN'), 1.0e-12, mol, bas)

        solver = ScfDriver(task.mpi_comm, task.ostream)
        solver.comp_2e_fock_split_comm(fock, dens, mol, bas, screening)

        self.assertEqual(fock.alpha_to_numpy(0).shape, dmat.shape)


if __name__ == "__main__":
    unittest.main()
