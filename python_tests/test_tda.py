from pathlib import Path

import numpy as np
import pytest

from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tdaeigensolver import TdaEigenSolver
from veloxchem.veloxchemlib import hartree_in_ev, is_mpi_master

from .addons import using_cppe


@pytest.mark.solvers
class TestTDA:

    def run_tda(self, inpfile, potfile, xcfun_label, data_lines):

        task = MpiTask([inpfile, None])
        task.input_dict['scf']['checkpoint_file'] = None

        if potfile is not None:
            task.input_dict['method_settings']['potfile'] = potfile

        if xcfun_label is not None:
            task.input_dict['method_settings']['xcfun'] = xcfun_label

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_results = scf_drv.compute(task.molecule, task.ao_basis,
                                      task.min_basis)

        ref_exc_ene = [float(line.split()[1]) for line in data_lines]
        ref_osc_str = [float(line.split()[3]) for line in data_lines]
        ref_rot_str = [float(line.split()[4]) for line in data_lines]

        tda_solver = TdaEigenSolver(task.mpi_comm, task.ostream)
        tda_solver.update_settings({'nstates': len(ref_exc_ene)},
                                   task.input_dict['method_settings'])
        tda_results = tda_solver.compute(task.molecule, task.ao_basis,
                                         scf_results)

        if is_mpi_master(task.mpi_comm):
            exc_ene = tda_results['eigenvalues'] * hartree_in_ev()
            osc_str = tda_results['oscillator_strengths']
            rot_str = tda_results['rotatory_strengths']

            assert np.max(np.abs(exc_ene - ref_exc_ene)) < 5.0e-4
            assert np.max(np.abs(osc_str - ref_osc_str)) < 5.0e-4
            assert np.max(np.abs(rot_str - ref_rot_str)) < 1.0e-2

    def test_tda_hf(self):

        # vlxtag: RHF, Absorption, ECD, CIS

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water.inp')

        potfile = None

        xcfun_label = None

        #  State Frequency   Oscillator Strength    Rotatory  Strength
        #          (eV)      Velocity     Length    Velocity    Length
        #  -----------------------------------------------------------
        raw_data = """
            1     8.7847     0.0695     0.0530     0.0000     0.0000
            2    10.4740     0.0000     0.0000    -0.0000    -0.0000
            3    11.1634     0.0869     0.1070    -0.0000    -0.0000
            4    12.1889     0.0014     0.0045    -0.0000    -0.0000
            5    12.8274     0.0117     0.0250     0.0000     0.0000
        """
        data_lines = raw_data.splitlines()[1:-1]

        self.run_tda(inpfile, potfile, xcfun_label, data_lines)

    def test_tda_dft(self):

        # vlxtag: RKS, Absorption, ECD, TDA

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water.inp')

        potfile = None

        xcfun_label = 'b3lyp'

        #  State Frequency   Oscillator Strength    Rotatory  Strength
        #          (eV)      Velocity     Length    Velocity    Length
        #  -----------------------------------------------------------
        raw_data = """
            1     6.9987     0.0571     0.0537     0.0000     0.0000
            2     8.4291     0.0000     0.0000    -0.0000    -0.0000
            3     9.2135     0.0611     0.0906     0.0000     0.0000
            4    10.3003     0.0005     0.0000    -0.0000    -0.0000
            5    10.6270     0.0044     0.0127    -0.0000    -0.0000
        """
        data_lines = raw_data.splitlines()[1:-1]

        self.run_tda(inpfile, potfile, xcfun_label, data_lines)

    def test_tda_dft_slda(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water.inp')

        potfile = None

        xcfun_label = 'slda'

        #  State Frequency   Oscillator Strength    Rotatory  Strength
        #          (eV)      Velocity     Length    Velocity    Length
        #  -----------------------------------------------------------
        raw_data = """
            1     6.7828     0.0588     0.0561    -0.0000     0.0000
            2     8.2221     0.0000     0.0000     0.0000     0.0000
            3     8.9101     0.0603     0.0901    -0.0000    -0.0000
            4    10.1323     0.0014     0.0003     0.0000     0.0000
            5    10.3444     0.0036     0.0115    -0.0000    -0.0000
        """
        data_lines = raw_data.splitlines()[1:-1]

        self.run_tda(inpfile, potfile, xcfun_label, data_lines)

    @using_cppe
    def test_tda_hf_pe(self):

        # vlxtag: RHF, Absorption, ECD, CIS, PE

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'pe_water.inp')
        potfile = str(here / 'inputs' / 'pe_water.pot')

        xcfun_label = None

        #  State Frequency   Oscillator Strength    Rotatory  Strength
        #          (eV)      Velocity     Length    Velocity    Length
        #  -----------------------------------------------------------
        raw_data = """
            1     9.1467     0.0758     0.0602     0.1440     0.3887
            2    11.1635     0.0113     0.0121    -0.4940    -0.2295
            3    11.4450     0.1062     0.1234    -0.1259     2.8161
            4    11.9792     0.0004     0.0013     0.0916     0.1963
            5    12.8941     0.0007     0.0006     0.0132    -0.4366
        """
        data_lines = raw_data.splitlines()[1:-1]

        self.run_tda(inpfile, potfile, xcfun_label, data_lines)

    @using_cppe
    def test_tda_dft_pe(self):

        # vlxtag: RKS, Absorption, ECD, TDA, PE

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'pe_water.inp')
        potfile = str(here / 'inputs' / 'pe_water.pot')

        xcfun_label = 'b3lyp'

        #  State Frequency   Oscillator Strength    Rotatory  Strength
        #          (eV)      Velocity     Length    Velocity    Length
        #  -----------------------------------------------------------
        raw_data = """
            1     7.3252     0.0609     0.0582     0.1219     0.3536
            2     9.0068     0.0067     0.0073    -0.0937     0.1897
            3     9.5211     0.0742     0.1000    -0.2837     2.3062
            4    10.1853     0.0003     0.0003    -0.1075    -0.2451
            5    11.0450     0.0076     0.0029     0.1461    -0.5130
        """
        data_lines = raw_data.splitlines()[1:-1]

        self.run_tda(inpfile, potfile, xcfun_label, data_lines)
