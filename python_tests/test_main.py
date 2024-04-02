from pathlib import Path
from unittest.mock import patch
import sys

from veloxchem.veloxchemlib import is_mpi_master, mpi_barrier
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.rspabsorption import Absorption
from veloxchem.rspc6 import C6
from veloxchem.rspcdspec import CircularDichroismSpectrum
from veloxchem.rsplinabscross import LinearAbsorptionCrossSection
from veloxchem.rsppolarizability import Polarizability
from veloxchem.rsptpa import TPA
from veloxchem.inputparser import get_random_string_parallel
from veloxchem.main import select_scf_driver, select_rsp_property, main


class TestMain:

    def create_input_file(self, lines, fname):

        if is_mpi_master():
            with open(fname, 'w') as f_inp:
                for line in lines.splitlines():
                    if line:
                        f_inp.write(f'{line}\n')
        mpi_barrier()

    def run_main_function(self, capsys, input_lines, input_fname, output):

        here = Path(__file__).parent
        random_string = get_random_string_parallel()
        input_file = here / 'inputs' / f'{input_fname}_{random_string}.inp'

        self.create_input_file(input_lines, input_file)
        with patch.object(sys, 'argv', ['vlx', str(input_file)]):
            main()
            captured = capsys.readouterr()
            if is_mpi_master():
                for line in output:
                    assert line in captured.out
            self.remove_input_and_h5_files(input_file)

    def remove_input_and_h5_files(self, input_file):

        if is_mpi_master():
            if input_file.is_file():
                input_file.unlink()

            scf_h5 = input_file.with_suffix('.scf.h5')
            if scf_h5.is_file():
                scf_h5.unlink()

            scf_final_h5 = scf_h5.with_suffix('.results.h5')
            if scf_final_h5.is_file():
                scf_final_h5.unlink()

            rsp_h5 = input_file.with_suffix('.rsp.h5')
            if rsp_h5.is_file():
                rsp_h5.unlink()

            rsp_solutions_h5 = rsp_h5.with_suffix('.solutions.h5')
            if rsp_solutions_h5.is_file():
                rsp_solutions_h5.unlink()

        mpi_barrier()

    def test_select_scf_driver(self):

        # restricted

        here = Path(__file__).parent
        inpfile = here / 'inputs' / 'water.inp'

        task = MpiTask([str(inpfile), None])
        scf_type = 'restricted'

        scf_drv = select_scf_driver(task, scf_type)
        assert isinstance(scf_drv, ScfRestrictedDriver)

        # unrestricted

        here = Path(__file__).parent
        inpfile = here / 'inputs' / 'water_triplet.inp'

        task = MpiTask([str(inpfile), None])
        scf_type = 'unrestricted'

        scf_drv = select_scf_driver(task, scf_type)
        assert isinstance(scf_drv, ScfUnrestrictedDriver)

    def test_select_rsp_driver(self):

        here = Path(__file__).parent
        inpfile = here / 'inputs' / 'water.inp'

        task = MpiTask([str(inpfile), None])
        task.input_dict['scf']['checkpoint_file'] = None

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        method_dict = {}

        rsp_dict = {'property': 'polarizability'}
        rsp_prop = select_rsp_property(task, scf_drv.mol_orbs, rsp_dict,
                                       method_dict)
        assert isinstance(rsp_prop, Polarizability)

        rsp_dict = {'property': 'absorption'}
        rsp_prop = select_rsp_property(task, scf_drv.mol_orbs, rsp_dict,
                                       method_dict)
        assert isinstance(rsp_prop, Absorption)

        rsp_dict = {'property': 'absorption (cpp)'}
        rsp_prop = select_rsp_property(task, scf_drv.mol_orbs, rsp_dict,
                                       method_dict)
        assert isinstance(rsp_prop, LinearAbsorptionCrossSection)

        rsp_dict = {'property': 'ecd (cpp)'}
        rsp_prop = select_rsp_property(task, scf_drv.mol_orbs, rsp_dict,
                                       method_dict)
        assert isinstance(rsp_prop, CircularDichroismSpectrum)

        rsp_dict = {'property': 'c6'}
        rsp_prop = select_rsp_property(task, scf_drv.mol_orbs, rsp_dict,
                                       method_dict)
        assert isinstance(rsp_prop, C6)

        rsp_dict = {'property': 'tpa'}
        rsp_prop = select_rsp_property(task, scf_drv.mol_orbs, rsp_dict,
                                       method_dict)
        assert isinstance(rsp_prop, TPA)

    def test_main(self, capsys):

        scf_job_input = """
            @jobs
            task: scf
            maximum hours: 1
            @end
        """

        rsp_job_input = """
            @jobs
            task: response
            maximum hours: 1
            @end
        """

        method_input = """
            @method settings
            basis: aug-cc-pvdz
            @end
        """

        rpa_input = """
            @response
            property: absorption
            nstates: 3
            @end
        """

        cpp_input = """
            @response
            property: absorption (cpp)
            frequencies: 0-0.25 (0.05)
            @end
        """

        molecule_input = """
            @molecule
            charge: 0
            multiplicity: 1
            units: au
            xyz:
            O   0.0   0.0   0.0
            H   0.0   1.4   1.1
            H   0.0  -1.4   1.1
            @end
        """

        scf_output = [
            'VeloxChem execution started',
            'SCF converged',
            'Total Energy                       :',
            'Electronic Energy                  :',
            'Nuclear Repulsion Energy           :',
            'VeloxChem execution completed',
        ]

        rpa_output = [
            'Linear Response EigenSolver Setup',
            'Linear response converged',
            'Electric Transition Dipole Moments (dipole length, a.u.)',
            'Electric Transition Dipole Moments (dipole velocity, a.u.)',
            'Magnetic Transition Dipole Moments (a.u.)',
            'One-Photon Absorption',
            'Electronic Circular Dichroism',
        ]

        cpp_output = [
            'Complex Response Solver Setup',
            'Complex response converged',
            'Response Functions at Given Frequencies',
            'Linear Absorption Cross-Section',
        ]

        # scf test

        scf_input_lines = (scf_job_input + method_input + molecule_input)

        self.run_main_function(capsys, scf_input_lines, 'vlx_scf_test',
                               scf_output)

        # rpa test

        rpa_input_lines = (rsp_job_input + method_input + rpa_input +
                           molecule_input)

        self.run_main_function(capsys, rpa_input_lines, 'vlx_rpa_test',
                               scf_output + rpa_output)

        # cpp test

        cpp_input_lines = (rsp_job_input + method_input + cpp_input +
                           molecule_input)

        self.run_main_function(capsys, cpp_input_lines, 'vlx_cpp_test',
                               scf_output + cpp_output)
