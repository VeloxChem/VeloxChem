from mpi4py import MPI
import time as tm
import sys
import os

from .veloxchemlib import mpi_initialized
from .veloxchemlib import mpi_master
from .mpitask import MpiTask
from .scfrestdriver import ScfRestrictedDriver
from .scfunrestdriver import ScfUnrestrictedDriver
from .gradientdriver import GradientDriver
from .optimizationdriver import OptimizationDriver
from .rsplinabscross import LinearAbsorptionCrossSection
from .rspcdspec import CircularDichroismSpectrum
from .rsppolarizability import Polarizability
from .rspabsorption import Absorption
from .rspc6 import C6
from .rspcustomproperty import CustomProperty
from .pulsedrsp import PulsedResponse
from .mp2driver import Mp2Driver
from .excitondriver import ExcitonModelDriver
from .visualizationdriver import VisualizationDriver
from .loprop import LoPropDriver
from .errorhandler import assert_msg_critical
from .slurminfo import get_slurm_maximum_hours
from .xtbdriver import XTBDriver



def main():

    program_start_time = tm.time()

    assert_msg_critical(mpi_initialized(), "MPI: Initialized")

    if len(sys.argv) <= 1 or sys.argv[1] in ['-h', '--help']:
        info_txt = [
            '',
            '=================   VeloxChem   =================',
            '',
            'Usage:',
            '    python3 -m veloxchem input_file [output_file]',
            '',
        ]
        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            print(os.linesep.join(info_txt), file=sys.stdout)
        sys.exit(0)

    # MPI task

    task = MpiTask(sys.argv[1:], MPI.COMM_WORLD)
    task_type = task.input_dict['jobs']['task'].lower()

    if 'method_settings' in task.input_dict:
        method_dict = task.input_dict['method_settings']
    else:
        method_dict = {}

    if 'pe' in task.input_dict:
        method_dict['pe_options'] = dict(task.input_dict['pe'])
    else:
        method_dict['pe_options'] = {}

    use_xtb = 'xtb' in method_dict

    # Timelimit in hours

    if 'maximum_hours' in task.input_dict['jobs']:
        try:
            maximum_hours = float(task.input_dict['jobs']['maximum_hours'])
        except ValueError:
            maximum_hours = None
    else:
        maximum_hours = get_slurm_maximum_hours()

    # Exciton model

    if task_type == 'exciton':
        if 'exciton' in task.input_dict:
            exciton_dict = task.input_dict['exciton']
        else:
            exciton_dict = {}

        exciton_dict['program_start_time'] = program_start_time
        exciton_dict['maximum_hours'] = maximum_hours

        exciton_drv = ExcitonModelDriver(task.mpi_comm, task.ostream)
        exciton_drv.update_settings(exciton_dict, method_dict)
        exciton_drv.compute(task.molecule, task.ao_basis, task.min_basis)

    # Self-consistent field

    run_scf = task_type in [
        'hf', 'rhf', 'uhf', 'scf', 'wavefunction', 'wave function', 'mp2',
        'gradient', 'optimize', 'response', 'pulses', 'visualization', 'loprop'
    ]

    if task_type == 'visualization' and 'visualization' in task.input_dict:
        run_scf = 'read_dalton' not in task.input_dict['visualization']['cubes']

    run_unrestricted = (task_type == 'uhf')

    if run_scf:
        if 'scf' in task.input_dict:
            scf_dict = task.input_dict['scf']
        else:
            scf_dict = {}

        scf_dict['program_start_time'] = program_start_time
        scf_dict['maximum_hours'] = maximum_hours
         
        if use_xtb:
            xtb_drv = XTBDriver(task.mpi_comm)
            xtb_drv.compute_energy(task.molecule, scf_dict, method_dict, 
                                   task.ostream)
        else:  
            nalpha = task.molecule.number_of_alpha_electrons()
            nbeta = task.molecule.number_of_beta_electrons()

            if nalpha == nbeta and not run_unrestricted:
                scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
            else:
                scf_drv = ScfUnrestrictedDriver(task.mpi_comm, task.ostream)
            scf_drv.update_settings(scf_dict, method_dict)
            scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

            mol_orbs = scf_drv.mol_orbs
            density = scf_drv.density
            scf_tensors = scf_drv.scf_tensors

            if not scf_drv.is_converged:
                return

    # Gradient

    if task_type == 'gradient' and scf_drv.restricted:
        grad_drv = GradientDriver(task.mpi_comm, task.ostream)
        grad_drv.update_settings(scf_dict, method_dict)
        grad_drv.compute(task.molecule, task.ao_basis, task.min_basis)

    # Geometry optimization

    if task_type == 'optimize' and scf_drv.restricted:
        if 'optimize' in task.input_dict:
            opt_dict = task.input_dict['optimize']
        else:
            opt_dict = {}
        opt_drv = OptimizationDriver(task.mpi_comm, task.ostream)
        opt_drv.update_settings(opt_dict, scf_dict, method_dict)
        opt_drv.compute(task.molecule, task.ao_basis, task.min_basis)

    # Response

    if task_type == 'response' and scf_drv.restricted:
        if 'response' in task.input_dict:
            rsp_dict = task.input_dict['response']
        else:
            rsp_dict = {}

        rsp_dict['program_start_time'] = program_start_time
        rsp_dict['maximum_hours'] = maximum_hours

        if 'eri_thresh' not in rsp_dict:
            rsp_dict['eri_thresh'] = scf_drv.eri_thresh
        if 'qq_type' not in rsp_dict:
            rsp_dict['qq_type'] = scf_drv.qq_type
        if not scf_drv.restart:
            rsp_dict['restart'] = 'no'

        if 'property' not in rsp_dict:
            rsp_dict['property'] = 'custom'
        prop_type = rsp_dict['property'].lower()

        if prop_type in ['polarizability', 'dipole polarizability']:
            rsp_prop = Polarizability(rsp_dict, method_dict)
        elif prop_type in ['absorption', 'uv-vis', 'ecd']:
            rsp_prop = Absorption(rsp_dict, method_dict)
        elif prop_type in [
                'linear absorption cross-section', 'linear absorption (cpp)',
                'absorption (cpp)'
        ]:
            rsp_prop = LinearAbsorptionCrossSection(rsp_dict, method_dict)
        elif prop_type in [
                'circular dichroism spectrum', 'circular dichroism (cpp)',
                'ecd (cpp)'
        ]:
            rsp_prop = CircularDichroismSpectrum(rsp_dict, method_dict)
        elif prop_type == 'c6':
            rsp_prop = C6(rsp_dict, method_dict)
        elif prop_type == 'custom':
            rsp_prop = CustomProperty(rsp_dict, method_dict)
        else:
            assert_msg_critical(False, 'input file: invalid response property')

        rsp_prop.init_driver(task.mpi_comm, task.ostream)
        rsp_prop.compute(task.molecule, task.ao_basis, scf_tensors)
        if not rsp_prop.converged():
            return

        if task.mpi_rank == mpi_master():
            rsp_prop.print_property(task.ostream)

    # Pulsed Linear Response Theory

    if ((task_type == 'pulses' or 'pulses' in task.input_dict) and
            scf_drv.restricted):
        if 'pulses' in task.input_dict:
            prt_dict = task.input_dict['pulses']
        else:
            prt_dict = {}
        cpp_dict = {}

        if 'eri_thresh' not in cpp_dict:
            cpp_dict['eri_thresh'] = scf_drv.eri_thresh
        if 'qq_type' not in cpp_dict:
            cpp_dict['qq_type'] = scf_drv.qq_type
        if not scf_drv.restart:
            cpp_dict['restart'] = 'no'

        pulsed_response = PulsedResponse(task.mpi_comm, task.ostream)
        pulsed_response.update_settings(prt_dict, cpp_dict)
        pulsed_response.compute(task.molecule, task.ao_basis, scf_tensors)

    # MP2 perturbation theory

    if task_type == 'mp2' and scf_drv.restricted:
        if 'mp2' in task.input_dict:
            mp2_dict = task.input_dict['mp2']
        else:
            mp2_dict = {}

        if 'eri_thresh' not in mp2_dict:
            mp2_dict['eri_thresh'] = scf_drv.eri_thresh
        if 'qq_type' not in mp2_dict:
            mp2_dict['qq_type'] = scf_drv.qq_type

        mp2_drv = Mp2Driver(task.mpi_comm, task.ostream)
        mp2_drv.update_settings(mp2_dict)
        mp2_drv.compute(task.molecule, task.ao_basis, mol_orbs)

    # Cube file

    if task_type == 'visualization':
        if 'visualization' in task.input_dict:
            cube_dict = task.input_dict['visualization']
        else:
            cube_dict = {}

        if 'read_dalton' not in task.input_dict['visualization']['cubes']:
            mol_orbs.broadcast(task.mpi_rank, task.mpi_comm)
            density.broadcast(task.mpi_rank, task.mpi_comm)
        else:
            mol_orbs = None
            density = None

        vis_drv = VisualizationDriver(task.mpi_comm)
        vis_drv.gen_cubes(cube_dict, task.molecule, task.ao_basis, mol_orbs,
                          density)

    # LoProp

    if task_type == 'loprop':
        task.input_dict['loprop']['density'] = density
        loprop_driver = LoPropDriver(task)
        loprop_driver.compute()

    # All done

    task.finish()
