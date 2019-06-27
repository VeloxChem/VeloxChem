from mpi4py import MPI
import sys
import os

from .veloxchemlib import mpi_initialized
from .veloxchemlib import mpi_master
from .mpitask import MpiTask
from .scfrestdriver import ScfRestrictedDriver
from .scfunrestdriver import ScfUnrestrictedDriver
from .rsplinabscross import LinearAbsorptionCrossSection
from .rspcdspec import CircularDichroismSpectrum
from .rsppolarizability import Polarizability
from .rspabsorption import Absorption
from .rspc6 import C6
from .pulsedrsp import PulsedResponse
from .mp2driver import Mp2Driver
from .excitondriver import ExcitonModelDriver
from .visualizationdriver import VisualizationDriver
from .errorhandler import assert_msg_critical


def main():

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

    # Exciton model

    if task_type == 'exciton':
        if 'exciton' in task.input_dict:
            exciton_dict = task.input_dict['exciton']
        else:
            exciton_dict = {}

        exciton_drv = ExcitonModelDriver(task.mpi_comm, task.ostream)
        exciton_drv.update_settings(exciton_dict, method_dict)
        exciton_drv.compute(task.molecule, task.ao_basis, task.min_basis)

    # Self-consistent field

    run_scf = task_type in [
        'hf', 'rhf', 'uhf', 'scf', 'wavefunction', 'wave function', 'mp2',
        'response', 'visualization'
    ]

    if task_type == 'visualization' and 'visualization' in task.input_dict:
        run_scf = 'read_dalton' not in task.input_dict['visualization']['cubes']

    run_unrestricted = (task_type == 'uhf')

    if run_scf:
        if 'scf' in task.input_dict:
            scf_dict = task.input_dict['scf']
        else:
            scf_dict = {}

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

    # Response

    if task_type == 'response' and scf_drv.restricted:
        if 'response' in task.input_dict:
            rsp_dict = task.input_dict['response']
        else:
            rsp_dict = {}

        if 'eri_thresh' not in rsp_dict:
            rsp_dict['eri_thresh'] = scf_drv.eri_thresh
        if 'qq_type' not in rsp_dict:
            rsp_dict['qq_type'] = scf_drv.qq_type

        assert_msg_critical('property' in rsp_dict,
                            'input file: response property not found')
        prop_type = rsp_dict['property'].lower()

        if prop_type == 'polarizability':
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
        else:
            assert_msg_critical(False, 'input file: invalid response property')

        rsp_prop.init_driver(task.mpi_comm, task.ostream)
        rsp_prop.compute(task.molecule, task.ao_basis, scf_tensors)
        if task.mpi_rank == mpi_master():
            rsp_prop.print_property(task.ostream)

    # Pulsed Linear Response Theory

    if 'pulses' in task.input_dict and scf_drv.restricted:
        prt_dict = task.input_dict['pulses']
        cpp_dict = {}

        pulsed_response = PulsedResponse(task.mpi_comm, task.ostream)
        pulsed_response.update_settings(prt_dict, cpp_dict)
        pulsed_response.compute(task.molecule, task.ao_basis, scf_tensors)

    # MP2 perturbation theory

    if task_type == 'mp2' and scf_drv.restricted:
        if 'mp2' in task.input_dict:
            mp2_dict = task.input_dict['mp2']
        else:
            mp2_dict = {}

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
        pass

    # All done

    task.finish()
