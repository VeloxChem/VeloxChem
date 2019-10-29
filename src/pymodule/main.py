from mpi4py import MPI
import sys
import os

from .veloxchemlib import mpi_initialized
from .veloxchemlib import mpi_master
from .veloxchemlib import CudaDevices
from .mpitask import MpiTask
from .scfrestdriver import ScfRestrictedDriver
from .scfunrestdriver import ScfUnrestrictedDriver
from .mointsdriver import MOIntegralsDriver
from .rsplinabscross import LinearAbsorptionCrossSection
from .rspcdspec import CircularDichroismSpectrum
from .rsppolarizability import Polarizability
from .rspabsorption import Absorption
from .cppsolver import ComplexResponse
from .pulsedrsp import PulsedResponse
from .mp2driver import Mp2Driver
from .adconedriver import AdcOneDriver
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
        print(os.linesep.join(info_txt), file=sys.stdout)
        sys.exit(0)

    # set up MPI task

    task = MpiTask(sys.argv[1:], MPI.COMM_WORLD)

    task_types = task.input_dict['jobs']['task'].lower().split(',')
    task_types = [x.strip() for x in task_types]

    if 'method_settings' in task.input_dict:
        method_dict = task.input_dict['method_settings']
    else:
        method_dict = {}

    # initialize CUDA capable devices

    gpu_devs = CudaDevices()
    if task.mpi_rank == mpi_master():
        if gpu_devs.get_number_devices() > 0:
            task.ostream.print_blank()
            task.ostream.print_block(str(gpu_devs))

    # Exciton model

    if 'exciton' in task_types:

        if 'exciton' in task.input_dict:
            exciton_dict = task.input_dict['exciton']
        else:
            exciton_dict = {}

        exciton_drv = ExcitonModelDriver(task.mpi_comm, task.ostream)
        exciton_drv.update_settings(exciton_dict, method_dict)
        exciton_drv.compute(task.molecule, task.ao_basis, task.min_basis)

    # Hartree-Fock

    run_scf = True in [
        x in [
            'hf', 'rhf', 'uhf', 'scf', 'wavefunction', 'wave function', 'mp2',
            'visualization', 'response', 'cpp', 'adc1'
        ] for x in task_types
    ]

    run_unrestricted = 'uhf' in task_types

    if run_scf:

        # initialize scf driver and run scf

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

        # molecular orbitals

        mol_orbs = scf_drv.mol_orbs
        density = scf_drv.density
        scf_tensors = scf_drv.scf_tensors

        # tranform integrals to MO basis
        if 'ao2mo' in scf_dict and scf_drv.restricted:

            moints_drv = MOIntegralsDriver(task.mpi_comm, task.ostream)

            grps = [p for p in range(task.mpi_comm.Get_size())]
            moints = moints_drv.compute(task.molecule, task.ao_basis, mol_orbs,
                                        scf_dict['ao2mo'].upper(), grps)

            # sketch for transforming MO integrals batches to antisymmetrized
            # integrals
            # Indexing scheme: occupied orbitals from 0..nocc
            #                  virtual orbitals from nocc..ntot
            #
            # Select here external indexes space (one need to generalize this)
            # if  mintstype == "OVOV":
            #     bra_idx = (0, nocc)
            #     ket_idx = (nocc, nocc + nvirt)
            # Assuming we add to tensor ten[i,j,k,l]
            # for idx, pair in enumerate(moints.get_gen_pairs()):
            #
            #   i = pair.first()
            #   j = pair.second()
            #
            #   fxy = moints.xy_to_numpy()
            #   fyx = moints.yx_to_numpy()
            #
            #    for k in bra_idx:
            #        for l in ket_idx:
            #           # aaaa, bbbb blocks
            #           ten_aaaa[i,j,k,l] = fxy[k,l] - fyx[l,k]
            #           # abab, baba
            #           ten_abab[i,j,k,l] = fxy[k,l]
            #           # abba, baba
            #           ten_abba[i,j,k,l] = -fyx[l,k]

    # Response

    if 'response' in task_types and scf_drv.restricted:

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
        else:
            assert_msg_critical(False, 'input file: invalid response property')

        rsp_prop.init_driver(task.mpi_comm, task.ostream)
        rsp_prop.compute(task.molecule, task.ao_basis, scf_tensors)
        if task.mpi_rank == mpi_master():
            rsp_prop.print_property(task.ostream)

    # Complex Response

    if 'cpp' in task_types and scf_drv.restricted:

        if 'cpp' in task.input_dict:
            cpp_dict = task.input_dict['cpp']
        else:
            cpp_dict = {}

        cpp_drv = ComplexResponse(task.mpi_comm, task.ostream)
        cpp_drv.update_settings(cpp_dict, method_dict)
        results = cpp_drv.compute(task.molecule, task.ao_basis, scf_tensors)
        if task.mpi_rank == mpi_master():
            cpp_drv.print_properties(results['properties'])

    # Pulsed Linear Response Theory

    if 'pulses' in task.input_dict and scf_drv.restricted:
        prt_dict = task.input_dict['pulses']
        cpp_dict = {}

        pulsed_response = PulsedResponse(task.mpi_comm, task.ostream)
        pulsed_response.update_settings(prt_dict, cpp_dict)
        pulsed_response.compute(task.molecule, task.ao_basis, scf_tensors)

    # ADC(1)

    if 'adc1' in task_types and scf_drv.restricted:

        if 'adc' in task.input_dict:
            adc_dict = task.input_dict['adc']
        else:
            adc_dict = {}

        adc1_drv = AdcOneDriver(task.mpi_comm, task.ostream)
        adc1_drv.update_settings(adc_dict)
        adc1_drv.compute(task.molecule, task.ao_basis, scf_tensors)

    # MP2 perturbation theory

    if 'mp2' in task_types and scf_drv.restricted:

        if 'mp2' in task.input_dict:
            mp2_dict = task.input_dict['mp2']
        else:
            mp2_dict = {}

        mp2_drv = Mp2Driver(task.mpi_comm, task.ostream)
        mp2_drv.update_settings(mp2_dict)
        mp2_drv.compute(task.molecule, task.ao_basis, mol_orbs)

    # Cube file

    if 'visualization' in task_types:

        if 'visualization' in task.input_dict:
            cube_dict = task.input_dict['visualization']
        else:
            cube_dict = {}

        mol_orbs.broadcast(task.mpi_rank, task.mpi_comm)
        density.broadcast(task.mpi_rank, task.mpi_comm)

        vis_drv = VisualizationDriver(task.mpi_comm)
        vis_drv.gen_cubes(cube_dict, task.molecule, task.ao_basis, mol_orbs,
                          density)

    # all done, print finish header to output stream

    task.finish()
