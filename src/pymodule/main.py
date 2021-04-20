#
#                           VELOXCHEM 1.0-RC
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
#  Contact: https://veloxchem.org/contact
#
#  SPDX-License-Identifier: LGPL-3.0-or-later
#
#  This file is part of VeloxChem.
#
#  VeloxChem is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Lesser General Public License as published by the Free
#  Software Foundation, either version 3 of the License, or (at your option)
#  any later version.
#
#  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

from mpi4py import MPI
import time as tm

from .veloxchemlib import mpi_initialized
from .veloxchemlib import mpi_master
from .cli import cli
from .errorhandler import assert_msg_critical
from .excitondriver import ExcitonModelDriver
from .loprop import LoPropDriver
from .mp2driver import Mp2Driver
from .mpitask import MpiTask
from .optimizationdriver import OptimizationDriver
from .pulsedrsp import PulsedResponse
from .rspabsorption import Absorption
from .rspc6 import C6
from .rspcdspec import CircularDichroismSpectrum
from .rspcustomproperty import CustomProperty
from .rsplinabscross import LinearAbsorptionCrossSection
from .rsppolarizability import Polarizability
from .rsptpa import TPA
from .scffirstorderprop import ScfFirstOrderProperties
from .scfgradientdriver import ScfGradientDriver
from .scfrestdriver import ScfRestrictedDriver
from .scfunrestdriver import ScfUnrestrictedDriver
from .slurminfo import get_slurm_maximum_hours
from .visualizationdriver import VisualizationDriver
from .xtbdriver import XTBDriver
from .xtbgradientdriver import XTBGradientDriver
from .tdhfgradientdriver import TdhfGradientDriver


def select_scf_driver(task, scf_type):
    """
    Selects SCF driver.

    :param task:
        The MPI task.
    :param scf_type:
        The type of SCF calculation (restricted or unrestricted).

    :return:
        The SCF driver object.
    """

    # check number of MPI nodes
    if task.mpi_rank == mpi_master():
        n_ao = task.ao_basis.get_dimensions_of_basis(task.molecule)
        assert_msg_critical(task.mpi_size == 1 or task.mpi_size <= n_ao,
                            'SCF: too many MPI processes')

    nalpha = task.molecule.number_of_alpha_electrons()
    nbeta = task.molecule.number_of_beta_electrons()

    if nalpha == nbeta and scf_type == 'restricted':
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
    else:
        scf_drv = ScfUnrestrictedDriver(task.mpi_comm, task.ostream)

    return scf_drv


def select_rsp_property(task, mol_orbs, rsp_dict, method_dict):
    """
    Selects response property.

    :param task:
        The MPI task.
    :param mol_orbs:
        The molecular orbitals.
    :param rsp_dict:
        The dictionary of response dict.
    :param method_dict:
        The dictionary of method settings.

    :return:
        The response property object.
    """

    # check number of MPI nodes
    if task.mpi_rank == mpi_master():
        nocc = task.molecule.number_of_alpha_electrons()
        n_ov = nocc * (mol_orbs.number_mos() - nocc)
        assert_msg_critical(task.mpi_size == 1 or task.mpi_size <= n_ov,
                            'Response: too many MPI processes')

    if 'property' in rsp_dict:
        prop_type = rsp_dict['property'].lower()
    else:
        prop_type = 'custom'

    if prop_type in [
            'polarizability',
            'dipole polarizability',
    ]:
        rsp_prop = Polarizability(rsp_dict, method_dict)

    elif prop_type in [
            'absorption',
            'uv-vis',
            'ecd',
    ]:
        rsp_prop = Absorption(rsp_dict, method_dict)

    elif prop_type in [
            'linear absorption cross-section',
            'linear absorption (cpp)',
            'absorption (cpp)',
    ]:
        rsp_prop = LinearAbsorptionCrossSection(rsp_dict, method_dict)

    elif prop_type in [
            'circular dichroism spectrum',
            'circular dichroism (cpp)',
            'ecd (cpp)',
    ]:
        rsp_prop = CircularDichroismSpectrum(rsp_dict, method_dict)

    elif prop_type == 'c6':
        rsp_prop = C6(rsp_dict, method_dict)

    elif prop_type == 'custom':
        rsp_prop = CustomProperty(rsp_dict, method_dict)

    elif prop_type == 'tpa':
        rsp_prop = TPA(rsp_dict, method_dict)

    else:
        assert_msg_critical(False, 'input file: invalid response property')

    return rsp_prop


def main():
    """
    Runs VeloxChem with command line arguments.
    """

    program_start_time = tm.time()

    assert_msg_critical(mpi_initialized(), "MPI not initialized")

    # Parse command line

    parser = cli()
    args = parser.parse_args()

    # MPI task

    task = MpiTask([args.input_file, args.output_file], MPI.COMM_WORLD)
    task_type = task.input_dict['jobs']['task'].lower()

    # Timelimit in hours

    if 'maximum_hours' in task.input_dict['jobs']:
        maximum_hours = float(task.input_dict['jobs']['maximum_hours'])
    else:
        maximum_hours = get_slurm_maximum_hours()

    # Method settings

    if 'method_settings' in task.input_dict:
        method_dict = dict(task.input_dict['method_settings'])
    else:
        method_dict = {}

    if 'pe' in task.input_dict:
        # add @pe group to method_dict as pe_options
        method_dict['pe_options'] = dict(task.input_dict['pe'])
    else:
        method_dict['pe_options'] = {}

    use_xtb = ('xtb' in method_dict)

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
        'hf', 'rhf', 'uhf', 'scf', 'uscf', 'wavefunction', 'wave function',
        'mp2', 'gradient', 'optimize', 'response', 'pulses', 'visualization',
        'loprop'
    ]

    if task_type == 'visualization' and 'visualization' in task.input_dict:
        run_scf = 'read_dalton' not in task.input_dict['visualization']['cubes']

    scf_type = 'unrestricted' if task_type in ['uhf', 'uscf'] else 'restricted'

    if run_scf:
        if 'scf' in task.input_dict:
            scf_dict = task.input_dict['scf']
        else:
            scf_dict = {}

        scf_dict['program_start_time'] = program_start_time
        scf_dict['maximum_hours'] = maximum_hours

        if use_xtb:
            if 'potfile' in method_dict:
                errmsg = 'XTBDriver: The \'potfile\' keyword is not supported '
                errmsg += 'in XTB calculation.'
                if task.mpi_rank == mpi_master():
                    assert_msg_critical(False, errmsg)
            xtb_drv = XTBDriver(task.mpi_comm)
            xtb_drv.set_method(method_dict['xtb'].lower())
            xtb_drv.compute(task.molecule, task.ostream)
        else:
            scf_drv = select_scf_driver(task, scf_type)
            scf_drv.update_settings(scf_dict, method_dict)
            scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

            mol_orbs = scf_drv.mol_orbs
            density = scf_drv.density
            scf_tensors = scf_drv.scf_tensors

            if not scf_drv.is_converged:
                return

            # SCF first-order properties
            scf_prop = ScfFirstOrderProperties(task.mpi_comm, task.ostream)
            scf_prop.compute(task.molecule, task.ao_basis, scf_tensors)
            if task.mpi_rank == mpi_master():
                scf_prop.print_properties(task.molecule)

            if (scf_drv.electric_field is not None and
                    task.molecule.get_charge() != 0):
                task.finish()
                return

    # Gradient

    if task_type == 'gradient':
        if use_xtb:
            grad_drv = XTBGradientDriver(xtb_drv, task.mpi_comm, task.ostream)
            grad_drv.compute(task.molecule)
        elif scf_drv.closed_shell:
            grad_drv = ScfGradientDriver(scf_drv, task.mpi_comm, task.ostream)
            grad_drv.compute(task.molecule, task.ao_basis, task.min_basis)

    # Geometry optimization

    if task_type == 'optimize':
        if 'potfile' in method_dict:
            errmsg = 'OptimizationDriver: The \'potfile\' keyword is not '
            errmsg += 'supported in geometry optimization.'
            if task.mpi_rank == mpi_master():
                assert_msg_critical(False, errmsg)

        if 'optimize' in task.input_dict:
            opt_dict = task.input_dict['optimize']
        else:
            opt_dict = {}

        if use_xtb:
            grad_drv = XTBGradientDriver(xtb_drv, task.mpi_comm, task.ostream)
            opt_drv = OptimizationDriver(task.input_dict['filename'], grad_drv,
                                         'XTB')
        elif scf_drv.closed_shell:
            grad_drv = ScfGradientDriver(scf_drv, task.mpi_comm, task.ostream)
            opt_drv = OptimizationDriver(task.input_dict['filename'], grad_drv,
                                         'SCF')

        opt_drv.update_settings(opt_dict)
        opt_drv.compute(task.molecule, task.ao_basis, task.min_basis)

    # Response

    if task_type == 'response' and scf_drv.closed_shell:
        if 'response' in task.input_dict:
            rsp_dict = dict(task.input_dict['response'])
        else:
            rsp_dict = {}

        rsp_dict['filename'] = task.input_dict['filename']

        rsp_dict['program_start_time'] = program_start_time
        rsp_dict['maximum_hours'] = maximum_hours

        if 'eri_thresh' not in rsp_dict:
            rsp_dict['eri_thresh'] = scf_drv.eri_thresh
        if 'qq_type' not in rsp_dict:
            rsp_dict['qq_type'] = scf_drv.qq_type
        if not scf_drv.restart:
            rsp_dict['restart'] = 'no'

        rsp_prop = select_rsp_property(task, mol_orbs, rsp_dict, method_dict)
        rsp_prop.init_driver(task.mpi_comm, task.ostream)
        rsp_prop.compute(task.molecule, task.ao_basis, scf_tensors)
        if not rsp_prop.converged():
            return

        if task.mpi_rank == mpi_master():
            rsp_prop.print_property(task.ostream)

        # Calculate the excited-state gradient (for now only
        # CIS relaxed dipole moments)
        if 'n_state_deriv' in rsp_dict:
            tdhfgrad_drv = TdhfGradientDriver(task.mpi_comm, task.ostream)
            tdhfgrad_drv.update_settings(rsp_dict, method_dict)
            tdhfgrad_drv.compute(task.molecule, task.ao_basis, scf_tensors,
                                rsp_prop.rsp_property)

    # Pulsed Linear Response Theory

    if ((task_type == 'pulses' or 'pulses' in task.input_dict) and
            scf_drv.closed_shell):
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
        pulsed_response.update_settings(prt_dict, cpp_dict, method_dict)
        pulsed_response.compute(task.molecule, task.ao_basis, scf_tensors)

    # MP2 perturbation theory

    if task_type == 'mp2' and scf_drv.closed_shell:
        if 'mp2' in task.input_dict:
            mp2_dict = task.input_dict['mp2']
        else:
            mp2_dict = {}

        if 'eri_thresh' not in mp2_dict:
            mp2_dict['eri_thresh'] = scf_drv.eri_thresh
        if 'qq_type' not in mp2_dict:
            mp2_dict['qq_type'] = scf_drv.qq_type

        mp2_drv = Mp2Driver(task.mpi_comm, task.ostream)
        mp2_drv.update_settings(mp2_dict, method_dict)
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
