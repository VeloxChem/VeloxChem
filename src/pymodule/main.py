#
#                           VELOXCHEM 1.0-RC2
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

from .veloxchemlib import mpi_initialized, mpi_master
from .mpitask import MpiTask
from .scfrestdriver import ScfRestrictedDriver
from .scfunrestdriver import ScfUnrestrictedDriver
from .scfrestopendriver import ScfRestrictedOpenDriver
from .scffirstorderprop import ScfFirstOrderProperties
from .respchargesdriver import RespChargesDriver
from .excitondriver import ExcitonModelDriver
from .mp2driver import Mp2Driver
from .loprop import LoPropDriver
from .trajectorydriver import TrajectoryDriver
from .scfgradientdriver import ScfGradientDriver
from .optimizationdriver import OptimizationDriver
from .pulsedrsp import PulsedResponse
from .rsppolarizability import Polarizability
from .rspabsorption import Absorption
from .rsplinabscross import LinearAbsorptionCrossSection
from .rspcdspec import CircularDichroismSpectrum
from .rspc6 import C6
from .rspshg import SHG
from .rsptpa import TPA
from .rspcustomproperty import CustomProperty
from .visualizationdriver import VisualizationDriver
from .xtbdriver import XTBDriver
from .xtbgradientdriver import XTBGradientDriver
from .cli import cli
from .errorhandler import assert_msg_critical
from .slurminfo import get_slurm_maximum_hours


def select_scf_driver(task, scf_type):
    """
    Selects SCF driver.

    :param task:
        The MPI task.
    :param scf_type:
        The type of SCF calculation (restricted, unrestricted, or restricted
        open-shell).

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

    if scf_type == 'restricted' and nalpha == nbeta:
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
    elif scf_type == 'restricted' and nalpha != nbeta:
        scf_drv = ScfUnrestrictedDriver(task.mpi_comm, task.ostream)
    elif scf_type == 'unrestricted':
        scf_drv = ScfUnrestrictedDriver(task.mpi_comm, task.ostream)
    elif scf_type == 'restricted_openshell':
        scf_drv = ScfRestrictedOpenDriver(task.mpi_comm, task.ostream)

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

    elif prop_type == 'tpa':
        rsp_prop = TPA(rsp_dict, method_dict)

    elif prop_type == 'shg':
        rsp_prop = SHG(rsp_dict, method_dict)

    elif prop_type == 'custom':
        rsp_prop = CustomProperty(rsp_dict, method_dict)

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
        if task.mpi_rank == mpi_master():
            maximum_hours = get_slurm_maximum_hours()
        else:
            maximum_hours = None
        maximum_hours = task.mpi_comm.bcast(maximum_hours, root=mpi_master())

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
        exciton_dict['filename'] = task.input_dict['filename']

        exciton_drv = ExcitonModelDriver(task.mpi_comm, task.ostream)
        exciton_drv.update_settings(exciton_dict, method_dict)
        exciton_drv.compute(task.molecule, task.ao_basis, task.min_basis)

    # Spectrum from trajectory

    if task_type == 'trajectory':
        if 'trajectory' in task.input_dict:
            traj_dict = dict(task.input_dict['trajectory'])
        else:
            traj_dict = {}
        if 'spectrum_settings' in task.input_dict:
            spect_dict = dict(task.input_dict['spectrum_settings'])
        else:
            spect_dict = {}
        if 'response' in task.input_dict:
            rsp_dict = dict(task.input_dict['response'])
        else:
            rsp_dict = {}

        traj_dict['filename'] = task.input_dict['filename']
        traj_dict['charges'] = task.input_dict['charges']
        traj_dict['polarizabilities'] = task.input_dict['polarizabilities']

        traj_drv = TrajectoryDriver(task.mpi_comm, task.ostream)
        traj_drv.update_settings(traj_dict, spect_dict, rsp_dict, method_dict)
        traj_drv.compute(task.molecule, task.ao_basis, task.min_basis)

    # Self-consistent field

    run_scf = task_type in [
        'hf', 'rhf', 'uhf', 'rohf', 'scf', 'uscf', 'roscf', 'wavefunction',
        'wave function', 'mp2', 'gradient', 'optimize', 'response', 'pulses',
        'visualization', 'loprop'
    ]

    if task_type == 'visualization' and 'visualization' in task.input_dict:
        run_scf = 'read_dalton' not in task.input_dict['visualization']['cubes']

    scf_type = 'restricted'
    if task_type in ['uhf', 'uscf']:
        scf_type = 'unrestricted'
    elif task_type in ['rohf', 'roscf']:
        scf_type = 'restricted_openshell'

    if run_scf:
        assert_msg_critical(task.molecule.number_of_atoms(),
                            'Molecule: no atoms found in molecule')

        if 'scf' in task.input_dict:
            scf_dict = task.input_dict['scf']
        else:
            scf_dict = {}

        scf_dict['program_start_time'] = program_start_time
        scf_dict['maximum_hours'] = maximum_hours
        scf_dict['filename'] = task.input_dict['filename']

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
        elif scf_drv.scf_type == 'restricted':
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

        opt_dict['filename'] = task.input_dict['filename']

        if use_xtb:
            grad_drv = XTBGradientDriver(xtb_drv, task.mpi_comm, task.ostream)
            opt_drv = OptimizationDriver(grad_drv, 'XTB')
        elif scf_drv.scf_type == 'restricted':
            grad_drv = ScfGradientDriver(scf_drv, task.mpi_comm, task.ostream)
            opt_drv = OptimizationDriver(grad_drv, 'SCF')

        opt_drv.update_settings(opt_dict)
        opt_drv.compute(task.molecule, task.ao_basis, task.min_basis)

    # Response

    if task_type == 'response' and scf_drv.scf_type == 'restricted':
        if 'response' in task.input_dict:
            rsp_dict = dict(task.input_dict['response'])
        else:
            rsp_dict = {}

        rsp_dict['program_start_time'] = program_start_time
        rsp_dict['maximum_hours'] = maximum_hours
        rsp_dict['filename'] = task.input_dict['filename']

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

    # Pulsed Linear Response Theory

    if ((task_type == 'pulses' or 'pulses' in task.input_dict) and
            scf_drv.scf_type == 'restricted'):
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

    if task_type == 'mp2' and scf_drv.scf_type == 'restricted':
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
        loprop_driver = LoPropDriver(task.mpi_comm, task.ostream)
        loprop_driver.compute(task.molecule, task.ao_basis, scf_tensors)

    # RESP and ESP charges

    if task_type in ['resp charges', 'esp charges']:
        if (task_type == 'resp charges' and 'resp_charges' in task.input_dict):
            charges_dict = task.input_dict['resp_charges']
        elif (task_type == 'esp charges' and 'esp_charges' in task.input_dict):
            charges_dict = task.input_dict['esp_charges']
        else:
            charges_dict = {}

        charges_dict['filename'] = task.input_dict['filename']

        chg_drv = RespChargesDriver(task.mpi_comm, task.ostream)
        chg_drv.update_settings(charges_dict, method_dict)

        if task_type == 'resp charges':
            chg_drv.compute(task.molecule, task.ao_basis, 'resp')
        elif task_type == 'esp charges':
            chg_drv.compute(task.molecule, task.ao_basis, 'esp')

    # All done

    task.finish()
