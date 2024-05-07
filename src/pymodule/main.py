#
#                              VELOXCHEM
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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
from datetime import datetime, timedelta

from .veloxchemlib import mpi_initialized, mpi_master
from .mpitask import MpiTask
from .scfrestdriver import ScfRestrictedDriver
from .cli import cli
from .errorhandler import assert_msg_critical


def select_scf_driver(task, scf_type):
    """
    Selects SCF driver.

    :param task:
        The MPI task.
    :param scf_type:
        The type of SCF calculation (restricted, unrestricted, or
        restricted_openshell).

    :return:
        The SCF driver object.
    """

    # check number of MPI nodes
    if task.mpi_rank == mpi_master():
        n_ao = task.ao_basis.get_dimensions_of_basis()
        assert_msg_critical(task.mpi_size == 1 or task.mpi_size <= n_ao,
                            'SCF: too many MPI processes')

    nalpha = task.molecule.number_of_alpha_electrons()
    nbeta = task.molecule.number_of_beta_electrons()

    if scf_type == 'restricted' and nalpha == nbeta:
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
    else:
        assert_msg_critical(False, f'SCF: invalide scf_type {scf_type}')

    return scf_drv


def main():
    """
    Runs VeloxChem with command line arguments.
    """

    program_start_time = datetime.now()

    assert_msg_critical(mpi_initialized(), "MPI not initialized")

    # Parse command line

    parser = cli()
    args = parser.parse_args()

    # MPI task

    task = MpiTask([args.input_file, args.output_file], MPI.COMM_WORLD)
    task_type = task.input_dict['jobs']['task'].lower()

    # Timelimit

    if 'maximum_hours' in task.input_dict['jobs']:
        maximum_hours = float(task.input_dict['jobs']['maximum_hours'])
        program_end_time = program_start_time + timedelta(hours=maximum_hours)
    else:
        program_end_time = None

    # Method settings

    method_dict = (dict(task.input_dict['method_settings'])
                   if 'method_settings' in task.input_dict else {})

    # Self-consistent field

    run_scf = task_type in ['hf', 'rhf', 'scf', 'wavefunction', 'wave function']

    scf_type = 'restricted'

    if run_scf:
        assert_msg_critical(task.molecule.number_of_atoms(),
                            'Molecule: no atoms found in molecule')

        scf_dict = (dict(task.input_dict['scf'])
                    if 'scf' in task.input_dict else {})
        scf_dict['program_end_time'] = program_end_time
        scf_dict['filename'] = task.input_dict['filename']

        scf_drv = select_scf_driver(task, scf_type)
        scf_drv.update_settings(scf_dict, method_dict)
        scf_results = scf_drv.compute(task.molecule, task.ao_basis,
                                      task.min_basis)

        if not scf_drv.is_converged:
            return

    # All done

    task.finish()
