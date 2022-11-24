#
#                           VELOXCHEM 1.0-RC3
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2022 by VeloxChem developers. All rights reserved.
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

from datetime import datetime
import subprocess
import os


def get_slurm_job_id():
    """
    Gets SLURM job ID as a string.
    """

    if 'SLURM_JOB_ID' in os.environ:
        return os.environ['SLURM_JOB_ID']

    elif 'SLURM_JOBID' in os.environ:
        return os.environ['SLURM_JOBID']

    else:
        return None


def executable_in_path(executable):
    """
    Checks if the executable is in PATH.

    :param executable:
        The executable.

    :return:
        True if the executable is in PATH, False otherwise.
    """

    for path in os.environ['PATH'].split(os.pathsep):
        fname = os.path.join(path, executable)
        if os.path.isfile(fname) and os.access(fname, os.X_OK):
            return True

    return False


def get_command_output(command):
    """
    Gets the output of a command.

    :param command:
        The command as a list of strings.

    :return:
        The output as a string.
    """

    try:
        output = subprocess.check_output(command).decode('utf-8')
    except subprocess.CalledProcessError:
        output = None

    return output


def get_slurm_end_time():
    """
    Gets SLURM end time.

    :return:
        The SLURM end time as a datetime object (None if not a SLURM job).
    """

    job_id = get_slurm_job_id()

    if job_id is not None and executable_in_path('scontrol'):
        scontrol_output = get_command_output(
            ['scontrol', 'show', 'job', job_id])

        if scontrol_output is not None:
            for line in scontrol_output.splitlines():
                if 'EndTime=' in line:
                    end_time_string = line.split('EndTime=')[1].split()[0]
                    return datetime_from_isoformat(end_time_string)

    return None


def datetime_from_isoformat(datetime_string):
    """
    Gets datetime object from string YYYY-MM-DDTHH:MM:SS.

    :return:
        The datetime object.
    """

    try:
        dt = datetime.fromisoformat(datetime_string)

    except AttributeError:
        date_string, time_string = datetime_string.split('T')
        year, month, day = date_string.split('-')
        hour, minute, second = time_string.split(':')
        dt = datetime(int(year), int(month), int(day), int(hour), int(minute),
                      int(second))

    return dt
