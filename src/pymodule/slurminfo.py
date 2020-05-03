import subprocess
import os
import re


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


def get_slurm_maximum_hours():
    """
    Gets SLURM timelimit in hours.

    :return:
        The SLURM timelimit in hours (None if not a SLURM job).
    """

    job_id = get_slurm_job_id()

    if job_id is not None and executable_in_path('scontrol'):
        scontrol_output = get_command_output(
            ['scontrol', 'show', 'job', job_id])

        if scontrol_output is not None:
            for line in scontrol_output.split(os.linesep):
                if 'TimeLimit=' in line:
                    match = re.search(
                        r'TimeLimit=(\d+-)?(\d{2}):(\d{2}):(\d{2})', line)

                    if match is not None:
                        if match.group(1) is not None:
                            days = match.group(1).replace('-', '')
                        else:
                            days = 0
                        hours = match.group(2)
                        minutes = match.group(3)
                        seconds = match.group(4)

                        try:
                            max_hours = (int(days) * 24 + int(hours) +
                                         int(minutes) / 60 +
                                         int(seconds) / 3600)
                        except ValueError:
                            max_hours = None

                        return max_hours

    return None
