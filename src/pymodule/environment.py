#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from os import environ, cpu_count
from pathlib import Path


def get_basis_path():
    """
    Returns location of basis files within module.

    :return:
        The location of basis files within module.
    """

    return Path(__file__).resolve().parent / "basis"


def set_vlxbasispath():
    """
    Sets location of basis files.
    """

    if 'VLXBASISPATH' not in environ:
        environ['VLXBASISPATH'] = str(get_basis_path())


def get_data_path():
    """
    Returns location of data files within module.

    :return:
        The location of data files within module.
    """

    return Path(__file__).resolve().parent / "database"


def set_vlxdatapath():
    """
    Sets location of data files.
    """

    if 'VLXDATAPATH' not in environ:
        environ['VLXDATAPATH'] = str(get_data_path())


def set_omp_num_threads(ncores=None):
    """
    Sets number of OpenMP threads.

    :param ncores:
        The number of cores available for OpenMP threads.
    """

    if ncores is not None:
        environ['OMP_NUM_THREADS'] = str(ncores)
    else:
        if 'OMP_NUM_THREADS' not in environ:
            ncores = cpu_count()
            if ncores is None:
                ncores = 1
            environ['OMP_NUM_THREADS'] = str(ncores)


def get_available_memory():
    """
    Gets the memory this process may still take, in bytes.

    :return:
        The free memory of the node, or less where a control group the process
        is in leaves less of its limit free.
    """

    try:
        import psutil
        available = psutil.virtual_memory().available
    except ImportError:
        available = 8 * 1024**3

    # NOTE: a batch system which shares a node between jobs confines each of
    # them to a control group, and the free memory of the node says nothing of
    # that limit: a job given a share of a node was promised the whole of it,
    # and the kernel killed it on reaching its share.
    confined = _get_cgroup_available_memory()

    if confined is not None:
        available = min(available, confined)

    return available


def get_ri_free_memory(comm, rank_reserve):
    """
    Gets the memory free to a rank for a resolution of the identity driver, less
    a reserve for the node and one for each rank on it, in bytes. Every rank of
    the communicator must call this.

    :param comm:
        The communicator, whose ranks sharing a node divide its memory.
    :param rank_reserve:
        The memory each rank keeps back for what the calculation allocates after
        this is measured, in bytes.

    :return:
        The memory in bytes, zero where the reserves leave nothing.
    """

    # NOTE: what is free to this process, which is less than what is free on
    # the node where a batch system confines the job to a share of it. What
    # the process already holds is not free, so it needs no reserve.
    available = get_available_memory()

    # NOTE: the ranks of one node see the same free memory and would each claim
    # the whole of it, so the node would be promised several times over. The
    # ranks sharing a node are counted by their host name and the memory is
    # divided between them.
    import platform
    here = platform.node()
    on_this_node = max(comm.allgather(here).count(here), 1)

    # NOTE: two reserves. The node keeps a few gigabytes for itself, and every
    # rank keeps what the rest of the calculation allocates after this is
    # measured.
    #
    # NOTE: the node's reserve is a guess at what is used outside the
    # calculation, and what the operating system holds is already not counted
    # as available, so on a small machine it is kept to a quarter of what is
    # free rather than a fixed amount which would leave nothing. The ranks'
    # reserves are what the calculation itself will allocate and are never
    # cut: where they leave nothing, the budget is zero and the build falls
    # back to the four-center integrals, rather than RI being handed memory
    # the rest of the calculation needs.
    node_reserve = min(4 * 1024**3, available // 4)

    free = available - node_reserve - on_this_node * rank_reserve

    return int(max(free, 0) / on_this_node)


def _get_cgroup_available_memory():
    """
    Gets the memory the control groups of this process leave free, in bytes.

    :return:
        The least any of them leaves free, or None where none of them has a
        memory limit.
    """

    try:
        lines = Path('/proc/self/cgroup').read_text().splitlines()
    except OSError:
        return None

    least = None

    for line in lines:
        parts = line.split(':', 2)
        if len(parts) != 3:
            continue
        hierarchy, controllers, path = parts

        # NOTE: version 2 is the one hierarchy numbered 0 with no controllers
        # named, and version 1 names memory among its own. A node may mount
        # both, and only the one which carries the memory controller has the
        # files read below; the other is passed over for want of them.
        if hierarchy == '0' and controllers == '':
            root = Path('/sys/fs/cgroup')
            names = ('memory.max', 'memory.current', 'inactive_file')
        elif 'memory' in controllers.split(','):
            root = Path('/sys/fs/cgroup/memory')
            names = ('memory.limit_in_bytes', 'memory.usage_in_bytes',
                     'total_inactive_file')
        else:
            continue

        # NOTE: the limit may be set on any group above this one, and Slurm sets
        # it on the job and the step rather than on the task the process is in,
        # so the groups are walked up to the root and the least is taken.
        group = root / path.lstrip('/')

        while True:
            free = _get_cgroup_free_memory(group, *names)
            if free is not None:
                least = free if least is None else min(least, free)
            if group == root or root not in group.parents:
                break
            group = group.parent

    return least


def _get_cgroup_free_memory(group, limit_name, usage_name, inactive_name):
    """
    Gets the memory one control group leaves free, in bytes.

    :param group:
        The directory of the group.
    :param limit_name:
        The file of its limit.
    :param usage_name:
        The file of its usage.
    :param inactive_name:
        The entry of memory.stat which counts its inactive file cache.

    :return:
        The memory free under its limit, or None where it has no limit.
    """

    try:
        limit = (group / limit_name).read_text().strip()
        usage = int((group / usage_name).read_text())
    except (OSError, ValueError):
        return None

    # NOTE: no limit is written as max by version 2, and as a number near the
    # largest a 64-bit signed integer holds by version 1.
    if limit == 'max' or int(limit) >= 2**62:
        return None

    # NOTE: the usage counts the file cache, of which the inactive part is
    # given back before the kernel kills anything, so it is not taken as used.
    inactive = 0
    try:
        for entry in (group / 'memory.stat').read_text().splitlines():
            name, value = entry.split()
            if name == inactive_name:
                inactive = int(value)
                break
    except (OSError, ValueError):
        pass

    return max(int(limit) - max(usage - inactive, 0), 0)
