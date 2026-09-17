"""What a rank has to do to the process before it measures anything.

Three things on the two-socket EPYC node, each of which fails in a way that looks
like something else:

  * numpy's pool is built from the affinity mask the first time it is used. The
    launcher pins the rank and importing veloxchem loads libgomp which pins it
    further, so left alone that pool is **one thread**: 105 Gflop/s where the
    machine gives 3500. The mask is widened first, then the pool is sized.
  * a BLAS built for fewer threads than the rank was given does not refuse them. It
    grows an auxiliary table, warns once per thread, and dies later in an unrelated
    allocation. The module's OpenBLAS is `MAX_THREADS=48`. That is checked here,
    before the calculation, rather than diagnosed afterwards from a core file.
  * `OPENBLAS_NUM_THREADS` must not be used to cap numpy. Every OpenBLAS in the
    process reads it, the driver's included, and setting it above a library's
    compiled ceiling is what causes the crash above. threadpoolctl selects by name
    and is the only safe way to cap one of them.

The mask a rank starts with is reported and never used as a core budget: a launcher
asked to bind a rank to a NUMA domain pins the process to a single core of it and
lets OpenMP spread over the rest from its own list of places, so the mask reads one
core where the rank has thirty two.
"""
import os
import platform
import socket
import subprocess
from pathlib import Path


def starting_mask():
    """The cores this process was pinned to, for the record and not for sizing."""
    if not hasattr(os, "sched_getaffinity"):
        return []
    return sorted(os.sched_getaffinity(0))


def _blas_libraries():
    """Every BLAS in the process, with the pool each of them actually has."""
    try:
        from threadpoolctl import ThreadpoolController
    except ImportError:
        return None

    return [{"lib": Path(lib.get("filepath", "?")).name,
             "path": str(lib.get("filepath", "?")),
             "api": str(lib.get("internal_api", "?")),
             "threads": lib.get("num_threads", 0)}
            for lib in ThreadpoolController().info()]


def widen_for_numpy(threads):
    """Unpins the process and sizes numpy's pool to this rank's threads.

    The driver's own library is left alone: it threads through OpenMP, which binds
    per region from OMP_PLACES rather than from the mask, and wants every thread it
    was given.

    :param threads:
        The cores this rank was given, which is OMP_NUM_THREADS and not the mask.

    :return:
        The pool numpy ended up with, or None when it could not be established.
        None rather than a guess: without threadpoolctl the pool is whatever the
        launcher's environment made it, and reporting the number that was asked
        for would be reporting an assumption as a measurement.
    """
    if hasattr(os, "sched_setaffinity"):
        try:
            os.sched_setaffinity(0, range(os.cpu_count()))
        except OSError:
            pass                                   # a cgroup may forbid it

    try:
        from threadpoolctl import ThreadpoolController
    except ImportError:
        return None

    controller = ThreadpoolController()
    controller.select(prefix="libscipy_openblas").limit(limits=threads)

    # NOTE: read back rather than assumed. The library numpy carries is built with
    # a ceiling of its own -- the scipy build is MAX_THREADS=64 -- so a limit above
    # it is silently clamped and a rank given more does not get them.
    pools = [lib.get("num_threads", 0) for lib in controller.info()
             if "openblas" in str(lib.get("internal_api", ""))]

    return max(pools) if pools else None


def check_blas_ceilings(threads):
    """The libraries which cannot serve the threads this rank was given.

    A library known to be below the count is fatal: it does not refuse the threads,
    it grows an auxiliary table, warns once per thread, and dies later somewhere
    else. Not being able to tell is a different thing and is only said out loud --
    refusing there would stop every machine which has no threadpoolctl, including
    the laptop, where there is no such library and nothing to trip over.

    :param threads:
        The cores this rank was given.

    :return:
        The complaints which should stop the run, and the ones which should not.
    """
    libraries = _blas_libraries()

    if libraries is None:
        return [], ["threadpoolctl is not installed, so no BLAS in this process "
                    "is verified and none is recorded; one built for fewer "
                    "threads than this rank was given will warn and may die in "
                    "its allocator. pip install threadpoolctl"]

    return ([f'{lib["lib"]} serves {lib["threads"]} threads and this rank was '
             f'given {threads} ({lib["path"]})'
             for lib in libraries if 0 < lib["threads"] < threads], [])


def topology(comm):
    """How the ranks of this job are spread over the machines of it.

    :param comm:
        The communicator.

    :return:
        The number of ranks, the number of distinct hosts, and the tasks on each.
    """
    hosts = comm.allgather(socket.gethostname())
    counts = sorted({hosts.count(h) for h in set(hosts)})

    return {
        "ranks": comm.Get_size(),
        "nodes": len(set(hosts)),
        # NOTE: a list, because a job whose last node is short is a job whose ranks
        # are not interchangeable, and a single number would hide that.
        "tasks_per_node": counts,
        "hosts": sorted(set(hosts)),
    }


def launcher():
    """What the launcher says it did, which is not always what it did.

    srun without `--mpi=pmix` starts one singleton per task rather than failing, so
    a world of one where the launcher was asked for many is the tell, and both
    halves of it are recorded so the comparison can be made afterwards.
    """
    keys = ("SLURM_JOB_ID", "SLURM_NTASKS", "SLURM_NNODES",
            "SLURM_NTASKS_PER_NODE", "SLURM_CPUS_PER_TASK", "SLURM_MPI_TYPE")
    slurm = {k: os.environ[k] for k in keys if k in os.environ}

    return {
        "slurm": slurm or None,
        "omp": {k: os.environ[k] for k in
                ("OMP_NUM_THREADS", "OMP_PROC_BIND", "OMP_PLACES")
                if k in os.environ} or None,
        # NOTE: the preload is the only way to reach the rebuilt OpenBLAS -- an
        # RPATH is baked into veloxchemlib.so and is read before LD_LIBRARY_PATH --
        # so whether it was set is part of what a number came from.
        "ld_preload": os.environ.get("LD_PRELOAD"),
    }


def _cpu_name():
    """The model rather than the architecture: platform.processor() says "arm"."""
    try:
        if platform.system() == "Darwin":
            return subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"],
                                  capture_output=True, text=True,
                                  check=True).stdout.strip()
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    except Exception:
        pass
    return platform.processor() or platform.machine()


def runtime(comm, threads):
    """Everything about this process which a measurement depends on."""
    return {
        "threads_per_rank": threads,
        "numpy_pool": widen_for_numpy(threads),
        "starting_mask": starting_mask(),
        "topology": topology(comm),
        "launcher": launcher(),
        "blas": _blas_libraries() or [],
    }
