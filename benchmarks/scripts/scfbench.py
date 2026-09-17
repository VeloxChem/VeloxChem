"""One SCF calculation, one record.

The record is what benchmarks/data holds: the calculation, its provenance, and the
timings split into the metric and B vectors the resolution of the identity forms
once, the two-electron build, the quadrature, and everything else.
Speedups are not recorded -- they are computed when a table is rendered, against
the `full` row of the same molecule, basis and functional, so that re-measuring one
method cannot leave a stale ratio behind.
"""
import json
import platform
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

from mpi4py import MPI

import veloxchem as vlx
from veloxchem.veloxchemlib import mpi_master
from veloxchem.outputstream import OutputStream
from veloxchem.scfrestdriver import ScfRestrictedDriver

GEOMETRIES = Path(__file__).resolve().parent.parent / "geometries"


def geometry(name):
    """The xyz of a molecule, named in the repo or given as a path.

    A run on a cluster should not have to copy its molecule into the repo to be
    measured, so a path is taken as well as a name.
    """
    candidate = Path(name)
    if candidate.suffix and candidate.is_file():
        return candidate
    return GEOMETRIES / f"{name}.xyz"


def fitting_set(basis_name):
    """The fitting set a basis is meant to be used with.

    The def2 sets take the universal jkfit; the correlation consistent ones take
    their own RIFIT. Kept here rather than in a runner so that two runners cannot
    come to disagree about which pair was measured.
    """
    name = basis_name.lower()

    if name.startswith("def2-"):
        return "def2-universal-jkfit"
    if name.endswith("-rifit") or name.endswith("-jkfit"):
        return basis_name
    if "cc-pv" in name:
        return f"{basis_name}-rifit"

    return None


# method -> (ri_jk, ri_jk_simd, ri_mode)
METHODS = {
    "full": (False, False, None),
    "ri_jk_conventional": (True, False, None),
    "ri_jk_simd": (True, True, "automatic"),
    "ri_jk_simd_direct": (True, True, "direct"),
}


def _git(*args):
    try:
        return subprocess.run(["git", *args], cwd=Path(__file__).resolve().parent,
                              capture_output=True, text=True,
                              check=True).stdout.strip()
    except Exception:
        return None


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


def _blas():
    try:
        from threadpoolctl import ThreadpoolController

        return [{"lib": Path(lib["filepath"]).name, "threads": lib["num_threads"]}
                for lib in ThreadpoolController().info()]
    except ImportError:
        return []


def provenance(name, ranks, threads, runtime=None):
    """What the numbers came from, which is half of what a benchmark is.

    :param runtime:
        What node.runtime() gathered -- the topology, the launcher, the pool numpy
        ended up with, and every BLAS in the process. Omitted on the laptop, where
        there is one rank and nothing to place.
    """
    return {
        "date": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
        "commit": _git("rev-parse", "--short", "HEAD"),
        "dirty": bool(_git("status", "--porcelain")),
        "veloxchem": vlx.__version__,
        "machine": {
            "name": name,
            "cpu": _cpu_name(),
            "cores": __import__("os").cpu_count(),
            "os": f"{platform.system()} {platform.release()}",
        },
        "config": {"mpi_ranks": ranks, "omp_threads": threads,
                   "blas": (runtime or {}).get("blas") or _blas(),
                   **({"runtime": runtime} if runtime else {})},
    }


def _fock_timings(driver):
    """The two-electron build and the quadrature, summed over the iterations.

    The profiler keys its dictionary by iteration and labels the two parts FockERI
    and FockXC. Iteration zero is the guess and is counted with the rest: it is
    work the calculation did.
    """
    totals = {"FockERI": 0.0, "FockXC": 0.0}
    builds = 0

    timings = getattr(getattr(driver, "_profiler", None), "timing_dict", None) or {}

    for key, labels in timings.items():
        if "FockERI" in labels:
            builds += 1
        for label in totals:
            totals[label] += labels.get(label, 0.0)

    return totals["FockERI"], totals["FockXC"], builds


def run(molecule_name, basis_name, aux_name, method, functional,
        conv_thresh=1.0e-8, max_iter=50, ostream_path=None, comm=None):
    """Runs one calculation and returns its record.

    :param ostream_path:
        Where to keep VeloxChem's own output. The iteration table and the timing
        breakdown are in there and nowhere else, so a run which is to be looked at
        afterwards wants it. None discards it, which is what the laptop suites do.
    :param comm:
        The communicator, needed to open an output which only master writes.
    """
    molecule = vlx.Molecule.read_xyz_file(str(geometry(molecule_name)))
    basis = vlx.MolecularBasis.read(molecule, basis_name.upper(), ostream=None)

    ri_jk, simd, ri_mode = METHODS[method]

    if ostream_path is None:
        ostream = OutputStream(None)
    elif comm is None:
        ostream = OutputStream(str(ostream_path))
    else:
        # NOTE: active on master and silent elsewhere, so the ranks of a job do not
        # race each other writing one file.
        ostream = OutputStream.create_mpi_ostream(comm, str(ostream_path))

    driver = ScfRestrictedDriver(comm=comm, ostream=ostream)
    driver.timing = True
    driver.conv_thresh = conv_thresh
    driver.max_iter = max_iter
    if functional.upper() != "HF":
        driver.xcfun = functional
    if ri_jk:
        driver.ri_jk = True
        driver.ri_jk_simd = simd
        driver.ri_auxiliary_basis = aux_name.upper()
        if ri_mode is not None:
            driver.ri_mode = ri_mode

    t0 = time.time()
    results = driver.compute(molecule, basis)
    wall = time.time() - t0

    # NOTE: flushed and closed before the record is returned, so the file is whole
    # even if the next calculation in the sweep takes the process down with it.
    ostream.flush()
    if ostream_path is not None:
        ostream.close()

    if driver.rank != mpi_master():
        return None

    eri, xc, builds = _fock_timings(driver)

    aux = None
    if ri_jk:
        aux = vlx.MolecularBasis.read(molecule, aux_name.upper(), ostream=None)

    mode = None
    if simd and getattr(driver, "_ri_drv", None) is not None:
        mode = str(driver._ri_drv.get_mode()).replace("rimode.", "")

    return {
        "molecule": molecule_name,
        "atoms": molecule.number_of_atoms(),
        "basis": basis_name,
        "nao": basis.get_dimensions_of_basis(),
        "aux_basis": aux_name if ri_jk else None,
        "naux": aux.get_dimensions_of_basis() if aux is not None else None,
        "occupied": molecule.number_of_alpha_electrons(),
        "method": method,
        "ri_mode": mode,
        "functional": functional,
        "conv_thresh": conv_thresh,
        "energy": results["scf_energy"] if results else None,
        "converged": results is not None,
        "iterations": driver.num_iter,
        "wall": round(wall, 3),
        "ri_setup": round(getattr(driver, "_ri_setup_time", 0.0), 3),
        "fock_2e_total": round(eri, 3),
        "fock_2e_mean": round(eri / builds, 4) if builds else None,
        "fock_xc_total": round(xc, 3),
        "fock_xc_mean": round(xc / builds, 4) if builds else None,
        "builds": builds,
        "remainder": round(wall - eri - xc
                           - getattr(driver, "_ri_setup_time", 0.0), 3),
    }


def write(path, suite, run_info, rows):
    """Writes one run to one file, which is never appended to.

    Only master writes. Every rank reaches here -- the suites loop over the whole
    grid on all of them, because compute() is collective -- and left unguarded the
    ranks of a node race each other over one path, each writing rows which are None
    off master.
    """
    if MPI.COMM_WORLD.Get_rank() != mpi_master():
        return None

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(
        {"schema": 1, "suite": suite, "run": run_info, "rows": rows},
        indent=2) + "\n")
    return path
