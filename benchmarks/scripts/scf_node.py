"""One molecule and one basis, four ways of building the Fock matrix, on a node.

    srun --mpi=pmix -N 2 --ntasks-per-node 8 --cpus-per-task 32 \\
         --export=ALL,OMP_NUM_THREADS=32,OMP_PROC_BIND=spread,OMP_PLACES=cores,\\
LD_PRELOAD=/nobackup/proj/disk/panor/personal/rinkevic/Libraries/openblas/lib/libopenblas.so \\
         python scf_node.py --machine epyc9755 --molecule tagrisso \\
                            --basis def2-tzvp

The laptop suite sweeps a grid of eight bases; this one measures what it is told to,
because on a node a grid is a day. The ranks and the threads are whatever the
launcher gave and are recorded rather than chosen here.

Three things about the launch line, each of which fails looking like something else:

  * **`--mpi=pmix` is not optional.** Without it srun offers PMI1/2, Open MPI wants
    PMIx, and rather than failing it starts one singleton per task: every rank
    reports a world of one and the job quietly measures the wrong thing. This
    script refuses to start when the launcher asked for more tasks than the
    communicator has.
  * **`mpirun` cannot reach a second node** on this machine. It dies in the pml
    framework with `fi_domain ... Function not implemented`. Use srun.
  * **Above 48 threads a rank, preload the rebuilt OpenBLAS.** The module's is
    compiled for 48; asked for more it warns once per thread and dies in its
    allocator, usually inside something unrelated. This script refuses to start
    when a BLAS in the process serves fewer threads than the rank was given.

One rank per NUMA domain is the design point: 8 x 32 on one node, 16 x 32 on two.

VeloxChem's own output is kept, one file per calculation, under --outdir. That is
where the iteration table and the timing breakdown live. It is written on the
machine the job ran on and is not tracked: a run of four methods is four outputs
and they are large.
"""
import argparse
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

# NOTE: the mask this rank started with, read before anything widens it. It is not
# the rank's core budget -- a launcher asked to bind a rank to a NUMA domain pins
# the process to a single core of it and lets OpenMP spread over the rest from its
# own list of places -- and is kept only to be reported.
_STARTING_MASK = (sorted(os.sched_getaffinity(0))
                  if hasattr(os, "sched_getaffinity") else [])

from mpi4py import MPI

import node
import scfbench
from veloxchem.veloxchemlib import mpi_master

METHODS = ["full", "ri_jk_conventional", "ri_jk_simd", "ri_jk_simd_direct"]


def _say(comm, text):
    """Printed once, by master, so a job of sixteen ranks says things once."""
    if comm.Get_rank() == mpi_master():
        print(text, flush=True)


def _check_launcher(comm, force):
    """Whether the world is the size the launcher was asked for.

    srun without --mpi=pmix starts one singleton per task rather than failing, and
    a job which measures a sixteenth of what it was asked to measure looks exactly
    like a job which is simply slow.
    """
    asked = os.environ.get("SLURM_NTASKS")

    if asked is None or comm.Get_size() == int(asked):
        return []

    return [f'the launcher asked for {asked} tasks and this communicator has '
            f'{comm.Get_size()}; srun needs --mpi=pmix, without which it starts '
            f'one singleton per task' + ('' if not force else ' (forced)')]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--machine", default="node")
    parser.add_argument("--molecule", default="caffeine",
                        help="a name in benchmarks/geometries, or a path to an xyz")
    parser.add_argument("--basis", default="def2-svp")
    parser.add_argument("--aux", default=None,
                        help="the fitting set; derived from the basis by default")
    parser.add_argument("--functional", default="B3LYP",
                        help="B3LYP by default; HF for Hartree-Fock")
    parser.add_argument("--methods", default=None,
                        help="a comma separated subset of "
                             + ",".join(METHODS))
    parser.add_argument("--out", default=None, help="the json record")
    parser.add_argument("--outdir", default=None,
                        help="where VeloxChem's own output is kept; beside the "
                             "json by default")
    parser.add_argument("--conv-thresh", type=float, default=1.0e-8)
    parser.add_argument("--max-iter", type=int, default=50)
    parser.add_argument("--force", action="store_true",
                        help="start even when the launcher or a BLAS looks wrong")
    args = parser.parse_args()

    comm = MPI.COMM_WORLD
    threads = int(os.environ.get("OMP_NUM_THREADS", os.cpu_count()))

    methods = ([m.strip() for m in args.methods.split(",")] if args.methods
               else METHODS)
    unknown = set(methods) - set(scfbench.METHODS)
    if unknown:
        raise SystemExit(f"no such method: {', '.join(sorted(unknown))}")

    aux = args.aux or scfbench.fitting_set(args.basis)
    if aux is None:
        raise SystemExit(f"no fitting set is known for {args.basis}; give --aux")

    geometry = scfbench.geometry(args.molecule)
    if not geometry.is_file():
        raise SystemExit(f"no such geometry: {geometry}")

    # NOTE: the pool is sized before anything touches a matrix, and after veloxchem
    # is imported, which is what pins the process in the first place.
    runtime = node.runtime(comm, threads)
    runtime["starting_mask"] = _STARTING_MASK

    fatal, warnings = node.check_blas_ceilings(threads)
    fatal = _check_launcher(comm, args.force) + fatal

    for note in warnings + fatal:
        _say(comm, f"# {note}")

    if fatal and not args.force:
        _say(comm, "# refusing to start; --force to measure anyway")
        return 1

    stem = (f"{time.strftime('%Y-%m-%d')}_{args.machine}_"
            f"{geometry.stem}_{args.basis}")
    out = Path(args.out) if args.out else (
        Path(__file__).resolve().parent.parent / "data" / "scf" / f"{stem}.json")
    outdir = Path(args.outdir) if args.outdir else out.parent / f"{stem}_out"

    if comm.Get_rank() == mpi_master():
        outdir.mkdir(parents=True, exist_ok=True)

    comm.barrier()

    info = scfbench.provenance(args.machine, ranks=comm.Get_size(),
                               threads=threads, runtime=runtime)

    topo = runtime["topology"]
    _say(comm, f"# {topo['ranks']} "
               f"{'rank' if topo['ranks'] == 1 else 'ranks'} over {topo['nodes']} "
               f"{'node' if topo['nodes'] == 1 else 'nodes'}, "
               f"{threads} threads a rank, numpy pool {runtime['numpy_pool']}")
    _say(comm, f"# {geometry.stem} / {args.basis} / {aux} / {args.functional}")
    _say(comm, f"# output kept in {outdir}")

    rows = []
    for method in methods:
        ostream_path = outdir / f"{args.functional}_{args.basis}_{method}.out"

        t0 = time.time()
        row = scfbench.run(args.molecule, args.basis, aux, method,
                           args.functional, conv_thresh=args.conv_thresh,
                           max_iter=args.max_iter, ostream_path=ostream_path,
                           comm=comm)
        rows.append(row)

        if row is not None:
            print(f"{args.functional:6s} {args.basis:12s} {method:20s} "
                  f"wall {row['wall']:9.2f}  ri {row['ri_setup']:8.2f}  "
                  f"2e {row['fock_2e_total']:9.2f}  "
                  f"xc {row['fock_xc_total']:8.2f}  E {row['energy']:.8f}  "
                  f"({time.time() - t0:.0f} s)", flush=True)

        scfbench.write(out, "scf", info, rows)

    _say(comm, f"\nwritten to {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
