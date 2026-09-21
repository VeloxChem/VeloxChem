"""The two-center Coulomb table: the plain driver against the range separated one.

    python t2c_laptop.py [--threads 14,1] [--omega 0.3]

The grid the existing table measures -- tagrisso, taxol, crambin and ubiquitin in
the def2 universal jfit and jkfit sets -- ordered so the largest case comes last.
That ordering is deliberate: the range separated driver holds two packed matrices
where the plain one holds a single one, so ubiquitin in jkfit asks for about
twenty four gigabytes of matrix alone and may not fit. Every row is written as it
finishes, and a case the system kills costs that row and not the run.

Each case runs as its own process, which is how the existing table was taken and
which is what lets a killed case be recorded rather than ending everything.
"""
import argparse
import json
import os
import subprocess
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from scfbench import provenance, write

GEOMETRIES = Path.home() / "Downloads"

# NOTE: smallest first, so the case which may not fit is the last thing attempted.
CASES = [
    ("tagrisso", "def2-universal-jfit"),
    ("tagrisso", "def2-universal-jkfit"),
    ("taxol", "def2-universal-jfit"),
    ("taxol", "def2-universal-jkfit"),
    ("crambin", "def2-universal-jfit"),
    ("crambin", "def2-universal-jkfit"),
    ("ubiquitin", "def2-universal-jfit"),
    ("ubiquitin", "def2-universal-jkfit"),
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--machine", default="m4max")
    parser.add_argument("--omega", type=float, default=0.3)
    parser.add_argument("--threads", default="14,1")
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    threads = [int(t) for t in args.threads.split(",")]

    out = Path(args.out) if args.out else (
        Path(__file__).resolve().parent.parent / "data" / "t2c" /
        f"{date.today():%Y-%m-%d}_{args.machine}_two_center.json")

    if out.exists():
        raise SystemExit(f"{out} exists: move it aside or name another with --out")

    rows = []
    for molecule, fitting in CASES:
        geometry = GEOMETRIES / f"{molecule}.xyz"
        if not geometry.is_file():
            print(f"  {molecule:10s} {fitting:22s} no geometry at {geometry}", flush=True)
            continue

        for nthreads in threads:
            env = dict(os.environ, OMP_NUM_THREADS=str(nthreads))

            # NOTE: the values are checked on the small cases alone. The check holds
            # a third and a fourth matrix, which the largest cases have no room for,
            # and what it tests does not depend on the size of the molecule.
            cmd = [sys.executable, str(Path(__file__).parent / "t2cbench.py"),
                   str(geometry), fitting, str(args.omega)]
            if molecule in ("tagrisso", "taxol"):
                cmd.append("--check")

            done = subprocess.run(cmd, env=env, capture_output=True, text=True)

            if done.returncode != 0 or not done.stdout.strip():
                tail = (done.stderr or "").strip().splitlines()[-1:] or ["no output"]
                print(f"  {molecule:10s} {fitting:22s} {nthreads:3d} thr  FAILED: {tail[0][:90]}",
                      flush=True)
                rows.append({"molecule": molecule, "basis": fitting,
                             "threads": nthreads, "failed": tail[0][:200]})
                write(out, "t2c", provenance(args.machine, 1, nthreads), rows)
                continue

            row = json.loads(done.stdout.strip().splitlines()[-1])
            rows.append(row)
            print(f'  {molecule:10s} {fitting:22s} {nthreads:3d} thr  nao {row["nao"]:6d}'
                  f'  plain {row["plain_wall"]:9.4f}  rs {row["rs_wall"]:9.4f}'
                  f'  {row["rs_wall"] / row["plain_wall"]:5.2f}x'
                  f'  packed {row["packed_gb"]:6.2f} GB', flush=True)
            write(out, "t2c", provenance(args.machine, 1, nthreads), rows)

    print(f"\nwrote {out}", flush=True)


if __name__ == "__main__":
    main()
