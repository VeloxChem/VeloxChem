"""The Coulomb only gradient on the laptop: the two fitted ways against each other.

    python rij_grad_laptop.py [--molecules caffeine,nitroxide] [--bases ...]

A pure functional's gradient is the Coulomb term of the resolution of the identity
and the quadrature, and nothing else of the two-electron part: there is no exchange
to differentiate. What this measures is the new simd driver against the
conventional one which fits the same approximation, **and the four-centre way
beside them**, because the ratio against a dense build is what says whether the
fitting is worth having at all.

The fitting set is `def2-universal-jfit`, the one a Coulomb only fitting takes, and
not the jkfit of the RI-JK tables.

Nitroxide is a doublet radical and is run unrestricted. Its rows are not comparable
with the closed shell ones: an open shell fits the total density where a closed
shell fits one spin's, and its quadrature runs over two densities.
"""
import argparse
import os
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from gradbench import run
from scfbench import provenance, write

CASES = [
    ("caffeine", 0, 1),
    ("nitroxide", 0, 2),
]

BASES = ["def2-svp", "def2-svpd", "def2-tzvp", "def2-tzvpd"]

WAYS = ["full", "ri_j_conventional", "ri_j_simd"]

AUX = "def2-universal-jfit"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--machine", default="m4max")
    parser.add_argument("--molecules", default=None)
    parser.add_argument("--bases", default=None)
    parser.add_argument("--ways", default=None)
    parser.add_argument("--functional", default="BLYP")
    parser.add_argument("--repeats", type=int, default=2)
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    cases = CASES
    if args.molecules:
        wanted = [m.strip().lower() for m in args.molecules.split(",")]
        cases = [c for c in CASES if c[0] in wanted]

    bases = args.bases.split(",") if args.bases else BASES
    ways = args.ways.split(",") if args.ways else WAYS

    threads = int(os.environ.get("OMP_NUM_THREADS", os.cpu_count()))
    info = provenance(args.machine, ranks=1, threads=threads)

    out = Path(args.out) if args.out else (
        Path(__file__).resolve().parent.parent / "data" / "gradient" /
        f"{date.today():%Y-%m-%d}_{args.machine}_rij_grad.json")

    rows = []
    for molecule, charge, multiplicity in cases:
        for basis_name in bases:
            for way in ways:
                row = run(molecule, basis_name, AUX, way, args.functional,
                          repeats=args.repeats, charge=charge,
                          multiplicity=multiplicity)
                rows.append(row)
                print(f"  {molecule:10s} {basis_name:11s} {way:18s} "
                      f"nao {row['nao']:5d} scf {row['scf_wall']:9.2f}  "
                      f"grad {row['grad_wall']:8.2f}",
                      flush=True)
                out.parent.mkdir(parents=True, exist_ok=True)
                write(out, "gradient", info, rows)

    print(f"\nwritten to {out}")


if __name__ == "__main__":
    main()
