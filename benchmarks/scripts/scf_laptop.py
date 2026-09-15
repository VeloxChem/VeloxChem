"""The laptop SCF suite: caffeine, three ways of building the Fock matrix.

    python scf_laptop.py [--machine m4max] [--out <path>] [--quick]

Every basis is paired with the fitting set it is meant to be used with. The def2
sets take the universal jkfit; the correlation consistent ones take their own
RIFIT. Both functionals are run over the same grid so that the quadrature's share
can be read off against the Hartree-Fock rows.
"""
import argparse
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import scfbench

PAIRS = [
    ("def2-svp", "def2-universal-jkfit"),
    ("def2-svpd", "def2-universal-jkfit"),
    ("def2-tzvp", "def2-universal-jkfit"),
    ("def2-tzvpd", "def2-universal-jkfit"),
    ("cc-pvdz", "cc-pvdz-rifit"),
    ("aug-cc-pvdz", "aug-cc-pvdz-rifit"),
    ("cc-pvtz", "cc-pvtz-rifit"),
    ("aug-cc-pvtz", "aug-cc-pvtz-rifit"),
]

# NOTE: the simd driver has two ways of doing the same thing and they are
# different calculations, not one: the held one forms the B vectors and keeps
# them, the direct one sweeps the integrals every build. They get a row each.
METHODS = ["full", "ri_jk_conventional", "ri_jk_simd", "ri_jk_simd_direct"]
FUNCTIONALS = ["HF", "B3LYP"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--machine", default="m4max")
    parser.add_argument("--molecule", default="caffeine")
    parser.add_argument("--out", default=None)
    parser.add_argument("--quick", action="store_true",
                        help="the two smallest pairs only, to check the wiring")
    args = parser.parse_args()

    pairs = PAIRS[:1] + PAIRS[4:5] if args.quick else PAIRS

    threads = int(os.environ.get("OMP_NUM_THREADS", os.cpu_count()))
    info = scfbench.provenance(args.machine, ranks=1, threads=threads)

    out = Path(args.out) if args.out else (
        Path(__file__).resolve().parent.parent / "data" / "scf" /
        f"{time.strftime('%Y-%m-%d')}_{args.machine}_{args.molecule}.json")

    rows = []
    for functional in FUNCTIONALS:
        for basis, aux in pairs:
            for method in METHODS:
                t0 = time.time()
                row = scfbench.run(args.molecule, basis, aux, method, functional)
                rows.append(row)
                print(f"{functional:6s} {basis:12s} {method:20s} "
                      f"wall {row['wall']:9.2f}  2e {row['fock_2e_total']:9.2f}  "
                      f"xc {row['fock_xc_total']:8.2f}  E {row['energy']:.8f}  "
                      f"({time.time()-t0:.0f} s)", flush=True)
                scfbench.write(out, "scf", info, rows)

    print(f"\nwritten to {out}")


if __name__ == "__main__":
    main()
