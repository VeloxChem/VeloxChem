"""The range separated SCF suite on the laptop: four centres against the simd way.

    python scf_rs_laptop.py [--molecule caffeine] [--charge 0] [--multiplicity 1]

A hybrid range separated functional splits its exchange between the plain operator
and the attenuated one, and both ways of building pay for that split twice over. The
four-centre way makes a second full pass of its own kernels every iteration, `kx_rs`
on top of `2jkx`. The simd resolution of the identity holds a second set of B vectors
and adds a second exchange inside the same sweep of the auxiliary basis.

**The conventional RI-JK column is not here.** That driver has no attenuated B
vectors and refuses a range separated functional, so it would be a column of
refusals rather than a column of numbers.

**A row of this suite is not comparable with a row of the plain suite.** What it
measures is a different Fock matrix, not the same one built differently, and the
interesting number is the ratio within a row rather than against `scf_laptop.py`.
The plain tables of the same molecule and bases are what say what the split costs,
and they are read against this one functional by functional, not row by row.
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
]

METHODS = ["full", "ri_jk_simd"]

FUNCTIONALS = ["CAM-B3LYP", "WB97X-D4"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--machine", default="m4max")
    parser.add_argument("--molecule", default="caffeine")
    parser.add_argument("--charge", type=int, default=0)
    parser.add_argument("--multiplicity", type=int, default=1)
    parser.add_argument("--out", default=None)
    parser.add_argument("--max-iter", type=int, default=100)
    parser.add_argument("--bases", default=None,
                        help="a comma separated subset of the orbital bases")
    parser.add_argument("--functionals", default=None,
                        help="a comma separated subset of the functionals")
    args = parser.parse_args()

    pairs = PAIRS
    if args.bases:
        wanted = [b.strip().lower() for b in args.bases.split(",")]
        pairs = [p for p in PAIRS if p[0] in wanted]
        missing = set(wanted) - {p[0] for p in pairs}
        if missing:
            raise SystemExit(f"no such pair: {', '.join(sorted(missing))}")

    functionals = FUNCTIONALS
    if args.functionals:
        functionals = [f.strip() for f in args.functionals.split(",")]

    threads = int(os.environ.get("OMP_NUM_THREADS", os.cpu_count()))
    info = scfbench.provenance(args.machine, ranks=1, threads=threads)

    # NOTE: the path is formed once, so a run which crosses midnight does not write
    # its last rows to a second file and leave the first a stale prefix of itself.
    # It names the spin state, as the same molecule is measured in both.
    spin = "closed" if args.multiplicity == 1 else f"m{args.multiplicity}"

    out = Path(args.out) if args.out else (
        Path(__file__).resolve().parent.parent / "data" / "scf" /
        f"{time.strftime('%Y-%m-%d')}_{args.machine}_{args.molecule}_rs_{spin}.json")

    if out.exists() and args.out is None:
        raise SystemExit(f"{out} exists; pass --out to write beside it")

    rows = []
    for functional in functionals:
        for basis, aux in pairs:
            for method in METHODS:
                t0 = time.time()
                row = scfbench.run(args.molecule, basis, aux, method, functional,
                                   charge=args.charge,
                                   multiplicity=args.multiplicity,
                                   max_iter=args.max_iter)
                rows.append(row)

                energy = ("did not converge" if row["energy"] is None
                          else f"E {row['energy']:.8f}")

                print(f"{functional:10s} {basis:11s} {method:12s} "
                      f"wall {row['wall']:9.2f}  ri {row['ri_setup']:7.2f}  "
                      f"2e {row['fock_2e_total']:9.2f}  "
                      f"xc {row['fock_xc_total']:8.2f}  it {row['iterations']:3d}  "
                      f"{energy}  ({time.time() - t0:.0f} s)", flush=True)

                # NOTE: written after every row, so a grid which is stopped or dies
                # leaves what it had rather than nothing.
                scfbench.write(out, "scf", info, rows)

    print(f"\nwritten to {out}")


if __name__ == "__main__":
    main()
