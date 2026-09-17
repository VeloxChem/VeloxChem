"""The open shell SCF suite on the laptop: three ways of building the Fock matrix.

    python scf_open_laptop.py [--molecule caffeine] [--bases def2-svp,def2-svpd]

A cation of the molecule, doublet, run unrestricted. Both functionals over the same
grid so that the quadrature's share can be read off against the Hartree-Fock rows,
and three methods: the four centre integrals, the conventional resolution of the
identity, and the simd one which holds its B vectors.

The way which forms the integrals again on every call is **not** here. It does not
serve an open shell -- its fitting is accumulated inside the sweep which builds the
exchange, and two spins there is a different piece of work -- so it would be a
column of refusals.

A row of this suite is not comparable with a row of the closed shell one however
alike the two read: an unrestricted iteration builds two Fock matrices and two
exchanges where a restricted one builds a single exchange and doubles the Coulomb.
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

METHODS = ["full", "ri_jk_conventional", "ri_jk_simd"]
FUNCTIONALS = ["HF", "B3LYP"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--machine", default="m4max")
    parser.add_argument("--molecule", default="caffeine")
    parser.add_argument("--charge", type=int, default=1)
    parser.add_argument("--multiplicity", type=int, default=2)
    parser.add_argument("--out", default=None)
    parser.add_argument("--max-iter", type=int, default=100,
                        help="a radical cation is slow: caffeine takes 45 "
                             "iterations by the four centre way and 51 by either "
                             "resolution of the identity, against 19 and 21 for "
                             "the closed shell neutral, so the default of 50 "
                             "stops two of the three ways one iteration short")
    parser.add_argument("--bases", default=None,
                        help="a comma separated subset of the orbital bases; each "
                             "keeps the fitting set it is paired with above")
    args = parser.parse_args()

    pairs = PAIRS

    if args.bases:
        wanted = [b.strip().lower() for b in args.bases.split(",")]
        pairs = [p for p in PAIRS if p[0] in wanted]
        missing = set(wanted) - {p[0] for p in pairs}
        if missing:
            raise SystemExit(f"no such pair: {', '.join(sorted(missing))}")

    threads = int(os.environ.get("OMP_NUM_THREADS", os.cpu_count()))
    info = scfbench.provenance(args.machine, ranks=1, threads=threads)

    # NOTE: the path is formed once, so a run which crosses midnight does not write
    # its last rows to a second file and leave the first a stale prefix of itself.
    ion = "cation" if args.charge > 0 else ("anion" if args.charge < 0 else "open")

    out = Path(args.out) if args.out else (
        Path(__file__).resolve().parent.parent / "data" / "scf" /
        f"{time.strftime('%Y-%m-%d')}_{args.machine}_{args.molecule}_{ion}.json")

    rows = []
    for functional in FUNCTIONALS:
        for basis, aux in pairs:
            for method in METHODS:
                t0 = time.time()
                row = scfbench.run(args.molecule, basis, aux, method, functional,
                                   charge=args.charge,
                                   multiplicity=args.multiplicity,
                                   max_iter=args.max_iter)
                rows.append(row)

                # NOTE: a row which did not converge has no energy, and is printed
                # and kept as what it is rather than stopping the run. The rest of
                # the grid is still worth having, and a gap which is recorded can
                # be looked at afterwards.
                energy = ("did not converge" if row["energy"] is None
                          else f"E {row['energy']:.8f}")

                print(f"{functional:6s} {basis:12s} {method:20s} "
                      f"wall {row['wall']:9.2f}  ri {row['ri_setup']:7.2f}  "
                      f"2e {row['fock_2e_total']:9.2f}  "
                      f"xc {row['fock_xc_total']:8.2f}  it {row['iterations']:3d}  "
                      f"{energy}  ({time.time() - t0:.0f} s)", flush=True)
                scfbench.write(out, "scf", info, rows)

    print(f"\nwritten to {out}")


if __name__ == "__main__":
    main()
