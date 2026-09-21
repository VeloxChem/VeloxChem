"""The RI-J suite on the laptop: four centres against the two fitted ways.

    python rij_laptop.py [--molecules caffeine,tagrisso,nitroxide] [--bases ...]

A pure functional has no exact exchange, so the Coulomb matrix is the whole of the
two-electron build and the fitting never has to be closed for an orbital. What the
suite measures is not only the ratio but **where the time goes afterwards**: once
the Coulomb build is a hundred times faster it stops being the cost, and the
quadrature is what remains.

Nitroxide is a doublet radical and is run unrestricted, which is two Fock matrices
an iteration. Its rows are not comparable with the closed shell ones.

**The fitting set is `def2-universal-jfit`**, a third the size of the jkfit the
RI-JK tables use. A Coulomb only fitting does not have to describe the products of
orbitals an exchange needs.
"""
import argparse
import os
from datetime import date
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

import rijbench
from scfbench import provenance, write

CASES = [
    ("caffeine", ["def2-svp", "def2-svpd", "def2-tzvp", "def2-tzvpd"], 0, 1),
    ("tagrisso", ["def2-svp", "def2-svpd"], 0, 1),
    ("nitroxide", ["def2-svp", "def2-svpd", "def2-tzvp", "def2-tzvpd"], 0, 2),
]

WAYS = ["full", "ri_j_conventional", "ri_j_simd", "ri_j_simd_direct"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--machine", default="m4max")
    parser.add_argument("--molecules", default=None,
                        help="a comma separated subset of the molecules")
    parser.add_argument("--bases", default=None,
                        help="a comma separated subset of the bases")
    parser.add_argument("--ways", default=None)
    parser.add_argument("--functional", default="BLYP")
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    cases = CASES
    if args.molecules:
        wanted = [m.strip().lower() for m in args.molecules.split(",")]
        cases = [c for c in CASES if c[0] in wanted]

    ways = args.ways.split(",") if args.ways else WAYS

    threads = int(os.environ.get("OMP_NUM_THREADS", os.cpu_count()))

    rows = []
    out = Path(args.out) if args.out else (
        Path(__file__).resolve().parent.parent / "data" / "scf" /
        f"{date.today():%Y-%m-%d}_{args.machine}_rij.json")

    for molecule, bases, charge, multiplicity in cases:
        if args.bases:
            subset = [b.strip().lower() for b in args.bases.split(",")]
            bases = [b for b in bases if b in subset]

        for basis_name in bases:
            for way in ways:
                row = rijbench.run(molecule, basis_name, way,
                                   functional=args.functional,
                                   charge=charge, multiplicity=multiplicity)
                rows.append(row)
                print(f"  {molecule:10s} {basis_name:11s} {way:18s} "
                      f"nao {row['nao']:5d} it {row['iterations']:3d} "
                      f"wall {row['wall']:9.2f}  J {row['coulomb']:9.2f}  "
                      f"XC {row['xc']:8.2f}  rest {row['rest']:7.2f}",
                      flush=True)
                out.parent.mkdir(parents=True, exist_ok=True)
                write(out, "rij", provenance(args.machine, 1, threads), rows)

    print(f"\nwrote {out}", flush=True)

    print(f"\n| molecule | basis | nao | naux | method | wall | J | XC | "
          f"speedup | J speedup |")
    print("| --- | --- | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: |")
    for molecule, bases, _, _ in cases:
        for basis_name in bases:
            picked = [r for r in rows
                      if r["molecule"] == molecule and r["basis"] == basis_name]
            if not picked:
                continue
            ref = next((r for r in picked if r["method"] == "full"), None)
            for r in picked:
                head = (f"| {molecule} | {basis_name} | {r['nao']} | "
                        f"{r['naux'] or '--'} " if r["method"] == "full"
                        else "| | | | ")
                speed = (f"{ref['wall'] / r['wall']:.2f}" if ref else "")
                jspeed = (f"{ref['coulomb'] / r['coulomb']:.0f}"
                          if ref and r["coulomb"] > 0 else "")
                print(f"{head}| {r['method']} | {r['wall']:.2f} | "
                      f"{r['coulomb']:.2f} | {r['xc']:.2f} | {speed} | {jspeed} |")


if __name__ == "__main__":
    sys.exit(main())
