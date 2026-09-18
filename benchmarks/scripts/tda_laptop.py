"""Caffeine TDA on the laptop: four-centre against RI-JK simd.

    python tda_laptop.py [molecule] [basis,basis] [tda|rpa] [functional,functional]

Five states, both ways in one process per case. The functionals default to HF and
B3LYP; naming range separated ones instead measures the path which splits the
exchange between the plain operator and the attenuated one, and the record goes to
a file of its own so the two are not mixed in one table.

Writes one record per calculation to ../data/tda/ and prints the table.
"""
import os
import sys
from datetime import date
from pathlib import Path

import numpy as np
from veloxchem.veloxchemlib import parse_xc_func

sys.path.insert(0, str(Path(__file__).resolve().parent))

from tdabench import run
from scfbench import provenance, write

MOLECULE = sys.argv[1] if len(sys.argv) > 1 else "caffeine"
AUX = "def2-universal-jkfit"
BASES = sys.argv[2].split(",") if len(sys.argv) > 2 else ["def2-svp"]
SOLVER = sys.argv[3] if len(sys.argv) > 3 else "tda"
FUNCTIONALS = sys.argv[4].split(",") if len(sys.argv) > 4 else ["HF", "B3LYP"]
METHODS = ["full", "ri_jk_simd"]

# NOTE: a range separated grid is kept apart from the plain one. They measure
# different Fock matrices, not the same one built differently, and a reader of one
# table should not have to tell its rows apart by the name of the functional.
# Whether a functional is one is asked of veloxchem rather than read off its name,
# which would have to be kept in step with the library by hand.
RANGE_SEPARATED = any(
    f.upper() != "HF" and parse_xc_func(f.upper()).is_range_separated()
    for f in FUNCTIONALS)

TAG = "_rs" if RANGE_SEPARATED else ""

THREADS = int(os.environ.get("OMP_NUM_THREADS", os.cpu_count()))

rows = []
for functional in FUNCTIONALS:
    for basis in BASES:
        for method in METHODS:
            row = run(MOLECULE, basis, AUX, method, functional,
                      solver=SOLVER)
            rows.append(row)
            print(f'  {functional:6s} {basis:10s} {method:12s}'
                  f'  scf {row["scf_wall"]:8.2f}  {SOLVER} {row["tda_wall"]:9.2f}'
                  f'  {row["iterations"]:3d} iter'
                  f'  converged {row["converged"]}', flush=True)
            out = (Path(__file__).resolve().parent.parent / "data" / "tda" /
                   f'{date.today():%Y-%m-%d}_m4max_{MOLECULE}{TAG}'
                   f'{"" if SOLVER == "tda" else "_" + SOLVER}.json')
            write(out, "tda", provenance("m4max", 1, THREADS), rows)

print(f'\nwrote {out}', flush=True)

print('\n| functional | basis | nao | method | TDA (s) | speedup | iter | SCF (s)'
      ' | max dE (a.u.) |')
print('| --- | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: |')
for functional in FUNCTIONALS:
    for basis in BASES:
        picked = [r for r in rows
                  if r["functional"] == functional and r["basis"] == basis]
        if not picked:
            continue
        ref = next(r for r in picked if r["method"] == "full")
        for r in picked:
            name = ('four-centre' if r["method"] == "full"
                    else f'RI-JK simd, {r["ri_mode"].replace("_", " ")}')
            speed = ref["tda_wall"] / r["tda_wall"]
            if r["method"] == "full":
                agree = ""
            else:
                d = (np.array(r["excitation_energies"]) -
                     np.array(ref["excitation_energies"]))
                agree = f'{np.abs(d).max():.1e}'
            head = (f'| {functional} | {basis} | {r["nao"]} '
                    if r["method"] == "full" else '| | | ')
            print(f'{head}| {name} | {r["tda_wall"]:.2f} | {speed:.2f} | '
                  f'{r["iterations"]} | {r["scf_wall"]:.2f} | {agree} |')
