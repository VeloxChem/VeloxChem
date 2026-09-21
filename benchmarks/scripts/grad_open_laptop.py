"""An open shell molecular gradient on the laptop: four-centre against RI-JK simd.

    python grad_open_laptop.py [molecule] [basis,basis] [charge] [multiplicity]

A cation of the molecule by default, doublet, run unrestricted. HF and B3LYP, both
ways in one process per case.

The conventional resolution of the identity is not a column: it has no open shell
gradient at all, so an unrestricted run asking for it is refused rather than
measured. The way which forms the integrals again on every call is not one either,
for the same reason it is not a column of the open shell SCF suite.

A row of this suite is not comparable with a row of the closed shell one: an open
shell gradient contracts an exchange for each spin where a closed shell contracts
one and doubles it, and the underlying SCF took its own number of iterations.
"""
import sys
from datetime import date
from pathlib import Path

import os

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from gradbench import run
from scfbench import provenance, write

# NOTE: the caffeine grid by default, which is the table this suite has always
# measured. A molecule and a comma separated list of bases may be named instead,
# so a second molecule does not need a second copy of this file.
MOLECULE = sys.argv[1] if len(sys.argv) > 1 else "caffeine"
CHARGE = int(sys.argv[3]) if len(sys.argv) > 3 else 1
MULTIPLICITY = int(sys.argv[4]) if len(sys.argv) > 4 else 2
AUX = "def2-universal-jkfit"
BASES = (sys.argv[2].split(",") if len(sys.argv) > 2
         else ["def2-svp", "def2-svpd", "def2-tzvp", "def2-tzvpd"])
FUNCTIONALS = ["HF", "B3LYP"]
METHODS = ["full", "ri_jk_simd"]

# NOTE: what the run actually used, so the table does not have to say
# "None threads" where the thread count belongs.
THREADS = int(os.environ.get("OMP_NUM_THREADS", os.cpu_count()))

ION = "cation" if CHARGE > 0 else ("anion" if CHARGE < 0 else "open")

# NOTE: formed once, so a run which crosses midnight does not split its output.
OUT = (Path(__file__).resolve().parent.parent / "data" / "gradient" /
       f'{date.today():%Y-%m-%d}_m4max_{MOLECULE}_{ION}.json')

rows = []
for functional in FUNCTIONALS:
    for basis in BASES:
        for method in METHODS:
            row = run(MOLECULE, basis, AUX, method, functional,
                      charge=CHARGE, multiplicity=MULTIPLICITY)
            rows.append(row)
            print(f'  {functional:6s} {basis:10s} {method:12s}'
                  f'  scf {row["scf_wall"]:8.2f}  grad {row["grad_wall"]:8.2f}'
                  f'  {row["grad_walls"]}', flush=True)
            write(OUT, "gradient", provenance("m4max", 1, THREADS), rows)

print(f'\nwrote {OUT}', flush=True)

# the table
print(f'\n| functional | basis | nao | method | gradient | speedup | '
      f'SCF | against four-centre |')
print('| --- | --- | ---: | --- | ---: | ---: | ---: | ---: |')
for functional in FUNCTIONALS:
    for basis in BASES:
        picked = [r for r in rows
                  if r["functional"] == functional and r["basis"] == basis]
        ref = next(r for r in picked if r["method"] == "full")
        for r in picked:
            name = ('four-centre' if r["method"] == "full"
                    else f'RI-JK simd, {r["ri_mode"].replace("_", " ")}')
            speed = ref["grad_wall"] / r["grad_wall"]
            if r["method"] == "full":
                agree = ""
            else:
                d = np.array(r["_gradient"]) - np.array(ref["_gradient"])
                agree = f'{np.abs(d).max():.1e}'
            head = (f'| {functional} | {basis} | {r["nao"]} '
                    if r["method"] == "full" else '| | | ')
            print(f'{head}| {name} | {r["grad_wall"]:.2f} | {speed:.2f} | '
                  f'{r["scf_wall"]:.2f} | {agree} |')
