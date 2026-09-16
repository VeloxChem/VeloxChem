"""Caffeine molecular gradient on the laptop: four-centre against RI-JK simd.

    python grad_laptop.py

HF and B3LYP, def2-svp and def2-svpd, both ways in one process per case. Writes
one record per calculation to ../data/gradient/ and prints the table.
"""
import sys
from datetime import date
from pathlib import Path

import os

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from gradbench import run
from scfbench import provenance, write

MOLECULE = "caffeine"
AUX = "def2-universal-jkfit"
BASES = ["def2-svp", "def2-svpd", "def2-tzvp", "def2-tzvpd"]
FUNCTIONALS = ["HF", "B3LYP"]
METHODS = ["full", "ri_jk_simd"]

# NOTE: what the run actually used, so the table does not have to say
# "None threads" where the thread count belongs.
THREADS = int(os.environ.get("OMP_NUM_THREADS", os.cpu_count()))

rows = []
for functional in FUNCTIONALS:
    for basis in BASES:
        for method in METHODS:
            row = run(MOLECULE, basis, AUX, method, functional)
            rows.append(row)
            print(f'  {functional:6s} {basis:10s} {method:12s}'
                  f'  scf {row["scf_wall"]:8.2f}  grad {row["grad_wall"]:8.2f}'
                  f'  {row["grad_walls"]}', flush=True)
            out = (Path(__file__).resolve().parent.parent / "data" / "gradient" /
                   f'{date.today():%Y-%m-%d}_m4max_{MOLECULE}.json')
            write(out, "gradient", provenance("m4max", 1, THREADS), rows)

print(f'\nwrote {out}', flush=True)

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
