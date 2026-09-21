"""Caffeine reduced two-photon absorption: four-centre against RI-JK simd.

    python redtpa_laptop.py [molecule] [basis,basis]

Five frequencies, HF and B3LYP, both ways in one process per case. Writes one
record per calculation to ../data/redtpa/ and prints the table.
"""
import os
import sys
from datetime import date
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from redtpabench import run
from scfbench import provenance, write

MOLECULE = sys.argv[1] if len(sys.argv) > 1 else "caffeine"
AUX = "def2-universal-jkfit"
BASES = sys.argv[2].split(",") if len(sys.argv) > 2 else ["def2-svp", "def2-svpd"]
FREQUENCIES = [0.050, 0.075, 0.100, 0.125, 0.150]
FUNCTIONALS = ["HF", "B3LYP"]
METHODS = ["full", "ri_jk_simd"]

THREADS = int(os.environ.get("OMP_NUM_THREADS", os.cpu_count()))

rows = []
for functional in FUNCTIONALS:
    for basis in BASES:
        for method in METHODS:
            row = run(MOLECULE, basis, AUX, method, functional, FREQUENCIES)
            rows.append(row)
            print(f'  {functional:6s} {basis:10s} {method:12s}'
                  f'  scf {row["scf_wall"]:8.2f}  tpa {row["tpa_wall"]:9.2f}',
                  flush=True)
            out = (Path(__file__).resolve().parent.parent / "data" / "redtpa" /
                   f'{date.today():%Y-%m-%d}_m4max_{MOLECULE}.json')
            write(out, "redtpa", provenance("m4max", 1, THREADS), rows)

print(f'\nwrote {out}', flush=True)

print('\n| functional | basis | nao | method | TPA (s) | speedup | SCF (s)'
      ' | gamma at 0.100 | max rel |')
print('| --- | --- | ---: | --- | ---: | ---: | ---: | --- | ---: |')
for functional in FUNCTIONALS:
    for basis in BASES:
        picked = [r for r in rows
                  if r["functional"] == functional and r["basis"] == basis]
        if not picked:
            continue
        ref = next(r for r in picked if r["method"] == "full")
        mid = FREQUENCIES.index(0.100)
        for r in picked:
            name = ('four-centre' if r["method"] == "full"
                    else f'RI-JK simd, {r["ri_mode"].replace("_", " ")}')
            speed = ref["tpa_wall"] / r["tpa_wall"]
            g = f'{r["gamma_real"][mid]:.3f}{r["gamma_imag"][mid]:+.3f}i'
            if r["method"] == "full":
                agree = ""
            else:
                a = np.array(ref["gamma_real"]) + 1j * np.array(ref["gamma_imag"])
                b = np.array(r["gamma_real"]) + 1j * np.array(r["gamma_imag"])
                agree = f'{np.abs((b - a) / a).max():.1e}'
            head = (f'| {functional} | {basis} | {r["nao"]} '
                    if r["method"] == "full" else '| | | ')
            print(f'{head}| {name} | {r["tpa_wall"]:.2f} | {speed:.2f} | '
                  f'{r["scf_wall"]:.2f} | {g} | {agree} |')
