"""Caffeine two-photon absorption on the laptop: four-centre against RI-JK simd.

    python tpa_laptop.py [molecule] [basis,basis] [nstates]

HF and B3LYP, both ways in one process per case. Writes one record per calculation
to ../data/tpa/ and prints the table.
"""
import os
import sys
from datetime import date
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from tpabench import run
from scfbench import provenance, write

MOLECULE = sys.argv[1] if len(sys.argv) > 1 else "caffeine"
AUX = "def2-universal-jkfit"
BASES = sys.argv[2].split(",") if len(sys.argv) > 2 else ["def2-svp"]
NSTATES = int(sys.argv[3]) if len(sys.argv) > 3 else 5
FUNCTIONALS = ["HF", "B3LYP"]
METHODS = ["full", "ri_jk_simd"]

THREADS = int(os.environ.get("OMP_NUM_THREADS", os.cpu_count()))

rows = []
for functional in FUNCTIONALS:
    for basis in BASES:
        for method in METHODS:
            row = run(MOLECULE, basis, AUX, method, functional, nstates=NSTATES)
            rows.append(row)
            print(f'  {functional:6s} {basis:10s} {method:12s}'
                  f'  scf {row["scf_wall"]:8.2f}  tpa {row["tpa_wall"]:9.2f}',
                  flush=True)
            out = (Path(__file__).resolve().parent.parent / "data" / "tpa" /
                   f'{date.today():%Y-%m-%d}_m4max_{MOLECULE}.json')
            write(out, "tpa", provenance("m4max", 1, THREADS), rows)

print(f'\nwrote {out}', flush=True)

print('\n| functional | basis | nao | method | TPA (s) | speedup | SCF (s)'
      ' | cross sections (a.u.) | max rel |')
print('| --- | --- | ---: | --- | ---: | ---: | ---: | --- | ---: |')
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
            speed = ref["tpa_wall"] / r["tpa_wall"]
            xs = r["cross_sections"] or []
            shown = ', '.join(f'{v:.4e}' for v in xs[:2])
            if r["method"] == "full" or not xs:
                agree = ""
            else:
                a = np.array(ref["cross_sections"]); b = np.array(xs)
                agree = f'{np.abs((b - a) / a).max():.1e}'
            head = (f'| {functional} | {basis} | {r["nao"]} '
                    if r["method"] == "full" else '| | | ')
            print(f'{head}| {name} | {r["tpa_wall"]:.2f} | {speed:.2f} | '
                  f'{r["scf_wall"]:.2f} | {shown} | {agree} |')
