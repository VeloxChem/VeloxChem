"""Caffeine geometry optimization on the laptop: four-centre against RI-JK simd.

    python opt_laptop.py

HF and B3LYP, def2-svp and def2-svpd, both ways, run to convergence with no cap
on the iterations. Writes one record per optimization to ../data/optimization/
and prints the table.
"""
import sys
from datetime import date
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from optbench import run
from scfbench import provenance, write

MOLECULE = "caffeine"
AUX = "def2-universal-jkfit"
BASES = ["def2-svp", "def2-svpd"]
FUNCTIONALS = ["HF", "B3LYP"]
METHODS = ["full", "ri_jk_simd"]

rows = []
for functional in FUNCTIONALS:
    for basis in BASES:
        for method in METHODS:
            row = run(MOLECULE, basis, AUX, method, functional)
            rows.append(row)
            print(f'  {functional:6s} {basis:10s} {method:12s}'
                  f'  {row["wall"]:9.2f} s  {row["steps"]:3d} steps'
                  f'  {row["wall_per_step"]:7.2f} s/step'
                  f'  E {row["energy"]:.8f}', flush=True)
            out = (Path(__file__).resolve().parent.parent / "data" /
                   "optimization" / f'{date.today():%Y-%m-%d}_m4max_{MOLECULE}.json')
            write(out, "optimization", provenance("m4max", 1, None), rows)

print(f'\nwrote {out}', flush=True)

print('\n| functional | basis | nao | method | total | speedup | steps | s/step'
      ' | energy | max displacement |')
print('| --- | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |')
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
            speed = ref["wall"] / r["wall"]
            if r["method"] == "full":
                moved = ""
            else:
                # NOTE: the largest distance an atom ended up from where the
                # four-centre optimization put it. The two minima are of slightly
                # different surfaces, so this is the fitting error as a geometry
                # rather than a disagreement to be explained away.
                d = np.array(r["geometry"]) - np.array(ref["geometry"])
                moved = f'{np.linalg.norm(d, axis=1).max():.2e}'
            head = (f'| {functional} | {basis} | {r["nao"]} '
                    if r["method"] == "full" else '| | | ')
            print(f'{head}| {name} | {r["wall"]:.1f} | {speed:.2f} | '
                  f'{r["steps"]} | {r["wall_per_step"]:.2f} | '
                  f'{r["energy"]:.8f} | {moved} |')
