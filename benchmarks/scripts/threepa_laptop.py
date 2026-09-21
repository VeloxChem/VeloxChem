"""Caffeine three-photon absorption on the laptop: four-centre against RI-JK simd.

    python threepa_laptop.py [molecule] [basis,basis] [nstates]

HF and B3LYP, both ways in one process per case. Writes one record per calculation
to ../data/3pa/ and prints the table.
"""
import os
import sys
from datetime import date
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from threepabench import run
from scfbench import provenance, write

MOLECULE = sys.argv[1] if len(sys.argv) > 1 else "caffeine"
AUX = "def2-universal-jkfit"
BASES = sys.argv[2].split(",") if len(sys.argv) > 2 else ["def2-svp"]
NSTATES = int(sys.argv[3]) if len(sys.argv) > 3 else 5
FUNCTIONALS = ["HF", "B3LYP"]
METHODS = ["full", "ri_jk_simd"]

# NOTE: everything compared is invariant to the phase of an excited state vector.
# The transition moments are not, and are recorded rather than compared.
COMPARED = ("strengths_circular", "strengths_linear", "oscillator_strengths",
            "photon_energies")

THREADS = int(os.environ.get("OMP_NUM_THREADS", os.cpu_count()))

# NOTE: the path is formed once. Formed inside the loop it follows the calendar, and
# a run which crosses midnight writes its last rows to a second file and leaves the
# first as a stale prefix of itself.
OUT = (Path(__file__).resolve().parent.parent / "data" / "3pa" /
       f'{date.today():%Y-%m-%d}_m4max_{MOLECULE}.json')
OUT.parent.mkdir(parents=True, exist_ok=True)

rows = []
for functional in FUNCTIONALS:
    for basis in BASES:
        for method in METHODS:
            row = run(MOLECULE, basis, AUX, method, functional, nstates=NSTATES)
            rows.append(row)
            print(f'  {functional:6s} {basis:10s} {method:12s}'
                  f'  scf {row["scf_wall"]:8.2f}  3pa {row["rsp_wall"]:9.2f}',
                  flush=True)
            write(OUT, "3pa", provenance("m4max", 1, THREADS), rows)

print(f'\nwrote {OUT}', flush=True)

print('\n| functional | basis | nao | method | 3PA (s) | speedup | SCF (s)'
      ' | circular strengths (a.u.) | max rel |')
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
            speed = ref["rsp_wall"] / r["rsp_wall"]
            xs = r["strengths_circular"] or []
            shown = ', '.join(f'{v:.4e}' for v in xs[:2])
            if r["method"] == "full" or not xs:
                agree = ""
            else:
                # NOTE: the worst of everything compared and not of the strengths
                # alone, so a column which agrees cannot hide one which does not.
                agree = max(
                    np.abs((np.array(r[k]) - np.array(ref[k])) /
                           np.array(ref[k])).max() for k in COMPARED)
                agree = f'{agree:.1e}'
            head = (f'| {functional} | {basis} | {r["nao"]} '
                    if r["method"] == "full" else '| | | ')
            print(f'{head}| {name} | {r["rsp_wall"]:.2f} | {speed:.2f} | '
                  f'{r["scf_wall"]:.2f} | {shown} | {agree} |')
