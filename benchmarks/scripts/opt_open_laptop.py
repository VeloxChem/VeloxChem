"""An open shell geometry optimization on the laptop: four-centre against RI-JK simd.

    python opt_open_laptop.py [molecule] [basis,basis] [charge] [multiplicity]

A radical by default, doublet, run unrestricted. HF and B3LYP, both ways, run to
convergence with no cap on the steps. Writes one record per optimization to
../data/optimization/ and prints the table.

The conventional resolution of the identity is not a column and neither is the way
which forms the integrals again on every call: the first has no open shell gradient
at all and the second has no gradient of any kind, so both are refused before the
first step rather than measured.

A row of this suite is not comparable with a row of the closed shell one. The steps
differ, the iterations inside each step differ, and an open shell builds two Fock
matrices where a closed shell builds one.
"""
import sys
from datetime import date
from pathlib import Path

import os

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from optbench import run
from scfbench import provenance, write

MOLECULE = sys.argv[1] if len(sys.argv) > 1 else "nitroxide"
AUX = "def2-universal-jkfit"
BASES = (sys.argv[2].split(",") if len(sys.argv) > 2 else ["def2-svp"])
CHARGE = int(sys.argv[3]) if len(sys.argv) > 3 else 0
MULTIPLICITY = int(sys.argv[4]) if len(sys.argv) > 4 else 2
FUNCTIONALS = ["HF", "B3LYP"]
METHODS = ["full", "ri_jk_simd"]

# NOTE: what the run actually used, so the table does not have to say
# "None threads" where the thread count belongs.
THREADS = int(os.environ.get("OMP_NUM_THREADS", os.cpu_count()))

ION = ("cation" if CHARGE > 0 else "anion" if CHARGE < 0
       else f"mult{MULTIPLICITY}")

# NOTE: the bases are in the name. Without them a second run of the same molecule
# on the same day writes to the same path and **replaces** the first run's record,
# which is what happened the first time this file was used: def2-svpd overwrote
# def2-svp and only the commit saved it.
STEM = f'{date.today():%Y-%m-%d}_m4max_{MOLECULE}_{ION}_{"_".join(BASES)}'

OUT = (Path(__file__).resolve().parent.parent / "data" / "optimization" /
       f'{STEM}.json')

# NOTE: and a run never overwrites another. One file per run is what makes every
# ratio inside a file comparable, and a record which is silently replaced is worse
# than one which is appended to.
if OUT.exists():
    raise SystemExit(
        f'{OUT} exists: a run of this molecule, spin state and bases was already '
        f'recorded today. Move it aside or name different bases.')

rows = []
for functional in FUNCTIONALS:
    for basis in BASES:
        for method in METHODS:
            row = run(MOLECULE, basis, AUX, method, functional,
                      charge=CHARGE, multiplicity=MULTIPLICITY)
            rows.append(row)
            print(f'  {functional:6s} {basis:10s} {method:12s}'
                  f'  {row["wall"]:9.2f} s  {row["steps"]:3d} steps'
                  f'  {row["wall_per_step"]:7.2f} s/step'
                  f'  E {row["energy"]:.8f}', flush=True)
            write(OUT, "optimization", provenance("m4max", 1, THREADS), rows)

print(f'\nwrote {OUT}', flush=True)

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
