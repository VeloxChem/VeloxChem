"""The range separated molecular gradient on the laptop: four centres against RI-JK simd.

    python grad_rs_laptop.py [molecule] [basis,basis,...] [functional,...]

The gradient of a hybrid range separated functional costs both ways a second
exchange. The four-centre way makes a whole further pass of its own derivative
kernels, `kx_rs` on top of `2jkx`, once per atom. The simd resolution of the
identity contracts a second set of B vectors against a second derivative tensor
which the same kernel wrote, on one sparsity pattern and in one sweep.

**The conventional RI-JK column is not here**, as in the range separated SCF suite:
that driver has no attenuated B vectors and refuses the functional.

**A row of this suite is not comparable with a row of `grad_laptop.py`.** It is the
gradient of a different functional, not the same one computed differently. What the
B3LYP rows of the plain suite are good for is the *ratio* -- what the split costs
each way -- and that is how the two tables are read against each other.
"""
import os
import sys
from datetime import date
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from gradbench import run
from scfbench import provenance, write

MOLECULE = sys.argv[1] if len(sys.argv) > 1 else "caffeine"
AUX = "def2-universal-jkfit"
BASES = (sys.argv[2].split(",") if len(sys.argv) > 2
         else ["def2-svp", "def2-svpd", "def2-tzvp", "def2-tzvpd"])
FUNCTIONALS = (sys.argv[3].split(",") if len(sys.argv) > 3 else ["CAM-B3LYP"])
METHODS = ["full", "ri_jk_simd"]

# NOTE: a tag on the records, so a control run of a plain hybrid does not overwrite
# the range separated one. The cost of the split is an absolute wall time divided by
# an absolute wall time, and that division is only meaningful when both were measured
# in one sitting on one tree -- so the control is run from this same file, minutes
# after, and lands beside it rather than on top of it.
TAG = sys.argv[4] if len(sys.argv) > 4 else "rs"

THREADS = int(os.environ.get("OMP_NUM_THREADS", os.cpu_count()))

rows = []
for functional in FUNCTIONALS:
    for basis in BASES:
        for method in METHODS:
            row = run(MOLECULE, basis, AUX, method, functional)
            rows.append(row)
            print(f'  {functional:10s} {basis:10s} {method:12s}'
                  f'  scf {row["scf_wall"]:8.2f}  grad {row["grad_wall"]:8.2f}'
                  f'  {row["grad_walls"]}', flush=True)
            out = (Path(__file__).resolve().parent.parent / "data" / "gradient" /
                   f'{date.today():%Y-%m-%d}_m4max_{MOLECULE}_{TAG}.json')
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
