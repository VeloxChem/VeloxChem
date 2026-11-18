"""HSO4- for the H2SO4 electrolyte: the same recipe as the solutes (make_solute.py).

    python make_hso4.py            # writes HSO4.pdb, HSO4.xml (residue HSO), HSO4_opt.xyz,
                                   # HSO4_provenance.json next to this file

1. B3LYP-D3/def2-SVP optimisation (no RI: VeloxChem has RI-J only for pure functionals), from the SMILES geometry, charge -1, singlet;
2. GAFF2 types/terms and RESP charges at HF/6-31G* (MMForceFieldGenerator.create_topology
   default) at the optimised geometry.
Not fitted to any H2SO4 solution property (status 'computed' in system_tools.SPECIES).
"""
import json
import os
import time
from pathlib import Path

import veloxchem as vlx

HERE = Path(__file__).resolve().parent
SMILES = 'S(=O)(=O)([O-])O'      # S first: it is the anion centre for slot counting
os.chdir(HERE)
t0 = time.time()
mol = vlx.Molecule.read_smiles(SMILES)
mol.set_charge(-1)
mol.set_multiplicity(1)
basis = vlx.MolecularBasis.read(mol, 'def2-svp')
scf = vlx.ScfRestrictedDriver()
scf.xcfun = 'b3lyp'
scf.dispersion = True
scf.conv_thresh = 1e-6
scf_res = scf.compute(mol, basis)
opt_res = vlx.OptimizationDriver(scf).compute(mol, basis, scf_res)
Path('HSO4_opt.xyz').write_text(opt_res['final_geometry'])
mol_opt = vlx.Molecule.read_xyz_string(opt_res['final_geometry'])
mol_opt.set_charge(-1)
mol_opt.set_multiplicity(1)
t1 = time.time()
ff = vlx.MMForceFieldGenerator()
ff.create_topology(mol_opt)                       # RESP HF/6-31G* + GAFF2
ff.write_openmm_files('HSO4', 'HSO')
n = mol_opt.number_of_atoms()
q = [ff.atoms[i]['charge'] for i in range(n)]
Path('HSO4_provenance.json').write_text(json.dumps(dict(
    smiles=SMILES, charge=-1, multiplicity=1,
    optimisation='B3LYP-D3/def2-SVP, VeloxChem OptimizationDriver',
    charges='RESP HF/6-31G* (MMForceFieldGenerator.create_topology default)', force_field='GAFF2',
    labels=mol_opt.get_labels(), resp_charges=q, gaff_types=[ff.atoms[i]['type'] for i in range(n)],
    total_charge=sum(q), veloxchem=getattr(vlx, '__version__', '?'),
    optimisation_s=t1 - t0, parametrisation_s=time.time() - t1), indent=1))
print(f'done: total charge {sum(q):+.4f}, {time.time() - t0:.0f} s')
