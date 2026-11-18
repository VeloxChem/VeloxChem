"""Standalone carbon structures: OpenMM system and files without solvent or solute.

carbon_system(model)        OpenMM System + Topology + positions for the carbon alone
write_carbon(model, outdir) structure.pdb (with CONECT), structure.xyz (angstrom),
                            forcefield.xml (nonbonded, one-atom residues), system.xml
                            (serialized, including the mobility terms), metadata.json

Slabs keep their periodic x/y cell and get vacuum along z; nanotubes keep their
periodic z length and get vacuum in x/y; flakes are non-periodic (NoCutoff).
Carbon-water NBFIX pairs of a carbon set are not in these files: they act only
between carbon and water and are added by assemble.build.
"""
from __future__ import annotations

import dataclasses
import io
import json
from pathlib import Path

import numpy as np
import openmm as mm
from openmm import app, unit

from . import params as P
from .forces import add_carbon_mobility, nonbonded_force
from .graphite import graphite_forcefield_xml, residue_names

NM = unit.nanometer


def carbon_system(model, *, vacuum_nm=3.0, cutoff_nm=1.0, switch_nm=0.9, minimize=False):
    """(System, Topology, positions nm, info) for the carbon model on its own.

    minimize=True first relaxes the structure with the flexible (GAFF2) terms and no
    restraint, then builds the requested mobility around the relaxed positions. Useful
    for flakes, whose ideal-lattice edges are not a force-field minimum; periodic
    slabs and nanotubes barely move (the ideal lattice is the minimum in its cell).
    """
    relaxed = None
    if minimize:
        flex = dataclasses.replace(model, mobility='flexible', restraint_k=0)
        s0, _, x0, _ = carbon_system(flex, vacuum_nm=vacuum_nm, cutoff_nm=cutoff_nm,
                                     switch_nm=switch_nm)
        ctx = mm.Context(s0, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('Reference'))
        ctx.setPositions(x0 * NM)
        mm.LocalEnergyMinimizer.minimize(ctx, 1.0, 20000)
        relaxed = ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(NM)
        relaxed = relaxed - (x0 - model.xyz)          # back to the model frame
        shift = float(np.abs(relaxed - model.xyz).max())
        model = dataclasses.replace(model, xyz=relaxed)
    xyz = np.array(model.xyz, float)
    box = None
    if model.kind == 'slab':
        lx, ly = model.box_xy
        box = np.array([lx, ly, model.thickness + vacuum_nm])
        xyz[:, 2] += box[2] / 2              # slab in the middle of the vacuum box
    elif model.kind == 'nanotube':
        side = 2 * model.info['radius_nm'] + vacuum_nm
        box = np.array([side, side, model.cell[2]])
        xyz[:, :2] += side / 2
    top = app.Topology()
    ch = top.addChain('G')
    elements = {'C': app.element.carbon, 'H': app.element.hydrogen, 'O': app.element.oxygen}
    atoms = [top.addAtom(e, elements[e], top.addResidue(r, ch))
             for e, r in zip(model.elements, residue_names(model))]
    if box is not None:
        top.setPeriodicBoxVectors([mm.Vec3(box[0], 0, 0), mm.Vec3(0, box[1], 0),
                                   mm.Vec3(0, 0, box[2])] * NM)
    ff = app.ForceField(io.StringIO(graphite_forcefield_xml(
        model.carbon_set, edge_charges=model.edge_charges)))
    method = app.CutoffPeriodic if box is not None else app.NoCutoff
    kw = dict(nonbondedCutoff=cutoff_nm * NM) if box is not None else {}
    system = ff.createSystem(top, nonbondedMethod=method, removeCMMotion=False,
                             residueTemplates={r: r.name for r in top.residues()}, **kw)
    nb = nonbonded_force(system)
    if box is not None and switch_nm:
        nb.setUseSwitchingFunction(True)
        nb.setSwitchingDistance(switch_nm * NM)
    nb.setUseDispersionCorrection(False)
    nb.setForceGroup(P.FORCE_GROUPS['nonbonded'])
    info = add_carbon_mobility(system, model, list(range(model.n_atoms)), xyz)
    for i, j in model.bonds():            # bonds only now: one-atom templates matched above
        top.addBond(atoms[i], atoms[j])
    return system, top, xyz, dict(box_nm=None if box is None else box.tolist(), mobility=info,
                                  minimized=dict(max_shift_nm=shift, tolerance_kj_mol_nm=1.0)
                                  if relaxed is not None else None)


def write_carbon(model, outdir, overwrite=False, **system_options):
    """Write the carbon model's structure and force field (see module docstring)."""
    out = Path(outdir)
    if out.exists() and any(out.iterdir()) and not overwrite:
        raise FileExistsError(f'{out} is not empty')
    out.mkdir(parents=True, exist_ok=True)
    system, top, xyz, info = carbon_system(model, **system_options)
    with open(out / 'structure.pdb', 'w') as fh:
        app.PDBFile.writeFile(top, xyz * NM, fh)
    lines = [str(model.n_atoms), f'{model.kind} {model.carbon_set} (angstrom)']
    lines += [f'{e} {10 * x:.6f} {10 * y:.6f} {10 * z:.6f}' for e, (x, y, z) in zip(model.elements, xyz)]
    (out / 'structure.xyz').write_text('\n'.join(lines) + '\n')
    (out / 'forcefield.xml').write_text(graphite_forcefield_xml(model.carbon_set,
                                                                edge_charges=model.edge_charges))
    (out / 'system.xml').write_text(mm.XmlSerializer.serialize(system))
    meta = dict(kind=model.kind, carbon_set=model.carbon_set,
                carbon_parameters=P.CARBON_SETS[model.carbon_set], n_atoms=model.n_atoms,
                n_carbon=model.n_carbon, n_bonds=len(model.bonds()), mass_da=model.mass,
                total_charge=float(model.charges().sum()), edge_charges=model.edge_charges,
                geometry=model.info, warning=model.warning, **info,
                note='Carbon-water NBFIX pairs (if any) are not in forcefield.xml/system.xml.')
    (out / 'metadata.json').write_text(json.dumps(meta, indent=1, default=str))
    return out
