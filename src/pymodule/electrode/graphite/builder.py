"""GraphiteBuilder: builds the carbon electrode and nothing else.

    gb = GraphiteBuilder('wu_aluru', layers=3, lateral_nm=(5, 5))     # Bernal slab, frozen
    gb.add_step(edge='zigzag').add_groups({'carboxylate': 2}, seed=1)  # optional surface features
    mol = gb.molecule()            # VeloxChem Molecule (angstrom), e.g. to pack solvent around
    xml = gb.forcefield_xml()      # OpenMM residue templates for the carbon atoms
    gb.apply_to_system(system, atoms, positions_nm)    # mobility, charges, group terms, exclusions
    gb.add_water_correction(system, carbon, water_O, water_H, cutoff_nm, switch_nm)   # NBFIX, if any
    gb.write('out/graphite')       # standalone files (no solvent)

Solvent, solute, field and FEP are added on top of this (veloxchem.electrode.build_system).
"""
from __future__ import annotations

import numpy as np

from . import surface
from .export import carbon_system, write_carbon
from .forces import add_carbon_mobility, add_nbfix_correction
from .graphite import (CarbonModel, graphite_forcefield_xml, make_flake, make_nanotube, make_slab,
                       residue_names)
from .params import CARBON_SETS


UC_CM2_TO_E_NM2 = 0.0624151      # 1 uC/cm2 in e/nm2


class GraphiteBuilder:
    """A carbon electrode model: geometry, surface features and force-field terms.

    The constructor builds a periodic slab (layers=1: graphene) centred on z = 0 (see
    graphite.make_slab for stacking, orientation, mobility and restraint_k). Use
    GraphiteBuilder.flake / .nanotube for finite flakes and tubes.
    """

    def __init__(self, carbon_set='steele', layers=3, lateral_nm=(5.0, 5.0), *, stacking='AB',
                 orientation='armchair_x', mobility='frozen', restraint_k=None,
                 allow_unverified=False, model=None):
        self.model = model if model is not None else make_slab(
            carbon_set, layers, tuple(lateral_nm), stacking=stacking, orientation=orientation,
            mobility=mobility, restraint_k=restraint_k, allow_unverified=allow_unverified)

    @classmethod
    def flake(cls, carbon_set='steele', size_nm=(2.0, 2.0), **kw):
        return cls(model=make_flake(carbon_set, size_nm, **kw))

    @classmethod
    def nanotube(cls, carbon_set, n, m, cells=4, **kw):
        return cls(model=make_nanotube(carbon_set, n, m, cells, **kw))

    @classmethod
    def from_model(cls, model):
        return model if isinstance(model, cls) else cls(model=model)

    def set_face_charge(self, sigma_top_uC_cm2, sigma_bottom_uC_cm2=None):
        """Constant-charge electrode (2026-11-16; the Elliott 2023 / Wang 2025 control for the
        uniform-field method): a uniform surface charge on the outer carbon layers, split
        evenly over their lattice C atoms. sigma_top on the top layer (the face the liquid sits
        on), sigma_bottom on the bottom layer (default -sigma_top: neutral slab). In
        microcoulomb/cm2; added to any existing charges. Use with the applied field off.
        For a slab: the top face is then a positive electrode at sigma_top > 0."""
        import numpy as np
        if sigma_bottom_uC_cm2 is None:
            sigma_bottom_uC_cm2 = -sigma_top_uC_cm2
        m = self.model
        if m.kind != 'slab':
            raise ValueError('set_face_charge: slab only')
        area = float(m.box_xy[0]) * float(m.box_xy[1])
        q = np.array(m.charges(), float)
        r = self.roles()
        per = {}
        for key, s in (('top_layer', sigma_top_uC_cm2), ('bottom_layer', sigma_bottom_uC_cm2)):
            idx = r[key]
            dq = s * UC_CM2_TO_E_NM2 * area / len(idx)
            q[idx] += dq
            per[key] = dq
        m.charge = q
        m.info['face_charge'] = dict(sigma_top_uC_cm2=float(sigma_top_uC_cm2),
                                     sigma_bottom_uC_cm2=float(sigma_bottom_uC_cm2),
                                     q_per_atom_top_e=per['top_layer'], q_per_atom_bottom_e=per['bottom_layer'],
                                     area_nm2=area)
        return self

    # ---- surface features (outer layer of `face`); each returns self, so calls chain
    def add_step(self, **kw):
        self.model = surface.add_step(self.model, **kw)
        return self

    def add_defects(self, **kw):
        self.model = surface.add_defects(self.model, **kw)
        return self

    def add_groups(self, counts, **kw):
        self.model = surface.add_groups(self.model, counts, **kw)
        return self

    # ---- description
    @property
    def kind(self):
        return self.model.kind

    @property
    def carbon_set(self):
        return self.model.carbon_set

    @property
    def n_atoms(self):
        return self.model.n_atoms

    @property
    def elements(self):
        return list(self.model.elements)

    @property
    def positions_nm(self):
        return np.array(self.model.xyz)

    @property
    def box_xy(self):
        return self.model.box_xy

    @property
    def thickness(self):
        return self.model.thickness

    @property
    def top_z(self):
        return self.model.top_z

    @property
    def layers(self):
        return self.model.layers

    @property
    def mass(self):
        return self.model.mass

    @property
    def charge(self):
        return float(round(self.model.charges().sum(), 6))

    @property
    def has_surface_features(self):
        m = self.model
        return bool(m.features) or not m.skeleton.all() or 'H' in m.elements

    def warnings(self):
        w = [self.model.warning] if self.model.warning else []
        return w + [x for x in self.model.info.get('warnings', []) if x not in w]

    def roles(self):
        """Index lists in model order: carbon (lattice C), carbon_h (edge H), surface_groups,
        top_layer / bottom_layer (lattice C of the outer layers)."""
        m = self.model
        skel = m.skeleton
        last = m.info.get('layers', 1) - 1
        out = {k: [] for k in ('carbon', 'carbon_h', 'surface_groups', 'top_layer', 'bottom_layer')}
        for i, (e, lay) in enumerate(zip(m.elements, m.layer)):
            if not skel[i]:
                out['surface_groups'].append(i)
            elif e != 'C':
                out['carbon_h'].append(i)
            else:
                out['carbon'].append(i)
                if lay == 0:
                    out['top_layer'].append(i)
                if lay == last:
                    out['bottom_layer'].append(i)
        return out

    def metadata(self):
        m = self.model
        groups = {}
        for f in m.features:
            if f['kind'] == 'group':
                groups[f['group']] = groups.get(f['group'], 0) + 1
        meta = dict(kind=m.kind, carbon_set=m.carbon_set, carbon_parameters=CARBON_SETS[m.carbon_set],
                    n_atoms=m.n_atoms, n_carbon=m.n_carbon, charge=self.charge, mass_da=m.mass,
                    mobility=m.mobility, restraint_k=m.restraint_k, geometry=dict(m.info),
                    features=m.features, groups=groups,
                    n_group_atoms=int((~m.skeleton).sum()),
                    n_edge_h=int(sum(1 for e, s in zip(m.elements, m.skeleton) if s and e == 'H')))
        if m.kind == 'slab':
            meta.update(layers=m.layers, thickness_nm=m.thickness, top_face_z=m.top_z,
                        box_xy_nm=list(m.box_xy))
        return meta

    # ---- outputs for code that builds the full system
    def molecule(self):
        """VeloxChem Molecule of the carbon model (angstrom, charge set)."""
        return self.model.to_vlx_molecule()

    def forcefield_xml(self):
        """OpenMM ForceField XML (one-atom GRA/GRH/GRO residues) for createSystem."""
        return graphite_forcefield_xml(self.model.carbon_set)

    def residue_names(self):
        return residue_names(self.model)

    def apply_to_system(self, system, atoms, positions_nm):
        """Carbon terms in a System created from forcefield_xml(): per-atom charges and group
        LJ, mobility (frozen / restrained / flexible), surface-group bonded terms and
        exclusions. atoms: system indices of the model atoms, in model order."""
        return add_carbon_mobility(system, self.model, atoms, positions_nm)

    def add_water_correction(self, system, carbon_atoms, water_O, water_H, cutoff_nm, switch_nm=None):
        """Carbon x water pair overrides of the carbon set (e.g. Wu-Aluru), as a correction
        on top of Lorentz-Berthelot. carbon_atoms: system indices of the lattice carbons.
        Returns the force, or None if the set has no overrides."""
        return add_nbfix_correction(system, carbon_atoms, water_O, water_H, self.model.carbon_set,
                                    cutoff_nm, switch_nm)

    # ---- standalone
    def system(self, **kw):
        """(System, Topology, positions nm, info) of the carbon alone (export.carbon_system)."""
        return carbon_system(self.model, **kw)

    def write(self, outdir, **kw):
        return write_carbon(self.model, outdir, **kw)


__all__ = ['GraphiteBuilder', 'CarbonModel']
