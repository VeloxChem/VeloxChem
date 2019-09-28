import numpy as np
import cppe

from .veloxchemlib import NuclearPotentialIntegralsDriver
from .veloxchemlib import ElectricFieldIntegralsDriver


class PolEmbed:

    def __init__(self, molecule, basis, comm, ostream, potfile, iso_pol=True):
        self.molecule = molecule
        self.basis = basis
        self.comm = comm
        self.ostream = ostream

        self.options = cppe.PeOptions()
        self.options.potfile = potfile
        self.options.iso_pol = iso_pol

        cppe_mol = cppe.Molecule()
        coords = np.vstack(
            (self.molecule.x_to_numpy(), self.molecule.y_to_numpy(),
             self.molecule.z_to_numpy())).T
        charges = self.molecule.elem_ids_to_numpy()

        for z, coord in zip(charges, coords):
            cppe_mol.append(cppe.Atom(int(z), *coord))
        self.cppe_state = cppe.CppeState(self.options, cppe_mol,
                                         self.print_header)
        self.cppe_state.calculate_static_energies_and_fields()
        self._enable_induction = False
        if self.cppe_state.get_polarizable_site_number():
            self._enable_induction = True
            coords = np.array([
                site.position
                for site in self.cppe_state.potentials
                if site.is_polarizable
            ])
            self.polarizable_coords = coords
        self.V_es = None

    def print_header(self, output):
        self.ostream.print_header(output)

    def get_pe_contribution(self, dm, elec_only=False):
        if self.V_es is None:
            self.V_es = self.compute_multipole_potential_integrals()

        if not elec_only:
            e_el = np.sum(self.V_es * dm)
            self.cppe_state.energies["Electrostatic"]["Electronic"] = e_el

        V_ind = np.zeros_like(self.V_es)
        if self._enable_induction:
            pass
            elec_fields = self.compute_electric_field_value(dm)
            # solve induced moments
            self.cppe_state.update_induced_moments(elec_fields.flatten(),
                                                   elec_only)
            induced_moments = np.array(
                self.cppe_state.get_induced_moments()).reshape(
                    self.polarizable_coords.shape)
            V_ind = self.compute_induction_operator(induced_moments)

        e = self.cppe_state.total_energy
        if not elec_only:
            vmat = self.V_es + V_ind
        else:
            vmat = V_ind
            e = self.cppe_state.energies["Polarization"]["Electronic"]
        return e, vmat

    def compute_multipole_potential_integrals(self):
        sites = np.empty((0, 3), dtype=np.float)
        dipole_sites = []
        charges = []
        dipoles = []
        for p in self.cppe_state.potentials:
            site = p.position
            for m in p.multipoles:
                # for now, we only do charges!
                if m.k == 0:
                    charges.append(m.values[0])
                    sites = np.vstack((sites, site))
                elif m.k == 1:
                    dipoles.append(m.values)
                    dipole_sites.append(site)
                else:
                    raise NotImplementedError(
                        "PE electrostatics only implemented through"
                        " first order.")
        # compute the 0th order operator (charges)
        np_charges = np.array(charges)
        np_dipoles = np.array(dipoles)
        npot_drv = NuclearPotentialIntegralsDriver(self.comm)
        V_es = -1.0 * npot_drv.compute(self.molecule, self.basis, np_charges,
                                       sites).to_numpy()
        if len(dipole_sites):
            ef_driver = ElectricFieldIntegralsDriver(self.comm)
            ret = ef_driver.compute(self.molecule, self.basis, np_dipoles,
                                    np.array(dipole_sites))
            V_es += -1.0 * (ret.x_to_numpy() + ret.y_to_numpy() +
                            ret.z_to_numpy())
        return V_es

    def compute_induction_operator(self, moments):
        ef_driver = ElectricFieldIntegralsDriver(self.comm)
        ret = ef_driver.compute(self.molecule, self.basis, moments,
                                self.polarizable_coords)
        V_ind = -1.0 * (ret.x_to_numpy() + ret.y_to_numpy() + ret.z_to_numpy())
        return V_ind

    def compute_electric_field_value(self, dm):
        ef_driver = ElectricFieldIntegralsDriver(self.comm)
        elec_field = np.zeros_like(self.polarizable_coords)
        for i, coord in enumerate(self.polarizable_coords):
            ret = ef_driver.compute(self.molecule, self.basis, *coord)
            elec_field[i] = [
                np.sum(dm * ret.x_to_numpy()),
                np.sum(dm * ret.y_to_numpy()),
                np.sum(dm * ret.z_to_numpy()),
            ]
        return elec_field
