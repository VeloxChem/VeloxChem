import numpy as np
try:
    import cppe
except ImportError:
    raise ImportError(
        'Unable to import cppe; '
        'please download and install from https://github.com/maxscheurer/cppe')

from .veloxchemlib import NuclearPotentialIntegralsDriver
from .veloxchemlib import ElectricFieldIntegralsDriver
from .subcommunicators import SubCommunicators
from .veloxchemlib import mpi_master


class PolEmbed:

    def __init__(self, molecule, basis, comm, potfile, iso_pol=True):
        self.molecule = molecule
        self.basis = basis
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()
        self.output = ''
        self.V_es = None

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
                                         self.print_callback)
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

    def print_callback(self, output):
        self.output += output

    def get_pe_contribution(self, dm, elec_only=False):
        if self.V_es is None:
            self.V_es = self.compute_multipole_potential_integrals()

        if not elec_only:
            if self.rank == mpi_master():
                e_el = np.sum(self.V_es * dm)
                self.cppe_state.energies["Electrostatic"]["Electronic"] = e_el

        V_ind = np.zeros(self.V_es.shape)

        if self._enable_induction:
            elec_fields = self.compute_electric_field_value(dm)
            # solve induced moments
            if self.rank == mpi_master():
                self.cppe_state.update_induced_moments(elec_fields.flatten(),
                                                       elec_only)
                induced_moments = np.array(
                    self.cppe_state.get_induced_moments()).reshape(
                        self.polarizable_coords.shape)
            else:
                induced_moments = None
            induced_moments = self.comm.bcast(induced_moments,
                                              root=mpi_master())
            V_ind = self.compute_induction_operator(induced_moments)

        if self.rank == mpi_master():
            if not elec_only:
                vmat = self.V_es + V_ind
                e = self.cppe_state.total_energy
            else:
                vmat = V_ind
                e = self.cppe_state.energies["Polarization"]["Electronic"]
            return e, vmat
        else:
            return 0.0, None

    def compute_multipole_potential_integrals(self):
        sites = []
        dipole_sites = []
        charges = []
        dipoles = []

        for p in self.cppe_state.potentials:
            site = p.position
            for m in p.multipoles:
                # zeroth order
                if m.k == 0:
                    charges.append(m.values[0])
                    sites.append(site)
                # first order
                elif m.k == 1:
                    dipoles.append(m.values)
                    dipole_sites.append(site)
                else:
                    raise NotImplementedError(
                        "PE electrostatics only implemented through first order"
                    )

        V_es = np.zeros(0)

        npot_drv = NuclearPotentialIntegralsDriver(self.comm)
        if self.rank == mpi_master():
            V_es = -1.0 * npot_drv.compute(self.molecule, self.basis,
                                           np.array(charges),
                                           np.array(sites)).to_numpy()

        if dipole_sites:
            ef_driver = ElectricFieldIntegralsDriver(self.comm)
            ret = ef_driver.compute(self.molecule, self.basis,
                                    np.array(dipoles), np.array(dipole_sites))
            if self.rank == mpi_master():
                V_es += -1.0 * (ret.x_to_numpy() + ret.y_to_numpy() +
                                ret.z_to_numpy())

        return V_es

    def compute_induction_operator(self, moments):

        node_grps = [p for p in range(self.nodes)]
        subcomm = SubCommunicators(self.comm, node_grps)
        local_comm = subcomm.local_comm
        cross_comm = subcomm.cross_comm

        ave, res = divmod(self.polarizable_coords.shape[0], self.nodes)
        counts = [ave + 1 if p < res else ave for p in range(self.nodes)]

        start = sum(counts[:self.rank])
        end = sum(counts[:self.rank + 1])

        ef_driver = ElectricFieldIntegralsDriver(local_comm)
        if start < end:
            ret = ef_driver.compute(self.molecule, self.basis,
                                    moments[start:end, :],
                                    self.polarizable_coords[start:end, :])
            V_size = ret.x_to_numpy().shape[0]
        else:
            V_size = 0

        if self.polarizable_coords.shape[0] < self.nodes:
            V_size = cross_comm.bcast(V_size, root=mpi_master())

        if local_comm.Get_rank() == mpi_master():
            if start < end:
                V_ind = -1.0 * (ret.x_to_numpy() + ret.y_to_numpy() +
                                ret.z_to_numpy())
            else:
                V_ind = np.zeros((V_size, V_size))
            V_ind = cross_comm.reduce(V_ind, root=mpi_master())
        else:
            V_ind = np.zeros(0)

        return V_ind

    def compute_electric_field_value(self, dm):

        node_grps = [p for p in range(self.nodes)]
        subcomm = SubCommunicators(self.comm, node_grps)
        local_comm = subcomm.local_comm
        cross_comm = subcomm.cross_comm

        ave, res = divmod(self.polarizable_coords.shape[0], self.nodes)
        counts = [ave + 1 if p < res else ave for p in range(self.nodes)]

        start = sum(counts[:self.rank])
        end = sum(counts[:self.rank + 1])

        ef_driver = ElectricFieldIntegralsDriver(local_comm)
        elec_field = []
        for i in range(start, end):
            coord = self.polarizable_coords[i]
            ret = ef_driver.compute(self.molecule, self.basis, *coord)
            if local_comm.Get_rank() == mpi_master():
                elec_field.append(np.sum(dm * ret.x_to_numpy()))
                elec_field.append(np.sum(dm * ret.y_to_numpy()))
                elec_field.append(np.sum(dm * ret.z_to_numpy()))

        if local_comm.Get_rank() == mpi_master():
            elec_field = cross_comm.gather(elec_field, root=mpi_master())
        else:
            elec_field = []

        if self.rank == mpi_master():
            elec_field = np.array([xyz for ef in elec_field for xyz in ef])
            elec_field = elec_field.reshape(self.polarizable_coords.shape)

        return elec_field
