import numpy as np
import time as tm
import math

from .veloxchemlib import KineticEnergyIntegralsDriver
from .veloxchemlib import NuclearPotentialIntegralsDriver
from .veloxchemlib import ElectricDipoleIntegralsDriver
from .veloxchemlib import LinearMomentumIntegralsDriver
from .veloxchemlib import AngularMomentumIntegralsDriver
from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import GridDriver
from .veloxchemlib import XCIntegrator
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .veloxchemlib import mpi_master
from .veloxchemlib import hartree_in_ev
from .veloxchemlib import rotatory_strength_in_cgs
from .veloxchemlib import get_dimer_ao_indices
from .veloxchemlib import parse_xc_func
from .molecule import Molecule
from .aodensitymatrix import AODensityMatrix
from .aofockmatrix import AOFockMatrix
from .scfrestdriver import ScfRestrictedDriver
from .rspabsorption import Absorption
from .errorhandler import assert_msg_critical
from .qqscheme import get_qq_scheme
from .lrmatvecdriver import read_rsp_hdf5
from .lrmatvecdriver import write_rsp_hdf5


class ExcitonModelDriver:
    """
    Implements the exciton model.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - H: The exciton model Hamiltonian matrix.
        - trans_dipoles: The diabatic electric transition dipole moments in
          length form.
        - velo_trans_dipoles: The diabatic electric transition dipole moments
          in velocity form.
        - magn_trans_dipoles: The diabatic magnetic transition dipole moments.
        - center_of_mass: The center of mass of the whole system.
        - state_info: Information of the diabatic excited states.
        - monomers: The monomer dictionaries.
        - natoms: The list containing number of atoms in each monomer.
        - charges: The list containing net charges of each monomer.
        - qq_type: The electron repulsion integrals screening scheme.
        - eri_thresh: The electron repulsion integrals screening threshold.
        - dft: The flag for running DFT.
        - grid_level: The accuracy level of DFT grid.
        - xcfun_label: The name of XC functional.
        - scf_conv_thresh: The convergence threshold for the SCF driver.
        - scf_max_iter: The maximum number of SCF iterations.
        - nstates: The number of locally excited states for each monomer.
        - ct_nocc: The number of occupied orbitals to be involved in
          charge-transfer excited states.
        - ct_nvir: The number of virtual orbitals to be involved in
          charge-transfer excited states.
        - tda_conv_thresh: The convergence threshold for the TDA driver.
        - tda_max_iter: The maximum number of TDA iterations.
        - comm: The MPI communicator.
        - rank: The MPI rank.
        - nodes: Number of MPI processes.
        - ostream: The output stream.
        - restart: The flag for restarting from checkpoint file.
        - checkpoint_file: Name of the exciton model checkpoint file.
    """

    def __init__(self, comm, ostream):
        """
        Initializes exciton model driver to default setup.
        """

        self.H = None
        self.trans_dipoles = None
        self.velo_trans_dipoles = None
        self.magn_trans_dipoles = None
        self.center_of_mass = None

        self.state_info = None

        # exciton monomers
        self.monomers = None
        self.natoms = None

        # eri settings
        self.qq_type = 'QQ_DEN'
        self.eri_thresh = 1.0e-12

        # dft settings
        self.dft = False
        self.grid_level = 4
        self.xcfun_label = 'Undefined'

        # scf settings
        self.scf_conv_thresh = 1.0e-6
        self.scf_max_iter = 150

        # tda settings
        self.nstates = 5
        self.ct_nocc = 1
        self.ct_nvir = 1
        self.tda_conv_thresh = 1.0e-4
        self.tda_max_iter = 100

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream

        # restart information
        self.restart = True
        self.checkpoint_file = None

    def update_settings(self, exciton_dict, method_dict={}):
        """
        Updates settings in exciton model driver.

        :param exciton_dict:
            The settings dictionary.
        """

        assert_msg_critical('fragments' in exciton_dict,
                            'ExcitonModel: fragments not defined')

        assert_msg_critical('atoms_per_fragment' in exciton_dict,
                            'ExcitonModel: atoms_per_fragment not defined')

        fragments = exciton_dict['fragments'].split(',')
        atoms_per_fragment = exciton_dict['atoms_per_fragment'].split(',')

        if 'charges' in exciton_dict:
            charges = exciton_dict['charges'].split(',')
        else:
            charges = ['0.0'] * len(fragments)

        assert_msg_critical(
            len(fragments) == len(atoms_per_fragment) and
            len(fragments) == len(charges),
            'ExcitonModel: mismatch in fragment input')

        self.natoms = []
        self.charges = []
        for n, x, q in zip(fragments, atoms_per_fragment, charges):
            self.natoms += [int(x)] * int(n)
            self.charges += [float(q)] * int(n)

        if 'nstates' in exciton_dict:
            self.nstates = int(exciton_dict['nstates'])

        if 'ct_nocc' in exciton_dict:
            self.ct_nocc = int(exciton_dict['ct_nocc'])
        if 'ct_nvir' in exciton_dict:
            self.ct_nvir = int(exciton_dict['ct_nvir'])

        if 'restart' in exciton_dict:
            key = exciton_dict['restart'].lower()
            self.restart = True if key == 'yes' else False
        if 'checkpoint_file' in exciton_dict:
            self.checkpoint_file = exciton_dict['checkpoint_file']

        if 'dft' in method_dict:
            key = method_dict['dft'].lower()
            self.dft = True if key == 'yes' else False
        if 'grid_level' in method_dict:
            self.grid_level = int(method_dict['grid_level'])
        if 'xcfun' in method_dict:
            if 'dft' not in method_dict:
                self.dft = True
            self.xcfun_label = method_dict['xcfun']

    def compute(self, molecule, basis, min_basis):
        """
        Executes exciton model calculation and writes checkpoint file.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis.
        :param min_basis:
            The minimal AO basis for generating initial guess.
        """

        # sanity check
        assert_msg_critical(
            sum(self.natoms) == molecule.number_of_atoms(),
            'ExcitonModelDriver: Inconsistent number of atoms')

        assert_msg_critical(self.nstates > 0,
                            'ExcitonModelDriver: Invalid number of LE states')

        assert_msg_critical(self.ct_nocc >= 0 and self.ct_nvir >= 0,
                            'ExcitonModelDriver: Invalid number of CT states')

        # exciton model setup
        nfragments = len(self.natoms)
        start_indices = [sum(self.natoms[:i]) for i in range(nfragments)]

        npairs = nfragments * (nfragments - 1) // 2
        ct_states = self.ct_nocc * self.ct_nvir

        total_LE_states = nfragments * self.nstates
        total_CT_states = npairs * ct_states * 2
        total_num_states = total_LE_states + total_CT_states

        self.H = np.zeros((total_num_states, total_num_states))
        self.trans_dipoles = np.zeros((total_num_states, 3))
        self.velo_trans_dipoles = np.zeros((total_num_states, 3))
        self.magn_trans_dipoles = np.zeros((total_num_states, 3))
        self.center_of_mass = molecule.center_of_mass()

        self.state_info = [{
            'type': '',
            'frag': '',
            'name': '',
        } for s in range(total_num_states)]

        self.monomers = [{} for i in range(nfragments)]

        if self.rank == mpi_master():
            self.print_title(total_LE_states, total_CT_states)

        # dimer indices
        dimer_id = np.zeros((nfragments, nfragments), dtype=np.int32)
        dimer_count = 0
        for ind_A in range(nfragments):
            dimer_id[ind_A, ind_A] = -1
            for ind_B in range(ind_A + 1, nfragments):
                dimer_id[ind_A, ind_B] = dimer_count
                dimer_id[ind_B, ind_A] = dimer_count
                dimer_count += 1

        # indices of diabatic excited states
        excitation_id = np.zeros((nfragments, nfragments), dtype=np.int32)
        for ind_A in range(nfragments):
            # LE
            excitation_id[ind_A, ind_A] = ind_A * self.nstates
            for ind_B in range(ind_A + 1, nfragments):
                ct_id = nfragments * self.nstates
                ct_id += dimer_id[ind_A][ind_B] * ct_states * 2
                # CT(A->B)
                excitation_id[ind_A, ind_B] = ct_id
                # CT(B->A)
                excitation_id[ind_B, ind_A] = ct_id + ct_states

        if self.dft:
            dft_func_label = self.xcfun_label
        else:
            dft_func_label = 'HF'

        rsp_vector_labels = [
            'dimer_indices',
            'num_states',
            'hamiltonian',
            'electric_transition_dipoles',
            'velocity_transition_dipoles',
            'magnetic_transition_dipoles',
            'state_info',
        ]

        # Monomers

        for ind in range(nfragments):

            monomer_start_time = tm.time()

            monomer_name = 'Monomer {}'.format(ind + 1)
            self.print_banner(monomer_name)

            # monomer molecule
            monomer = molecule.get_sub_molecule(start_indices[ind],
                                                self.natoms[ind])
            monomer.set_charge(self.charges[ind])
            monomer.check_multiplicity()

            if self.rank == mpi_master():
                self.ostream.print_block(monomer.get_string())

                valstr = 'Molecular charge            : {:.0f}'.format(
                    monomer.get_charge())
                self.ostream.print_header(valstr.ljust(70))
                valstr = 'Spin multiplicity           : {:d}'.format(
                    monomer.get_multiplicity())
                self.ostream.print_header(valstr.ljust(70))
                valstr = 'Number of alpha electrons   : {:d}'.format(
                    monomer.number_of_alpha_electrons())
                self.ostream.print_header(valstr.ljust(70))
                valstr = 'Number of beta  electrons   : {:d}'.format(
                    monomer.number_of_beta_electrons())
                self.ostream.print_header(valstr.ljust(70))
                self.ostream.print_blank()
                self.ostream.print_blank()

                self.ostream.print_block(
                    basis.get_string('Atomic Basis', monomer))
                self.ostream.flush()

            # 1e integral
            dipole_drv = ElectricDipoleIntegralsDriver(self.comm)
            dipole_drv.set_origin(*self.center_of_mass)
            dipole_mats = dipole_drv.compute(monomer, basis)

            linmom_drv = LinearMomentumIntegralsDriver(self.comm)
            linmom_mats = linmom_drv.compute(monomer, basis)

            angmom_drv = AngularMomentumIntegralsDriver(self.comm)
            angmom_drv.set_origin(*self.center_of_mass)
            angmom_mats = angmom_drv.compute(monomer, basis)

            if self.rank == mpi_master():
                dipole_ints = (dipole_mats.x_to_numpy(),
                               dipole_mats.y_to_numpy(),
                               dipole_mats.z_to_numpy())

                linmom_ints = (linmom_mats.x_to_numpy(),
                               linmom_mats.y_to_numpy(),
                               linmom_mats.z_to_numpy())

                angmom_ints = (angmom_mats.x_to_numpy(),
                               angmom_mats.y_to_numpy(),
                               angmom_mats.z_to_numpy())

            if self.dft:
                method_dict = {
                    'dft': 'yes',
                    'grid_level': self.grid_level,
                    'xcfun': self.xcfun_label
                }
            else:
                method_dict = {'dft': 'no'}

            # SCF calculation
            scf_drv = ScfRestrictedDriver(self.comm, self.ostream)
            scf_drv.update_settings(
                {
                    'qq_type': self.qq_type,
                    'eri_thresh': self.eri_thresh,
                    'conv_thresh': self.scf_conv_thresh,
                    'max_iter': self.scf_max_iter,
                    'restart': 'yes' if self.restart else 'no',
                    'checkpoint_file': 'monomer_{:d}.scf.h5'.format(ind + 1),
                }, method_dict)
            scf_drv.compute(monomer, basis, min_basis)

            if self.rank == mpi_master():
                self.monomers[ind]['C'] = scf_drv.scf_tensors['C'].copy()

            # update ERI screening threshold
            if self.eri_thresh > scf_drv.eri_thresh and scf_drv.eri_thresh > 0:
                self.eri_thresh = scf_drv.eri_thresh

            # TDA calculation
            abs_spec = Absorption(
                {
                    'tamm_dancoff': 'yes',
                    'nstates': self.nstates,
                    'qq_type': self.qq_type,
                    'eri_thresh': self.eri_thresh,
                    'conv_thresh': self.tda_conv_thresh,
                    'max_iter': self.tda_max_iter,
                    'restart': 'yes' if self.restart else 'no',
                    'checkpoint_file': 'monomer_{:d}.rsp.h5'.format(ind + 1),
                }, method_dict)
            abs_spec.init_driver(self.comm, self.ostream)
            abs_spec.compute(monomer, basis, scf_drv.scf_tensors)

            if self.rank == mpi_master():
                abs_spec.print_property(self.ostream)
                self.ostream.flush()

                self.monomers[ind]['exc_energies'] = abs_spec.get_property(
                    'eigenvalues').copy()
                self.monomers[ind]['exc_vectors'] = abs_spec.get_property(
                    'eigenvectors').copy()

                nocc = monomer.number_of_alpha_electrons()
                nvir = self.monomers[ind]['C'].shape[1] - nocc
                mo_occ = self.monomers[ind]['C'][:, :nocc]
                mo_vir = self.monomers[ind]['C'][:, nocc:]

                for s in range(self.nstates):
                    # LE excitation energy
                    h = excitation_id[ind, ind] + s
                    self.H[h, h] = self.monomers[ind]['exc_energies'][s]

                    # LE transition dipole
                    vec = self.monomers[ind]['exc_vectors'][:, s]
                    tdens = math.sqrt(2.0) * np.matmul(
                        mo_occ, np.matmul(vec.reshape(nocc, nvir), mo_vir.T))
                    self.trans_dipoles[h, :] = np.array(
                        [np.sum(tdens * dipole_ints[d]) for d in range(3)])
                    self.velo_trans_dipoles[h, :] = np.array(
                        [np.sum(tdens * linmom_ints[d]) for d in range(3)])
                    self.magn_trans_dipoles[h, :] = np.array(
                        [np.sum(tdens * angmom_ints[d]) for d in range(3)])

            valstr = '*** Time used in monomer calculation:'
            valstr += ' {:.2f} sec'.format(tm.time() - monomer_start_time)
            self.ostream.print_block(valstr.ljust(92))
            self.ostream.print_blank()

        # Dimers

        dimer_indices = None
        dimer_index_A = 0

        potfile_text = ''

        # read checkpoint file
        if self.restart:
            if self.rank == mpi_master():
                (dimer_indices, num_states, H, tdip, vdip, mdip,
                 state_info) = read_rsp_hdf5(
                     self.checkpoint_file, rsp_vector_labels,
                     molecule.nuclear_repulsion_energy(),
                     molecule.elem_ids_to_numpy(), basis.get_label(),
                     dft_func_label, potfile_text, self.ostream)
                read_success = (dimer_indices is not None and
                                num_states is not None and H is not None and
                                tdip is not None and vdip is not None and
                                mdip is not None and state_info is not None)
                if read_success:
                    shape_match = (
                        H.shape == self.H.shape and
                        num_states[0] == self.nstates and
                        num_states[1] == self.ct_nocc and
                        num_states[2] == self.ct_nvir and
                        tdip.shape == self.trans_dipoles.shape and
                        vdip.shape == self.velo_trans_dipoles.shape and
                        mdip.shape == self.magn_trans_dipoles.shape and
                        len(state_info) == len(self.state_info))
                    read_success = read_success and shape_match
                self.restart = read_success
            self.restart = self.comm.bcast(self.restart, root=mpi_master())

        if self.restart:
            dimer_indices = self.comm.bcast(dimer_indices, root=mpi_master())
            dimer_index_A = dimer_indices[0]

            if self.rank == mpi_master():
                self.H = H.copy()
                self.trans_dipoles = tdip.copy()
                self.velo_trans_dipoles = vdip.copy()
                self.magn_trans_dipoles = mdip.copy()

                for info_string in state_info:
                    info = info_string.decode('utf-8').split()
                    state_id = int(info[0])
                    if len(info) > 1:
                        self.state_info[state_id]['type'] = info[1]
                        self.state_info[state_id]['frag'] = info[2]
                        self.state_info[state_id]['name'] = info[3]

        for ind_A in range(dimer_index_A, nfragments):
            monomer_a = molecule.get_sub_molecule(start_indices[ind_A],
                                                  self.natoms[ind_A])
            monomer_a.set_charge(self.charges[ind_A])
            monomer_a.check_multiplicity()

            if self.restart and ind_A == dimer_index_A:
                dimer_index_B = dimer_indices[1] + 1
            else:
                dimer_index_B = ind_A + 1

            for ind_B in range(dimer_index_B, nfragments):
                monomer_b = molecule.get_sub_molecule(start_indices[ind_B],
                                                      self.natoms[ind_B])
                monomer_b.set_charge(self.charges[ind_B])
                monomer_b.check_multiplicity()

                dimer_start_time = tm.time()

                dimer_name = 'Dimer {} {}'.format(ind_A + 1, ind_B + 1)
                self.print_banner(dimer_name)

                # dimer molecule
                dimer = Molecule(monomer_a, monomer_b)
                dimer.check_multiplicity()

                if self.rank == mpi_master():
                    self.ostream.print_block(dimer.get_string())

                    valstr = 'Molecular charge            : {:.0f}'.format(
                        dimer.get_charge())
                    self.ostream.print_header(valstr.ljust(70))
                    valstr = 'Spin multiplicity           : {:d}'.format(
                        dimer.get_multiplicity())
                    self.ostream.print_header(valstr.ljust(70))
                    valstr = 'Number of alpha electrons   : {:d}'.format(
                        dimer.number_of_alpha_electrons())
                    self.ostream.print_header(valstr.ljust(70))
                    valstr = 'Number of beta  electrons   : {:d}'.format(
                        dimer.number_of_beta_electrons())
                    self.ostream.print_header(valstr.ljust(70))
                    self.ostream.print_blank()
                    self.ostream.print_blank()

                    self.ostream.print_block(
                        basis.get_string('Atomic Basis', dimer))
                    self.ostream.flush()

                # 1e integrals
                kin_drv = KineticEnergyIntegralsDriver(self.comm)
                kin_mat = kin_drv.compute(dimer, basis)

                npot_drv = NuclearPotentialIntegralsDriver(self.comm)
                npot_mat = npot_drv.compute(dimer, basis)

                dipole_drv = ElectricDipoleIntegralsDriver(self.comm)
                dipole_drv.set_origin(*self.center_of_mass)
                dipole_mats = dipole_drv.compute(dimer, basis)

                linmom_drv = LinearMomentumIntegralsDriver(self.comm)
                linmom_mats = linmom_drv.compute(dimer, basis)

                angmom_drv = AngularMomentumIntegralsDriver(self.comm)
                angmom_drv.set_origin(*self.center_of_mass)
                angmom_mats = angmom_drv.compute(dimer, basis)

                # dft grid
                if self.dft:
                    grid_drv = GridDriver(self.comm)
                    grid_drv.set_level(self.grid_level)

                    grid_t0 = tm.time()
                    dimer_molgrid = grid_drv.generate(dimer)
                    n_grid_points = dimer_molgrid.number_of_points()
                    dimer_molgrid.distribute(self.rank, self.nodes, self.comm)
                    self.ostream.print_info(
                        'Molecular grid with {} points generated in {:.2f} sec.'
                        .format(n_grid_points,
                                tm.time() - grid_t0))
                    self.ostream.print_blank()

                if self.rank == mpi_master():
                    dipole_ints = (dipole_mats.x_to_numpy(),
                                   dipole_mats.y_to_numpy(),
                                   dipole_mats.z_to_numpy())

                    linmom_ints = (linmom_mats.x_to_numpy(),
                                   linmom_mats.y_to_numpy(),
                                   linmom_mats.z_to_numpy())

                    angmom_ints = (angmom_mats.x_to_numpy(),
                                   angmom_mats.y_to_numpy(),
                                   angmom_mats.z_to_numpy())

                    # get indices of monomer AOs in dimer
                    ao_inds_A, ao_inds_B = get_dimer_ao_indices(
                        monomer_a, monomer_b, basis, basis)

                    # assemble MO coefficient matrix

                    #         occ.A     occ.B     vir.A     vir.B
                    #      +---------+---------+---------+---------+
                    #      |         |         |         |         |
                    # ao.A |  CoccA  |         |  CvirA  |         |
                    #      |         |         |         |         |
                    #      +---------+---------+---------+---------+
                    #      |         |         |         |         |
                    # ao.B |         |  CoccB  |         |  CvirB  |
                    #      |         |         |         |         |
                    #      +---------+---------+---------+---------+

                    CA = self.monomers[ind_A]['C']
                    CB = self.monomers[ind_B]['C']

                    nao_A = CA.shape[0]
                    nao_B = CB.shape[0]
                    nmo_A = CA.shape[1]
                    nmo_B = CB.shape[1]

                    nocc_A = monomer_a.number_of_alpha_electrons()
                    nocc_B = monomer_b.number_of_alpha_electrons()
                    nvir_A = nmo_A - nocc_A
                    nvir_B = nmo_B - nocc_B

                    nocc = nocc_A + nocc_B
                    nvir = nvir_A + nvir_B

                    mo = np.zeros((nao_A + nao_B, nmo_A + nmo_B))

                    for row in range(nao_A):
                        mo[ao_inds_A[row], :nocc_A] = CA[row, :nocc_A]
                        mo[ao_inds_A[row], nocc:nocc +
                           nvir_A] = CA[row, nocc_A:]

                    for row in range(nao_B):
                        mo[ao_inds_B[row], nocc_A:nocc] = CB[row, :nocc_B]
                        mo[ao_inds_B[row], nocc + nvir_A:] = CB[row, nocc_B:]

                    mo_occ = mo[:, :nocc]
                    mo_vir = mo[:, nocc:]

                    # compute density matrix
                    dens = np.matmul(mo_occ, mo_occ.T)
                    dens_mat = AODensityMatrix([dens], denmat.rest)
                else:
                    dens_mat = AODensityMatrix()

                dens_mat.broadcast(self.rank, self.comm)

                # compute Fock matrix
                fock_mat = AOFockMatrix(dens_mat)

                if self.dft:
                    xcfun = parse_xc_func(self.xcfun_label.upper())
                    if xcfun.is_hybrid():
                        fock_mat.set_fock_type(fockmat.restjkx, 0)
                        fock_mat.set_scale_factor(
                            xcfun.get_frac_exact_exchange(), 0)
                    else:
                        fock_mat.set_fock_type(fockmat.restj, 0)
                else:
                    fock_mat.set_fock_type(fockmat.restjk, 0)

                eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
                screening = eri_drv.compute(get_qq_scheme(self.qq_type),
                                            self.eri_thresh, dimer, basis)
                eri_drv.compute(fock_mat, dens_mat, dimer, basis, screening)
                fock_mat.reduce_sum(self.rank, self.nodes, self.comm)

                if self.dft:
                    if not xcfun.is_hybrid():
                        fock_mat.scale(2.0, 0)

                    xc_drv = XCIntegrator(self.comm)
                    vxc_mat = xc_drv.integrate(dens_mat, dimer, basis,
                                               dimer_molgrid, self.xcfun_label)
                    vxc_mat.reduce_sum(self.rank, self.nodes, self.comm)

                if self.rank == mpi_master():

                    # compute dimer energy

                    hcore = kin_mat.to_numpy() - npot_mat.to_numpy()
                    fock = hcore + fock_mat.to_numpy(0)

                    dimer_energy = dimer.nuclear_repulsion_energy()
                    dimer_energy += np.sum(dens * (hcore + fock))
                    if self.dft:
                        dimer_energy += vxc_mat.get_energy()

                    self.ostream.print_blank()

                    valstr = 'Excitonic Couplings'
                    self.ostream.print_header(valstr)
                    self.ostream.print_header('=' * (len(valstr) + 2))
                    self.ostream.print_blank()

                    valstr = 'Dimer Energy:{:20.10f} au'.format(dimer_energy)
                    self.ostream.print_header(valstr.ljust(72))
                    self.ostream.print_blank()
                    self.ostream.flush()

                    # compute Fock in MO basis

                    if self.dft:
                        fock += vxc_mat.get_matrix().to_numpy()

                    fock_mo = np.matmul(mo.T, np.matmul(fock, mo))
                    fock_occ = fock_mo[:nocc, :nocc]
                    fock_vir = fock_mo[nocc:, nocc:]

                    # assemble TDA CI vectors

                    vectors_A = self.monomers[ind_A]['exc_vectors']
                    vectors_B = self.monomers[ind_B]['exc_vectors']

                    CI_vectors = []

                    #          vir.A     vir.B                vir.A     vir.B
                    #       +---------+---------+          +---------+---------+
                    # occ.A |  LE(A)  |         |    occ.A |         |         |
                    #       +---------+---------+          +---------+---------+
                    # occ.B |         |         |    occ.B |         |  LE(B)  |
                    #       +---------+---------+          +---------+---------+

                    for sA in range(self.nstates):
                        vec_A = vectors_A[:, sA].reshape(nocc_A, nvir_A)
                        ci_vec = np.zeros((nocc, nvir))
                        ci_vec[:nocc_A, :nvir_A] = vec_A[:, :]
                        CI_vectors.append({
                            'vec': ci_vec,
                            'frag': 'A',
                            'type': 'LE',
                            'index': excitation_id[ind_A, ind_A] + sA,
                            'name': '{}e({})'.format(ind_A + 1, sA + 1),
                        })

                    for sB in range(self.nstates):
                        vec_B = vectors_B[:, sB].reshape(nocc_B, nvir_B)
                        ci_vec = np.zeros((nocc, nvir))
                        ci_vec[nocc_A:, nvir_A:] = vec_B[:, :]
                        CI_vectors.append({
                            'vec': ci_vec,
                            'frag': 'B',
                            'type': 'LE',
                            'index': excitation_id[ind_B, ind_B] + sB,
                            'name': '{}e({})'.format(ind_B + 1, sB + 1),
                        })

                    #          vir.A     vir.B                vir.A     vir.B
                    #       +---------+---------+          +---------+---------+
                    # occ.A |         |  CT(AB) |    occ.A |         |         |
                    #       +---------+---------+          +---------+---------+
                    # occ.B |         |         |    occ.B |  CT(BA) |         |
                    #       +---------+---------+          +---------+---------+

                    for oA in range(self.ct_nocc):
                        for vB in range(self.ct_nvir):
                            ctAB = oA * self.ct_nvir + vB
                            ci_vec = np.zeros((nocc, nvir))
                            ci_vec[nocc_A - 1 - oA][nvir_A + vB] = 1.0
                            CI_vectors.append({
                                'vec': ci_vec,
                                'frag': 'AB',
                                'type': 'CT',
                                'index': excitation_id[ind_A, ind_B] + ctAB,
                                'name': '{}+(H{}){}-(L{})'.format(
                                    ind_A + 1, oA, ind_B + 1, vB),
                            })

                    for oB in range(self.ct_nocc):
                        for vA in range(self.ct_nvir):
                            ctBA = oB * self.ct_nvir + vA
                            ci_vec = np.zeros((nocc, nvir))
                            ci_vec[nocc - 1 - oB][vA] = 1.0
                            CI_vectors.append({
                                'vec': ci_vec,
                                'frag': 'BA',
                                'type': 'CT',
                                'index': excitation_id[ind_B, ind_A] + ctBA,
                                'name': '{}-(L{}){}+(H{})'.format(
                                    ind_A + 1, vA, ind_B + 1, oB),
                            })

                    # update excited state information in self.state_info
                    for vec in CI_vectors:
                        state_id = vec['index']
                        self.state_info[state_id]['frag'] = vec['frag']
                        self.state_info[state_id]['type'] = vec['type']
                        self.state_info[state_id]['name'] = vec['name']

                    # create masked CI_vectors containing LE(A), CT(AB), CT(BA)
                    mask_ci_vectors = CI_vectors[:self.nstates] + CI_vectors[
                        self.nstates * 2:]

                # compute sigma vectors for LE(A), CT(AB), CT(BA)
                if self.rank == mpi_master():
                    tdens = []
                    for vec in mask_ci_vectors:
                        tdens.append(
                            np.matmul(mo_occ, np.matmul(vec['vec'], mo_vir.T)))
                    tdens_mat = AODensityMatrix(tdens, denmat.rest)
                else:
                    tdens_mat = AODensityMatrix()

                tdens_mat.broadcast(self.rank, self.comm)

                tfock_mat = AOFockMatrix(tdens_mat)

                if self.dft:
                    if xcfun.is_hybrid():
                        fock_flag = fockmat.rgenjkx
                        fact_xc = xcfun.get_frac_exact_exchange()
                        for s in range(tfock_mat.number_of_fock_matrices()):
                            tfock_mat.set_scale_factor(fact_xc, s)
                    else:
                        fock_flag = fockmat.rgenj
                else:
                    fock_flag = fockmat.rgenjk

                for s in range(tfock_mat.number_of_fock_matrices()):
                    tfock_mat.set_fock_type(fock_flag, s)

                eri_drv.compute(tfock_mat, tdens_mat, dimer, basis, screening)
                if self.dft:
                    if not xcfun.is_hybrid():
                        for s in range(tfock_mat.number_of_fock_matrices()):
                            tfock_mat.scale(2.0, s)
                    xc_drv.integrate(tfock_mat, tdens_mat, dens_mat, dimer,
                                     basis, dimer_molgrid, self.xcfun_label)
                tfock_mat.reduce_sum(self.rank, self.nodes, self.comm)

                if self.rank == mpi_master():
                    sigma_vectors = []

                    for s in range(tfock_mat.number_of_fock_matrices()):
                        tfock = tfock_mat.to_numpy(s)
                        sigma_vec = np.matmul(mo_occ.T,
                                              np.matmul(tfock, mo_vir))

                        sigma_vec -= np.matmul(fock_occ,
                                               mask_ci_vectors[s]['vec'])
                        sigma_vec += np.matmul(mask_ci_vectors[s]['vec'],
                                               fock_vir)

                        sigma_vectors.append({
                            'vec': sigma_vec,
                            'frag': mask_ci_vectors[s]['frag'],
                            'type': mask_ci_vectors[s]['type'],
                            'index': mask_ci_vectors[s]['index'],
                            'name': mask_ci_vectors[s]['name'],
                        })

                    # compute couplings
                    # sigma_vectors contains LE(A), CT(AB), CT(BA)
                    # CI_vectors[self.nstates:] contains LE(B), CT(AB), CT(BA)

                    for svec in sigma_vectors:
                        for cvec in CI_vectors[self.nstates:]:
                            if svec['index'] == cvec['index']:
                                continue

                            coupling = np.sum(svec['vec'] * cvec['vec'])

                            self.H[svec['index'], cvec['index']] = coupling
                            self.H[cvec['index'], svec['index']] = coupling

                            if svec['type'] == 'CT' and cvec['type'] == 'LE':
                                valstr = '{}-{} coupling:'.format(
                                    cvec['type'], svec['type'])
                                valstr += '  {:>15s}  {:>15s}'.format(
                                    cvec['name'], svec['name'])
                            else:
                                valstr = '{}-{} coupling:'.format(
                                    svec['type'], cvec['type'])
                                valstr += '  {:>15s}  {:>15s}'.format(
                                    svec['name'], cvec['name'])

                            valstr += '  {:20.12f}'.format(coupling)

                            if not (svec['type'] == 'CT' and
                                    cvec['type'] == 'CT' and
                                    svec['index'] > cvec['index']):
                                self.ostream.print_header(valstr.ljust(72))

                        self.ostream.print_blank()

                    # compute CT excitation energies and transition dipoles
                    # sigma_vectors[self.nstates:] contains CT(AB), CT(BA)

                    for ivec, svec in enumerate(sigma_vectors[self.nstates:]):

                        # the CI vector that corresponds to svec
                        cvec = CI_vectors[ivec + self.nstates * 2]

                        # CT excitation energy
                        energy = np.sum(svec['vec'] * cvec['vec'])

                        self.H[svec['index'], svec['index']] = energy

                        valstr = '{} excitation energy:'.format(svec['type'])
                        valstr += '  {:>26s}'.format(svec['name'])

                        valstr += '  {:20.12f}'.format(energy)
                        self.ostream.print_header(valstr.ljust(72))

                        # CT transition dipole
                        tdens = math.sqrt(2.0) * np.matmul(
                            mo_occ, np.matmul(cvec['vec'], mo_vir.T))
                        self.trans_dipoles[svec['index'], :] = np.array(
                            [np.sum(tdens * dipole_ints[d]) for d in range(3)])
                        self.velo_trans_dipoles[svec['index'], :] = np.array(
                            [np.sum(tdens * linmom_ints[d]) for d in range(3)])
                        self.magn_trans_dipoles[svec['index'], :] = np.array(
                            [np.sum(tdens * angmom_ints[d]) for d in range(3)])

                    self.ostream.print_blank()

                    # three-body CT-CT couplings

                    for ind_C in range(nfragments):
                        if ind_C == ind_A or ind_C == ind_B:
                            continue

                        # C(i)->A(a) vs C(i)->B(b): f_ab
                        for oC in range(self.ct_nocc):
                            for vA in range(self.ct_nvir):
                                for vB in range(self.ct_nvir):
                                    ctCA = oC * self.ct_nvir + vA
                                    ctCB = oC * self.ct_nvir + vB
                                    ctCA += excitation_id[ind_C, ind_A]
                                    ctCB += excitation_id[ind_C, ind_B]

                                    coupling = fock_vir[vA, nvir_A + vB]

                                    self.H[ctCA, ctCB] = coupling
                                    self.H[ctCB, ctCA] = coupling

                                    valstr = 'CT-CT coupling:'
                                    name_CA = '{}+(H{}){}-(L{})'.format(
                                        ind_C + 1, oC, ind_A + 1, vA)
                                    name_CB = '{}+(H{}){}-(L{})'.format(
                                        ind_C + 1, oC, ind_B + 1, vB)
                                    valstr += '  {:>15s}  {:>15s}'.format(
                                        name_CA, name_CB)
                                    valstr += '  {:20.12f}'.format(coupling)
                                    self.ostream.print_header(valstr.ljust(72))

                        # A(i)->C(a) vs B(j)->C(a): -f_ij
                        for vC in range(self.ct_nvir):
                            for oA in range(self.ct_nocc):
                                for oB in range(self.ct_nocc):
                                    ctAC = oA * self.ct_nvir + vC
                                    ctBC = oB * self.ct_nvir + vC
                                    ctAC += excitation_id[ind_A, ind_C]
                                    ctBC += excitation_id[ind_B, ind_C]

                                    coupling = -fock_occ[nocc_A - 1 - oA, nocc -
                                                         1 - oB]

                                    self.H[ctAC, ctBC] = coupling
                                    self.H[ctBC, ctAC] = coupling

                                    valstr = 'CT-CT coupling:'
                                    name_AC = '{}+(H{}){}-(L{})'.format(
                                        ind_A + 1, oA, ind_C + 1, vC)
                                    name_BC = '{}+(H{}){}-(L{})'.format(
                                        ind_B + 1, oB, ind_C + 1, vC)
                                    valstr += '  {:>15s}  {:>15s}'.format(
                                        name_AC, name_BC)
                                    valstr += '  {:20.12f}'.format(coupling)
                                    self.ostream.print_header(valstr.ljust(72))

                        if self.ct_nocc * self.ct_nvir > 0:
                            self.ostream.print_blank()

                valstr = '*** Time used in dimer calculation:'
                valstr += ' {:.2f} sec'.format(tm.time() - dimer_start_time)
                self.ostream.print_block(valstr.ljust(72))
                self.ostream.print_blank()

                state_info_list = [
                    '{} {} {} {}'.format(state_id, info['type'], info['frag'],
                                         info['name'])
                    for state_id, info in enumerate(self.state_info)
                ]

                rsp_vector_list = [
                    np.array([ind_A, ind_B]),
                    np.array([self.nstates, self.ct_nocc,
                              self.ct_nvir]), self.H, self.trans_dipoles,
                    self.velo_trans_dipoles, self.magn_trans_dipoles,
                    np.string_(state_info_list)
                ]

                if self.rank == mpi_master():
                    write_rsp_hdf5(self.checkpoint_file, rsp_vector_list,
                                   rsp_vector_labels,
                                   molecule.nuclear_repulsion_energy(),
                                   molecule.elem_ids_to_numpy(),
                                   basis.get_label(), dft_func_label,
                                   potfile_text, self.ostream)

        if self.rank == mpi_master():
            self.print_banner('Summary')

            eigvals, eigvecs = np.linalg.eigh(self.H)
            adia_trans_dipoles = np.matmul(eigvecs.T, self.trans_dipoles)
            adia_velo_trans_dipoles = np.matmul(eigvecs.T,
                                                self.velo_trans_dipoles)
            adia_magn_trans_dipoles = np.matmul(eigvecs.T,
                                                self.magn_trans_dipoles)

            for s in range(adia_velo_trans_dipoles.shape[0]):
                adia_velo_trans_dipoles[s, :] /= -eigvals[s]
                adia_magn_trans_dipoles[s, :] *= 0.5

            valstr = 'Adiabatic excited states:'
            self.ostream.print_header(valstr.ljust(80))
            self.ostream.print_blank()

            osc_str = []
            rot_str = []

            for i, e in enumerate(eigvals):
                dip_strength = np.sum(adia_trans_dipoles[i, :]**2)
                f = 2.0 / 3.0 * dip_strength * e
                osc_str.append(f)

                R = np.vdot(adia_velo_trans_dipoles[i, :],
                            adia_magn_trans_dipoles[i, :])
                R *= rotatory_strength_in_cgs()
                rot_str.append(R)

            self.print_absorption('One-Photon Absorption', eigvals, osc_str)
            self.print_ecd('Electronic Circular Dichroism', rot_str)

            valstr = 'Characters of excited states:'
            self.ostream.print_header(valstr.ljust(80))
            self.ostream.print_blank()

            for s in range(eigvecs.shape[1]):
                valstr = 'Excited state {}:'.format(s + 1)
                self.ostream.print_header(valstr.ljust(80))
                self.ostream.print_header(('-' * len(valstr)).ljust(80))

                for k in range(eigvecs.shape[0]):
                    composition = eigvecs[k, s]**2
                    if (composition > 0.05):
                        valstr = '{:<7s} {:<15s} {:>6.1f}%'.format(
                            self.state_info[k]['type'],
                            self.state_info[k]['name'],
                            composition * 100,
                        )
                        self.ostream.print_header(valstr.ljust(80))
                self.ostream.print_blank()

            self.ostream.flush()

    def print_banner(self, title):
        """
        Prints header.

        :param title:
            The text to be shown in the header.
        """

        valstr = '|' + ' ' * 10 + title + ' ' * 10 + '|'
        line = '+' + '-' * (len(valstr) - 2) + '+'
        self.ostream.print_header(line)
        self.ostream.print_header(valstr)
        self.ostream.print_header(line)
        self.ostream.print_blank()

    def print_title(self, num_LE, num_CT):
        """
        Prints exciton model title.

        :param num_LE:
            The total number of locally excited states.
        :param num_CT:
            The total number of charge-transfer excited states.
        """

        self.print_banner('Exciton Model')

        num_atoms = sum(self.natoms)
        num_frags = len(self.natoms)

        valstr = 'Total number of atoms:        {}'.format(num_atoms)
        self.ostream.print_header(valstr.ljust(54))
        valstr = 'Total number of monomers:     {}'.format(num_frags)
        self.ostream.print_header(valstr.ljust(54))
        valstr = 'Total number of LE states:    {}'.format(num_LE)
        self.ostream.print_header(valstr.ljust(54))
        valstr = 'Total number of CT states:    {}'.format(num_CT)
        self.ostream.print_header(valstr.ljust(54))
        self.ostream.print_blank()

        for i, n in enumerate(self.natoms):
            valstr = 'Monomer  {}  has  {}  atoms'.format(i + 1, n)
            self.ostream.print_header(valstr.ljust(54))
        self.ostream.print_blank()

    def print_absorption(self, title, eigvals, osc_str):
        """
        Prints absorption spectrum.

        :param title:
            The title.
        :param eigvals:
            The eigenvalues (adiabatic excitation energies).
        :param osc_str:
            The oscillator strengths.
        """

        valstr = title
        self.ostream.print_header(valstr.ljust(80))
        self.ostream.print_header(('-' * len(valstr)).ljust(80))
        for s, (e, f) in enumerate(zip(eigvals, osc_str)):
            valstr = 'Excited State {:>5s}: '.format('S' + str(s + 1))
            valstr += '{:15.8f} a.u. '.format(e)
            valstr += '{:12.5f} eV'.format(e * hartree_in_ev())
            valstr += '    Osc.Str. {:9.4f}'.format(f)
            self.ostream.print_header(valstr.ljust(80))
        self.ostream.print_blank()

    def print_ecd(self, title, rot_str):
        """
        Prints electronic circular dichroism spectrum.

        :param title:
            The title.
        :param rot_str:
            The rotatory strengths in 10**(-40) (esu**2)*(cm**2).
        """

        valstr = title
        self.ostream.print_header(valstr.ljust(80))
        self.ostream.print_header(('-' * len(valstr)).ljust(80))
        for s, R in enumerate(rot_str):
            valstr = 'Excited State {:>5s}: '.format('S' + str(s + 1))
            valstr += '    Rot.Str. {:11.4f}'.format(R)
            valstr += '    [10**(-40) (esu**2)*(cm**2)]'
            self.ostream.print_header(valstr.ljust(80))
        self.ostream.print_blank()
