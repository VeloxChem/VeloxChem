import numpy as np
import math

from .veloxchemlib import KineticEnergyIntegralsDriver
from .veloxchemlib import NuclearPotentialIntegralsDriver
from .veloxchemlib import ElectricDipoleIntegralsDriver
from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .veloxchemlib import mpi_master
from .veloxchemlib import hartree_in_ev
from .veloxchemlib import get_dimer_ao_indices
from .molecule import Molecule
from .aodensitymatrix import AODensityMatrix
from .aofockmatrix import AOFockMatrix
from .scfrestdriver import ScfRestrictedDriver
from .rspabsorption import Absorption
from .errorhandler import assert_msg_critical
from .qqscheme import get_qq_scheme


class ExcitonModelDriver:

    def __init__(self, comm, ostream):

        self.H = None
        self.trans_dipoles = None
        self.nuc_chg_center = None

        # exciton monomers
        self.monomers = None
        self.natoms = None

        # eri settings
        self.qq_type = 'QQ_DEN'
        self.eri_thresh = 1.0e-12

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

    def update_settings(self, exciton_dict):

        assert_msg_critical('fragments' in exciton_dict,
                            'ExcitonModelDriver: fragments not defined')

        self.natoms = []
        natoms_list = exciton_dict['fragments'].replace(' ', '').split(',')
        for x in natoms_list:
            if '*' in x:
                content = x.split('*')
                self.natoms += [int(content[0])] * int(content[1])
            elif x:
                self.natoms.append(int(x))

        if 'nstates' in exciton_dict:
            self.nstates = int(exciton_dict['nstates'])

        if 'ct_nocc' in exciton_dict:
            self.ct_nocc = int(exciton_dict['ct_nocc'])
        if 'ct_nvir' in exciton_dict:
            self.ct_nvir = int(exciton_dict['ct_nvir'])

    def compute(self, molecule, basis, min_basis):

        assert_msg_critical(
            sum(self.natoms) == molecule.number_of_atoms(),
            'ExcitonModelDriver: Inconsistent number of atoms')

        assert_msg_critical(self.nstates > 0,
                            'ExcitonModelDriver: Invalid number of LE states')

        if self.rank == mpi_master():
            self.print_title(self.natoms)

        nfragments = len(self.natoms)
        start_indices = [sum(self.natoms[:i]) for i in range(nfragments)]

        npairs = nfragments * (nfragments - 1) // 2
        ct_states = self.ct_nocc * self.ct_nvir

        self.H = np.zeros((nfragments * self.nstates + npairs * ct_states * 2,
                           nfragments * self.nstates + npairs * ct_states * 2))

        self.trans_dipoles = np.zeros(
            (nfragments * self.nstates + npairs * ct_states * 2, 3))

        self.nuc_chg_center = molecule.center_of_nuclear_charge()

        self.monomers = [{} for i in range(nfragments)]

        # Monomers

        for ind in range(nfragments):

            monomer_name = 'Monomer {}'.format(ind + 1)
            self.print_banner(monomer_name)

            # monomer molecule
            monomer = molecule.get_sub_molecule(start_indices[ind],
                                                self.natoms[ind])
            if self.rank == mpi_master():
                self.ostream.print_block(monomer.get_string())
                self.ostream.flush()

            # 1e integral
            dipole_drv = ElectricDipoleIntegralsDriver(self.comm)
            dipole_drv.set_origin(*self.nuc_chg_center)
            dipole_mats = dipole_drv.compute(monomer, basis)

            if self.rank == mpi_master():
                dipole_ints = (dipole_mats.x_to_numpy(),
                               dipole_mats.y_to_numpy(),
                               dipole_mats.z_to_numpy())

            # SCF calculation
            scf_drv = ScfRestrictedDriver(self.comm, self.ostream)
            scf_drv.update_settings({
                'qq_type': self.qq_type,
                'eri_thresh': self.eri_thresh,
                'conv_thresh': self.scf_conv_thresh,
                'max_iter': self.scf_max_iter,
            })
            scf_drv.compute(monomer, basis, min_basis)

            if self.rank == mpi_master():
                self.monomers[ind]['C'] = scf_drv.scf_tensors['C'].copy()
                self.monomers[ind]['F'] = scf_drv.scf_tensors['F'][0].copy()

            # TDA calculation
            abs_spec = Absorption({
                'nstates': self.nstates,
                'qq_type': self.qq_type,
                'eri_thresh': self.eri_thresh,
                'conv_thresh': self.tda_conv_thresh,
                'max_iter': self.tda_max_iter,
            })
            abs_spec.init_driver(self.comm, self.ostream)
            abs_spec.compute(monomer, basis, scf_drv.scf_tensors)

            if self.rank == mpi_master():
                abs_spec.print_property(self.ostream)
                self.ostream.flush()

                self.monomers[ind]['exc_energies'] = abs_spec.get_property(
                    'eigenvalues')
                self.monomers[ind]['exc_vectors'] = abs_spec.get_property(
                    'eigenvectors')

                nocc = monomer.number_of_alpha_electrons()
                nvir = self.monomers[ind]['C'].shape[1] - nocc
                mo_occ = self.monomers[ind]['C'][:, :nocc]
                mo_vir = self.monomers[ind]['C'][:, nocc:]

                for s in range(self.nstates):
                    h = s + ind * self.nstates

                    # LE excitation energy
                    self.H[h, h] = self.monomers[ind]['exc_energies'][s]

                    # LE transition dipole
                    vec = self.monomers[ind]['exc_vectors'][:, s]
                    tdens = math.sqrt(2.0) * np.matmul(
                        mo_occ, np.matmul(vec.reshape(nocc, nvir), mo_vir.T))
                    self.trans_dipoles[h, :] = np.array(
                        [np.sum(tdens * dipole_ints[d]) for d in range(3)])

        # Dimers

        dimer_id = -1

        for ind_A in range(nfragments):
            monomer_a = molecule.get_sub_molecule(start_indices[ind_A],
                                                  self.natoms[ind_A])

            for ind_B in range(ind_A + 1, nfragments):
                monomer_b = molecule.get_sub_molecule(start_indices[ind_B],
                                                      self.natoms[ind_B])

                dimer_id += 1

                dimer_name = 'Dimer {} {}'.format(ind_A + 1, ind_B + 1)
                self.print_banner(dimer_name)

                # dimer molecule
                dimer = Molecule(monomer_a, monomer_b)
                if self.rank == mpi_master():
                    self.ostream.print_block(dimer.get_string())
                    self.ostream.flush()

                # 1e integrals
                kin_drv = KineticEnergyIntegralsDriver(self.comm)
                npot_drv = NuclearPotentialIntegralsDriver(self.comm)
                dipole_drv = ElectricDipoleIntegralsDriver(self.comm)
                dipole_drv.set_origin(*self.nuc_chg_center)

                kin_mat = kin_drv.compute(dimer, basis)
                npot_mat = npot_drv.compute(dimer, basis)
                dipole_mats = dipole_drv.compute(dimer, basis)

                if self.rank == mpi_master():
                    dipole_ints = (dipole_mats.x_to_numpy(),
                                   dipole_mats.y_to_numpy(),
                                   dipole_mats.z_to_numpy())

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
                fock_mat.set_fock_type(fockmat.restjk, 0)

                eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
                screening = eri_drv.compute(get_qq_scheme(self.qq_type),
                                            self.eri_thresh, dimer, basis)
                eri_drv.compute(fock_mat, dens_mat, dimer, basis, screening)
                fock_mat.reduce_sum(self.rank, self.nodes, self.comm)

                if self.rank == mpi_master():
                    hcore = kin_mat.to_numpy() - npot_mat.to_numpy()
                    fock = hcore + fock_mat.to_numpy(0)

                    fock_mo = np.matmul(mo.T, np.matmul(fock, mo))
                    fock_occ = fock_mo[:nocc, :nocc]
                    fock_vir = fock_mo[nocc:, nocc:]

                    # compute dimer energy
                    dimer_energy = dimer.nuclear_repulsion_energy()
                    dimer_energy += np.sum(dens * (hcore + fock))

                    valstr = 'Dimer Energy:{:20.10f} au'.format(dimer_energy)
                    self.ostream.print_header(valstr.ljust(92))
                    self.ostream.print_blank()
                    self.ostream.flush()

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
                            'index': ind_A * self.nstates + sA,
                            'name': '{}e({}){}g'.format(ind_A + 1, sA + 1,
                                                        ind_B + 1),
                        })

                    for sB in range(self.nstates):
                        vec_B = vectors_B[:, sB].reshape(nocc_B, nvir_B)
                        ci_vec = np.zeros((nocc, nvir))
                        ci_vec[nocc_A:, nvir_A:] = vec_B[:, :]
                        CI_vectors.append({
                            'vec': ci_vec,
                            'frag': 'B',
                            'type': 'LE',
                            'index': ind_B * self.nstates + sB,
                            'name': '{}g{}e({})'.format(ind_A + 1, ind_B + 1,
                                                        sB + 1),
                        })

                    #          vir.A     vir.B                vir.A     vir.B
                    #       +---------+---------+          +---------+---------+
                    # occ.A |         |  CT(AB) |    occ.A |         |         |
                    #       +---------+---------+          +---------+---------+
                    # occ.B |         |         |    occ.B |  CT(BA) |         |
                    #       +---------+---------+          +---------+---------+

                    ct_ind = nfragments * self.nstates
                    ct_ind += dimer_id * ct_states * 2

                    for oA in range(self.ct_nocc):
                        for vB in range(self.ct_nvir):
                            ci_vec = np.zeros((nocc, nvir))
                            ci_vec[nocc_A - 1 - oA][nvir_A + vB] = 1.0
                            CI_vectors.append({
                                'vec': ci_vec,
                                'frag': 'AB',
                                'type': 'CT',
                                'index': ct_ind,
                                'name': '{}+(H{}){}-(L{})'.format(
                                    ind_A + 1, oA, ind_B + 1, vB),
                            })
                            ct_ind += 1

                    for oB in range(self.ct_nocc):
                        for vA in range(self.ct_nvir):
                            ci_vec = np.zeros((nocc, nvir))
                            ci_vec[nocc - 1 - oB][vA] = 1.0
                            CI_vectors.append({
                                'vec': ci_vec,
                                'frag': 'BA',
                                'type': 'CT',
                                'index': ct_ind,
                                'name': '{}-(L{}){}+(H{})'.format(
                                    ind_A + 1, vA, ind_B + 1, oB),
                            })
                            ct_ind += 1

                # compute sigma vectors for LE(A), CT(AB), CT(BA)
                if self.rank == mpi_master():
                    tdens = []
                    for vec in CI_vectors:
                        if vec['frag'] == 'B':
                            continue
                        tdens.append(
                            np.matmul(mo_occ, np.matmul(vec['vec'], mo_vir.T)))
                    tdens_mat = AODensityMatrix(tdens, denmat.rest)
                else:
                    tdens_mat = AODensityMatrix()

                tdens_mat.broadcast(self.rank, self.comm)

                tfock_mat = AOFockMatrix(tdens_mat)
                for s in range(tfock_mat.number_of_fock_matrices()):
                    tfock_mat.set_fock_type(fockmat.rgenjk, s)

                eri_drv.compute(tfock_mat, tdens_mat, dimer, basis, screening)
                tfock_mat.reduce_sum(self.rank, self.nodes, self.comm)

                if self.rank == mpi_master():
                    sigma_vectors = []

                    for s in range(tfock_mat.number_of_fock_matrices()):
                        tfock = tfock_mat.to_numpy(s)
                        sigma_vec = np.matmul(mo_occ.T,
                                              np.matmul(tfock, mo_vir))

                        # skip LE(B) in CI_vectors
                        # since sigma_vectors does not contain LE(B)
                        s2 = s
                        if s >= self.nstates:
                            s2 += self.nstates

                        sigma_vec -= np.matmul(fock_occ, CI_vectors[s2]['vec'])
                        sigma_vec += np.matmul(CI_vectors[s2]['vec'], fock_vir)

                        sigma_vectors.append({
                            'vec': sigma_vec,
                            'frag': CI_vectors[s2]['frag'],
                            'type': CI_vectors[s2]['type'],
                            'index': CI_vectors[s2]['index'],
                            'name': CI_vectors[s2]['name'],
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
                            self.ostream.print_header(valstr.ljust(92))

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
                        self.ostream.print_header(valstr.ljust(92))

                        # CT transition dipole
                        tdens = math.sqrt(2.0) * np.matmul(
                            mo_occ, np.matmul(cvec['vec'], mo_vir.T))
                        self.trans_dipoles[svec['index'], :] = np.array(
                            [np.sum(tdens * dipole_ints[d]) for d in range(3)])

                    self.ostream.print_blank()

        if self.rank == mpi_master():
            self.print_banner('Summary')

            eigvals, eigvecs = np.linalg.eigh(self.H)
            adia_trans_dipoles = np.matmul(eigvecs.T, self.trans_dipoles)

            valstr = 'Adiabatic excited states:'
            self.ostream.print_header(valstr.ljust(92))
            self.ostream.print_blank()

            for i, e in enumerate(eigvals):
                enestr = 'E[{}]='.format(i + 1)
                dip_strength = np.sum(adia_trans_dipoles[i, :]**2)
                osc_strength = 2.0 / 3.0 * dip_strength * e

                valstr = '{:<8s} {:12.6f} a.u. {:12.5f} eV'.format(
                    enestr, e, e * hartree_in_ev())
                valstr += '    osc.str. {:12.5f}'.format(osc_strength)
                self.ostream.print_header(valstr.ljust(92))

            self.ostream.print_blank()
            self.ostream.flush()

    def print_banner(self, title):

        valstr = '|' + ' ' * 10 + title + ' ' * 10 + '|'
        line = '+' + '-' * (len(valstr) - 2) + '+'
        self.ostream.print_header(line)
        self.ostream.print_header(valstr)
        self.ostream.print_header(line)
        self.ostream.print_blank()

    def print_title(self, natoms):

        self.print_banner('Exciton Model')

        valstr = 'Total number of atoms:        {}'.format(sum(natoms))
        self.ostream.print_header(valstr.ljust(72))
        valstr = 'Total number of monomers:     {}'.format(len(natoms))
        self.ostream.print_header(valstr.ljust(72))
        valstr = 'Total number of LE states:    {}'.format(
            len(natoms) * self.nstates)
        self.ostream.print_header(valstr.ljust(72))
        self.ostream.print_blank()

        for i, n in enumerate(natoms):
            valstr = 'Monomer  {}  has  {}  atoms'.format(i + 1, n)
            self.ostream.print_header(valstr.ljust(72))
        self.ostream.print_blank()
