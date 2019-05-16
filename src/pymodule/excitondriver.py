import numpy as np

from .veloxchemlib import KineticEnergyIntegralsDriver
from .veloxchemlib import NuclearPotentialIntegralsDriver
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

        # exciton monomers
        self.monomers = None

        # eri settings
        self.qq_type = 'QQ_DEN'
        self.eri_thresh = 1.0e-12

        # scf settings
        self.scf_conv_thresh = 1.0e-6
        self.scf_max_iter = 150

        # tda settings
        self.nstates = 2
        self.tda_conv_thresh = 1.0e-4
        self.tda_max_iter = 100

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream

    def compute(self, molecule, basis, min_basis, fragments_input):

        natoms = [int(n) for n in fragments_input.split(',')]
        assert_msg_critical(
            sum(natoms) == molecule.number_of_atoms(),
            'ExcitonModelDriver: Inconsistent number of atoms')

        if self.rank == mpi_master():
            self.print_title(natoms)

        nfragments = len(natoms)
        start_indices = [sum(natoms[:i]) for i in range(nfragments)]

        self.H = np.zeros(
            (nfragments * self.nstates, nfragments * self.nstates))

        self.monomers = [{} for i in range(nfragments)]

        for ind in range(nfragments):

            monomer_name = 'Monomer {}'.format(ind + 1)
            self.print_banner(monomer_name)

            # monomer molecule
            monomer = molecule.get_sub_molecule(start_indices[ind], natoms[ind])
            if self.rank == mpi_master():
                self.ostream.print_block(monomer.get_string())
                self.ostream.flush()

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

                for s in range(self.nstates):
                    h = s + ind * self.nstates
                    self.H[h, h] = self.monomers[ind]['exc_energies'][s]

        for ind_A in range(nfragments):
            monomer_a = molecule.get_sub_molecule(start_indices[ind_A],
                                                  natoms[ind_A])

            for ind_B in range(ind_A + 1, nfragments):
                monomer_b = molecule.get_sub_molecule(start_indices[ind_B],
                                                      natoms[ind_B])

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

                if self.rank == mpi_master():
                    kin_mat = kin_drv.compute(dimer, basis)
                    npot_mat = npot_drv.compute(dimer, basis)

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
                fock_mat.set_fock_type(fockmat.rgenjk, 0)

                eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
                screening = eri_drv.compute(get_qq_scheme(self.qq_type),
                                            self.eri_thresh, dimer, basis)
                eri_drv.compute(fock_mat, dens_mat, dimer, basis, screening)
                fock_mat.reduce_sum(self.rank, self.nodes, self.comm)

                if self.rank == mpi_master():
                    fock_mo = np.matmul(mo.T, np.matmul(fock_mat.to_numpy(0),
                                                        mo))
                    fock_occ = fock_mo[:nocc, :nocc]
                    fock_vir = fock_mo[nocc:, nocc:]

                    # compute dimer energy
                    dimer_energy = dimer.nuclear_repulsion_energy()
                    dimer_energy += fock_mat.get_energy(0, dens_mat, 0)
                    dimer_energy += 2.0 * kin_mat.get_energy(dens_mat, 0)
                    dimer_energy -= 2.0 * npot_mat.get_energy(dens_mat, 0)

                    valstr = 'Dimer Energy:{:20.10f} au'.format(dimer_energy)
                    self.ostream.print_header(valstr.ljust(92))
                    self.ostream.print_blank()
                    self.ostream.flush()

                    # assemble TDA CI vectors

                    #          vir.A     vir.B                vir.A     vir.B
                    #       +---------+---------+          +---------+---------+
                    # occ.A |  LE(A)  |         |    occ.A |         |         |
                    #       +---------+---------+          +---------+---------+
                    # occ.B |         |         |    occ.B |         |  LE(B)  |
                    #       +---------+---------+          +---------+---------+

                    #          vir.A     vir.B                vir.A     vir.B
                    #       +---------+---------+          +---------+---------+
                    # occ.A |         |  CT(AB) |    occ.A |         |         |
                    #       +---------+---------+          +---------+---------+
                    # occ.B |         |         |    occ.B |  CT(BA) |         |
                    #       +---------+---------+          +---------+---------+

                    vectors_A = self.monomers[ind_A]['exc_vectors']
                    vectors_B = self.monomers[ind_B]['exc_vectors']

                    CI_vectors_A = []
                    CI_vectors_B = []

                    for sA in range(self.nstates):
                        vec_A = vectors_A[:, sA].reshape(nocc_A, nvir_A)
                        CI_vec_A = np.zeros((nocc, nvir))
                        CI_vec_A[:nocc_A, :nvir_A] = vec_A[:, :]
                        CI_vectors_A.append(CI_vec_A)

                    for sB in range(self.nstates):
                        vec_B = vectors_B[:, sB].reshape(nocc_B, nvir_B)
                        CI_vec_B = np.zeros((nocc, nvir))
                        CI_vec_B[nocc_A:, nvir_A:] = vec_B[:, :]
                        CI_vectors_B.append(CI_vec_B)

                # compute LE-LE couplings
                if self.rank == mpi_master():
                    tdens_A = []
                    for sA in range(self.nstates):
                        tdens_A.append(
                            np.matmul(mo_occ,
                                      np.matmul(CI_vectors_A[sA], mo_vir.T)))
                    tdens_mat = AODensityMatrix(tdens_A, denmat.rest)
                else:
                    tdens_mat = AODensityMatrix()

                tdens_mat.broadcast(self.rank, self.comm)

                tfock_mat = AOFockMatrix(tdens_mat)
                for sA in range(self.nstates):
                    tfock_mat.set_fock_type(fockmat.rgenjk, sA)

                eri_drv.compute(tfock_mat, tdens_mat, dimer, basis, screening)
                tfock_mat.reduce_sum(self.rank, self.nodes, self.comm)

                if self.rank == mpi_master():
                    for sA in range(self.nstates):
                        tfock = tfock_mat.to_numpy(sA)
                        sigma_vec_A = np.matmul(mo_occ.T,
                                                np.matmul(tfock, mo_vir))

                        sigma_vec_A -= np.matmul(fock_occ, CI_vectors_A[sA])
                        sigma_vec_A += np.matmul(CI_vectors_A[sA], fock_vir)

                        for sB in range(self.nstates):
                            coupling = np.sum(sigma_vec_A * CI_vectors_B[sB])

                            hA = sA + ind_A * self.nstates
                            hB = sB + ind_B * self.nstates
                            self.H[hA, hB] = coupling
                            self.H[hB, hA] = coupling

                            valstr = 'LE-LE coupling:'
                            valstr += '  {}e({}){}g  {}g{}e({})'.format(
                                ind_A + 1, sA + 1, ind_B + 1, ind_A + 1,
                                ind_B + 1, sB + 1)
                            valstr += '  {:20.12f}'.format(coupling)
                            self.ostream.print_header(valstr.ljust(92))
                        self.ostream.print_blank()
                        self.ostream.flush()

        if self.rank == mpi_master():
            self.print_banner('Summary')

            eigvals, eigvecs = np.linalg.eigh(self.H)

            valstr = 'Eigenvalues:'
            self.ostream.print_header(valstr.ljust(92))
            for i, e in enumerate(eigvals):
                valstr = 'E[{}]= {:12.6f} a.u. {:12.5f} eV'.format(
                    i + 1, e, e * hartree_in_ev())
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
