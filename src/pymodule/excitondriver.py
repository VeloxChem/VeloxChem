import numpy as np

from .veloxchemlib import KineticEnergyIntegralsDriver
from .veloxchemlib import NuclearPotentialIntegralsDriver
from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .veloxchemlib import mpi_master
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

        nfragments = len(natoms)
        start_indices = [sum(natoms[:i]) for i in range(nfragments)]

        self.monomers = [{} for i in range(nfragments)]

        for ind in range(nfragments):

            # monomer molecule
            monomer = molecule.get_sub_molecule(start_indices[ind], natoms[ind])
            self.ostream.print_block(monomer.get_string())

            # SCF calculation
            scf_drv = ScfRestrictedDriver(self.comm, self.ostream)
            scf_drv.update_settings({
                'qq_type': self.qq_type,
                'eri_thresh': self.eri_thresh,
                'conv_thresh': self.scf_conv_thresh,
                'max_iter': self.scf_max_iter,
            })
            scf_drv.compute(monomer, basis, min_basis)

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

            self.monomers[ind]['exc_energies'] = abs_spec.get_property(
                'eigenvalues')
            self.monomers[ind]['exc_vectors'] = abs_spec.get_property(
                'eigenvectors')

        for ind_a in range(nfragments):
            monomer_a = molecule.get_sub_molecule(start_indices[ind_a],
                                                  natoms[ind_a])

            for ind_b in range(ind_a + 1, nfragments):
                monomer_b = molecule.get_sub_molecule(start_indices[ind_b],
                                                      natoms[ind_b])

                # dimer molecule
                dimer = Molecule(monomer_a, monomer_b)
                self.ostream.print_block(dimer.get_string())

                # 1e integrals
                kin_drv = KineticEnergyIntegralsDriver(self.comm)
                kin_mat = kin_drv.compute(dimer, basis)

                npot_drv = NuclearPotentialIntegralsDriver(self.comm)
                npot_mat = npot_drv.compute(dimer, basis)

                # get indices of monomer AOs in dimer
                ao_inds_1, ao_inds_2 = get_dimer_ao_indices(
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

                CA = self.monomers[ind_a]['C']
                CB = self.monomers[ind_b]['C']

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
                    mo[ao_inds_1[row], :nocc_A] = CA[row, :nocc_A]
                    mo[ao_inds_1[row], nocc:nocc + nvir_A] = CA[row, nocc_A:]

                for row in range(nao_B):
                    mo[ao_inds_2[row], nocc_A:nocc] = CB[row, :nocc_B]
                    mo[ao_inds_2[row], nocc + nvir_A:] = CB[row, nocc_B:]

                mo_occ = mo[:, :nocc]
                mo_vir = mo[:, nocc:]

                # compute density matrix
                dens = np.matmul(mo_occ, mo_occ.T)
                dens_mat = AODensityMatrix([dens], denmat.rest)
                dens_mat.broadcast(self.rank, self.comm)

                # compute Fock matrix
                fock_mat = AOFockMatrix(dens_mat)
                fock_mat.set_fock_type(fockmat.rgenjk, 0)

                eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
                screening = eri_drv.compute(get_qq_scheme(self.qq_type),
                                            self.eri_thresh, dimer, basis)
                eri_drv.compute(fock_mat, dens_mat, dimer, basis, screening)
                fock_mat.reduce_sum(self.rank, self.nodes, self.comm)

                fock_mo = np.matmul(mo.T, np.matmul(fock_mat.to_numpy(0), mo))
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

                vecs_A = self.monomers[ind_a]['exc_vectors']
                vecs_B = self.monomers[ind_b]['exc_vectors']

                for sA in range(vecs_A.shape[1]):
                    vec_A = vecs_A[:, sA].reshape(nocc_A, nvir_A)

                    CI_vec_A = np.zeros((nocc, nvir))
                    CI_vec_A[:nocc_A, :nvir_A] = vec_A[:, :]

                    tdens = np.matmul(mo_occ, np.matmul(CI_vec_A, mo_vir.T))

                    tdens_mat = AODensityMatrix([tdens], denmat.rest)
                    tdens_mat.broadcast(self.rank, self.comm)

                    tfock_mat = AOFockMatrix(tdens_mat)
                    tfock_mat.set_fock_type(fockmat.rgenjk, 0)

                    eri_drv.compute(tfock_mat, tdens_mat, dimer, basis,
                                    screening)
                    tfock_mat.reduce_sum(self.rank, self.nodes, self.comm)
                    tfock = tfock_mat.to_numpy(0)

                    sigma_vec_A = np.matmul(mo_occ.T, np.matmul(tfock, mo_vir))

                    sigma_vec_A -= np.matmul(fock_occ, CI_vec_A)
                    sigma_vec_A += np.matmul(CI_vec_A, fock_vir)

                    for sB in range(vecs_B.shape[1]):
                        vec_B = vecs_B[:, sB].reshape(nocc_B, nvir_B)

                        CI_vec_B = np.zeros((nocc, nvir))
                        CI_vec_B[nocc_A:, nvir_A:] = vec_B[:, :]

                        coupling = np.sum(sigma_vec_A * CI_vec_B)
                        valstr = 'LE-LE coupling:'
                        valstr += '  {}e({}){}g  {}g{}e({}) {:20.12f}'.format(
                            ind_a + 1, sA + 1, ind_b + 1, ind_a + 1, ind_b + 1,
                            sB + 1, coupling)
                        self.ostream.print_header(valstr.ljust(92))
                    self.ostream.print_blank()
