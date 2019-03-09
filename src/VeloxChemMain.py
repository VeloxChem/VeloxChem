#!/usr/bin/env python3

from mpi4py import MPI
import veloxchem as vlx

import sys
import os.path
import time
import numpy as np


def main():

    vlx.assert_msg_critical(vlx.mpi_initialized(), "MPI: Initialized")

    # set up MPI task

    task = vlx.MpiTask(sys.argv[1:], MPI.COMM_WORLD)

    task_type = task.input_dict['jobs']['task'].lower()

    # Hartree-Fock

    if task_type in ['hf', 'mp2', 'cube', 'response']:

        # initialize scf driver and run scf

        scf_drv = vlx.ScfRestrictedDriver()

        scf_drv.compute_task(task)

        # write hdf5 files after convergence

        if task.mpi_rank == vlx.mpi_master() and scf_drv.is_converged:
            scf_drv.mol_orbs.write_hdf5("mol_orbs.h5")
            scf_drv.density.write_hdf5("density.h5")

    # Response

    if task_type == 'response':

        # molecular orbitals

        if task.mpi_rank == vlx.mpi_master():
            mol_orbs = scf_drv.mol_orbs
        else:
            mol_orbs = vlx.MolecularOrbitals()

        rsp_drv = vlx.ResponseDriver()

        rsp_drv.compute_task(mol_orbs, task)

    # MP2 perturbation theory

    if task_type == 'mp2':

        # molecular orbitals

        if task.mpi_rank == vlx.mpi_master():
            mol_orbs = scf_drv.mol_orbs
        else:
            mol_orbs = vlx.MolecularOrbitals()

        # MO integrals

        moints_drv = vlx.MOIntegralsDriver()
        oovv = moints_drv.compute_task(task, mol_orbs, "OOVV")

        # MP2 energy

        nocc = task.molecule.number_of_alpha_electrons()
        e_mp2 = moints_drv.compute_mp2_energy(mol_orbs, nocc, oovv)

        if task.mpi_rank == vlx.mpi_master():
            mp2_str = "*** MP2 correlation energy: %20.12f a.u." % e_mp2
            task.ostream.print_header(mp2_str.ljust(92))
            task.ostream.print_blank()

    # Cube

    if task_type == 'cube':

        # generate cube file

        vis_drv = vlx.VisualizationDriver()

        if task.mpi_rank == vlx.mpi_master():
            mol_orbs = vlx.MolecularOrbitals.read_hdf5("mol_orbs.h5")
            density = vlx.AODensityMatrix.read_hdf5("density.h5")

            nelec = task.molecule.number_of_electrons()
            homo = nelec // 2 - 1

            vis_grid = vis_drv.gen_grid(task.molecule)
            vis_drv.write_cube('homo.cube', task.molecule, task.ao_basis,
                               mol_orbs, homo, "alpha", vis_grid)
            vis_drv.write_cube('density.cube', task.molecule, task.ao_basis,
                               density, 0, "alpha", vis_grid)

    # all done, print finish header to output stream

    task.finish()


if __name__ == "__main__":
    main()
