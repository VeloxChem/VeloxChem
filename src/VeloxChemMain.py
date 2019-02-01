#!/usr/bin/env python3

from mpi4py import MPI
import veloxchem as vlx

import sys
import os.path
import time


def main():

    vlx.assert_msg_critical(vlx.mpi_initialized(), "MPI: Initialized")

    # set up MPI task

    task = vlx.MpiTask(sys.argv[1:], MPI.COMM_WORLD)

    task_type = task.input_dict['jobs']['task'].lower()

    # Hartree-Fock

    if task_type == 'hf':

        # initialize scf driver and run scf

        scf_drv = vlx.ScfRestrictedDriver()

        scf_drv.compute_task(task)

        # write hdf5 files after convergence

        if task.mpi_rank == vlx.mpi_master() and scf_drv.is_converged:
            scf_drv.mol_orbs.write_hdf5("mol_orbs.h5")
            scf_drv.density.write_hdf5("density.h5")

        # test multiple AO to MO matrices computation

        moints_drv = vlx.MOIntegralsDriver()

        moints_drv.compute_task(task, scf_drv.mol_orbs)

    # Cube

    elif task_type == 'cube':

        # generate cube file

        vis_drv = vlx.VisualizationDriver()

        if task.mpi_rank == vlx.mpi_master():
            mol_orbs = vlx.MolecularOrbitals.read_hdf5("mol_orbs.h5")
            density = vlx.AODensityMatrix.read_hdf5("density.h5")

            nelec = task.molecule.number_of_electrons()
            homo = nelec // 2 - 1

            vis_grid = vis_drv.gen_grid(task.molecule)
            vis_drv.write_cube(task.molecule, task.ao_basis, mol_orbs, homo,
                               "alpha", vis_grid)
            vis_drv.write_cube_dens(task.molecule, task.ao_basis, density, 0,
                                    "alpha", vis_grid)

    # all done, print finish header to output stream

    task.finish()


if __name__ == "__main__":
    main()
