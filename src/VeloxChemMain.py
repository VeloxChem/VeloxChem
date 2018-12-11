#!/usr/bin/env python3

from mpi4py import MPI
import veloxchem as vlx

import sys
import os.path
import time


def main():

    vlx.assert_msg_critical(vlx.mpi_initialized(), "MPI: Initialized")

    # set up MPI communicator data

    comm = MPI.COMM_WORLD
    glob_rank = comm.Get_rank()
    glob_nodes = comm.Get_size()

    # read command line parameters on master node

    input_fname = ""
    output_fname = ""

    if glob_rank == vlx.mpi_master():

        # input syntax error

        if len(sys.argv) <= 2:
            print("Usage: %s [input file] [output file]\n" % sys.argv[0])
            print("[input  file] - name of the input  file.")
            print("[output file] - name of the output file.")
            return

        # set up names of input and output files

        input_fname = sys.argv[1]
        output_fname = sys.argv[2]

        # critical error: input file does not exist

        if not os.path.isfile(input_fname):
            print("**** Error: %s does not exist ****" % input_fname)
            return

        # critical error: input file and output files have same name

        if input_fname == output_fname:
            print("**** Error: input/output files cannot be the same ****")
            return

    # set up MPI task

    task = vlx.MpiTask(input_fname, output_fname, comm)

    # initialize scf driver

    scf_drv = vlx.ScfRestrictedDriver()

    # read minimal basis if needed

    scf_drv.compute_task(task)

    # write hdf5 files for MOs and density after SCF convergence

    if glob_rank == vlx.mpi_master():
        if scf_drv.is_converged:
            scf_drv.mol_orbs.write_hdf5("mol_orbs.h5")
            scf_drv.density.write_hdf5("density.h5")

    # generate cube file

    """
    vis_drv = vlx.VisualizationDriver()

    if glob_rank == vlx.mpi_master():
        mol_orbs = vlx.MolecularOrbitals.read_hdf5("mol_orbs.h5")
        density = vlx.AODensityMatrix.read_hdf5("density.h5")

        nelec = task.molecule.number_of_electrons()
        homo = nelec // 2 - 1

        vis_grid = vis_drv.gen_grid(task.molecule)
        vis_drv.write_cube(task.molecule, task.ao_basis, scf_drv.mol_orbs,
                           homo, "alpha", vis_grid)
        vis_drv.write_cube_dens(task.molecule, task.ao_basis,
                                scf_drv.density, 0, "alpha", vis_grid)
    """

    # all done, print finish header to output stream

    task.finish()


if __name__ == "__main__":
    main()
