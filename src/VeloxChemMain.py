#!/usr/bin/env python3

from mpi4py import MPI
from VeloxChemLib import AppManager
from VeloxChemLib import assert_msg_critical
from VeloxChemLib import mpi_initialized
from VeloxChemLib import mpi_master

import sys
import os.path


def main():

    assert_msg_critical(mpi_initialized(), "MPI: Uninitialized")

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if len(sys.argv) <= 2:
        if rank == mpi_master():
            print("Usage: VeloxChemMP.x [input file] [output file]\n")
            print("[input  file] - name of the input  file.")
            print("[output file] - name of the output file.")
        return

    input_fname = sys.argv[1]
    output_fname = sys.argv[2]

    if not os.path.isfile(input_fname):
        if rank == mpi_master():
            print("**** Error: %s does not exist ****" % input_fname)
        return

    if input_fname == output_fname:
        if rank == mpi_master():
            print("**** Error: input/output files cannot be the same ****")
        return

    app = AppManager(input_fname, output_fname)
    assert_msg_critical(app.get_state(), "AppManager: Initialization failed")

    app.execute()
    assert_msg_critical(app.get_state(), "AppManager: Execution failed")


if __name__ == "__main__":
    main()
