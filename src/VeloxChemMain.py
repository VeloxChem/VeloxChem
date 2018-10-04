#!/usr/bin/env python3

from mpi4py import MPI
import veloxchem as vlx

import sys
import os.path


def main():

    vlx.assert_msg_critical(vlx.mpi_initialized(), "MPI: Initialized")

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if len(sys.argv) <= 2:
        if rank == vlx.mpi_master():
            print("Usage: %s [input file] [output file]\n" % sys.argv[0])
            print("[input  file] - name of the input  file.")
            print("[output file] - name of the output file.")
        return

    input_fname = sys.argv[1]
    output_fname = sys.argv[2]

    if not os.path.isfile(input_fname):
        if rank == vlx.mpi_master():
            print("**** Error: %s does not exist ****" % input_fname)
        return

    if input_fname == output_fname:
        if rank == vlx.mpi_master():
            print("**** Error: input/output files cannot be the same ****")
        return

    app = vlx.AppManager(input_fname, output_fname)
    vlx.assert_msg_critical(app.get_state(), "AppManager: Initialization")

    app.execute()
    vlx.assert_msg_critical(app.get_state(), "AppManager: Execution")


if __name__ == "__main__":
    main()
