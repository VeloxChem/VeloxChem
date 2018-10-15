#!/usr/bin/env python3

from mpi4py import MPI
import veloxchem as vlx

import sys
import os.path
import time

def print_start_header(glob_nodes, ostream):
    """
    Prints start header to output stream.
    """
    # get start time
    
    start_time = time.time()
    
    # print start header
    
    ostream.put_separator()
    ostream.put_title("")
    ostream.put_title("VELOX CHEM MP (V.0.0 2018)")
    ostream.put_title("AN ELECTRONIC STRUCTURE CODE FOR NANOSCALE")
    ostream.put_title("")
    ostream.put_title("Copyright (c) 2018 Velox Chem MP developers.")
    ostream.put_title("All rights reserved.")
    ostream.put_separator()
    exec_str = "VeloxChem MP execution started"
    if glob_nodes > 1:
        exec_str += " on " + str(glob_nodes) + " compute nodes"
    exec_str += " at " + time.asctime( time.localtime(start_time) ) + "."
    ostream.put_title(exec_str)
    ostream.put_separator()
    ostream.new_line()

    # return start time

    return start_time

def print_finish_header(start_time, ostream):
    """
    Prints start header to output stream.
    """
    # get end time
    
    end_time = time.time()
    
    # print finish header
    
    ostream.put_separator()
    exec_str = "VeloxChem MP execution completed at ";
    exec_str += time.asctime( time.localtime(end_time) ) + "."
    ostream.put_title(exec_str)
    ostream.put_separator()
    exec_str = "Total execution time is "
    exec_str += str(int(end_time - start_time)) + " sec."
    ostream.put_title(exec_str)
    ostream.put_separator()

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

    # set up input/output streams
    
    ostream = vlx.OutputStream(output_fname)
    istream = vlx.InputStream(input_fname, ostream)

    # print start header and set up start time
    
    start_time = None
    
    if glob_rank == vlx.mpi_master():
        start_time = print_start_header(glob_nodes, ostream)

    # read input data from input file on master node
    
    input_data = vlx.InputData()
    
    if glob_rank == vlx.mpi_master():
        istream.read(input_data, ostream)
    
    # read molecular geometry
    
    mol_geom = vlx.Molecule()
    
    if glob_rank == vlx.mpi_master():
        molxyz_reader = vlx.MolXYZReader()
        molxyz_reader.parse(mol_geom, input_data, ostream)
        mol_geom.print_geometry(ostream)
        ostream.flush()
    
    mol_geom.broadcast(glob_rank, comm)

    # read environment variables on master node

    path_to_basis_lib = None

    if glob_rank == vlx.mpi_master():
        env_reader = vlx.EnvironmentReader()
        env_reader.parse(input_data, ostream)
        path_to_basis_lib = env_reader.get_path_to_basis_sets()

    # set up molecular basis reader on master node

    molbas_reader = None

    if glob_rank == vlx.mpi_master():
        molbas_reader = vlx.BasisReader()
        molbas_reader.parse(input_data, ostream)

    # read molecular basis

    mol_basis = vlx.MolecularBasis()

    if glob_rank == vlx.mpi_master():
        mol_basis = molbas_reader.get_ao_basis(path_to_basis_lib, mol_geom,
                                               ostream)
        mol_basis.print_basis("Atomic Basis", mol_geom, ostream)
        ostream.flush()

    mol_basis.broadcast(glob_rank, comm)

    # compute SCF energy, molecular orbitals

    scf_drv = vlx.ScfDriver()

    #scf_drv.acc_type = "DIIS"

    scf_drv.compute(comm, ostream)

    # all done, print finish header to output stream

    if glob_rank == vlx.mpi_master():
        print_finish_header(start_time, ostream)
        ostream.flush()

if __name__ == "__main__":
    main()
