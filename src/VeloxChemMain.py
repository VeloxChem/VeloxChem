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

    # set up input/output streams
    
    ostream = vlx.OutputStream(output_fname)
    istream = vlx.InputStream(input_fname, ostream)

    # print start header and set up start time
    
    if glob_rank == vlx.mpi_master():
        start_time = ostream.print_start_header(glob_nodes)

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

    if glob_rank == vlx.mpi_master():
        env_reader = vlx.EnvironmentReader()
        env_reader.parse(input_data, ostream)
        path_to_basis_lib = env_reader.get_path_to_basis_sets()

    # set up molecular basis reader on master node

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

    # initialize scf driver

    scf_drv = vlx.ScfRestrictedDriver()

    # read minimal basis if needed

    min_basis = vlx.MolecularBasis()

    if scf_drv.need_min_basis():
        if glob_rank == vlx.mpi_master():
            min_basis = molbas_reader.get_min_basis(path_to_basis_lib, mol_geom,
                                                    ostream)
        min_basis.broadcast(glob_rank, comm)

    scf_drv.compute(mol_geom, mol_basis, min_basis, comm, ostream)

    # write hdf5 files for MOs and density after SCF convergence

    if scf_drv.is_converged:
        scf_drv.mol_orbs.write_hdf5("mol_orbs.h5")
        scf_drv.density.write_hdf5("density.h5")

    # initialize visualization driver

    vis_drv = vlx.VisualizationDriver()

    # TODO: generate grid and write cube file
    # TODO: should also be able to compute psi/density along a line

    if scf_drv.is_converged:
        vis_grid = vis_drv.gen_grid(mol_geom)
        vis_drv.write_cube(mol_geom, mol_basis, scf_drv.mol_orbs, "HOMO",
                           "alpha")

    # all done, print finish header to output stream

    if glob_rank == vlx.mpi_master():
        ostream.print_finish_header(start_time)

if __name__ == "__main__":
    main()
