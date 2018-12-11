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

    # set up output stream

    ostream = vlx.OutputStream(output_fname)

    # print start header and set up start time

    if glob_rank == vlx.mpi_master():
        start_time = ostream.print_start_header(glob_nodes)

    # initialize molecule and basis set

    molecule = vlx.Molecule()
    mol_basis = vlx.MolecularBasis()

    # read molecule and basis set

    if glob_rank == vlx.mpi_master():
        ostream.put_info("Reading input file %s..." % input_fname)

        input_dict = vlx.InputParser(input_fname).parse()

        ostream.put_info("Found %d control groups." % len(input_dict.keys()))
        ostream.put_info("...done.")
        ostream.new_line()

        ostream.put_info("Parsing @molecule group...")
        ostream.put_info("...done.")
        ostream.new_line()

        molecule = vlx.Molecule.from_xyz(input_dict["molecule"]["atom_labels"],
                                         input_dict["molecule"]["x_coords"],
                                         input_dict["molecule"]["y_coords"],
                                         input_dict["molecule"]["z_coords"])

        if "charge" in input_dict["molecule"].keys():
            molecule.set_charge(int(input_dict["molecule"]["charge"]))

        if "multiplicity" in input_dict["molecule"].keys():
            molecule.set_multiplicity(
                int(input_dict["molecule"]["multiplicity"]))

        molecule.check_multiplicity()
        molecule.check_proximity(0.1, ostream)
        molecule.print_geometry(ostream)

        ostream.put_info("Parsing @method settings group...")
        ostream.put_info("...done.")
        ostream.new_line()

        mol_basis = vlx.MolecularBasis.from_lib(
            input_dict["method settings"]["basis"],
            input_dict["method settings"]["basis path"], molecule, ostream)
        mol_basis.print_basis("Atomic Basis", molecule, ostream)

        ostream.flush()

    # broadcast molecule and basis set

    molecule.broadcast(glob_rank, comm)
    mol_basis.broadcast(glob_rank, comm)

    # initialize scf driver

    scf_drv = vlx.ScfRestrictedDriver()

    # read minimal basis if needed

    min_basis = vlx.MolecularBasis()

    if scf_drv.need_min_basis():
        if glob_rank == vlx.mpi_master():
            min_basis = vlx.MolecularBasis.from_lib(
                "MIN-CC-PVDZ", input_dict["method settings"]["basis path"],
                molecule, ostream)
        min_basis.broadcast(glob_rank, comm)

    scf_drv.compute(molecule, mol_basis, min_basis, comm, ostream)

    # write hdf5 files for MOs and density after SCF convergence

    if glob_rank == vlx.mpi_master():
        if scf_drv.is_converged:
            scf_drv.mol_orbs.write_hdf5("mol_orbs.h5")
            scf_drv.density.write_hdf5("density.h5")
        """
        # initialize visualization driver

        vis_drv = vlx.VisualizationDriver()

        # TODO: generate grid and write cube file
        # TODO: should also be able to compute psi/density along a line

        if scf_drv.is_converged:
            vis_grid = vis_drv.gen_grid(molecule)
            nelec = molecule.number_of_electrons()
            homo = nelec // 2 -1
            vis_drv.write_cube(molecule, mol_basis, scf_drv.mol_orbs, homo,
                               "alpha", vis_grid)
            vis_drv.write_cube_dens(molecule, mol_basis, scf_drv.density, 0,
                                    "alpha", vis_grid)
        """

    # all done, print finish header to output stream

    if glob_rank == vlx.mpi_master():
        ostream.print_finish_header(start_time)


if __name__ == "__main__":
    main()
