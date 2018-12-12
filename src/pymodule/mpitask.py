from .veloxchemlib import Molecule
from .veloxchemlib import MolecularBasis
from .veloxchemlib import AtomBasis
from .veloxchemlib import BasisFunction
from .veloxchemlib import ChemicalElement
from .veloxchemlib import OutputStream
from .veloxchemlib import mpi_master
from .veloxchemlib import assert_msg_critical
from .veloxchemlib import to_angular_momentum
from .inputparser import InputParser


class MpiTask:

    def __init__(self, input_fname, output_fname, mpi_comm):

        # mpi settings

        self.mpi_comm = mpi_comm
        self.mpi_rank = mpi_comm.Get_rank()
        self.mpi_size = mpi_comm.Get_size()

        # initialize molecule, basis set and output stream

        self.molecule = Molecule()
        self.ao_basis = MolecularBasis()
        self.min_basis = MolecularBasis()
        self.ostream = OutputStream(output_fname)

        # process input file on master node

        if (self.mpi_rank == mpi_master()):

            self.start_time = self.ostream.print_start_header(self.mpi_size)

            self.ostream.put_info("Reading input file %s..." % input_fname)

            # read input file

            input_dict = InputParser(input_fname).parse()

            self.ostream.put_info(
                "Found %d control groups." % len(input_dict.keys()))
            self.ostream.put_info("...done.")
            self.ostream.new_line()

            # create molecule

            self.ostream.put_info("Parsing @self.molecule group...")
            self.ostream.put_info("...done.")
            self.ostream.new_line()

            self.molecule = Molecule.from_xyz(
                input_dict["molecule"]["atom_labels"],
                input_dict["molecule"]["x_coords"],
                input_dict["molecule"]["y_coords"],
                input_dict["molecule"]["z_coords"])

            if "charge" in input_dict["molecule"].keys():
                self.molecule.set_charge(int(input_dict["molecule"]["charge"]))

            if "multiplicity" in input_dict["molecule"].keys():
                self.molecule.set_multiplicity(
                    int(input_dict["molecule"]["multiplicity"]))

            self.molecule.check_multiplicity()
            self.molecule.check_proximity(0.1, self.ostream)
            self.molecule.print_geometry(self.ostream)

            # create basis set (new code)

            basis_path = input_dict["method_settings"]["basis_path"]
            basis_label = input_dict["method_settings"]["basis"].upper()
            basis_fname = basis_path + '/' + basis_label
            basis_dict = InputParser(basis_fname).parse()

            assert_msg_critical(
                basis_label == basis_dict['basis_set_name'].upper(),
                "basis set name")

            mol_basis = MolecularBasis()

            elem_comp = self.molecule.get_elemental_composition()

            for elem_id in elem_comp:

                elem = ChemicalElement()
                err = elem.set_atom_type(elem_id)
                assert_msg_critical(err, "ChemicalElement.set_atom_type")

                basis_key = 'atombasis_%s' % elem.get_name().lower()
                basis_list = [entry for entry in basis_dict[basis_key]]

                atom_basis = AtomBasis()

                while basis_list:
                    shell_title = basis_list.pop(0).split()
                    assert_msg_critical(
                        len(shell_title) == 3,
                        "Basis set parser (shell): %s" % ' '.join(shell_title))

                    angl = to_angular_momentum(shell_title[0])
                    npgto = int(shell_title[1])
                    ncgto = int(shell_title[2])

                    expons = [0.0] * npgto
                    coeffs = [0.0] * npgto * ncgto

                    for i in range(npgto):
                        prims = basis_list.pop(0).split()
                        assert_msg_critical(
                            len(prims) == ncgto + 1,
                            "Basis set parser (primitive): %s" %
                            ' '.join(prims))

                        expons[i] = float(prims[0])
                        for k in range(ncgto):
                            coeffs[k * npgto + i] = float(prims[k + 1])

                    bf = BasisFunction.from_list(expons, coeffs, ncgto, angl)
                    bf.normalize()

                    atom_basis.add_basis_function(bf)

                atom_basis.set_elemental_id(elem_id)

                mol_basis.add_atom_basis(atom_basis)

            mol_basis.set_label(basis_label)

            mol_basis.print_basis("Atomic Basis", self.molecule, self.ostream)

            # create basis set (old code)

            self.ostream.put_info("Parsing @method settings group...")
            self.ostream.put_info("...done.")
            self.ostream.new_line()

            self.ao_basis = MolecularBasis.from_lib(
                input_dict["method_settings"]["basis"],
                input_dict["method_settings"]["basis_path"], self.molecule,
                self.ostream)
            self.ao_basis.print_basis("Atomic Basis", self.molecule,
                                      self.ostream)

            assert (self.ao_basis == mol_basis)

            self.min_basis = MolecularBasis.from_lib(
                "MIN-CC-PVDZ", input_dict["method_settings"]["basis_path"],
                self.molecule, self.ostream)

            self.ostream.flush()

        # broadcast molecule and basis set

        self.molecule.broadcast(self.mpi_rank, self.mpi_comm)
        self.ao_basis.broadcast(self.mpi_rank, self.mpi_comm)
        self.min_basis.broadcast(self.mpi_rank, self.mpi_comm)

    def finish(self):

        if (self.mpi_rank == mpi_master()):
            self.ostream.print_finish_header(self.start_time)
