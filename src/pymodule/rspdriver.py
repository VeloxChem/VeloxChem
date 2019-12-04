from .cppsolver import ComplexResponse
from .lrsolver import LinearResponseSolver
from .lreigensolver import LinearResponseEigenSolver
from .c6solver import C6Solver
from .tdaexcidriver import TDAExciDriver


class ResponseDriver:
    """
    Implements response driver for molecular property calculations using
    conventional Hartree-Fock/Kohn-Sham response theory.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - prop_type: The type of the property to be calculated.
        - tamm_dancoff: The flag for using Tamm-Dancoff approximation.
        - rsp_dict: The dictionary of response input.
        - method_dict: The dictionary of method settings.
        - comm: The MPI communicator.
        - rank: The rank of MPI process.
        - nodes: The number of MPI processes.
    """

    def __init__(self, comm, ostream):
        """
        Initializes response driver to default setup.
        """

        # default calculation type
        self.prop_type = 'generic'
        self.tamm_dancoff = False
        self.rsp_dict = {}
        self.method_dict = {}

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream

    def update_settings(self, rsp_dict, method_dict={}):
        """
        Updates settings in response solver.

        :param rsp_dict:
            The dictionary of response input.
        :param method_dict:
            The dictionary of method settings.
        """

        self.rsp_dict = dict(rsp_dict)
        self.method_dict = dict(method_dict)

        if 'property' in rsp_dict:
            self.prop_type = rsp_dict['property'].lower()

        if 'tamm_dancoff' in rsp_dict:
            key = rsp_dict['tamm_dancoff'].lower()
            self.tamm_dancoff = True if key in ['yes', 'y'] else False

    def compute(self, molecule, ao_basis, scf_tensors):
        """
        Performs molecular property calculation using molecular data

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.

        :return:
            The results from the actual response solver.
        """

        # Linear response eigensolver

        if (self.rsp_dict['response'] == 'linear' and
                self.rsp_dict['residue'] == 'single' and
                self.rsp_dict['complex'] == 'no'):

            if self.tamm_dancoff:
                solver = TDAExciDriver(self.comm, self.ostream)
            else:
                solver = LinearResponseEigenSolver(self.comm, self.ostream)

            solver.update_settings(self.rsp_dict, self.method_dict)

            return solver.compute(molecule, ao_basis, scf_tensors)

        # Linear response solver

        if (self.rsp_dict['response'] == 'linear' and
                self.rsp_dict['residue'] == 'none' and
                self.rsp_dict['complex'] == 'no'):

            solver = LinearResponseSolver(self.comm, self.ostream)

            solver.update_settings(self.rsp_dict, self.method_dict)

            return solver.compute(molecule, ao_basis, scf_tensors)

        # Complex linear response solver

        if (self.rsp_dict['response'] == 'linear' and
                self.rsp_dict['residue'] == 'none' and
                self.rsp_dict['onlystatic'] == 'no' and
                self.rsp_dict['complex'] == 'yes'):

            clr_solver = ComplexResponse(self.comm, self.ostream)

            clr_solver.update_settings(self.rsp_dict, self.method_dict)

            return clr_solver.compute(molecule, ao_basis, scf_tensors)

        # C6 linear response solver

        if (self.rsp_dict['response'] == 'linear' and
                self.rsp_dict['residue'] == 'none' and
                self.rsp_dict['onlystatic'] == 'yes' and
                self.rsp_dict['complex'] == 'yes'):

            c6_solver = C6Solver(self.comm, self.ostream)

            c6_solver.update_settings(self.rsp_dict, self.method_dict)

            return c6_solver.compute(molecule, ao_basis, scf_tensors)

    def prop_str(self):
        """
        Gets string with type of molecular property calculation (Excited
        states, linear and non-linear spectroscopies).

        :return:
            The string with type of molecular property calculation.
        """

        if self.prop_type == 'polarizability':
            return 'Polarizability'

        if self.prop_type == 'absorption':
            return 'Singlet Excited States'

        if self.prop_type == 'linear absorption cross-section':
            return 'Linear Absorption Cross-Section'

        if self.prop_type == 'circular dichroism spectrum':
            return 'Circular Dichroism Spectrum'

        if self.prop_type == 'c6':
            return 'C6 values'

        return 'Undefined'
