#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from pathlib import Path
from mpi4py import MPI
import sys

from .veloxchemlib import mpi_master, bohr_in_angstrom
from .molecule import Molecule
from .scfrestdriver import ScfRestrictedDriver
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical


class PEDriver:
    """
    Driver for performing polarizable embedding (PE) calculations.
    
    This class automates the creation of PE potentials from solvated systems
    and executes PE calculations.

    : param comm:
        The MPI communicator.
    : param ostream:
        The output stream.

    Instance variables:
        - qm_molecule: The QM molecule.
        - solvated system: The complete solvates system from SolvationBuilder.
        - pe_parameters: Dictionary containing charges and polarizabilities.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initialize the PE driver.
        """
        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()

        # output stream
        self.ostream = ostream
      
        self.qm_molecule = None
        self.solvated_system = None
        self.pe_parameters = self._initialize_pe_parameters()

    def _initialize_pe_parameters(self):
        """
        Initialize parameters commonly used in PE.

        These parameters have been benchmarked and validated for accuracy,
        and can be divided into several categories:
          - SEP for water and ions of solvation
          - CP3 for amino acids in proteins
          - ALEP for lipids in a lipid bilayer

        More information here: 
        doi:10.26434/chemrxiv-2025-vzz1r

        : return: 
             A dictionary containing the charges and polarizabilities.

        Parameters list to be progressively populated.
        """
        parameters ={
            "water": {
                "H": {
                    "charge": 0.33722206,
                    "polarizability": [2.30839051, 0.0, 0.0, 2.30839051, 0.0, 2.30839051]
                },
                "O": {
                    "charge": -0.67444408,
                    "polarizability": [5.73935090, 0.0, 0.0, 5.73935090, 0.0, 5.73935090]
                }
            }
        }

        # SolvationBuilder works with 'spce' or 'tip3p'
        # The SEP parameters for water are the same regardless of 'spce' or 'tip3p'
        parameters["spce"] = parameters["water"]
        parameters["tip3p"] = parameters["water"]

        return parameters



    def write_pe_potential(self, solvation_builder, qm_molecule, potfile='pe.pot',
                           solvent_name="water", residue_name="water"):
        """
        Create a PE potential from a solvated system.

        :param solvation_builder:
            The SolvationBuilder object containing the solvated system.
        :param qm_molecule:
            The QM molecule (typically the solute).
        :param potfile:
            The name of the output potential file.
        :param solvent_name:
            The name of the solvent in the parameters dictionary.
        :param residue_name:
            The name of the residue to use in the PE potential.

        : return:
            The path to the created PE potential.
        """
        self.qm_molecule = qm_molecule
        self.solvated_system = solvation_builder.system_molecule

        # Get the solvent name from the solvation builder if not provided
        if hasattr(solvation_builder, 'solvent_name') and solvation_builder.solvent_name:
            solvent_name = solvation_builder.solvent_name

        # Get the residue name from the solvation builder if not provided
        if hasattr(solvation_builder, 'residue_name') and solvation_builder.residue_name:
            residue_name = solvation_builder.residue_name

        # Check if we have parameters for the specified solvent
        # if solvent_name not in self.pe_parameters:
        #     raise ValueError(f"No PE parameters found for solvent: {solvent_name}. "
        #                      f"Available solvents: {list(self.pe_parameters.keys())}")
        
        # Extract MM enviroment atoms (all atoms except QM region)
        qm_natoms = qm_molecule.number_of_atoms()
        total_natoms = self.solvated_system.number_of_atoms()

        assert_msg_critical(
            total_natoms > qm_natoms,
            "PEDriver: The solvated system must contain more atoms than the QM molecule."
        )

        # Get coordinates and labels for all atoms
        all_coords = self.solvated_system.get_coordinates_in_angstrom()
        all_labels = self.solvated_system.get_labels()

        # MM region: atoms after the QM region
        mm_coords = all_coords[qm_natoms:]
        mm_labels = all_labels[qm_natoms:]

        # Get number of atoms per solvent molecule:
        atoms_per_molecule = len(solvation_builder.solvent_labels[0])

        # Write PE potential
        self._write_pe_file(mm_coords, mm_labels, potfile, solvent_name, residue_name, atoms_per_molecule)

        return potfile

    def _write_pe_file(self, mm_coords, mm_labels, potfile, solvent_name, residue_name, atoms_per_molecule):
        """
        Write the PE potential.

        :param mm_coords:
            Coordinates of the MM atoms in Angstrom.
        :param mm_labels:
            Element labels of the MM atoms.
        :param potfile:
            The name of the output potential file.
        :param solvent_name:
            The name of the solvent in the parameters dictionary.
        :param residue_name:
            The name of the residue to use in the PE potential.
        """
        params = self.pe_parameters[solvent_name]

        with open(potfile, 'w') as f:
            # Write environment section
            f.write("@environment\n")
            f.write("units: angstrom\n")
            f.write("xyz:\n")

            # Group atoms by molecules
            n_molecules = len(mm_coords) // atoms_per_molecule

            for mol_id in range(n_molecules):
                start_idx = mol_id * atoms_per_molecule
                end_idx = start_idx + atoms_per_molecule

                for i in range(start_idx, end_idx):
                    label = mm_labels[i]
                    x, y, z  = mm_coords[i]
                    f.write(f"{label}  {x:12.7f}  {y:12.7f}  {z:12.7f}  {residue_name}  {mol_id + 1}\n")

            f.write("@end\n")

            # Write charges section
            f.write("@charges\n")
            for i in range(atoms_per_molecule):
                label = mm_labels[i] # This will be O, H, H for water
                if label in params:
                    charge = params[label]['charge']
                    f.write(f"{label} {charge:11.8f}  {residue_name}\n")
                else:
                    raise ValueError(f"No charge parameter found for element: {label}")
            f.write("@end\n\n")

            # Write polarizabilities section
            f.write("@polarizabilities\n")
            for i in range(atoms_per_molecule):
                label = mm_labels[i]
                if label in params:
                    pols = params[label]['polarizability']
                    pols_str = "  ".join([f"{p:12.8f}" for p in pols])
                    f.write(f"{label}   {pols_str}  {residue_name}\n")
                else:
                    raise ValueError(f"No polarizability parameter found for element: {label}")
            f.write("@end\n")

            if self.rank == mpi_master():
                self.ostream.print_info(f"PE potential written to {potfile}")
                self.ostream.print_info(f"Number of MM atoms: {len(mm_coords)}")
                self.ostream.print_info(f"Number of molecules: {n_molecules}")
                self.ostream.flush()
        
    def compute(self, qm_molecule, basis, solvation_builder=None,
                potfile='pe.pot', solvent_name="water", **scf_kwargs):
        """
        Perform a PE calculation.

        :param qm_molecule:
            The QM molecule.
        :param basis:
            The basis set for the QM calculation.
        :param solvation_builder:
            Optional SolvationBuilder object. If provided, a PE potential file
            will be created automatically from the solvated system.
        :param potfile:
            The name of the PE potential file. If solvation_builder is provided,
            this file will be created automatically. Otherwise, it should point to an
            existing PE potential file.
        :param solvent_name:
            The name of the solvent in the parameters dictionary.
        :param scf_kwargs:
            Additional keyword arguments for the SCF calculation.

        : return:
            SCF results dictionary.
        """
        # Create PE potential if solvation builder is provided
        if solvation_builder is not None:
            potfile = self.write_pe_potential(
                solvation_builder, qm_molecule, potfile, solvent_name
            )
        else:
            # Check if potfile exists
            if not Path(potfile).exists():
                raise FileNotFoundError(
                    f"PE potential file not found: {potfile}."
                    "Provide a solvation_builder or ensure the potfile exists."
                    )
        # Set up and run SCF calculation with PE
        scf_drv = ScfRestrictedDriver(comm=self.comm, ostream=self.ostream)
        scf_drv.potfile = potfile

        # Apply any additional SCF options
        for key, value in scf_kwargs.items():
            setattr(scf_drv, key, value)

        if self.rank == mpi_master():
            self.ostream.print_blank()
            self.ostream.print_info("Starting PE calculation...")
            #self.ostream.print_info(f"Using PE potential file: {potfile}")
            self.ostream.flush()

        # Run the calculation
        scf_results = scf_drv.compute(qm_molecule, basis)

        return scf_results