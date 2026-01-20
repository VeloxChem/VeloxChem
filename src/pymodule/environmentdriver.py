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
from .veloxchemlib import mpi_master
from .outputstream import OutputStream

class EnvironmentDriver:
    """
    Handles the environment from the snapshots of a trajectory.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables:
        - solvent_models: Dictionary containing solvent model parameters.
        - resname_to_model: Mapping from residue names to solvent models.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initialize the environment driver.
        See this reference, Figure 4, for different "cp3 approach" parameters:
        https://doi.org/10.1021/acs.jctc.5c01719
        """
        self.solvent_pe_models = {
            "SEP": {
                    "pattern": ["O", "H", "H"],
                    "charges": {"O": -0.67444000, "H": 0.33722000},
                    "polarizabilities": {
                        "O": [0.0, 0.0, 5.73935000, 0.0, 0.0, 0.0],
                        "H": [0.0, 0.0, 2.30839000, 0.0, 0.0, 0.0],
                    },
            }
        }

        self.resname_to_pe_model = {
            "SOL": "SEP",
            "WAT": "SEP",
            "HOH": "SEP",
        }

        self.solvent_npe_models = {
            "tip3p": {
                    "pattern": ["O", "H", "H"],
                    "charges": {"O": -0.83400000, "H": 0.41700000},
            }
        }

        self.resname_to_npe_model = {
            "SOL": "tip3p",
            "WAT": "tip3p",
            "HOH": "tip3p",
        }

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

    def write_pot_files(self, snapshots, outdir: str | Path):
        """
        Write environment snapshots to .pot files.

        Generates one .pot file per snapshot.

        :param snapshots:
            A list of snapshot dictionaries.
        :param outdir:
            Output directory where the .pot files will be written.
        :return:
            None.
        """
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        for snap in snapshots:
            self.write_single_pot(snap, outdir)

    def write_single_pot(self, snapshot, outdir: Path):
        """
        Write a single snapshot to a .pot file.

        The file contains the @environment,
        @charges, and @polarizabilities sections.

        :param snapshot:
            A snapshot dictionary containing:

            - frame (int):
                Frame number used in the output file name.
            - mm_coords (numpy.ndarray):
                MM region Cartesian coordinates, shape (N, 3), in Angstrom.
            - mm_elements (numpy.ndarray):
                Element symbols for each MM atom, length N.
            - mm_resids (numpy.ndarray):
                Residue id for each MM atom, length N.
            - mm_resnames (numpy.ndarray):
                Residue name for each MM atom, length N.

        :param outdir:
            Directory where the .pot file will be written.
        :return:
            None.
        """
        frame = snapshot['frame']
        mm_coords = snapshot['mm_coords']
        mm_elements = snapshot['mm_elements']
        mm_resids = snapshot['mm_resids']
        mm_resnames = snapshot['mm_resnames']
        present_resnames = sorted(set(str(r) for r in mm_resnames))

        model_to_resnames = {}
        for resn in present_resnames:
            model = self.resname_to_model.get(resn)
            if model is None:
                raise KeyError(f"No model registered for residue name '{resn}'")
            model_to_resnames.setdefault(model, set()).add(resn)

        pot_path = outdir / f"frame_{frame:06d}.pot"
        with pot_path.open('w') as fh:
            fh.write("@environment\n")
            fh.write("units: angstrom\n")
            fh.write("xyz:\n")
            for (x, y, z), elem, resn, resid in zip(mm_coords, mm_elements, mm_resnames, mm_resids):
                fh.write(f"{elem:<2} {x:12.6f} {y:12.6f} {z:12.6f}  {resn:>3}  {resid}\n")
            fh.write("@end\n\n")
            fh.write("@charges\n")

            for model, resnames in model_to_resnames.items():
                pattern = self.solvent_models[model]["pattern"]
                charges = self.solvent_models[model]["charges"]
                for resn in sorted(resnames):
                    for atom in pattern:
                        fh.write(f"{atom:<2} {charges[atom]:12.8f}  {resn}\n")
            fh.write("@end\n\n")
            fh.write("@polarizabilities\n")

            for model, resnames in model_to_resnames.items():
                pattern = self.solvent_models[model]["pattern"]
                polar = self.solvent_models[model]["polarizabilities"]
                for resn in sorted(resnames):
                    for atom in pattern:
                        vals = polar[atom]
                        fh.write(
                            f"{atom:<2} {vals[0]:12.8f} {vals[1]:12.8f} {vals[2]:12.8f} "
                            f"{vals[3]:12.8f} {vals[4]:12.8f} {vals[5]:12.8f}  {resn}\n"
                        )           
            fh.write("@end\n")

