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
import numpy as np
from .outputstream import OutputStream
from .molecule import Molecule
from .molecularbasis import MolecularBasis

from .veloxchemlib import (mpi_master, bohr_in_angstrom)

class EnvironmentDriver:
    """
    Handles the environment from the snapshots from EnsembleParser.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables:
        - solvent_pe_models: Dictionary containing solvent PE model parameters.
        - solvent_npe_models: Dictionary containing solvent NPE model parameters.
        - pe_model: PE model name.
        - npe_model: NPE model name.
    """

    # def __init__(self, comm=None, ostream=None, *, scf_drv=None, rsp_drv=None):
    def __init__(self, comm=None, ostream=None):
        """
        Initialize the environment driver.
        See this reference, Figure 4, for a summary of 
        different "cp3 approach" parameters:
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

        self.solvent_npe_models = {
            "tip3p": {
                    "pattern": ["O", "H", "H"],
                    "charges": {"O": -0.83400000, "H": 0.41700000},
            }
        }
        # set them to None:
        self.pe_model = None
        self.npe_model = None
        self.set_env_models()

        # fetch them from .csv file 

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()
        self.ostream = ostream

    def set_env_models(self,
                       pe_model="SEP",
                       npe_model="tip3p"):
        """
        Set the parameter PE and NPE models

        :param pe_model:
            Name of the PE parameter model.
        :param npe_model:
            Name of the NPE parameter model.
        :return:
            None.
        :raises KeyError:
            If an unknown model name is requested.
        """
        if pe_model not in self.solvent_pe_models:
            raise KeyError(
                f"Unknown PE model '{pe_model}'. Available: {sorted(self.solvent_pe_models)}"
            )
        if npe_model not in self.solvent_npe_models:
            raise KeyError(
                f"Unknown NPE model '{npe_model}'. Available: {sorted(self.solvent_npe_models)}"
            )
        
        self.pe_model = self.solvent_pe_models[pe_model]
        self.npe_model = self.solvent_npe_models[npe_model]
        return self.pe_model, self.npe_model

    def _build_point_charges(self, coords_ang, elements):
        """
        Build the point charges array expected by ScfDriver.
    
        : param coords_ang:
            Cartesian coordinates in Angstrom
        : param elements:
            List of element symbols
        :return:
            pe (numpy.ndarray):
            point charges array, shape (6, N).
        :raises KeyError:
            If an elements has no charge in the selected NPE model.
        """
        coords_ang = np.asarray(coords_ang, dtype=float)
        if coords_ang.size == 0:
            return None
    
        elements = np.asarray(elements, dtype=object)
        charges_map = self.npe_model["charges"]

        charges = []
        for elem in elements:
            key = str(elem)
            if key not in charges_map:
                available = ", ".join(sorted(charges_map))
                raise KeyError(
                    f"NPE model missing charge for element '{key}'. "
                    f"Available keys: {available}"
                )
            charges.append(charges_map[key])
        q = np.array(charges, dtype=float)

        pc = np.zeros((6, coords_ang.shape[0]), dtype=float)
        pc[0:3, :] = coords_ang.T / bohr_in_angstrom()
        pc[3, :] = q
        return pc
        
    def write_pot_files(self,
                        snapshots, 
                        outdir: str | Path):
        """
        Write pe environment snapshots to .pot files.

        Generates one .pot file per snapshot.
        The file contains the @environment,
        @charges, and @polarizabilities sections.

        :param snapshots:
            A list of snapshot dictionaries returned from TrajectoryDriver.
        :param outdir:
            Output directory where the .pot files will be written.
        :param pe_model:
            PE model name
        :return:
            None.
        """
        # Normalize input (single snapshot or list)
        if isinstance(snapshots, dict):
            snapshots = [snapshots]

        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)

        pattern = self.pe_model["pattern"]
        charges = self.pe_model["charges"]
        polar = self.pe_model["polarizabilities"]

        for snap in snapshots:
            frame = int(snap["frame"])

            pe_coords = np.asarray(snap["pe_coords"], dtype=float)
            pe_elements = np.asarray(snap["pe_elements"], dtype=object)
            pe_resids = np.asarray(snap["pe_resids"], dtype=int)
            pe_resnames = np.asarray(snap["pe_resnames"], dtype=object)

            if pe_coords.size == 0:
                continue

            resname_set = sorted({str(r) for r in pe_resnames}) if pe_resnames.size else []

            pe_pot_path = outdir / f"pe_frame_{frame:06d}.pot"
            with pe_pot_path.open("w") as fh:

                fh.write("@environment\n")
                fh.write("units: angstrom\n")
                fh.write("xyz:\n")
                for (x, y, z), elem, resn, resid in zip(pe_coords, pe_elements, pe_resnames, pe_resids):
                    fh.write(
                        f"{str(elem):<2} {x:12.6f} {y:12.6f} {z:12.6f}  {str(resn):>3}  {int(resid)}\n"
                    )
                fh.write("@end\n\n")

                fh.write("@charges\n")
                for resn in resname_set:
                    for atom in pattern:
                        fh.write(f"{atom:<2} {charges[atom]:12.8f}  {resn}\n")
                fh.write("@end\n\n")

                fh.write("@polarizabilities\n")
                for resn in resname_set:
                    for atom in pattern:
                        vals = polar[atom]
                        fh.write(
                            f"{atom:<2} {vals[0]:12.8f} {vals[1]:12.8f} {vals[2]:12.8f} "
                            f"{vals[3]:12.8f} {vals[4]:12.8f} {vals[5]:12.8f}  {resn}\n"
                        )
                fh.write("@end\n")

    def compute(
        self,
        snapshots,
        *,
        basis_label: str,
        scf_drv,
        rsp_drv,
        potdir: str | Path = "pot_frames",
        write_pe_potfiles: bool = True,
    ):
        """
        Drives the computation over the ensemble of snapshots.

        :param snapshots:
            A list of snapshot dictionaries (or a single dict).
        :param basis_label: (str)
            Basis set label.
        :param scf_drv:
            An initialized SCF driver instance.
        :param rsp_drv:
            An initialized LR eigen solver.
        :param potdir : (str or Path)
            Directory to store/read PE potfiles.
        :param write_pe_potfiles: (bool)
            If True, PE potfiles are (re)generated before the loop when needed.
        :return:
            Dictionary with:
              - scf_all: list of (frame, scf_results)
              - rsp_all: list of (frame, rsp_results)
        """
        if isinstance(snapshots, dict):
            snapshots = [snapshots]

        potdir = Path(potdir)

        # Write PE potfiles only if requested and if any snapshot has PE atoms
        if write_pe_potfiles:
            has_any_pe = False
            for snap in snapshots:
                pe_coords = snap.get("pe_coords", [])
                if np.asarray(pe_coords).size > 0:
                    has_any_pe = True
                    break
        if has_any_pe:
            self.write_pot_files(snapshots, outdir=potdir)

        scf_all = []
        rsp_all = []

        for snap in snapshots:
            frame = int(snap["frame"])

            labels = [str(x) for x in snap["qm_elements"]]
            coords = np.asarray(snap["qm_coords"], dtype=float)
            molecule = Molecule(labels, coords)
            basis = MolecularBasis.read(molecule, basis_label)

            pe_coords = np.asarray(snap.get("pe_coords", []), dtype=float)
            npe_coords = np.asarray(snap.get("npe_coords", []), dtype=float)

            has_pe = pe_coords.size > 0
            has_npe = npe_coords.size > 0
            
            if has_pe:
                scf_drv.potfile = str(potdir / f"pe_frame_{frame:06d}.pot")

            if has_npe:
                scf_drv.point_charges = self._build_point_charges(
                    snap.get("npe_coords", []), 
                    snap.get("npe_elements", []),
                    )

            scf_results = scf_drv.compute(molecule, basis)
            rsp_results = rsp_drv.compute(molecule, basis, scf_results)

            scf_all.append((frame, scf_results))
            rsp_all.append((frame, rsp_results))

        return {"scf_all": scf_all, "rsp_all": rsp_all}
