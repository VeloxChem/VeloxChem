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
import pyframe
import numpy as np
from .molecularbasis import MolecularBasis
from .outputstream import OutputStream

class SpectrumDriver:
    """
    Driver for the sppectra in polarizable embedding (PE) calculations.
    
    This class automates the creation of PE potentials from solvated systems
    and executes PE calculations.

    : param comm:
        The MPI communicator.
    : param ostream:
        The output stream.

    Instance variables:
        - qm_molecule: The QM molecule.
        - solvated system: The complete solvates system from SolvationBuilder.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initialize the SpectrumDriver.
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

    @staticmethod
    def trajectory_parser(topology_file, trajectory_file, 
                          qm_region,  
                          num_snapshots):
        """"
        Trajectory parser for managing snapshots with pyframe.
        """
        
        traj = pyframe.Trajectory(topology_file=topology_file, 
                                  trajectory_file=trajectory_file,
                                  center_selection=qm_region,
                                  num_snapshots=num_snapshots)
        traj.set_core_region(qm_region)
        traj.add_region(
            name='water and/or ions',
            selection='resname SOL or resname Na+',
            use_standard_potentials=True,
            standard_potential_model='SEP'
        )
        
        traj.create_potential()

        return traj
        # traj.write_core()
        # traj.write_potential()

    def write_potential(self, trajectory=None, filename_pot: str = None) -> None:
        """Write environment potential files in custom format for each snapshot.
        
        Format includes coordinates, charges, and polarizabilities organized by fragment type.
        """
        if trajectory is None:
            raise ValueError("Trajectory must be provided to write potential files.")
        
        for frame_idx, snapshot in zip(trajectory._frame_indices, trajectory.snapshots):
            if filename_pot is None:
                base_filename = trajectory.name
            else:
                base_filename = filename_pot
            
            output_filename = f'{base_filename}_frame-{frame_idx}'
            
            print(f"==== Processing frame {frame_idx} ====")
            print(f"Number of MM regions: {len(snapshot.regions)}")
            print(f"Number of MM sites: {len(snapshot.potential)}")
            
            # Collect data organized by fragment type
            fragments_data = {}
            potential_sites = list(snapshot.potential.values())
            site_idx = 0

            def _site_pol_tensor(site):
                pol = np.zeros((3, 3))
                if hasattr(site, 'polarizability') and site.polarizability is not None:
                    p = np.array(site.polarizability)
                    if p.shape == (3, 3):
                        return p
                if hasattr(site, 'P') and site.P is not None:
                    p = np.array(site.P)
                    if p.shape == (3, 3):
                        return p
                if hasattr(site, 'P00') and site.P00:
                    pol[0, 0] = site.P00[0]
                if hasattr(site, 'P10') and site.P10:
                    pol[0, 1] = pol[1, 0] = site.P10[0]
                if hasattr(site, 'P11') and site.P11:
                    pol[0, 2] = pol[2, 0] = site.P11[0]
                if hasattr(site, 'P20') and site.P20:
                    pol[1, 1] = site.P20[0]
                if hasattr(site, 'P21') and site.P21:
                    pol[1, 2] = pol[2, 1] = site.P21[0]
                if hasattr(site, 'P22') and site.P22:
                    pol[2, 2] = site.P22[0]
                return pol

            # First pass: organize data by fragment
            for region_name, region in snapshot.regions.items():
                for fragment in region.fragments.values():
                    frag_name = fragment.name
                    if frag_name not in fragments_data:
                        fragments_data[frag_name] = {
                            'atoms': [],
                            'elements': [],
                            'coords': [],
                            'charges': [],
                            'polarizabilities': []
                        }

                    for atom in fragment.atoms:
                        site = potential_sites[site_idx]
                        site_idx += 1

                        fragments_data[frag_name]['atoms'].append(atom.name)      # use atom name
                        fragments_data[frag_name]['elements'].append(atom.element)
                        fragments_data[frag_name]['coords'].append(site.coordinate)

                        charge = site.M0[0] if site.M0 else 0.0
                        fragments_data[frag_name]['charges'].append(charge)

                        pol_tensor = _site_pol_tensor(site)
                        pol6 = [
                            pol_tensor[0, 0],  # xx
                            pol_tensor[0, 1],  # xy
                            pol_tensor[0, 2],  # xz
                            pol_tensor[1, 1],  # yy
                            pol_tensor[1, 2],  # yz
                            pol_tensor[2, 2],  # zz
                        ]
                        fragments_data[frag_name]['polarizabilities'].append(pol6)

            
            # Write output file
            with open(f'{output_filename}.pot', 'w') as f:
                # Write coordinates section
                f.write('@environment\n')
                f.write('units: angstrom\n')
                f.write('xyz:\n')
                
                fragment_counter = {}
                for frag_name, data in fragments_data.items():
                    if frag_name not in fragment_counter:
                        fragment_counter[frag_name] = 0

                    # Use atom names to calculate atoms per fragment
                    atoms_per_fragment = len(set(data['atoms']))
                    
                    for i, (element, coord) in enumerate(zip(data['elements'], data['coords'])):
                        if i % atoms_per_fragment == 0:  # New fragment every N atoms
                            fragment_counter[frag_name] += 1
                        
                        frag_num = fragment_counter[frag_name]                    
                        
                        f.write(f'{element:2} {coord[0]:14.7f} {coord[1]:14.7f} {coord[2]:14.7f}  '
                               f'{frag_name}  {frag_num}\n')                
                f.write('@end\n\n')
                
                # Write charges section (one line per atom name; keeps both H entries)
                f.write('@charges\n')
                for frag_name, data in fragments_data.items():
                    seen = set()
                    for atom_name, elem, q in zip(data['atoms'], data['elements'], data['charges']):
                        if atom_name in seen:
                            continue
                        seen.add(atom_name)
                        f.write(f'{elem:2} {q:12.8f}  {frag_name}\n')
                f.write('@end\n\n')

                # Write polarizabilities section (one line per atom name; keeps both H entries)
                f.write('@polarizabilities\n')
                for frag_name, data in fragments_data.items():
                    seen = set()
                    for atom_name, elem, pol in zip(data['atoms'], data['elements'], data['polarizabilities']):
                        if atom_name in seen:
                            continue
                        seen.add(atom_name)
                        pol_str = '    '.join(f'{p:12.8f}' for p in pol)
                        f.write(f'{elem:2} {pol_str}  {frag_name}\n')
                f.write('@end\n')


    def _snapshot_qm_region_to_molecule(self, snapshot):
        """
        Build a Molecule from the snapshot QM region atoms.
        """
        assert_msg_critical(
            hasattr(snapshot, "core_region") and hasattr(snapshot.core_region, "fragments"),
            "SpectrumDriver.compute: snapshot.core_region.fragments not found."
        )

        atoms = []
        for frag in snapshot.core_region.fragments.values():
            atoms.extend(frag.atoms)

        assert_msg_critical(
            len(atoms) > 0,
            "SpectrumDriver.compute: QM atom list is empty."
        )

        lines = []
        for atom in atoms:
            x, y, z = atom.coordinate  # Angstrom
            lines.append(f"{atom.element} {x:16.8f} {y:16.8f} {z:16.8f}")
        mol_str = "\n".join(lines)
        return Molecule.read_molecule_string(mol_str, units="angstrom")

    def compute(self,
                trajectory,
                basis,
                functional,
                filename_pot=None,
                basis_path=".",
                workdir="."):
        """
        Run SCF for all snapshots in a trajectory using PE pot files.

        Args:
            trajectory: pyframe.Trajectory returned by trajectory_parser().
            basis: e.g. 'def2-svp'.
            functional: e.g. 'CAM-B3LYP' or 'HF'.
            filename_pot: base name used when writing .pot (defaults to trajectory.name).
            basis_path: path to basis files (if not in default).
            workdir: directory where .pot files are located.

        Returns:
            List of dicts with keys: frame, molecule, scf_results.
        """
        results = []
        base_filename = trajectory.name if filename_pot is None else filename_pot

        for frame_idx, snapshot in zip(trajectory._frame_indices, trajectory.snapshots):
            potfile = Path(workdir) / f"{base_filename}_frame-{frame_idx}.pot"
            assert_msg_critical(
                potfile.is_file(),
                f"SpectrumDriver.compute: pot file not found: {potfile}"
            )

            # Build QM molecule from the QM region
            mol = self._snapshot_qm_region_to_molecule(snapshot)

            # Basis
            mol_basis = MolecularBasis.read(mol, basis, basis_path=basis_path, ostream=None)

            # SCF driver
            scf_drv = ScfRestrictedDriver(self.comm, self.ostream)
            scf_drv.restart = False
            scf_drv.filename = f"{base_filename}_frame-{frame_idx}"
            scf_drv.potfile = str(potfile)
            scf_drv.xcfun = functional

            scf_res = scf_drv.compute(mol, mol_basis)
            results.append({
                "frame": frame_idx,
                "molecule": mol,
                "scf_results": scf_res,
            })

        return results