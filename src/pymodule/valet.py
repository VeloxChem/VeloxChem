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

import numpy as np
import sys

from .molecule import Molecule
from .molecularbasis import MolecularBasis
from .visualizationdriver import VisualizationDriver
from .densityviewer import DensityViewer
from .errorhandler import assert_msg_critical

try:
    import matplotlib.pyplot as plt
    import matplotlib.path as mpath
    import matplotlib.patches as mpatches
except ImportError:
    pass

try:
    import py3Dmol
except ImportError:
    pass


class ValetAnalyzer:
    """
    Visual Analysis of Electronic Transitionss (VALET).
    
    Implements NTO (Natural Transition Orbital) analysis, transition matrix
    calculation, and interactive visualization of excited state properties.
    Based on the VALET method by Talha Bin Masood et al. [https://doi.org/10.1111/cgf.14307]
    
    Instance variables:
        - density_isovalue: Isovalue for density isosurfaces (default: 0.002)
        - density_opacity: Opacity of density isosurfaces (default: 0.7)
    """

    def __init__(self):
        """
        Initialize VALET analyzer.
        """
        self.density_isovalue = 0.002
        self.density_opacity = 0.7
        self.detachment_color = "#FFFFFF"  # White
        self.attachment_color = "#FFFFFF"  # White
        
        # Results and molecular data
        self._scf_results = None
        self._rsp_results = None
        self._molecule = None
        self._basis = None
        self._num_states = None
        self._subgroups = None
        
        # Density viewer for computing densities on grid
        self._density_viewer = None
        
        # Interactive viewer state
        self._mol_grid_viewer = None
    
    def initialize(self, scf_results, rsp_results, subgroups=None):
        """
        Initialize VALET analyzer with SCF and response results.
        Extracts and stores molecule, basis, and number of states.
        
        :param scf_results:
            The dictionary of results from converged SCF wavefunction.
        :param rsp_results:
            The dictionary containing the rsp results.
        :param subgroups:
            Optional list of subgroup tuples for transition diagrams.
            Each tuple contains (name, atom_index or list_of_atom_indices).
        """
        self._scf_results = scf_results
        self._rsp_results = rsp_results
        self._molecule = self._extract_molecule_from_scf(scf_results)
        self._basis = self._extract_basis_from_scf(scf_results, self._molecule)
        self._num_states = len(rsp_results['eigenvalues'])
        self._subgroups = subgroups
        
        # Initialize density viewer for computing densities
        self._density_viewer = DensityViewer()
        self._density_viewer.initialize(self._molecule, self._basis)
        
        print(f"VALET initialized with {self._num_states} excited states")
        print(f"Molecule: {self._molecule.number_of_atoms()} atoms")
        print(f"Basis: {self._basis.get_label()}")
    
    def set_subgroups(self, subgroups):
        """
        Set or update subgroups for transition diagrams.
        
        :param subgroups:
            List of subgroup tuples for transition diagrams.
            Each tuple contains (name, atom_index or list_of_atom_indices).
        """
        self._subgroups = subgroups
        print(f"Subgroups updated: {[name for name, _ in subgroups]}")

    @staticmethod
    def _extract_molecule_from_scf(scf_results):
        """
        Extracts molecule from SCF results.
        
        :param scf_results:
            The dictionary of results from converged SCF wavefunction.
            
        :return:
            The Molecule object.
        """
        coordinates = scf_results['atom_coordinates']
        nuclear_charges = np.array(scf_results['nuclear_charges'])
        nuclear_charges = nuclear_charges.astype(int)
        molecule = Molecule(nuclear_charges, coordinates, units="au")
        return molecule

    @staticmethod
    def _extract_basis_from_scf(scf_results, molecule):
        """
        Extracts basis from SCF results.
        
        :param scf_results:
            The dictionary of results from converged SCF wavefunction.
        :param molecule:
            The Molecule object.
            
        :return:
            The MolecularBasis object.
        """
        basis_set_label = scf_results['basis_set'][0].decode("utf-8")
        basis = MolecularBasis.read(molecule, basis_set_label)
        return basis

    def compute_nto_densities(self, scf_results, rsp_results, state_index=1):
        """
        Computes NTO attachment and detachment densities for a given excited state.
        
        :param scf_results:
            The dictionary of results from converged SCF wavefunction.
        :param rsp_results:
            The dictionary containing the rsp results.
        :param state_index:
            The excited state for which the NTO densities are computed (1-based).
            
        :return:
            A dictionary containing detachment and attachment density matrices.
        """
        
        # Extract molecule from SCF results
        molecule = self._extract_molecule_from_scf(scf_results)
        
        # Get MO coefficients and orbital information
        mo = scf_results["C_alpha"]
        nocc = molecule.number_of_alpha_electrons()
        norb = mo.shape[1]
        nvirt = norb - nocc
        mo_occ = mo[:, :nocc].copy()
        mo_vir = mo[:, nocc:].copy()

        # Get eigenvector for the specified state
        excstate_ident = "S" + str(state_index)
        
        if excstate_ident not in rsp_results:
            raise ValueError(f"State {excstate_ident} not found in rsp_results. "
                           f"Available keys: {list(rsp_results.keys())}")
        
        eigvec = rsp_results[excstate_ident]

        # Construct transition density matrix
        nexc = nocc * nvirt
        z_mat = eigvec[:nexc]
        if eigvec.shape[0] == nexc:
            tdens_mo = np.reshape(z_mat, (nocc, nvirt))
        else:
            y_mat = eigvec[nexc:]
            tdens_mo = np.reshape(z_mat - y_mat, (nocc, nvirt))

        # Compute hole (detachment) density matrix in MO basis
        hole_dens_mo = np.matmul(tdens_mo, tdens_mo.T)
        hole_dens_ao = np.linalg.multi_dot([mo_occ, hole_dens_mo, mo_occ.T])

        # Compute particle (attachment) density matrix in MO basis
        part_dens_mo = np.matmul(tdens_mo.T, tdens_mo)
        part_dens_ao = np.linalg.multi_dot([mo_vir, part_dens_mo, mo_vir.T])

        return {
            'detachment_density_matrix_AO': hole_dens_ao,
            'attachment_density_matrix_AO': part_dens_ao,
        }

    def compute_atom_charges(self, scf_results, rsp_results, state_index=1):
        """
        Computes Mulliken-like atomic charges for detachment and attachment densities.
        
        Uses Mulliken population analysis: q_A = sum_{mu in A} (P*S)_{mu,mu}
        where P is the density matrix and S is the overlap matrix.
        
        :param scf_results:
            The dictionary of results from converged SCF wavefunction.
        :param rsp_results:
            The dictionary containing the rsp results.
        :param state_index:
            The excited state for which the charges are computed (1-based).
            
        :return:
            A dictionary containing detachment and attachment charges per atom.
        """
        
        # Extract molecule and basis from SCF results
        molecule = self._extract_molecule_from_scf(scf_results)
        basis = self._extract_basis_from_scf(scf_results, molecule)
        
        # Get density matrices
        densities = self.compute_nto_densities(scf_results, rsp_results, state_index)
        
        # Get overlap matrix
        S = scf_results["S"]
        
        # Map atoms to AOs
        vis_drv = VisualizationDriver()
        atom_to_aos = vis_drv.map_atom_to_atomic_orbitals(molecule, basis)
        
        # Compute detachment charges using Mulliken analysis
        detach_dens_ao = densities['detachment_density_matrix_AO']
        # Note: detachment density is negative, so we take absolute values
        detach_dens_ao = np.abs(detach_dens_ao)
        PS_detach = np.matmul(detach_dens_ao, S)
        
        detachment_charges = []
        for atom_aos in atom_to_aos:
            # Sum diagonal elements of PS for this atom's AOs
            charge = np.sum([PS_detach[i, i] for i in atom_aos])
            detachment_charges.append(charge)
        
        # Compute attachment charges using Mulliken analysis
        attach_dens_ao = densities['attachment_density_matrix_AO']
        PS_attach = np.matmul(attach_dens_ao, S)
        
        attachment_charges = []
        for atom_aos in atom_to_aos:
            # Sum diagonal elements of PS for this atom's AOs
            charge = np.sum([PS_attach[i, i] for i in atom_aos])
            attachment_charges.append(charge)
        
        # Normalize charges to sum to 1.0
        detach_sum = sum(detachment_charges)
        attach_sum = sum(attachment_charges)
        
        if detach_sum > 0:
            detachment_charges = [c / detach_sum for c in detachment_charges]
        if attach_sum > 0:
            attachment_charges = [c / attach_sum for c in attachment_charges]
        
        return {
            'detachment_charges': detachment_charges,
            'attachment_charges': attachment_charges,
        }

    @staticmethod
    def compute_transition_matrix(detachment_charges, attachment_charges):
        """
        Computes the transition matrix from detachment and attachment charges.
        
        Ported from VALET (https://github.com/tbmasood/VALET) with permission
        from the author, Talha Bin Masood.
        
        :param detachment_charges:
            List of detachment charges per atom/fragment.
        :param attachment_charges:
            List of attachment charges per atom/fragment.
            
        :return:
            The transition matrix as a 2D list.
        """
        T_complete = []
        for i in range(len(detachment_charges)):
            T_complete.append([0] * len(detachment_charges))
        donors = []
        acceptors = []
        charge_diff = []
        for i in range(len(detachment_charges)):
            gs_charge = detachment_charges[i]
            es_charge = attachment_charges[i]
            if gs_charge > es_charge:
                donors.append(i)
            else:
                acceptors.append(i)
            diff = es_charge - gs_charge
            T_complete[i][i] = min(gs_charge, es_charge)
            charge_diff.append(diff)

        total_acceptor_charge = 0
        for acceptor in acceptors:
            total_acceptor_charge = total_acceptor_charge + charge_diff[acceptor]

        # heuristic
        for donor in donors:
            charge_deficit = -charge_diff[donor]
            for acceptor in acceptors:
                contrib = charge_deficit * \
                    (charge_diff[acceptor]) / total_acceptor_charge
                T_complete[acceptor][donor] = contrib

        return T_complete

    def compute_transition_data(self,
                                scf_results,
                                rsp_results,
                                subgroups,
                                state_index=1):
        """
        Computes all data needed for transition diagram visualization.
        
        :param scf_results:
            The dictionary of results from converged SCF wavefunction.
        :param rsp_results:
            The dictionary containing the rsp results.
        :param subgroups:
            A list of subgroup tuples, where each tuple contains the subgroup
            name and an atom index (or list of atom indices).
        :param state_index:
            The index of the excited state (1-based).
            
        :return:
            A dictionary containing all data needed for plotting the transition diagram.
        """
        
        # Compute NPA charges
        charges = self.compute_atom_charges(scf_results, rsp_results, state_index)
        src_detachment_charges = charges['detachment_charges']
        src_attachment_charges = charges['attachment_charges']

        num_atoms = len(src_detachment_charges)

        assert_msg_critical(
            num_atoms == len(src_attachment_charges),
            'Detachment and attachment charges must have the same length.')

        num_subgroups = len(subgroups) + 1  # +1 for 'Unassigned' subgroup
        # Initialize subgroup map to 0, which corresponds to 'Unassigned'
        sg_atom_map = {atom: 0 for atom in range(num_atoms)}
        sg_names = ['Unassigned']

        # Assign atoms to their respective subgroups
        # Note: subgroups use 1-based atom indices, convert to 0-based for internal use
        for i, (name, atom_indices) in enumerate(subgroups):
            # Handle both single atom and list of atoms
            if isinstance(atom_indices, (list, tuple)):
                for atom in atom_indices:
                    sg_atom_map[atom - 1] = i + 1  # Convert 1-based to 0-based
            else:
                sg_atom_map[atom_indices - 1] = i + 1  # Convert 1-based to 0-based
            sg_names.append(name)

        # Assign charges to subgroups
        sg_detachment_charges = [0.0] * num_subgroups
        sg_attachment_charges = [0.0] * num_subgroups

        for atom_idx, sg_idx in sg_atom_map.items():
            sg_detachment_charges[sg_idx] += src_detachment_charges[atom_idx]
            sg_attachment_charges[sg_idx] += src_attachment_charges[atom_idx]

        # Compute transition matrix
        sg_matrix = self.compute_transition_matrix(sg_detachment_charges, sg_attachment_charges)
        
        # Define colors
        colors = ["#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
                  "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3"]
        
        # Assign subgroup colors (first is gray for unassigned)
        sg_colors = ["#727272"] + [colors[i % len(colors)] for i in range(1, num_subgroups)]
        
        return {
            'subgroup_names': sg_names,
            'subgroup_colors': sg_colors,
            'detachment_charges': sg_detachment_charges,
            'attachment_charges': sg_attachment_charges,
            'transition_matrix': sg_matrix,
            'atom_to_subgroup_map': sg_atom_map,
        }

    def plot_transition_diagram(self,
                                transition_data,
                                title=None,
                                width=700,
                                height=500,
                                dpi=150):
        """
        Renders a transition diagram from pre-computed transition data.
        
        Ported from VALET (https://github.com/tbmasood/VALET) with permission
        from the author, Talha Bin Masood.
        
        :param transition_data:
            A dictionary containing:
                - 'subgroup_names': List of subgroup names
                - 'subgroup_colors': List of colors for each subgroup
                - 'detachment_charges': List of detachment charges per subgroup
                - 'attachment_charges': List of attachment charges per subgroup
                - 'transition_matrix': 2D list representing charge flow
                - 'atom_to_subgroup_map': Dictionary mapping atoms to subgroups
        :param title:
            The title of the plot.
        :param width:
            Width of the plot.
        :param height:
            Height of the plot.
            
        :return:
            A tuple containing (fig, ax) for the matplotlib figure.
        """

        assert_msg_critical('matplotlib' in sys.modules,
                            'matplotlib is required for transition diagrams.')
        
        # Extract data from dictionary
        sg_names = transition_data['subgroup_names']
        sg_colors = transition_data['subgroup_colors']
        sg_detachment_charges = transition_data['detachment_charges']
        sg_attachment_charges = transition_data['attachment_charges']
        sg_matrix = transition_data['transition_matrix']
        
        # Skip first 'unassigned' group if it has no charge
        skip_first = (len(sg_names) > 0 and 
                     sg_names[0].lower() == 'unassigned' and
                     sg_detachment_charges[0] == 0.0 and 
                     sg_attachment_charges[0] == 0.0)
        
        if skip_first:
            sg_names = sg_names[1:]
            sg_colors = sg_colors[1:]
            sg_detachment_charges = sg_detachment_charges[1:]
            sg_attachment_charges = sg_attachment_charges[1:]
            # Remove first row and column from matrix
            sg_matrix = [row[1:] for row in sg_matrix[1:]]
        
        num_subgroups = len(sg_names)

        WIDTH = width
        HEIGHT = height
        
        # Create figure and get renderer for text measurements
        fig, ax = plt.subplots(figsize=(WIDTH/100, HEIGHT/100), dpi=dpi, constrained_layout=False)
        
        # Calculate actual text heights in data coordinates
        temp_text_10pt = ax.text(0, 0, "Test", family="sans-serif", size=10)
        temp_text_16pt = ax.text(0, 0, "Test", family="serif", size=16)
        
        bbox_10pt = temp_text_10pt.get_window_extent(renderer=fig.canvas.get_renderer())
        bbox_16pt = temp_text_16pt.get_window_extent(renderer=fig.canvas.get_renderer())
        
        temp_text_10pt.remove()
        temp_text_16pt.remove()
        
        # Convert heights from display to data coordinates (pixels)
        text_height_10pt = bbox_10pt.height
        text_height_16pt = bbox_16pt.height
        
        # Layout parameters
        PADDING_LEFT = 0.08
        BAR_THICKNESS = 0.04
        GAP = 0.035
        FLOW_GAP = 0.0
        SPACING = 3  # pixels between text rows
        EDGE_PADDING = 10  # pixels at very top and bottom
        
        # Calculate heights in pixels
        leftPadding = WIDTH * PADDING_LEFT
        barHeight = HEIGHT * BAR_THICKNESS
        
        # Title height (16pt text + spacing)
        titleHeight = (text_height_16pt + SPACING) if title is not None else 0
        
        # 2 rows of 10pt text (percentage + name) + spacing
        textHeight = 2 * text_height_10pt + 2 * SPACING

        availableWidth = (1.0 - 2 * PADDING_LEFT - (num_subgroups - 1) * GAP) * WIDTH

        # Calculate vertical positions
        # Top bar: HEIGHT - edge - title - top text - bar
        topY = HEIGHT - EDGE_PADDING - titleHeight - textHeight - barHeight
        
        # Bottom bar: edge + bottom text (below bar) + bar
        bottomY = EDGE_PADDING + textHeight * 2 + barHeight

        # draw ES bars (attachment at top)
        xCoordsES = [0] * num_subgroups
        total = 0.0

        if title is not None:
            plt.text(WIDTH / 2, HEIGHT - EDGE_PADDING - titleHeight/2, title,
                     ha="center", va="center", family="serif", size=16)

        for i in range(num_subgroups):
            x = leftPadding + total * availableWidth + i * GAP * WIDTH
            barWidth = sg_attachment_charges[i] * availableWidth

            percent = "%.2f%%" % (sg_attachment_charges[i] * 100)
            plt.text(x + barWidth / 2, topY + barHeight + SPACING + text_height_10pt/2,
                    sg_names[i], ha="center", va="center", family="sans-serif", size=10)
            plt.text(x + barWidth / 2, topY + barHeight + SPACING + text_height_10pt + SPACING + text_height_10pt/2,
                    percent, ha="center", va="center", family="monospace", size=10)
            rect = mpatches.Rectangle(
                [x, topY], barWidth, barHeight, facecolor=sg_colors[i], ec="black", lw=0.5)
            ax.add_patch(rect)

            xCoordsES[i] = x
            total = total + sg_attachment_charges[i]

        # draw GS bars (detachment at bottom)
        xCoordsGS = [0] * num_subgroups
        total = 0.0

        for i in range(num_subgroups):
            x = leftPadding + total * availableWidth + i * GAP * WIDTH
            barWidth = sg_detachment_charges[i] * availableWidth

            percent = "%.2f%%" % (sg_detachment_charges[i] * 100)
            plt.text(x + barWidth / 2, bottomY - barHeight - SPACING - text_height_10pt/2,
                    sg_names[i], ha="center", va="center", family="sans-serif", size=10)
            plt.text(x + barWidth / 2, bottomY - barHeight - SPACING - text_height_10pt - SPACING - text_height_10pt/2,
                    percent, ha="center", va="center", family="monospace", size=10)
            rect = mpatches.Rectangle(
                [x, bottomY - barHeight], barWidth, barHeight, facecolor=sg_colors[i], ec="black", lw=0.5)
            ax.add_patch(rect)

            xCoordsGS[i] = x
            total = total + sg_detachment_charges[i]

        # draw connectors
        propoFilledES = [0] * num_subgroups
        propoFilledGS = [0] * num_subgroups

        for i in range(num_subgroups):
            xGS = xCoordsGS[i]
            for j in range(num_subgroups):
                if sg_matrix[j][i] == 0.0:
                    continue

                xES = xCoordsES[j]
                flow = sg_matrix[j][i]
                propGS = propoFilledGS[i]
                propES = propoFilledES[j]

                bottomXleft = xGS + availableWidth * propGS
                bottomXright = xGS + availableWidth * (propGS + flow)
                yb = bottomY + FLOW_GAP * HEIGHT
                topXleft = xES + availableWidth * propES
                topXright = xES + availableWidth * (propES + flow)
                yt = topY - FLOW_GAP * HEIGHT
                yMiddle = (yt + yb) / 2

                # add a path patch
                Path = mpath.Path
                path_data = [
                    (Path.MOVETO, [bottomXleft, yb]),
                    (Path.CURVE4, [bottomXleft, yMiddle]),
                    (Path.CURVE4, [topXleft, yMiddle]),
                    (Path.CURVE4, [topXleft, yt]),
                    (Path.LINETO, [topXright, yt]),
                    (Path.CURVE4, [topXright, yMiddle]),
                    (Path.CURVE4, [bottomXright, yMiddle]),
                    (Path.CURVE4, [bottomXright, yb]),
                    (Path.CLOSEPOLY, [bottomXleft, yb])
                ]
                codes, verts = zip(*path_data)
                path = mpath.Path(verts, codes)
                patch = mpatches.PathPatch(
                    path, facecolor=sg_colors[i], alpha=0.4, ec="none")
                ax.add_patch(patch)

                if flow > 0.05:
                    percent = "%.2f%%" % (flow * 100)
                    plt.text((topXleft + topXright + bottomXleft + bottomXright)/ 4,
                        yMiddle, percent, ha="center", family="monospace", size=8)

                propoFilledGS[i] = propoFilledGS[i] + flow
                propoFilledES[j] = propoFilledES[j] + flow

        ax.set_xlim(0, WIDTH)
        ax.set_ylim(0, HEIGHT)
        plt.axis('off')

        return fig, ax
    
    def show(self, initial_state=1, aspect_ratio=0.6):
        """
        Display interactive viewer with transition diagram and density views.
        Diagrams are generated on-demand when changing states.
        
        Must call initialize() first with scf_results, rsp_results, and subgroups.
        
        :param initial_state:
            Initial excited state to display (default: 1).
        :param aspect_ratio:
            Height to width ratio for the viewer (default: 0.6).
            Height will be calculated as width * aspect_ratio.
        """
        
        # Check if initialized
        assert_msg_critical(
            self._scf_results is not None,
            'VALET.show: Must call initialize() first with scf_results and rsp_results.'
        )
        assert_msg_critical(
            self._subgroups is not None,
            'VALET.show: Must provide subgroups via initialize() or set_subgroups().'
        )

        
        assert_msg_critical('py3Dmol' in sys.modules,
                    'py3Dmol is required for interactive viewer.')
        
        try:
            from IPython.display import display, HTML
            import ipywidgets as widgets
            import base64
            from io import BytesIO
            import json
        except ImportError:
            raise ImportError('IPython.display and ipywidgets are required for interactive viewer.')
        
        # Calculate component sizes based on layout
        # Width: Use 100% of available space
        # Height: Calculated from aspect ratio
        # Layout: [Left: Diagram | Right: [Top: Attachment | Bottom: Detachment]]
        # Note: We'll use percentage-based widths, but need pixel heights for matplotlib/py3Dmol
        
        # Assume a reference width for calculations (will use % in layout)
        reference_width = 1000  # Reference for ratio calculations
        height = int(reference_width * aspect_ratio)
        
        left_width_pct = 50  # 50% for diagram
        right_width_pct = 50  # 50% for viewers
        
        # Calculate pixel sizes for matplotlib figure and py3Dmol viewers
        diagram_width = int(reference_width * left_width_pct * 0.01)
        diagram_height = height
        viewer_width = int(reference_width * right_width_pct * 0.01)
        viewer_height = height
        
        num_states = self._num_states
        
        # Validate initial state
        if initial_state < 1 or initial_state > num_states:
            print(f"Warning: initial_state {initial_state} out of range. Using state 1.")
            initial_state = 1
        
        # Compute initial transition data
        transition_data = self.compute_transition_data(
            self._scf_results, self._rsp_results,
            self._subgroups, state_index=initial_state
        )
        
        # Create grid viewer for molecule structure with subgroup coloring
        self._mol_grid_viewer = py3Dmol.view(viewergrid=(2,1), width=viewer_width, height=viewer_height, linked=True)
        
        # Add molecule structure
        self._mol_grid_viewer.addModel(self._molecule.get_xyz_string(), "xyz")

        # Apply colors based on subgroups
        sg_colors = transition_data['subgroup_colors']
        sg_atom_map = transition_data['atom_to_subgroup_map']
        num_atoms = self._molecule.number_of_atoms()
        
        for atom_idx in range(num_atoms):
            sg_idx = sg_atom_map.get(atom_idx, 0)
            color = sg_colors[sg_idx]
            # py3Dmol uses 1-based indexing
            self._mol_grid_viewer.setStyle({'serial': atom_idx}, {
                'stick': {'color': color},
                'sphere': {'scale': 0.25, 'color': color},
            })

        self._mol_grid_viewer.zoomTo()
        
        # Store calculated sizes for use in on_state_change
        self._viewer_width = viewer_width
        self._viewer_height = viewer_height
        self._diagram_width = diagram_width
        self._diagram_height = diagram_height
        
        # Use persistent output widgets for dynamic updates
        diagram_output = widgets.Output()
        viewer_output  = widgets.Output()
        
        # Helper function to display transition diagram
        def display_diagram(state_idx):
            transition_data = self.compute_transition_data(
                self._scf_results, self._rsp_results,
                self._subgroups, state_index=state_idx
            )
            fig, ax = self.plot_transition_diagram(
                transition_data,
                width=self._diagram_width, height=self._diagram_height
            )
            # Display the figure directly in the output widget
            plt.show()
            return transition_data
        
        # Display initial state
        with diagram_output:
            transition_data = display_diagram(initial_state)
        
        with viewer_output:
            # Display initial viewer with molecule only
            display(HTML(self._mol_grid_viewer._make_html()))
        
        # Store outputs for updates
        self._diagram_output = diagram_output
        self._viewer_output  = viewer_output
        
        def on_state_change(change):
            state_idx = change['new']
            
            # Update diagram
            self._diagram_output.clear_output(wait=True)
            with self._diagram_output:
                trans_data = display_diagram(state_idx)
            
            # Recreate viewer with new state
            self._viewer_output.clear_output(wait=True)
            with self._viewer_output:
                # Create new viewer with densities
                new_viewer = self._create_viewer_with_densities(state_idx, trans_data)
                self._mol_grid_viewer = new_viewer
                display(HTML(new_viewer._make_html()))
        
        # Create dropdown
        dropdown = widgets.Dropdown(
            options=list(range(1, num_states + 1)),
            value=initial_state,
            description='State:'
        )
        dropdown.observe(on_state_change, names='value')
        
        # Layout with percentage-based widths and calculated heights
        diagram_box = widgets.VBox([
            diagram_output
        ], layout=widgets.Layout(width=f'{left_width_pct}%', height=f'{height}px', overflow='hidden'))
        
        viewer_box = widgets.VBox([
            viewer_output
        ], layout=widgets.Layout(width=f'{right_width_pct}%', height=f'{height}px', overflow='hidden'))
        
        main_layout = widgets.HBox([diagram_box, viewer_box],
                                   layout=widgets.Layout(width='100%', height=f'{height}px', overflow='hidden'))
        
        display(dropdown, main_layout)
        
        # Invoke initial density addition
        on_state_change({'new': initial_state})
        
        return None
    
    def _create_viewer_with_densities(self, state_idx, transition_data):
        """
        Create a new grid viewer with molecule, subgroup coloring, and density isosurfaces.
        
        :param state_idx:
            The excited state index (1-based).
        :param transition_data:
            Dictionary containing subgroup colors and atom mapping.
            
        :return:
            py3Dmol viewer object with all visual elements.
        """
        # Create grid viewer
        viewer = py3Dmol.view(
            viewergrid=(2,1), 
            width=self._viewer_width, 
            height=self._viewer_height, 
            linked=True
        )
        
        # Add molecule - addModel without viewer param adds to all grid positions
        xyz_str = self._molecule.get_xyz_string()
        viewer.addModel(xyz_str, "xyz")
        
        # Apply subgroup coloring to all viewers
        sg_colors = transition_data['subgroup_colors']
        sg_atom_map = transition_data['atom_to_subgroup_map']
        num_atoms = self._molecule.number_of_atoms()
        
        for atom_idx in range(num_atoms):
            sg_idx = sg_atom_map.get(atom_idx, 0)
            color = sg_colors[sg_idx]
            style = {
                'stick': {'color': color},
                'sphere': {'scale': 0.25, 'color': color}
            }
            viewer.setStyle({'serial': atom_idx}, style)
        
        # Zoom before adding densities
        viewer.zoomTo()
        
        # Compute densities
        densities = self.compute_nto_densities(
            self._scf_results, self._rsp_results, state_idx
        )
        
        # Compute attachment density
        attachment_data = self._density_viewer.compute_density(
            densities['attachment_density_matrix_AO']
        )
        attachment_cube_str = self._density_viewer._get_density_cube_str(
            attachment_data, f"S{state_idx}_attachment"
        )
        
        # Compute detachment density
        detachment_data = self._density_viewer.compute_density(
            densities['detachment_density_matrix_AO']
        )
        detachment_cube_str = self._density_viewer._get_density_cube_str(
            detachment_data, f"S{state_idx}_detachment"
        )
        
        # Add volumetric data directly using viewer parameter
        # Attachment density (green) to top viewer (0,0)
        viewer.addVolumetricData(
            attachment_cube_str, "cube",
            {'isoval': self.density_isovalue, 'color': self.attachment_color, 'opacity': self.density_opacity},
            viewer=(0,0)
        )
        
        # Detachment density (red) to bottom viewer (1,0)  
        viewer.addVolumetricData(
            detachment_cube_str, "cube",
            {'isoval': self.density_isovalue, 'color': self.detachment_color, 'opacity': self.density_opacity},
            viewer=(1,0)
        )
        
        return viewer
