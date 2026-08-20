"""Driver for computing atom- and molecule-level QM/empirical descriptors."""
 
from collections import defaultdict
from typing import Any, Dict, Optional
 
import numpy as np
import veloxchem as vlx
 
 
# Unit conversion / physical constants used below.
HARTREE_TO_EV = 27.211407953
HARTREE_TO_J_PER_MOL = 2625500.2
GAS_CONSTANT_J_PER_MOL_K = 8.3144626  # R
ROOM_TEMPERATURE_K = 298.15
 
 
class DescriptorDriver:
    """
    Computes a set of atom- and molecule-level descriptors from a converged
    SCF calculation:
 
    - RESP atomic partial charges
    - Local surface ionization energy (IE) and electron affinity (EA),
      evaluated on a CPCM-derived molecular surface and averaged per atom
    - GAFF atom types
    - logP (water / 1-octanol) and total SASA
 
    Typical usage
    -------------
        driver = DescriptorDriver()
        results = driver.compute_descriptors(molecule)
    """
 
    def __init__(self):
        self.molecule = None       # molecule for descriptor calculation
        self.basis = None
        self.scf_drv = None
        self.scf_results = None
 
        self.molgrid = None
        self.surface_points = None
        self.point_atom_indices = None
        self.occ_mo_points_amplitude = []      # MO amplitude for occupied orbitals, per point
        self.unocc_mo_points_amplitude = []    # MO amplitude for unoccupied orbitals, per point
        self.ea_values = []
        self.ie_values = []
        self.atomtypes = []                    # GAFF atom types, indexed like self.molecule's atoms
        self.sasa = None                       # Solvent-accessible surface area (Å^2)
        self.log_p = None
        self.resp_charges = None
 
        self.ea_ie_results = None
        self.results_dict = {}
 
        # Level of theory for the auxiliary logP calculation. Deliberately
        # decoupled from self.basis/self.scf_drv: logP is estimated from a
        # water/1-octanol solvation free-energy difference, which is
        # conventionally computed at a fixed reference level of theory
        # rather than whatever level the caller chose for the main
        # calculation. Override on the instance before calling
        # compute_descriptors() if you need something else.
        self.logp_basis_label = "def2-svp"
        self.logp_xcfun = "b3lyp"
 
        # Basis used for the RESP charge fit. Same rationale as above.
        self.resp_basis_label = "6-31G*"
 
    def compute_respcharges(self):
        """Compute RESP atomic partial charges and store them in self.resp_charges."""
        resp_drv = vlx.RespChargesDriver()
        resp_basis = vlx.MolecularBasis.read(self.molecule, self.resp_basis_label)
        self.resp_charges = resp_drv.compute(self.molecule, resp_basis)
 
    def compute_gto_values_at_custom_points(self, custom_points: np.ndarray) -> np.ndarray:
        """
        Evaluate the GTO basis functions at an arbitrary set of points.
 
        Work-around for vlx.XCIntegrator, which normally operates on its own
        internally generated grid — this lets us reuse it on the CPCM
        surface grid instead (which showed less instability with respect to 
        grid density compared to molgridriver).
 
        Parameters
        ----------
        custom_points
            Array of shape (N, 3), point coordinates in bohr.
 
        Returns
        -------
        np.ndarray
            GTO values at each point, shape (n_basis_functions, N).
        """
        from veloxchem.veloxchemlib import DenseMatrix, MolecularGrid
 
        weights = np.ones(len(custom_points))
        molgrid_points = np.vstack((
            custom_points[:, 0],
            custom_points[:, 1],
            custom_points[:, 2],
            weights,
        ))
 
        coords_mat = DenseMatrix(molgrid_points)
        mol_grid = MolecularGrid(coords_mat)
        mol_grid.partition_grid_points()
        mol_grid.distribute_counts_and_displacements(0, 1)
 
        xc_drv = vlx.XCIntegrator()
        chi_g = xc_drv.compute_gto_values(self.molecule, self.basis, mol_grid)
        return chi_g
 
    def create_surface_points(self):
        """Generate a CPCM molecular surface grid and store it in self.surface_points."""
        cpcm_drv = vlx.CpcmDriver()
        raw_cpcm_points = cpcm_drv.generate_cpcm_grid(self.molecule)
        self.molgrid = raw_cpcm_points[0][:, :3]
 
        chi_g = self.compute_gto_values_at_custom_points(self.molgrid)
        D = self.scf_results["D_alpha"] + self.scf_results["D_beta"]  # total electron density matrix
        G = np.einsum("ab,bg->ag", D, chi_g)
        n_g = np.einsum("ag,ag->g", chi_g, G)  
 
        self.surface_points = self.molgrid.copy()
        print(f"No. of surface points created: {len(self.surface_points)}")
 
    def assign_points_to_atoms(self):
        """
        Assign each surface point to its nearest atom, using van der Waals
        radius-adjusted distances, and store the result in
        self.point_atom_indices (same length as self.surface_points).
        """
        atom_positions = self.molecule.get_coordinates_in_bohr()      # (n_atoms, 3)
        atom_vdw_radii = self.molecule.vdw_radii_to_numpy()           # (n_atoms,)
 
        point_atom_indices = []
        for point in self.surface_points:
            distances = np.linalg.norm(atom_positions - point, axis=1)
            adjusted_distances = distances - atom_vdw_radii
            closest_atom_index = np.argmin(adjusted_distances)
            point_atom_indices.append(closest_atom_index)
        self.point_atom_indices = point_atom_indices
 
        assert len(self.surface_points) == len(self.point_atom_indices), (
            f"Mismatch! Created {len(self.surface_points)} surface points "
            f"but generated {len(self.point_atom_indices)} point-to-atom indices"
        )
 
    def molecular_orbital_amplitude_per_point(self):
        """
        Compute the amplitude of every occupied and unoccupied molecular
        orbital at each surface point, storing them in
        self.occ_mo_points_amplitude / self.unocc_mo_points_amplitude.
        """
        chi_g = self.compute_gto_values_at_custom_points(self.molgrid)
        mol_orb_coeff = self.scf_results["C_alpha"]
        occ_orb = sum(self.scf_results["occ_alpha"])
        tot_orb = len(self.scf_results["occ_alpha"])
 
        for i in range(int(occ_orb)):
            mo_val = np.dot(chi_g.T, mol_orb_coeff[:, i])
            self.occ_mo_points_amplitude.append(mo_val)
 
        for i in range(int(occ_orb), int(tot_orb)):
            mo_val = np.dot(chi_g.T, mol_orb_coeff[:, i])
            self.unocc_mo_points_amplitude.append(mo_val)
 
    def compute_IE_and_EA(self):
        """
        Compute local surface ionization energy (IE) and electron affinity
        (EA) at each surface point, in eV, from the amplitude-weighted
        average of orbital energies (analogous to the average local
        ionization energy formalism). Results are stored in self.ie_values
        and self.ea_values, one entry per surface point.
        """
        n_occ = int(sum(self.scf_results["occ_alpha"]))
        energy_occ = self.scf_results["E_alpha"][:n_occ]      # Hartree
        energy_unocc = self.scf_results["E_alpha"][n_occ:]    # Hartree
 
        for idx, _point in enumerate(self.surface_points):
            ea_numerator, ea_denominator = 0.0, 0.0
            for j in range(len(self.unocc_mo_points_amplitude)):
                amp_sq = self.unocc_mo_points_amplitude[j][idx] ** 2
                ea_numerator += amp_sq * energy_unocc[j]
                ea_denominator += amp_sq
            ea = (-ea_numerator / ea_denominator) * HARTREE_TO_EV
            self.ea_values.append(ea)
 
            ie_numerator, ie_denominator = 0.0, 0.0
            for j in range(len(self.occ_mo_points_amplitude)):
                amp_sq = self.occ_mo_points_amplitude[j][idx] ** 2
                ie_numerator += amp_sq * abs(energy_occ[j])
                ie_denominator += amp_sq
            ie = (ie_numerator / ie_denominator) * HARTREE_TO_EV
            self.ie_values.append(ie)
 
        print(
            f"{len(self.ea_values)}, {len(self.ie_values)} values of electron "
            f"affinity and ionization energy calculated for "
            f"{len(self.surface_points)} surface points"
        )
 
    def compute_atom_statistics(self):
        """
        Aggregate per-point IE/EA values into per-atom statistics
        (min/max/mean/std/n_points), stored in self.ea_ie_results, keyed by
        atom index.
        """
        point_data = []
        for i in range(len(self.surface_points)):
            point_data.append({
                "point_index": i + 1,
                "atom_label": self.molecule.get_label(self.point_atom_indices[i]),
                "atom_index": int(self.point_atom_indices[i]),
                "IE": self.ie_values[i],
                "EA": self.ea_values[i],
            })
 
        ea_data = defaultdict(list)
        ie_data = defaultdict(list)
        for point in point_data:
            atom_idx = point["atom_index"]
            ea_data[atom_idx].append(point["EA"])
            ie_data[atom_idx].append(point["IE"])
 
        atom_statistics = {}
        for atom in range(len(self.molecule.get_labels())):
            if not ea_data[atom] or not ie_data[atom]:
                atom_statistics[atom] = {
                    "min_EA": np.nan, "max_EA": np.nan, "mean_EA": np.nan, "std_EA": np.nan,
                    "n_points": 0,
                    "min_IE": np.nan, "max_IE": np.nan, "mean_IE": np.nan, "std_IE": np.nan,
                }
                continue
            atom_statistics[atom] = {
                "min_EA": min(ea_data[atom]),
                "max_EA": max(ea_data[atom]),
                "mean_EA": np.mean(ea_data[atom]),
                "std_EA": np.std(ea_data[atom]),
                "n_points": len(ea_data[atom]),
                "min_IE": min(ie_data[atom]),
                "max_IE": max(ie_data[atom]),
                "mean_IE": np.mean(ie_data[atom]),
                "std_IE": np.std(ie_data[atom]),
            }
 
        self.ea_ie_results = atom_statistics
 
    def define_atom_types(self):
        """Assign a GAFF atom type to every atom via vlx.AtomTypeIdentifier."""
        atomtypeidentifier = vlx.AtomTypeIdentifier()
        self.atomtypes = atomtypeidentifier.generate_gaff_atomtypes(self.molecule)
 
    def compute_log_p_and_SASA(self):
        """
        Estimate logP from the water/1-octanol SMD solvation free-energy
        difference, and compute total SASA (Å^2).
 
        Runs two additional SCF calculations (water, 1-octanol solvation)
        at self.logp_basis_label / self.logp_xcfun, independent of the
        basis/driver used for the main descriptor calculation — see the
        note on those attributes in __init__.
        """
        smd = vlx.SmdDriver()
        smd.solute = self.molecule
        sasa_list = smd._get_SASA()
        self.sasa = sum(sasa_list)  # Å^2
 
        logp_basis = vlx.MolecularBasis.read(self.molecule, self.logp_basis_label)
 
        scf_drv = vlx.ScfRestrictedDriver()
        scf_drv.xcfun = self.logp_xcfun
        scf_drv.solvation_model = "smd"
        scf_drv.smd_solvent = "water"
        scf_results_water = scf_drv.compute(self.molecule, logp_basis)
 
        scf_drv = vlx.ScfRestrictedDriver()
        scf_drv.xcfun = self.logp_xcfun
        scf_drv.solvation_model = "smd"
        scf_drv.smd_solvent = "1-octanol"
        scf_results_octanol = scf_drv.compute(self.molecule, logp_basis)
 
        ddg_solv_hartree = scf_results_octanol["scf_energy"] - scf_results_water["scf_energy"]
        ddg_solv_j_mol = ddg_solv_hartree * HARTREE_TO_J_PER_MOL
        self.log_p = -ddg_solv_j_mol / (
            np.log(10) * GAS_CONSTANT_J_PER_MOL_K * ROOM_TEMPERATURE_K
        )
 
    def compile_results(self):
        """Assemble every computed descriptor into self.results_dict."""
        self.results_dict["IE_EA"] = self.ea_ie_results
        self.results_dict["atomtypes"] = self.atomtypes
        self.results_dict["log_p"] = self.log_p
        self.results_dict["sasa"] = self.sasa
        self.results_dict["RESP_charges"] = self.resp_charges
        self.results_dict["ie_surface_average"] = np.mean(self.ie_values)
        self.results_dict["ea_surface_average"] = np.mean(self.ea_values)
        self.results_dict["scf_results"] = self.scf_results
 
    def compute_descriptors(
        self,
        molecule: "vlx.Molecule",
        basis: Optional["vlx.MolecularBasis"] = None,
        scf_drv: Optional["vlx.ScfRestrictedDriver"] = None,
        scf_results: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        """
        Compute the full set of descriptors for `molecule` and return them
        as a dict.
 
        Parameters
        ----------
        molecule
            Molecule to compute descriptors for.
        basis
            Basis set to use for the SCF/surface calculations. Defaults to
            def2-svpd if not provided.
        scf_drv
            SCF driver to use. Defaults to a restricted B3LYP driver if not
            provided.
        scf_results
            Pre-computed SCF results for molecule/basis/scf_drv, to avoid
            recomputation if already available.
 
        Returns
        -------
        dict
            See compile_results() for the keys included.
        """
        self.molecule = molecule
 
        if basis is None:
            basis = vlx.MolecularBasis.read(molecule, "def2-svpd")
        self.basis = basis
 
        if scf_drv is None:
            scf_drv = vlx.ScfRestrictedDriver()
            scf_drv.xcfun = "b3lyp"
        self.scf_drv = scf_drv
 
        if scf_results is None:
            scf_results = scf_drv.compute(self.molecule, self.basis)
        self.scf_results = scf_results
 
        print("computing RESP-charges...")
        self.compute_respcharges()
 
        print("Creating surface grid...")
        self.create_surface_points()
 
        print("Assigning grid points to atoms...")
        self.assign_points_to_atoms()
 
        print("Computing atomic IE and EA energies...")
        self.molecular_orbital_amplitude_per_point()
        self.compute_IE_and_EA()
 
        self.compute_atom_statistics()
        self.define_atom_types()
        self.compute_log_p_and_SASA()
 
        self.compile_results()
        return self.results_dict


