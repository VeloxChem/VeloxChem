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

from mpi4py import MPI
import numpy as np
from typing import List, Tuple, Dict, Optional
from sys import stdout
import os
import weakref
import h5py
from contextlib import redirect_stderr
from io import StringIO

from .interpolationdatapoint import (InterpolationDatapoint)
from .outputstream import OutputStream
from .veloxchemlib import mpi_master, bohr_in_angstrom
from .errorhandler import assert_msg_critical
from .inputparser import parse_input


with redirect_stderr(StringIO()) as fg_err:
    import geometric


class InterpolationDriver:
    """
    Implements the modified Shepard interpolation (MSI) of a potential
    energy surface from a database of second-order Taylor expansions in
    internal coordinates.

    Every datapoint i (InterpolationDatapoint) provides the local model

        T_i = E_i + g_i^T s + 1/2 s^T H_i s,

    where s is the internal-coordinate displacement of the current geometry
    from the datapoint: Delta q for bonds and angles, 2 sin(Delta/2) for
    proper dihedrals and 2 tan(Delta/2) for impropers, with Delta wrapped to
    [-pi, pi). The interpolated energy and gradient are

        E = sum_i W_i T_i,   grad E = sum_i (W_i grad T_i + T_i grad W_i),
        W_i = w_i / sum_j w_j,   w_i = 1 / (u^p + u^q),   u = (d_i / R_i)^2,

    with d_i the distance to datapoint i (Kabsch-aligned Cartesian distance
    over the atoms in the weights for weightfunction_type='cartesian', an
    internal-coordinate distance for 'internal'), R_i its confidence radius
    and p, q = exponent_p, exponent_q, so that the effective distance powers
    are 2p and 2q. Raw weights with u > 1e6 are zero. With use_tc_weights the
    weights also depend on the important coordinates of each datapoint.

    :param z_matrix:
        The z-matrix dictionary {'bonds', 'angles', 'dihedrals', 'impropers'}
        of 0-based atom-index tuples. It must equal the z-matrix of the
        datapoints (checked at the first evaluation of every datapoint).
    :param comm:
        The MPI communicator (used for the output stream only; every rank
        computes the same result).
    :param ostream:
        The output stream.

    Settings (keys of update_settings)
        - interpolation_type: 'shepard', the only implemented scheme.
        - weightfunction_type: 'cartesian' (default) or 'internal'.
        - exponent_p: the exponent p of u (int); required before compute.
        - exponent_q: the exponent q of u (int). Default p/2; the default
          follows later changes of exponent_p until q is given.
        - use_inverse_bond_length, use_eq_bond_length: the bond coordinate,
          1/r in 1/bohr, r_eq - r_eq^2/r in bohr, or r in bohr if both are
          False. use_inverse_bond_length defaults to True; a settings
          dictionary that gives only one of the two keys as True sets the
          other one to False, and both True raises. The mode selects the HDF5
          label suffix ('_rinv', '_eq', '_r') and must be the mode the
          database was built with.
        - use_mass_weight: the construction setting. The Taylor gradients are
          converted from mass-weighted coordinates whenever
          impes_coordinate.inv_sqrt_masses is set (the B-matrix is then
          mass-weighted), independent of this flag; set inv_sqrt_masses on
          impes_coordinate and on the datapoints exactly when the database
          was built with use_mass_weight=True.
        - use_tc_weights: use the target-customised weights.
        - tc_weight_mode: 'multiplicative' (default): w_i exp(-lambda A_i),
          A_i the mean squared displacement of the important coordinates of
          datapoint i in units of their tolerances; 'additive_rhee':
          u = (d_i^2 + D_tc) / R_i^2 (Cartesian weights only).
        - imforcefield_file: the HDF5 database, read by compute if
          qm_data_points is None.
        - labels: the datapoint labels to read (default: all).
        - confidence_radius: unused; the radii are stored per datapoint.

    Parameters that are attributes only (set them directly)
        - tc_imp_gate_lambda (0.7), tc_imp_bond_sigma_angstrom (0.04),
          tc_imp_angle_sigma_degrees (6), tc_imp_dihedral_sigma_degrees (12),
          tc_imp_improper_sigma_degrees (6), tc_imp_gate_exp_clip (50),
          tc_imp_min_sigma: the multiplicative tc weights.
        - alpha: the decay coefficient of the tc factor exp(-alpha D_imp)
          of the internal weights.
        - weight_eps (1e-12): if the sum of the raw weights is not finite or
          not above weight_eps, the weights are divided by their maximum
          before the normalisation (mathematically identical).
        - calc_optim_trust_radius: also compute dw/dR and d(grad w)/dR for
          IMTrustRadiusOptimizer.

    Input objects
        - qm_data_points: the list of InterpolationDatapoint objects, built
          with the same z-matrix, bond mode and inv_sqrt_masses.
        - impes_coordinate: the InterpolationDatapoint of the current
          geometry; set impes_coordinate.inv_sqrt_masses and, in eq mode,
          impes_coordinate.eq_bond_lengths before compute.
        - symmetry_information: list built by IMForceFieldGenerator or
          OpenMMIMDynamics: [0] all atom indices, [1] rotor groups,
          [2] regrouped symmetry groups, [3] the atoms in the weights (the
          sorted complement of [4]), [4] the atoms excluded from the weights
          (CX3 hydrogens with OpenMMIMDynamics.use_symmetry; empty during
          construction), [5] the 0-based rotatable bonds, [6]-[8] symmetry
          dihedral bookkeeping, [-1] (dihedral_start, dihedral_end) of the
          proper dihedrals in the flat z-matrix. The driver reads [3] and
          [4]; IMTrustRadiusOptimizer reads [3], the collector reads [5].

    Outputs after compute (read by the radius optimiser, the collector and
    analysis scripts)
        - impes_coordinate.energy (Hartree) and impes_coordinate.gradient
          (Hartree/bohr, (natoms, 3)); see get_energy and get_gradient.
        - weights: dictionary {index in qm_data_points: W_i}.
        - potentials, gradients: T_i and grad T_i of every datapoint.
        - raw_weights_array, normalized_weights_array: w_i and W_i.
        - sum_of_weights and raw_sum_of_weights (the same value),
          sum_of_weights_grad, effective_sum_of_weights(_grad) (after the
          rescaling), weight_rescale_factor, used_weight_rescale,
          used_weight_fallback (no finite positive weight: the nearest
          datapoint gets W = 1 and grad W = 0).
        - dw_dalpha_list, dw_dX_dalpha_list: dw_i/dR_i and d(grad w_i)/dR_i
          (only with calc_optim_trust_radius; "alpha" is the radius R).
        - pair_term_records: None by default. If it is set to a list before
          compute, cartesian_distance appends (d, v, A, grad A) for every
          datapoint: the distance, the aligned displacement of the atoms in
          the weights, and the multiplicative tc gate metric with its
          gradient (None, None in the other weight modes). These terms do
          not depend on the radius; the cached confidence-radius fit of
          IMDatabasePointCollecter reads them.
        - bond_rmsd, angle_rmsd, dihedral_rmsd: per datapoint the L2 norms
          (not RMS values) of the bond, angle and proper-dihedral parts of s.
    """

    def __init__(self, z_matrix=None, comm=None, ostream=None):
        """
        Initializes the interpolation driver.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(stdout)
            else:
                ostream = OutputStream(None)

        self.comm = comm
        self.ostream = ostream
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()
        self.impes_coordinate = InterpolationDatapoint(z_matrix, ostream=self.ostream)

        # only Shepard interpolation is implemented
        self.interpolation_type = 'shepard'
        self.weightfunction_type = 'cartesian'
        self.exponent_p = None
        self.exponent_q = None
        self.confidence_radius = 0.5
        self.alpha = 1.0

        self.z_matrix = z_matrix
        self.impes_dict = None

        self.sum_of_weights = None
        self.sum_of_weights_grad = None

        # Raw and normalized Shepard-weight bookkeeping.
        self.raw_weights_array = None
        self.normalized_weights_array = None
        self.raw_sum_of_weights = None
        self.used_weight_fallback = False
        self.weight_eps = 1.0e-12

        # symmetry metadata produced during database-point generation
        # (layout in the class docstring)
        self.symmetry_information = None

        self.weights = {}
        self.potentials = []
        self.gradients = []

        # HDF5 database with the QM datapoints
        self.imforcefield_file = None
        self.qm_data_points = None

        # Confidence radii optimization information
        self.dw_dalpha_list = []
        self.dw_dX_dalpha_list = []
        self.calc_optim_trust_radius = False

        # radius-independent terms of every datapoint (see the class docstring)
        self.pair_term_records = None

        self.bond_rmsd = None
        self.angle_rmsd = None
        self.dihedral_rmsd = None
        # labels of the QM datapoints
        self.labels = None
        self.use_inverse_bond_length = True
        self.use_eq_bond_length = False
        self.use_mass_weight = False

        # General target customized schemes
        self.use_tc_weights = False
        # 'multiplicative': w * exp(-lambda A) with the coordinate-specific
        # tolerances below; 'additive_rhee': u = (d^2 + D_tc) / R^2
        self.tc_weight_mode = "multiplicative"

        # multiplicative weight function scheme
        self.tc_imp_gate_lambda = 0.7
        # Coordinate tolerances defining the sphere of influence.
        self.tc_imp_bond_sigma_angstrom = 0.04
        self.tc_imp_angle_sigma_degrees = 6.0
        self.tc_imp_dihedral_sigma_degrees = 12.0
        self.tc_imp_improper_sigma_degrees = 6.0
        # Numerical guards.
        self.tc_imp_gate_exp_clip = 50.0
        self.tc_imp_min_sigma = 1.0e-12

        # datapoints whose z-matrix was checked against impes_coordinate,
        # mapped to the (datapoint z-matrix, driver z-matrix) objects checked
        self._z_matrix_checked_datapoints = weakref.WeakKeyDictionary()

        # exponent_q derived from exponent_p (None if not derived)
        self._derived_exponent_q = None

        self._input_keywords = {
            'im_settings': {
                'interpolation_type':
                    ('str', 'type of interpolation (Shepard)'),
                'weightfunction_type':
                    ('str', 'type of interpolation (cartesian/internal)'),
                'exponent_p': ('int', 'the main exponent'),
                'exponent_q': ('int', 'the additional exponent (Shepard IM)'),
                'confidence_radius': ('float', 'unused: the radii are stored per datapoint'),
                'imforcefield_file':
                    ('str', 'the name of the chk file with QM data'),
                    'use_inverse_bond_length': ('bool', 'whether to use inverse bond lengths in the Z-matrix'),
                    'use_eq_bond_length': ('bool', 'whether to use eq bond lengths in the Z-matrix'),
                    'use_tc_weights':('bool', 'whether to use target-customised weights'),
                    'tc_weight_mode':('str', 'the mode for the target customized weights (multiplicative/additive)'),
                    'use_mass_weight':('bool', 'construction setting; the metric follows impes_coordinate.inv_sqrt_masses'),
                'labels': ('seq_fixed_str', 'the list of QM data point labels'),
            }
        }

    def update_settings(self, impes_dict=None):
        """
        Updates the settings of the driver and of impes_coordinate.

        Raises for interpolation_type other than 'shepard', for both bond
        flags True, and for unknown tc_weight_mode values. exponent_q is set
        to exponent_p / 2 unless it is given.

        :param impes_dict:
            The dictionary of interpolation settings (see the class
            docstring).
        """

        if impes_dict is None:
            impes_dict = {}

        self.impes_dict = dict(impes_dict)

        self.impes_coordinate.update_settings(impes_dict)

        im_keywords = {
            key: val[0]
            for key, val in self._input_keywords['im_settings'].items()
        }

        parse_input(self, im_keywords, impes_dict)

        if self.use_eq_bond_length and not self.use_inverse_bond_length:
            self.use_inverse_bond_length = False
        elif self.use_inverse_bond_length and not self.use_eq_bond_length:
            self.use_eq_bond_length = False
        elif self.use_eq_bond_length and self.use_inverse_bond_length:
            assert_msg_critical(
                False, 'use_inverse_bond_length and use_eq_bond_length '
                'cannot both be True. Set exactly one of them (or neither, for '
                'plain bond lengths in bohr).')

        if self.interpolation_type != 'shepard':
            raise ValueError(
                "Unrecognized interpolation type: "
                f"{self.interpolation_type}. Only 'shepard' is supported.")

        # exponent_q defaults to exponent_p / 2. The derived value follows
        # exponent_p on later updates, a value given by the user is kept.
        q_given = impes_dict.get('exponent_q', None) is not None
        q_is_derived = (self.exponent_q is None or
                        self.exponent_q is self._derived_exponent_q)
        if not q_given and q_is_derived:
            if self.exponent_p is not None:
                self.exponent_q = self.exponent_p / 2.0
                self._derived_exponent_q = self.exponent_q
            else:
                self.exponent_q = None
                self._derived_exponent_q = None
        else:
            self._derived_exponent_q = None

        valid_tc_modes = {"multiplicative", "additive_rhee"}

        if self.tc_weight_mode not in valid_tc_modes:
            raise ValueError(
                "Unknown tc_weight_mode: "
                f"{self.tc_weight_mode}. "
                "Allowed values are 'multiplicative' and 'additive_rhee'."
            )

        if (
            self.use_tc_weights
            and self.tc_weight_mode == "additive_rhee"
            and self.weightfunction_type != "cartesian"
        ):
            raise NotImplementedError(
                "The additive Rhee target-customized formalism is currently "
                "implemented only for Cartesian base weighting coordinates."
            )

    def mark_runtime_data_cache_dirty(self):
        """
        Does nothing; kept for backwards compatibility.

        The driver keeps no runtime cache of datapoint data: every compute()
        reads qm_data_points (radii, coordinates, gradients) as they are, so
        changes to the datapoints take effect without calling this method.
        """

        return None

    def _get_internal_coordinate_partitions(self):
        """
        Returns the index ranges of the coordinate types in the flat
        z-matrix of impes_coordinate.

        Layout:
            bonds:            [0, bond_end)
            angles:           [bond_end, angle_end)
            proper dihedrals: [dihedral_start, dihedral_end)
            impropers/tail:   [dihedral_end, n_ic)
        """
        zmat_dict = getattr(self.impes_coordinate, 'z_matrix_dict', None)
        assert_msg_critical(
            isinstance(zmat_dict, dict),
            'InterpolationDriver: impes_coordinate has no z-matrix dictionary.')

        n_bonds = len(zmat_dict.get('bonds', []))
        n_angles = len(zmat_dict.get('angles', []))
        n_dihedrals = len(zmat_dict.get('dihedrals', []))
        dihedral_start = n_bonds + n_angles
        dihedral_end = dihedral_start + n_dihedrals
        return {
            'bond_end': int(n_bonds),
            'angle_end': int(dihedral_start),
            'dihedral_start': int(dihedral_start),
            'dihedral_end': int(dihedral_end),
        }

    def read_labels(self):
        """
        Reads the datapoint labels and the z-matrix from imforcefield_file.

        Only datapoints of the bond-coordinate mode of impes_coordinate
        (label suffix '_rinv', '_eq' or '_r') are listed.

        :return:
            A tuple (labels, z_matrix): the labels without suffix in HDF5 key
            order, which is lexicographic (point_1, point_10, point_100, ...),
            and the z-matrix dictionary of the first complete datapoint.
            Callers that build qm_data_points in this order get the weights
            indexed in this order.
        """

        assert_msg_critical(
            self.imforcefield_file is not None,
            'ImpesDriver: Please provide a chekpoint file name.')

        with h5py.File(os.fspath(self.imforcefield_file), 'r') as h5f:

            keys = h5f.keys()

            remove_from_label = ""
            z_matrix_label = ''
            if self.impes_coordinate.use_inverse_bond_length:
                remove_from_label += "_rinv"
                z_matrix_label += '_rinv'

            elif self.impes_coordinate.use_eq_bond_length:
                remove_from_label += "_eq"
                z_matrix_label += '_eq'

            else:
                remove_from_label += "_r"
                z_matrix_label += '_r'

            remove_from_label += "_dihedral"
            z_matrix_label += '_dihedral'

            remove_from_label += "_energy"

            z_matrix_bonds = z_matrix_label + '_bonds'
            z_matrix_angles = z_matrix_label + '_angles'
            z_matrix_dihedrals = z_matrix_label + '_dihedrals'
            z_matrix_impropers = z_matrix_label + '_impropers'

            z_matrix_labels = [z_matrix_bonds, z_matrix_angles, z_matrix_dihedrals, z_matrix_impropers]
            z_matrix = {"bonds": [], "angles": [], "dihedrals": [], "impropers": []}

            labels = []
            counter = 0

            for key in keys:

                if remove_from_label in key:

                    label = key.replace(remove_from_label, "")
                    if label not in labels:
                        labels.append(label)

                    # the z-matrix is taken from the first point that has all
                    # four z-matrix datasets (a damaged point is skipped)
                    if counter == 0 and all(
                            (label + label_obj) in h5f for label_obj in z_matrix_labels):
                        for label_obj in z_matrix_labels:

                            z_label = label + label_obj
                            current_z_list = [tuple(z_list.tolist()) for z_list in list(h5f.get(z_label))]
                            if 'bonds' in z_label:
                                z_matrix['bonds'] = current_z_list
                            elif 'angles' in z_label:
                                z_matrix['angles'] = current_z_list
                            elif 'dihedrals' in z_label:
                                z_matrix['dihedrals'] = current_z_list
                            elif 'impropers' in z_label:
                                z_matrix['impropers'] = current_z_list
                        counter = 1

        return labels, z_matrix

    def define_impes_coordinate(self, coordinates):
        """
        Sets the current geometry of impes_coordinate (B-matrix and
        internal coordinates).

        :param coordinates:
            The Cartesian coordinates in bohr, (natoms, 3).
        """

        self.impes_coordinate.reset_coordinates_impes_driver(coordinates)

    def compute(self, molecule):
        """
        Computes the energy and gradient at the geometry of molecule by
        Shepard interpolation over qm_data_points (read from
        imforcefield_file if qm_data_points is None).

        The results are stored in impes_coordinate.energy (Hartree) and
        impes_coordinate.gradient (Hartree/bohr), see get_energy and
        get_gradient, together with the outputs listed in the class
        docstring.

        :param molecule:
            The molecule; its coordinates in bohr are used.
        """

        assert_msg_critical(
            self.exponent_p is not None and self.exponent_q is not None,
            'InterpolationDriver: exponent_p is not set; pass it to '
            'update_settings before compute.')

        self.bond_rmsd = []
        self.angle_rmsd = []
        self.dihedral_rmsd = []

        self.molecule = molecule

        self.define_impes_coordinate(molecule.get_coordinates_in_bohr())

        if self.qm_data_points is None:
            self.qm_data_points = self.read_qm_data_points()

        if self.interpolation_type == 'shepard':
            self.shepard_interpolation()
        else:
            errtxt = "Unrecognized interpolation type: "
            errtxt += self.interpolation_type
            raise ValueError(errtxt)

    def shepard_interpolation(self):
        """
        Performs the Shepard interpolation over qm_data_points at the
        current impes_coordinate (set by compute).

        E = sum_i W_i T_i and grad E = sum_i (W_i grad T_i + T_i grad W_i)
        with W_i = w_i / sum_j w_j. If the sum of the raw weights is not
        finite or not above weight_eps, the weights are first divided by
        their maximum (mathematically identical). If no finite positive
        weight exists, the nearest datapoint gets W = 1 and grad W = 0.
        """

        natms = self.impes_coordinate.cartesian_coordinates.shape[0]

        potentials = []
        gradients = []
        self.weights = {}
        weights_cart = []
        weight_gradients_cart = []
        used_labels = []
        self.potentials = []
        self.gradients = []
        self.dw_dalpha_list = []
        self.dw_dX_dalpha_list = []

        weight_fallback_metrics = []
        masses = self.molecule.get_masses().copy()
        masses_cart = np.repeat(masses, 3)
        sqrt_masses = np.sqrt(masses_cart)

        # The B-matrix of impes_coordinate is mass-weighted whenever
        # inv_sqrt_masses is set (independent of use_mass_weight), so the
        # Taylor gradients are converted back from mass-weighted coordinates
        # exactly in that case, with the masses used for the B-matrix.
        inv_sqrt_masses = self.impes_coordinate.inv_sqrt_masses
        mass_weighted_b_matrix = (inv_sqrt_masses is not None)
        if mass_weighted_b_matrix and not np.allclose(
                sqrt_masses * inv_sqrt_masses, 1.0, rtol=0.0, atol=1.0e-12):
            sqrt_masses = 1.0 / np.asarray(inv_sqrt_masses, dtype=np.float64)

        distances_and_gradients = []

        if self.weightfunction_type == 'cartesian':
            for dp_index, data_point in enumerate(self.qm_data_points[:]):

                self._check_datapoint_z_matrix(data_point)

                distance, dihedral_dist, denominator, weight_gradient, distance_vec, _, _, _ = self.cartesian_distance(data_point)

                distances_and_gradients.append((distance, dihedral_dist, dp_index, denominator, weight_gradient, distance_vec))

        elif self.weightfunction_type == 'internal':
            for dp_index, data_point in enumerate(self.qm_data_points[:]):

                self._check_datapoint_z_matrix(data_point)

                distance, dihedral_dist, denominator, weight_gradient, distance_vec, _, _, _ = self.internal_distance(data_point)

                distances_and_gradients.append((distance, dihedral_dist, dp_index, denominator, weight_gradient, distance_vec))

        else:
            errtxt = "Unrecognized weight function type: "
            errtxt += self.weightfunction_type
            raise ValueError(errtxt)

        close_distances = [
            (self.qm_data_points[index], distance, dihedral_dist, denom, wg, distance_vec, index)
            for distance, dihedral_dist, index, denom, wg, distance_vec in distances_and_gradients]

        for (qm_data_point, distance, dihedral_dist, denominator_cart, weight_grad_cart, _, label_idx) in close_distances:

            weight_cart = 1.0 / (denominator_cart)

            potential, gradient_mw, _ = self.compute_potential(
                qm_data_point,
                self.impes_coordinate.internal_coordinates_values,
            )

            gradient = gradient_mw
            if mass_weighted_b_matrix:
                gradient = sqrt_masses * gradient_mw.reshape(-1)
                gradient = gradient.reshape(gradient_mw.shape)
            gradients.append(gradient)

            potentials.append(potential)
            used_labels.append(label_idx)
            weights_cart.append(weight_cart)
            self.potentials.append(potential)
            self.gradients.append(gradient)
            weight_gradients_cart.append(weight_grad_cart)
            weight_fallback_metrics.append(float(abs(distance)) if np.isfinite(distance) else np.inf)

        # --- initialise accumulators -------------------------------------------------
        self.impes_coordinate.energy = 0.0
        self.impes_coordinate.gradient = np.zeros((natms, 3))
        self.impes_coordinate.NAC = np.zeros((natms, 3))       # if you need it

        # --- 1.  raw (unnormalised) weights and their gradients ----------------------
        w_i = np.array(weights_cart, dtype=np.float64)
        grad_w_i = np.array(weight_gradients_cart, dtype=np.float64)   # shape (n_pts, natms, 3)

        # --- 2. Normalize raw weights -----------------------------------------------
        S_raw = float(np.sum(w_i))
        sum_grad_w_raw = grad_w_i.sum(axis=0)

        eps = float(getattr(self, "weight_eps", 1.0e-12))
        finite_positive = np.isfinite(w_i) & (w_i > 0.0)

        used_weight_fallback = False
        used_weight_rescale = False
        weight_rescale_factor = 1.0

        if np.any(finite_positive):
            # Scale-invariant Shepard normalization.
            # This preserves relative locality even when all raw weights are tiny.
            max_w = float(np.max(w_i[finite_positive]))

            if np.isfinite(S_raw) and S_raw > eps:
                w_norm = w_i
                grad_w_norm = grad_w_i
            else:
                used_weight_rescale = True
                weight_rescale_factor = max_w
                w_norm = np.where(finite_positive, w_i / max_w, 0.0)
                grad_w_norm = np.where(np.isfinite(w_i)[:, None, None],
                                       grad_w_i / max_w, 0.0)

            S_norm = float(np.sum(w_norm))
            if (not np.isfinite(S_norm)) or S_norm <= 0.0:
                raise RuntimeError(
                    "InterpolationDriver.shepard_interpolation: failed to normalize "
                    "finite positive Shepard weights."
                )

            sum_grad_w_norm = grad_w_norm.sum(axis=0)
            W_i = w_norm / S_norm

            grad_W_i = (
                grad_w_norm * S_norm
                - w_norm[:, None, None] * sum_grad_w_norm[None, :, :]
            ) / (S_norm * S_norm)

        else:
            # No finite positive weight (every u > 1e6, or NaN weights): the
            # nearest datapoint is used with W = 1 and a zero weight gradient.
            # E is discontinuous here; this is reached only far outside the
            # database.
            used_weight_fallback = True

            metric = np.asarray(weight_fallback_metrics, dtype=np.float64)
            finite_metric = np.isfinite(metric)

            if not np.any(finite_metric):
                raise RuntimeError(
                    "InterpolationDriver.shepard_interpolation: no finite positive "
                    "Shepard weights and no finite distances. The structure is outside "
                    "the interpolation domain."
                )

            nearest = int(np.argmin(np.where(finite_metric, metric, np.inf)))

            W_i = np.zeros_like(w_i, dtype=np.float64)
            W_i[nearest] = 1.0
            grad_W_i = np.zeros_like(grad_w_i)

            S_norm = 1.0
            sum_grad_w_norm = np.zeros_like(sum_grad_w_raw)

        self.raw_sum_of_weights = S_raw
        self.sum_of_weights = S_raw
        self.sum_of_weights_grad = sum_grad_w_raw
        self.effective_sum_of_weights = S_norm
        self.effective_sum_of_weights_grad = sum_grad_w_norm
        self.weight_rescale_factor = weight_rescale_factor
        self.used_weight_rescale = bool(used_weight_rescale)
        self.raw_weights_array = w_i.copy()
        self.normalized_weights_array = W_i.copy()
        self.used_weight_fallback = bool(used_weight_fallback)
        self.weight_eps = eps

        # --- 3.  accumulate energy and gradient --------------------------------------
        potentials = np.array(potentials, dtype=np.float64)        # Uᵢ
        gradients = np.array(gradients, dtype=np.float64)        # ∇Uᵢ  shape (n_pts, natms, 3)

        self.impes_coordinate.energy = np.dot(W_i, potentials)     # Σ Wᵢ Uᵢ

        self.impes_coordinate.gradient = (np.tensordot(W_i, gradients, axes=1) + np.tensordot(potentials, grad_W_i, axes=1))

        # --- 4.  book-keeping (optional) ---------------------------------------------
        for lbl, Wi in zip(used_labels, W_i):
            self.weights[lbl] = Wi

    def _check_datapoint_z_matrix(self, data_point):
        """
        Raises if the z-matrix of a datapoint differs from the z-matrix of
        impes_coordinate. The internal-coordinate differences are taken
        element by element, so both must define the same coordinates in the
        same order. The result is cached per datapoint object (and per
        z-matrix objects), so repeated calls cost a dictionary lookup.

        :param data_point:
            The InterpolationDatapoint.
        """

        driver_z_matrix = self.impes_coordinate.z_matrix
        dp_z_matrix = getattr(data_point, 'z_matrix', None)

        try:
            checked = self._z_matrix_checked_datapoints.get(data_point)
        except TypeError:
            checked = None
        if (checked is not None and checked[0] is dp_z_matrix and
                checked[1] is driver_z_matrix):
            return

        assert_msg_critical(
            dp_z_matrix is not None and driver_z_matrix is not None,
            'InterpolationDriver: the datapoint or the driver has no z-matrix.')

        ref = [tuple(int(x) for x in z) for z in driver_z_matrix]
        other = [tuple(int(x) for x in z) for z in dp_z_matrix]

        if other != ref:
            label = getattr(data_point, 'point_label', None)
            msg = 'InterpolationDriver: the z-matrix of datapoint '
            msg += f'{label} does not match the z-matrix of the driver '
            msg += f'({len(other)} vs {len(ref)} internal coordinates'
            first = next((i for i, (a, b) in enumerate(zip(other, ref))
                          if a != b), None)
            if first is not None:
                msg += f'; first difference at index {first}: '
                msg += f'{other[first]} vs {ref[first]}'
            msg += '). The database was probably built with other z-matrix '
            msg += 'rules; use the z-matrix stored in the database.'
            assert_msg_critical(False, msg)

        try:
            self._z_matrix_checked_datapoints[data_point] = (
                dp_z_matrix, driver_z_matrix)
        except TypeError:
            pass

    def read_qm_data_points(self):
        """
        Reads the datapoints of labels (all labels of imforcefield_file if
        labels is not set) with the settings of the driver and the
        inv_sqrt_masses of impes_coordinate.

        :return:
            The list of InterpolationDatapoint objects (also stored in
            qm_data_points).
        """

        qm_data_points = []

        assert_msg_critical(
            self.imforcefield_file is not None,
            'ImpesDriver: Please provide a chekpoint file name.')

        if not self.labels:
            labels, z_matrix = self.read_labels()
            self.labels = labels
            self.z_matrix = z_matrix
        # InterpolationDatapoint expects the z-matrix dictionary, not the
        # flattened list stored in impes_coordinate.z_matrix.
        z_matrix = self.impes_coordinate.z_matrix_dict

        for label in self.labels:
            data_point = InterpolationDatapoint(z_matrix)
            data_point.update_settings(self.impes_dict)
            data_point.read_hdf5(self.imforcefield_file, label)
            data_point.inv_sqrt_masses = self.impes_coordinate.inv_sqrt_masses
            qm_data_points.append(data_point)

        self.qm_data_points = qm_data_points

        return qm_data_points

    def compute_potential(self, data_point, org_int_coords):
        """
        Evaluates the second-order Taylor model T_i of a datapoint.

        :param data_point:
            The InterpolationDatapoint.
        :param org_int_coords:
            The internal coordinates at which the model is evaluated
            (normally impes_coordinate.internal_coordinates_values).

        :return:
            A tuple (T_i, grad T_i, c * (g_i + H_i s)): the energy in
            Hartree; its gradient B^T [c * (g_i + H_i s)] with the B-matrix of
            impes_coordinate, (natoms, 3), with respect to mass-weighted
            coordinates if impes_coordinate.inv_sqrt_masses is set; and the
            chain-scaled internal gradient. s is the torsion-safe
            displacement of the class docstring and c its chain factor
            (cos(Delta/2) for proper dihedrals, 1/cos^2(Delta/2) for
            impropers, 1 otherwise). The L2 norms of the bond, angle and
            proper-dihedral parts of s are appended to bond_rmsd,
            angle_rmsd and dihedral_rmsd.
        """

        natm = data_point.cartesian_coordinates.shape[0]

        bounds = self._get_internal_coordinate_partitions()
        bond_end = int(bounds['bond_end'])
        angle_end = int(bounds['angle_end'])
        dihedral_start = int(bounds['dihedral_start'])
        dihedral_end = int(bounds['dihedral_end'])

        energy = data_point.energy
        grad = data_point.internal_gradient.copy()
        hessian = data_point.internal_hessian.copy()
        dist_org = (org_int_coords.copy() - data_point.internal_coordinates_values)
        dist_check = (org_int_coords.copy() - data_point.internal_coordinates_values)
        chain = np.ones_like(dist_org)

        d_prop = self._principal_torsion_delta(dist_org[dihedral_start:dihedral_end])
        dist_check[dihedral_start:dihedral_end] = 2.0 * np.sin(0.5 * d_prop)
        d_imp = self._principal_torsion_delta(
            dist_org[dihedral_end:]
        )
        chain[dihedral_start:dihedral_end] = np.cos(0.5 * d_prop)
        imp_slice = slice(dihedral_end, len(dist_org))
        dist_check[imp_slice] = 2.0 * np.tan(0.5 * d_imp)
        chain[imp_slice] = 1.0 / np.maximum(np.cos(0.5 * d_imp)**2, 1.0e-12)

        self.bond_rmsd.append(np.sqrt(np.mean(np.sum((dist_org[:bond_end])**2))))
        self.angle_rmsd.append(np.sqrt(np.mean(np.sum(dist_org[bond_end:angle_end]**2))))
        self.dihedral_rmsd.append(np.sqrt(np.mean(np.sum(dist_check[dihedral_start:dihedral_end]**2))))

        pes = (
            energy
            + np.matmul(dist_check.T, grad)
            + 0.5 * np.linalg.multi_dot([dist_check.T, hessian, dist_check])
        )

        dist_hessian_eff = np.matmul(dist_check.T, hessian)
        grad *= chain
        dist_hessian_eff *= chain

        pes_prime = np.matmul(self.impes_coordinate.b_matrix.T, (grad + dist_hessian_eff)).reshape(natm, 3)
        return pes, pes_prime, (grad + dist_hessian_eff)

    def _iter_imp_internal_coordinate_rows(self, data_point):
        """
        Yields (section, coordinate, z-matrix index) for the important
        coordinates of data_point that are in the z-matrix of
        impes_coordinate; the others are skipped.
        """

        imp = getattr(data_point, "imp_int_coordinates", None)
        if not isinstance(imp, dict):
            return

        zmat = [tuple(int(x) for x in z) for z in self.impes_coordinate.z_matrix]

        for section in ("bonds", "angles", "dihedrals", "impropers"):
            for coord_raw in imp.get(section, []):
                coord = tuple(int(x) for x in coord_raw)

                try:
                    idx = zmat.index(coord)
                except ValueError:
                    # Important coordinate is absent from the active z-matrix.
                    # This should not happen for a consistent database, but skip
                    # defensively rather than killing interpolation.
                    continue

                yield section, coord, idx

    def _imp_coordinate_sigma(self, section, idx, q_ref):
        """
        Returns the tolerance sigma of important coordinate idx in the
        representation of the internal coordinate row idx.

        - bonds: sigma_r = tc_imp_bond_sigma_angstrom in bohr. For inverse
          bonds (q = 1/r) it is mapped to q-space, sigma_q = sigma_r / r_ref^2.
          For eq bonds (q = r_eq - r_eq^2/r) and plain r it is used
          unchanged, which is exact for r and for eq bonds at r = r_eq
          (dq/dr = r_eq^2/r^2).
        - angles: tc_imp_angle_sigma_degrees in radians.
        - dihedrals: 2 sin(sigma/2); impropers: 2 tan(sigma/2), matching the
          displacements 2 sin(Delta/2) and 2 tan(Delta/2).

        :param section:
            'bonds', 'angles', 'dihedrals' or 'impropers'.
        :param idx:
            The index of the coordinate in the flat z-matrix.
        :param q_ref:
            The internal coordinates of the datapoint.
        """

        sigma_bond_bohr = (
            float(self.tc_imp_bond_sigma_angstrom) / bohr_in_angstrom()
        )

        if section == "bonds":
            if self.use_inverse_bond_length:
                # q_ref is 1/r in bohr^-1.
                r_ref = 1.0 / max(abs(float(q_ref[idx])), self.tc_imp_min_sigma)
                sigma = sigma_bond_bohr / max(r_ref * r_ref, self.tc_imp_min_sigma)
            else:
                sigma = sigma_bond_bohr

        elif section == "angles":
            sigma = np.deg2rad(float(self.tc_imp_angle_sigma_degrees))

        elif section == "dihedrals":
            sigma = 2.0 * np.sin(np.deg2rad(0.5 * float(self.tc_imp_dihedral_sigma_degrees)))

        elif section == "impropers":
            sigma = 2.0 * np.tan(0.5 * np.deg2rad(float(self.tc_imp_improper_sigma_degrees)))

        else:
            raise ValueError(f"Unknown important-coordinate section: {section}")

        return max(float(abs(sigma)), self.tc_imp_min_sigma)

    def _important_coordinate_gate_metric(
        self,
        data_point,
        active_dofs,
    ):
        """
        Returns (A, grad A) of the multiplicative tc gate: A = mean_k z_k^2
        over the important coordinates of data_point that are in the
        z-matrix, with z_k = 2 sin(Delta/2)/sigma (dihedrals),
        2 tan(Delta/2)/sigma (impropers) or Delta q/sigma (bonds, angles),
        and grad A its gradient with respect to the Cartesian coordinates
        (not mass-weighted) of active_dofs. (0, 0) without important
        coordinates.
        """

        q_cur = np.asarray(self.impes_coordinate.internal_coordinates_values, dtype=np.float64)
        q_ref = np.asarray(data_point.internal_coordinates_values, dtype=np.float64)
        B = np.asarray(self.impes_coordinate.b_matrix, dtype=np.float64)

        inv_sqrt_masses = getattr(self.impes_coordinate, "inv_sqrt_masses", None)
        if inv_sqrt_masses is not None:
            inv_sqrt_active = np.asarray(inv_sqrt_masses, dtype=np.float64)[active_dofs]
        else:
            inv_sqrt_active = None

        selected_idx = []
        sigmas = []
        z_values = []
        dzdx_rows = []

        for section, coord, idx in self._iter_imp_internal_coordinate_rows(data_point):
            dq_raw = float(q_cur[idx] - q_ref[idx])
            sigma = self._imp_coordinate_sigma(section, idx, q_ref)

            b_row = B[idx, active_dofs].copy()
            if inv_sqrt_active is not None:
                b_row = b_row / inv_sqrt_active

            if section == "dihedrals":
                d = self._principal_torsion_delta(dq_raw)
                z = 2.0 * np.sin(0.5 * d) / sigma
                dz_dx = np.cos(0.5 * d) * b_row / sigma

            elif section == "impropers":
                d = self._principal_torsion_delta(dq_raw)
                eta = 2.0 * np.tan(0.5 * d)
                cos_half = np.cos(0.5 * d)
                chain_eta = 1.0 / np.maximum(cos_half * cos_half, 1.0e-12)
                z = eta / sigma
                dz_dx = chain_eta * b_row / sigma

            elif section in ("angles", "bonds"):
                z = dq_raw / sigma
                dz_dx = b_row / sigma

            else:
                continue

            selected_idx.append(int(idx))
            sigmas.append(float(sigma))
            z_values.append(float(z))
            dzdx_rows.append(np.asarray(dz_dx, dtype=np.float64).reshape(-1))

        nsel = len(z_values)
        if nsel == 0:
            return 0.0, np.zeros(active_dofs.size, dtype=np.float64)

        z = np.asarray(z_values, dtype=np.float64)
        J = np.vstack(dzdx_rows)  # shape: nsel x n_active_dofs

        A_iso = float(np.dot(z, z) / nsel)
        grad_iso = (2.0 / nsel) * (J.T @ z)

        return A_iso, grad_iso

    def _apply_imp_coordinate_penalty_gate(
            self,
            denominator_base,
            raw_weight_gradient_base,
            A_imp,
            grad_A_imp,
    ):
        """
        Multiplies the raw weight 1/denominator_base by the gate
        exp(-lambda A), lambda = tc_imp_gate_lambda. lambda A is clipped at
        tc_imp_gate_exp_clip, where the derivative of the gate is set to
        zero.

        :return:
            (1/w_mod, grad w_mod); (inf, 0) for a non-finite or non-positive
            base denominator or w_mod <= 1e-300.
        """

        if not np.isfinite(denominator_base) or denominator_base <= 0.0:
            return np.inf, np.zeros_like(raw_weight_gradient_base)

        lam = float(self.tc_imp_gate_lambda)

        raw_exponent = lam * float(A_imp)
        clip = float(self.tc_imp_gate_exp_clip)

        if raw_exponent >= clip:
            exponent = clip
            lambda_eff = 0.0
        else:
            exponent = raw_exponent
            lambda_eff = lam

        gate = np.exp(-exponent)
        w_base = 1.0 / denominator_base

        grad_base_flat = np.asarray(raw_weight_gradient_base, dtype=np.float64).reshape(-1)
        grad_A_flat = np.asarray(grad_A_imp, dtype=np.float64).reshape(-1)

        grad_gate = -lambda_eff * gate * grad_A_flat

        w_mod = w_base * gate
        grad_w_mod = gate * grad_base_flat + w_base * grad_gate

        if w_mod <= 1.0e-300:
            return np.inf, np.zeros_like(raw_weight_gradient_base)

        denominator_mod = 1.0 / w_mod

        return denominator_mod, grad_w_mod.reshape(raw_weight_gradient_base.shape)

    def _target_customized_additive_distance_term(self, data_point, active_dofs):
        """
        Computes the additive target-customized weighting-coordinate term
        following the Kim/Rhee formalism.

        The modified squared distance is D_total = D_cart + D_tc with
        D_tc = sum_k delta_k^2 over the important coordinates of data_point
        (delta_k = Delta q / sigma for bonds and angles,
        2 sin(Delta/2) / sigma for dihedrals and impropers).

        :param data_point:
            The InterpolationDatapoint.
        :param active_dofs:
            The Cartesian components of the atoms in the weights.

        :return:
            (D_tc, grad_D_tc): the dimensionless contribution and its
            derivative d(D_tc)/dX over active_dofs (not mass-weighted).
        """

        q_cur = np.asarray(
            self.impes_coordinate.internal_coordinates_values,
            dtype=np.float64,
        )

        q_ref = np.asarray(
            data_point.internal_coordinates_values,
            dtype=np.float64,
        )

        B = np.asarray(
            self.impes_coordinate.b_matrix,
            dtype=np.float64,
        )

        inv_sqrt_masses = getattr(
            self.impes_coordinate,
            "inv_sqrt_masses",
            None,
        )

        if inv_sqrt_masses is not None:
            inv_sqrt_active = np.asarray(
                inv_sqrt_masses,
                dtype=np.float64,
            )[active_dofs]
        else:
            inv_sqrt_active = None

        D_tc = 0.0
        grad_D_tc = np.zeros(active_dofs.size, dtype=np.float64)

        for section, coord, idx in self._iter_imp_internal_coordinate_rows(
            data_point
        ):
            dq_raw = float(q_cur[idx] - q_ref[idx])

            b_row = B[idx, active_dofs].copy()

            # Convert dq/dQ back to dq/dX if the B-matrix is expressed
            # with respect to mass-weighted Cartesian coordinates.
            if inv_sqrt_active is not None:
                b_row = b_row / inv_sqrt_active

            if section == "bonds":
                sigma = self._rhee_coordinate_sigma(section, idx, q_ref)

                delta = dq_raw / sigma
                ddelta_dx = b_row / sigma

            elif section == "angles":
                sigma = self._rhee_coordinate_sigma(section, idx, q_ref)

                delta = dq_raw / sigma
                ddelta_dx = b_row / sigma

            elif section in ("dihedrals", "impropers"):
                sigma = self._rhee_coordinate_sigma(section, idx, q_ref)

                d = self._principal_torsion_delta(dq_raw)

                delta = 2.0 * np.sin(0.5 * d) / sigma

                ddelta_dx = (
                    np.cos(0.5 * d) * b_row / sigma
                )
            else:
                continue

            D_tc += delta * delta
            grad_D_tc += 2.0 * delta * ddelta_dx

        return float(D_tc), grad_D_tc

    def _rhee_coordinate_sigma(self, section, idx, q_ref):
        """
        Returns the scaling parameter for the additive Kim/Rhee
        target-customized weighting-coordinate formalism.

        For bonds:
            sigma_r = tc_imp_bond_sigma_angstrom in bohr. For inverse bonds
            (q = 1/r) it is mapped to q-space, sigma_q = sigma_r / r_ref^2;
            for eq bonds and plain r it is used unchanged (exact for r, and
            for eq bonds at r = r_eq).

        For angles:
            sigma is in radians.

        For proper and improper dihedrals:
            sigma = 2 sin(tau / 2), so that

                delta = 2 sin(Delta / 2) / sigma
                    = sin(Delta / 2) / sin(tau / 2).

            This is the weighting-coordinate convention of Kim and Rhee.
        """

        sigma_bond_bohr = (
            float(self.tc_imp_bond_sigma_angstrom) / bohr_in_angstrom()
        )

        if section == "bonds":
            if self.use_inverse_bond_length:
                r_ref = 1.0 / max(
                    abs(float(q_ref[idx])),
                    self.tc_imp_min_sigma,
                )

                sigma = sigma_bond_bohr / max(
                    r_ref * r_ref,
                    self.tc_imp_min_sigma,
                )
            else:
                sigma = sigma_bond_bohr

        elif section == "angles":
            sigma = np.deg2rad(
                float(self.tc_imp_angle_sigma_degrees)
            )

        elif section == "dihedrals":
            sigma = 2.0 * np.sin(
                0.5 * np.deg2rad(
                    float(self.tc_imp_dihedral_sigma_degrees)
                )
            )

        elif section == "impropers":
            sigma = 2.0 * np.sin(
                0.5 * np.deg2rad(
                    float(self.tc_imp_improper_sigma_degrees)
                )
            )

        else:
            raise ValueError(
                f"Unknown important-coordinate section: {section}"
            )

        return max(float(abs(sigma)), self.tc_imp_min_sigma)

    def trust_radius_weight_gradient(self, confidence_radius, distance):
        """
        Returns dw/dR of the raw weight w = 1/(u^p + u^q), u = d^2/R^2:

            dw/dR = (2/R) (p u^p + q u^q) / (u^p + u^q)^2   (0 for u > 1e6).

        :param confidence_radius:
            The confidence radius R in bohr.
        :param distance:
            The distance d in bohr.
        """

        R = confidence_radius
        D = distance**2

        p = float(self.exponent_p)
        q = float(self.exponent_q)

        # Dimensionless ratio u = (d/R)^2
        u = D / (R**2)

        # Early Exit: If the distance is massively outside the confidence radius,
        # the trust weight gradient flattens to absolute zero.
        if np.isscalar(u) and u > 1e6:
            return 0.0

        denom_base = u**p + u**q
        denom_base = np.clip(denom_base, a_min=1e-300, a_max=None)

        # Algebraically simplified numerator
        numerator = p * u**p + q * u**q

        trust_radius_weight_gradient = (2.0 / R) * (numerator / denom_base**2)

        return trust_radius_weight_gradient

    def trust_radius_weight_gradient_gradient(self, confidence_radius, distance, distance_vector, scale2):
        """
        Returns d(grad w)/dR of the raw weight (with d^2 = scale2 |v|^2):

            (4 scale2 v / R^3) [(p^2 u^(p-1) + q^2 u^(q-1)) / den^2
                                - 2 (p u^(p-1) + q u^(q-1)) (p u^p + q u^q) / den^3],

        den = u^p + u^q; zero for u > 1e6.

        :param confidence_radius:
            The confidence radius R in bohr.
        :param distance:
            The distance d in bohr.
        :param distance_vector:
            The aligned displacement v of the atoms in the weights.
        :param scale2:
            The distance scaling factor (1.0).
        """

        R = confidence_radius
        D = distance**2
        v = distance_vector

        p = float(self.exponent_p)
        q = float(self.exponent_q)

        u = D / (R**2)

        # Early Exit
        if np.isscalar(u) and u > 1e6:
            return np.zeros_like(v, dtype=np.float64)

        denom_base = u**p + u**q
        denom_base = np.clip(denom_base, a_min=1e-300, a_max=None)

        # Factored numerator components (the constants and R factors are moved to the multiplier)
        term1 = (p**2 * u**(p - 1.0) + q**2 * u**(q - 1.0)) / (denom_base**2)
        term2 = 2.0 * ((p * u**(p - 1.0) + q * u**(q - 1.0)) * (p * u**p + q * u**q)) / (denom_base**3)

        # Common scaling factor analytically pulled out of the original expressions
        multiplier = (4.0 * scale2 * v) / (R**3)

        trust_radius_weight_gradient_gradient = multiplier * (term1 - term2)

        return trust_radius_weight_gradient_gradient

    def shepard_weight_gradient(self, distance_vector, distance, scale2, confidence_radius):
        """
        Returns (den, grad w) of the raw weight w = 1/den, den = u^p + u^q,
        u = d^2/R^2 (den clipped at 1e-300):

            grad w = -2 scale2 (p u^(p-1) + q u^(q-1)) / (R^2 den^2) v;

        (inf, 0) for u > 1e6.

        :param distance_vector:
            The aligned displacement v of the atoms in the weights.
        :param distance:
            The distance d in bohr.
        :param scale2:
            The distance scaling factor (1.0).
        :param confidence_radius:
            The confidence radius R in bohr.
        """

        R = confidence_radius
        D = distance**2
        v = distance_vector

        p = float(self.exponent_p)
        q = float(self.exponent_q)

        u = D / (R**2)

        # Early Exit: Return np.inf for the denominator and an array of zeros for the gradient
        if np.isscalar(u) and u > 1e6:
            return np.inf, np.zeros_like(v, dtype=np.float64)

        denom_base = u**p + u**q
        denom_base = np.clip(denom_base, a_min=1e-300, a_max=None)

        numerator = p * u**(p - 1.0) + q * u**(q - 1.0)

        weight_gradient = -2.0 * scale2 * (numerator / (R**2 * denom_base**2)) * v

        return denom_base, weight_gradient

    def trust_radius_tc_weight_gradient(self, confidence_radius, distance, imp_int_coords_distance):
        """
        Returns dw/dR for u = (d^2 + D_extra)/R^2 with the formula of
        trust_radius_weight_gradient. Used for the additive Rhee weights
        (D_extra = D_tc) and for the internal weights. For internal weights
        with use_tc_weights the factor exp(-alpha D_imp) of the weight is not
        included (known issue K12a).

        :param confidence_radius:
            The confidence radius R in bohr.
        :param distance:
            The distance d in bohr.
        :param imp_int_coords_distance:
            The additional squared distance D_extra.
        """

        R = confidence_radius
        # Incorporating the internal coordinates distance
        D = distance**2 + imp_int_coords_distance

        p = float(self.exponent_p)
        q = float(self.exponent_q)

        u = D / (R**2)

        # Early Exit
        if np.isscalar(u) and u > 1e6:
            return 0.0

        denom_base = u**p + u**q
        denom_base = np.clip(denom_base, a_min=1e-300, a_max=None)

        numerator = p * u**p + q * u**q

        trust_radius_weight_gradient = (2.0 / R) * (numerator / denom_base**2)

        return trust_radius_weight_gradient

    def trust_radius_tc_weight_gradient_gradient(self, confidence_radius, distance, distance_vector, imp_int_coords_distance, imp_int_coords_dist_derivative, scale2=1.0,):
        """
        Returns d(grad w)/dR for u = (d^2 + D_extra)/R^2 with
        grad(d^2 + D_extra) = 2 scale2 v + grad D_extra, analogous to
        trust_radius_weight_gradient_gradient; zero for u > 1e6.

        :param confidence_radius:
            The confidence radius R in bohr.
        :param distance:
            The distance d in bohr.
        :param distance_vector:
            The aligned displacement v of the atoms in the weights.
        :param imp_int_coords_distance:
            The additional squared distance D_extra.
        :param imp_int_coords_dist_derivative:
            The Cartesian gradient of D_extra (flat).
        :param scale2:
            The distance scaling factor (1.0).
        """

        # 1. Define base variables
        R = confidence_radius
        D = distance**2 + imp_int_coords_distance
        v = imp_int_coords_dist_derivative + 2.0 * scale2 * distance_vector.reshape(-1)

        p = float(self.exponent_p)
        q = float(self.exponent_q)

        # 2. Compute the dimensionless ratio
        u = D / (R**2)

        # Early exit: for u > 1e6 the weight is zero, which also avoids the
        # float64 overflow of u**(3 * max(p, q)).
        if np.isscalar(u) and u > 1e6:
            return np.zeros((distance_vector.shape[0], 3), dtype=np.float64)

        # 3. Base denominator
        denom_base = u**p + u**q
        denom_base = np.clip(denom_base, a_min=1e-300, a_max=None)

        # 4. Compute the grouped terms
        term1 = (p**2 * u**(p - 1.0) + q**2 * u**(q - 1.0)) / (denom_base**2)
        term2 = 2.0 * ((p * u**(p - 1.0) + q * u**(q - 1.0)) * (p * u**p + q * u**q)) / (denom_base**3)

        # 5. Apply the common multiplier
        multiplier = (2.0 * v) / (R**3)

        trust_radius_weight_gradient_gradient = multiplier * (term1 - term2)

        return trust_radius_weight_gradient_gradient.reshape(distance_vector.shape[0], 3)

    def additive_rhee_shepard_weight_gradient(
        self,
        distance_vector,
        distance,
        scale2,
        confidence_radius,
        D_tc,
        grad_D_tc,
    ):
        """
        Computes the raw two-part Shepard weight denominator and its
        Cartesian derivative for the additive Kim/Rhee target-customized
        weighting-distance formalism:

            D_total = d_cart^2 + D_tc,   u = D_total / R^2,
            w = 1 / (u^p + u^q).

        :param distance_vector:
            The aligned Cartesian displacement, (n_active_atoms, 3).
        :param distance:
            The Cartesian distance, distance^2 = scale2 |distance_vector|^2.
        :param scale2:
            The Cartesian distance scaling factor.
        :param confidence_radius:
            The confidence radius R of the datapoint in bohr.
        :param D_tc:
            The additive customized-coordinate distance contribution.
        :param grad_D_tc:
            The Cartesian derivative of D_tc, (n_active_dofs,).

        :return:
            (den, grad w): the denominator u^p + u^q (inf for u > 1e6) and
            the gradient of w, (n_active_atoms, 3).
        """

        R = float(confidence_radius)
        p = float(self.exponent_p)
        q = float(self.exponent_q)

        D_cart = float(distance * distance)
        D_total = D_cart + float(D_tc)

        grad_D_cart = (
            2.0 * float(scale2) * distance_vector.reshape(-1)
        )

        grad_D_total = (
            grad_D_cart
            + np.asarray(grad_D_tc, dtype=np.float64).reshape(-1)
        )

        u = D_total / (R * R)

        if np.isscalar(u) and u > 1.0e6:
            return (
                np.inf,
                np.zeros_like(distance_vector, dtype=np.float64),
            )

        denominator = u**p + u**q
        denominator = np.clip(
            denominator,
            a_min=1.0e-300,
            a_max=None,
        )

        numerator = (
            p * u**(p - 1.0)
            + q * u**(q - 1.0)
        )

        prefactor = -numerator / (
            (R * R) * denominator**2
        )

        weight_gradient = prefactor * grad_D_total

        return denominator, weight_gradient.reshape(distance_vector.shape)

    def calculate_translation_coordinates(self, given_coordinates):
        """
        Returns the coordinates translated to their geometric centre.

        :param given_coordinates:
            The coordinates, (natoms, 3).
        """

        center = np.mean(given_coordinates, axis=0)
        translated_coordinates = given_coordinates - center

        return translated_coordinates

    def internal_distance(self, data_point, store_alpha_gradients=True):
        """
        Computes the raw Shepard weight of a datapoint from an
        internal-coordinate distance (weightfunction_type='internal').

        D_int = |s|^2 with s = Delta q for bonds and angles, sin(Delta) for
        proper dihedrals and 2 tan(Delta/2) for impropers (floored at 1e-8);
        u = D_int / R^2. With use_tc_weights the weight is multiplied by
        exp(-alpha D_imp), D_imp built from the important coordinates scaled
        by the diagonal curvature of the datapoint Hessian.

        :param data_point:
            The InterpolationDatapoint.
        :param store_alpha_gradients:
            If calc_optim_trust_radius is True, append dw/dR and
            d(grad w)/dR to dw_dalpha_list and dw_dX_dalpha_list.

        :return:
            An 8-tuple (sqrt(D_int), 0.0, 1/w, grad w, zeros, zeros, dw/dR,
            d(grad w)/dR): the distance, an unused placeholder, the
            denominator of the raw weight, its Cartesian gradient
            (natoms, 3; zero for excluded atoms), two unused placeholders
            (arrays of the active atoms) and the radius derivatives (0
            unless calc_optim_trust_radius).
        """

        reference_coordinates = self.impes_coordinate.cartesian_coordinates.copy()
        natms = reference_coordinates.shape[0]

        active_atoms = np.delete(np.arange(natms), self.symmetry_information[4])
        active_dofs = np.array([3 * a + i for a in active_atoms for i in range(3)], dtype=np.int64)

        q_current = np.asarray(self.impes_coordinate.internal_coordinates_values, dtype=np.float64)
        q_ref = np.asarray(data_point.internal_coordinates_values, dtype=np.float64)

        dq_raw = q_current - q_ref
        dq_eff = dq_raw.copy()
        chain = np.ones_like(dq_raw)

        bounds = self._get_internal_coordinate_partitions()
        dstart = int(bounds["dihedral_start"])
        dend = int(bounds["dihedral_end"])

        for idx in range(len(dq_raw)):
            if dstart <= idx < dend:  # proper dihedrals only
                d = self._principal_torsion_delta(dq_raw[idx])
                dq_eff[idx] = np.sin(d)
                chain[idx] = np.cos(d)
            elif idx >= dend:  # impropers
                d = self._principal_torsion_delta(dq_raw[idx])
                dq_eff[idx] = self._improper_dihedral_displacement(d)
                chain[idx] = self._improper_dihedral_chain(d)
            else:  # bonds and angles
                dq_eff[idx] = dq_raw[idx]
                chain[idx] = 1.0

        grad_Dint = np.zeros(active_dofs.shape[0], dtype=np.float64)
        if self.impes_coordinate.inv_sqrt_masses is not None:
            inv_sqrt_masses = self.impes_coordinate.inv_sqrt_masses[active_dofs]
        else:
            inv_sqrt_masses = np.ones(active_dofs.shape[0], dtype=np.float64)

        for idx in range(len(dq_eff)):
            b_row = self.impes_coordinate.b_matrix[idx, active_dofs] / inv_sqrt_masses
            grad_Dint += 2.0 * dq_eff[idx] * chain[idx] * b_row

        Dint_sq = float(np.dot(dq_eff, dq_eff))
        if Dint_sq < 1e-8:
            Dint_sq = 1e-8
            grad_Dint[:] = 0.0

        distance = np.sqrt(Dint_sq)
        dihedral_dist = 0.0
        distance_vector_sub = np.zeros((active_atoms.shape[0], 3), dtype=np.float64)
        grad_s = np.zeros_like(distance_vector_sub)

        Dimp_sq = 0.0
        dq_dx = np.zeros_like(grad_Dint)

        if self.use_tc_weights and isinstance(data_point.imp_int_coordinates, dict):
            H_eff = data_point.internal_hessian
            eigval, V = np.linalg.eigh(H_eff)
            coord_curvature = np.sum(eigval * (V**2), axis=1)

            for _, entries in data_point.imp_int_coordinates.items():
                for element in entries:
                    idx = self.impes_coordinate.z_matrix.index(element)
                    delta_q = (
                        self.impes_coordinate.internal_coordinates_values[idx] -
                        data_point.internal_coordinates_values[idx]
                    )
                    b_row = self.impes_coordinate.b_matrix[idx, active_dofs] / inv_sqrt_masses

                    if len(element) == 4:
                        arg = coord_curvature[idx] * delta_q
                        dq_dx += (
                            coord_curvature[idx] *
                            2.0 *
                            b_row *
                            np.cos(arg) *
                            np.sin(arg)
                        )
                        Dimp_sq += np.sin(arg)**2
                    else:
                        dq_dx += coord_curvature[idx]**2 * 2.0 * b_row * delta_q
                        Dimp_sq += (coord_curvature[idx] * delta_q)**2

            if Dimp_sq < 1e-8:
                Dimp_sq = 1e-8
                dq_dx[:] = 0.0

        R = data_point.confidence_radius
        p = float(self.exponent_p)
        q = float(self.exponent_q)
        u = Dint_sq / (R**2)

        if np.isscalar(u) and u > 1e6:
            denom_base = np.inf
            base_weight_grad = np.zeros_like(grad_Dint, dtype=np.float64)
        else:
            denom_base = u**p + u**q
            denom_base = np.clip(denom_base, a_min=1e-300, a_max=None)
            numerator = p * u**(p - 1.0) + q * u**(q - 1.0)
            base_weight_grad = -1.0 * (numerator / (R**2 * denom_base**2)) * grad_Dint

        if self.use_tc_weights:
            exp_decay = np.exp(-self.alpha * Dimp_sq)
            denominator = denom_base * np.exp(self.alpha * Dimp_sq)
            weight_grad_sub = (
                base_weight_grad * exp_decay -
                (self.alpha * exp_decay / denom_base) * dq_dx
            )
        else:
            denominator = denom_base
            weight_grad_sub = base_weight_grad

        weight_gradient_sub = weight_grad_sub.reshape(active_atoms.shape[0], 3)
        weight_gradient = np.zeros_like(reference_coordinates)
        weight_gradient[self.symmetry_information[3]] = weight_gradient_sub

        dw_dalpha_i = 0.0
        dw_dX_dalpha_i = np.zeros_like(weight_gradient_sub)

        if self.calc_optim_trust_radius:
            if self.use_tc_weights:
                dw_dalpha_i = self.trust_radius_tc_weight_gradient(
                    data_point.confidence_radius,
                    distance,
                    Dimp_sq)
                dw_dX_dalpha_i = self.trust_radius_tc_weight_gradient_gradient(
                    data_point.confidence_radius,
                    distance,
                    distance_vector_sub,
                    Dimp_sq,
                    dq_dx)
            else:
                dw_dalpha_i = self.trust_radius_tc_weight_gradient(
                    data_point.confidence_radius,
                    0.0,
                    Dint_sq)
                dw_dX_dalpha_i = self.trust_radius_tc_weight_gradient_gradient(
                    data_point.confidence_radius,
                    0.0,
                    distance_vector_sub,
                    Dint_sq,
                    grad_Dint)

            if store_alpha_gradients:
                self.dw_dalpha_list.append(dw_dalpha_i)
                dw_dX_dalpha_i_full = np.zeros_like(reference_coordinates)
                dw_dX_dalpha_i_full[self.symmetry_information[3]] = dw_dX_dalpha_i
                self.dw_dX_dalpha_list.append(dw_dX_dalpha_i_full)

        return (
            distance,
            dihedral_dist,
            denominator,
            weight_gradient,
            distance_vector_sub,
            grad_s,
            dw_dalpha_i,
            dw_dX_dalpha_i,
        )

    def cartesian_distance(self, data_point, store_alpha_gradients=True):
        """
        Computes the Cartesian distance between the current geometry and a
        datapoint and the raw Shepard weight of the datapoint.

        The atoms in the weights (all atoms except symmetry_information[4])
        of both geometries are centred and the datapoint is rotated onto the
        current geometry (Kabsch); d = |v| with v = X_c - R Y_c in bohr (not
        mass-weighted, floored at 1e-8 with v = 0).

        :param data_point:
            The InterpolationDatapoint.
        :param store_alpha_gradients:
            If calc_optim_trust_radius is True, append dw/dR and
            d(grad w)/dR to dw_dalpha_list and dw_dX_dalpha_list.

        :return:
            An 8-tuple (d, 0.0, 1/w, grad w, v, 1.0, dw/dR, d(grad w)/dR):
            the distance, an unused placeholder, the denominator of the raw
            weight w (inf for w = 0), its Cartesian gradient (natoms, 3; zero
            for excluded atoms), the aligned displacement of the active
            atoms, an unused placeholder and the radius derivatives (the
            latter for the active atoms; 0 unless calc_optim_trust_radius).
        """

        # First, translate the cartesian coordinates to zero
        target_coordinates = data_point.cartesian_coordinates.copy()
        reference_coordinates = self.impes_coordinate.cartesian_coordinates.copy()

        active_atoms = np.delete(np.arange(reference_coordinates.shape[0]), self.symmetry_information[4])
        active_dofs = np.array([3*a + i for a in active_atoms for i in range(3)])

        target_coordinates_core = target_coordinates[active_atoms]
        reference_coordinates_core = reference_coordinates[active_atoms]

        center_target_coordinates_core = self.calculate_translation_coordinates(target_coordinates_core)
        center_reference_coordinates_core = self.calculate_translation_coordinates(reference_coordinates_core)

        grad_s = 1.0
        Xc = center_reference_coordinates_core     # variable (current geometry)
        Yc = center_target_coordinates_core        # fixed (datapoint)
        R_x = geometric.rotate.get_rot(Yc, Xc)

        rotated_current = (R_x @ Yc.T).T
        distance_vector_sub = Xc - rotated_current

        scale2 = 1.0
        distance = np.sqrt(scale2 * np.dot(distance_vector_sub.ravel(), distance_vector_sub.ravel()))

        dihedral_dist = 0.0
        if distance < 1e-8:
            distance = 1e-8
            distance_vector_sub[:] = 0

        dw_dalpha_i = 0
        dw_dX_dalpha_i = 0
        D_tc = None
        grad_D_tc = None
        A_imp = None
        grad_A_imp = None

        if not self.use_tc_weights:
            denominator_imp_coord, weight_gradient_sub_imp_coord = (
                self.shepard_weight_gradient(
                    distance_vector_sub,
                    distance,
                    scale2,
                    data_point.confidence_radius,
                )
            )

        elif self.tc_weight_mode == "multiplicative":
            denominator_base, weight_gradient_sub_base = (
                self.shepard_weight_gradient(
                    distance_vector_sub,
                    distance,
                    scale2,
                    data_point.confidence_radius,
                )
            )

            A_imp, grad_A_imp = self._important_coordinate_gate_metric(
                data_point,
                active_dofs,
            )

            denominator_imp_coord, weight_gradient_sub_imp_coord = (
                self._apply_imp_coordinate_penalty_gate(
                    denominator_base,
                    weight_gradient_sub_base,
                    A_imp,
                    grad_A_imp,
                )
            )

        elif self.tc_weight_mode == "additive_rhee":
            D_tc, grad_D_tc = self._target_customized_additive_distance_term(
                data_point,
                active_dofs,
            )

            denominator_imp_coord, weight_gradient_sub_imp_coord = (
                self.additive_rhee_shepard_weight_gradient(
                    distance_vector=distance_vector_sub,
                    distance=distance,
                    scale2=scale2,
                    confidence_radius=data_point.confidence_radius,
                    D_tc=D_tc,
                    grad_D_tc=grad_D_tc,
                )
            )

        else:
            raise ValueError(
                f"Unsupported tc_weight_mode: {self.tc_weight_mode}"
            )

        if self.calc_optim_trust_radius:

            if not self.use_tc_weights:
                dw_dalpha_i = self.trust_radius_weight_gradient(
                    data_point.confidence_radius,
                    distance,
                )

                dw_dX_dalpha_i = self.trust_radius_weight_gradient_gradient(
                    data_point.confidence_radius,
                    distance,
                    distance_vector_sub,
                    scale2,
                )

            elif self.tc_weight_mode == "multiplicative":
                # the gate metric of this datapoint was computed above
                lam = float(self.tc_imp_gate_lambda)
                gate = np.exp(
                    -min(lam * A_imp, self.tc_imp_gate_exp_clip)
                )

                grad_gate = -lam * gate * grad_A_imp

                dw_base_dR = self.trust_radius_weight_gradient(
                    data_point.confidence_radius,
                    distance,
                )

                dgradw_base_dR = self.trust_radius_weight_gradient_gradient(
                    data_point.confidence_radius,
                    distance,
                    distance_vector_sub,
                    scale2,
                )

                dw_dalpha_i = gate * dw_base_dR

                dw_dX_dalpha_i = (
                    gate * dgradw_base_dR.reshape(-1)
                    + grad_gate * dw_base_dR
                ).reshape(distance_vector_sub.shape)

            elif self.tc_weight_mode == "additive_rhee":
                if D_tc is None or grad_D_tc is None:
                    D_tc, grad_D_tc = self._target_customized_additive_distance_term(
                        data_point,
                        active_dofs,
                    )

                dw_dalpha_i = self.trust_radius_tc_weight_gradient(
                    data_point.confidence_radius,
                    distance,
                    D_tc,
                )

                dw_dX_dalpha_i = (
                    self.trust_radius_tc_weight_gradient_gradient(
                        confidence_radius=data_point.confidence_radius,
                        distance=distance,
                        distance_vector=distance_vector_sub,
                        imp_int_coords_distance=D_tc,
                        imp_int_coords_dist_derivative=grad_D_tc,
                        scale2=scale2,
                    )
                )

            else:
                raise ValueError(
                    f"Unsupported tc_weight_mode: {self.tc_weight_mode}"
                )

            if store_alpha_gradients:
                self.dw_dalpha_list.append(dw_dalpha_i)

                dw_dX_dalpha_i_full = np.zeros_like(reference_coordinates)
                dw_dX_dalpha_i_full[active_atoms] = dw_dX_dalpha_i

                self.dw_dX_dalpha_list.append(dw_dX_dalpha_i_full)

        weight_gradient = np.zeros_like(reference_coordinates)   # (natms,3)
        weight_gradient[self.symmetry_information[3]] = weight_gradient_sub_imp_coord

        if self.pair_term_records is not None:
            self.pair_term_records.append((
                float(distance),
                distance_vector_sub.copy(),
                None if A_imp is None else float(A_imp),
                None if grad_A_imp is None else np.array(grad_A_imp, dtype=np.float64),
            ))

        return distance, dihedral_dist, denominator_imp_coord, weight_gradient, distance_vector_sub, grad_s, dw_dalpha_i, dw_dX_dalpha_i

    def _effective_dimension(self, scores, eps=1e-12):
        """
        Returns the participation number 1 / sum_i p_i^2 of the normalised
        non-negative scores (0 if all scores vanish).
        """
        scores = np.asarray(scores, dtype=float)
        scores = np.maximum(scores, 0.0)
        s = scores.sum()
        if s < eps:
            return 0.0
        p = scores / s
        return 1.0 / (np.sum(p * p) + eps)

    def determine_important_internal_coordinates(
        self,
        qm_gradient,
        molecule,
        z_matrix,
        datapoints,
    ):
        """
        Selects the internal coordinates that explain the error of the
        interpolation at the current geometry, for the constrained
        optimisation and the important coordinates of a new datapoint.

        1) Where does the local model fail? Response score |g_qm - g_model|
           per internal coordinate.
        2) Which displaced coordinates drive that failure? Source-blame score
           from a leave-one-source-out re-evaluation of each datapoint model.

        compute(molecule) must have been called before: the B-matrix and
        the internal coordinates of impes_coordinate are used. The QM
        gradient is projected with that B-matrix after conversion to its
        metric (multiplied by impes_coordinate.inv_sqrt_masses if set).

        :param qm_gradient:
            The Cartesian QM gradient at the current geometry in
            Hartree/bohr, (natoms, 3).
        :param molecule:
            The molecule at the current geometry.
        :param z_matrix:
            The z-matrix dictionary.
        :param datapoints:
            A list of (InterpolationDatapoint, weight) pairs, normally the
            datapoints that carry most of the interpolation weight.

        :return:
            A 4-tuple (primary_constraint, candidate_constraints,
            block_coords, fallback_constraints) of lists of 0-based z-matrix
            coordinate tuples: the dominant coordinate, the selected
            constraints, all coordinates of the best block, and further
            local coordinates. Empty lists if nothing is found.
        """

        w_source = 0.60
        w_response = 0.30
        w_support = 0.10
        helpful_penalty = 0.15
        support_top_frac = 0.20

        # block construction / ranking
        pair_weight = 0.35

        component_score_z = 2.0
        component_pair_z = 1.5
        component_coverage_mass = 0.60
        component_size_penalty = 0.01
        max_component_size = None

        # final selection from best block
        max_constraints_to_return = None
        singleton_dominance = 1.35
        singleton_source_frac = 0.85
        singleton_response_frac = 0.85
        secondary_score_frac = 0.60
        min_secondary_pair_frac = 0.35

        eps = 1e-12

        coords_flat, kinds_flat = self._flatten_internal_coordinates(z_matrix)
        N = len(coords_flat)
        if max_component_size is None:
            max_component_size = min(20, max(5, int(np.ceil(np.sqrt(N)))))

        if N == 0:
            return [], [], [], []

        constraints_to_exclude = []
        per_dp_results = []

        # Exclude near-linear angles from direct constraint selection

        for element in z_matrix.get("angles", []):
            current_angle = molecule.get_angle_in_degrees(
                (element[0] + 1, element[1] + 1, element[2] + 1)
            )
            if abs(current_angle) < 2.0 or abs(current_angle) > 178.0:
                constraints_to_exclude.append(tuple(int(x) for x in element))

        # The QM gradient is projected with impes_coordinate.b_matrix, which
        # is mass-weighted exactly when inv_sqrt_masses is set; transform the
        # gradient into the same metric (with the same factors).
        inv_sqrt_masses = self.impes_coordinate.inv_sqrt_masses
        if inv_sqrt_masses is not None:
            qm_gradient_mw = qm_gradient.reshape(-1) * inv_sqrt_masses
        else:
            qm_gradient_mw = qm_gradient.reshape(-1).copy()

        q_current = np.asarray(self.impes_coordinate.internal_coordinates_values, dtype=float)

        # The projected QM gradient does not depend on the datapoint.
        pred_qm_G_int = self.transform_gradient_to_internal_coordinates(
            molecule, qm_gradient_mw.reshape(qm_gradient.shape), self.impes_coordinate.b_matrix
        )

        # ------------------------------------------------------------------
        # Per-datapoint diagnostic pass
        # ------------------------------------------------------------------
        for datapoint, total_weight_contribution in datapoints:
            q_dp = np.asarray(datapoint.internal_coordinates_values, dtype=float)

            dq_raw, dq_eff, chain = self._compute_internal_differences(
                q_current=q_current,
                q_ref=q_dp,
                z_matrix=z_matrix,
            )

            diag_result = self.compute_internal_gradient_diagnostics_runtime_model(
                dq_raw=dq_raw,
                dq_eff=dq_eff,
                chain=chain,
                datapoint=datapoint,
                molecule=molecule,
                z_matrix=z_matrix,
                q_current=q_current,
                q_ref=q_dp,
                g_qm_int=np.asarray(pred_qm_G_int, dtype=float),
            )

            response_score = np.asarray(diag_result["response_score"], dtype=float)
            source_blame_score = np.asarray(diag_result["source_blame_score"], dtype=float)
            helpful_source_score = np.asarray(diag_result["helpful_source_score"], dtype=float)
            pair_score = np.asarray(diag_result["pair_score"], dtype=float)
            rho = np.asarray(diag_result["rho"], dtype=float)

            # local ranking used only for support counting
            local_source_for_support = self._normalize_nonnegative(source_blame_score)
            local_response_for_support = self._normalize_nonnegative(response_score)

            local_coord_score = self._normalize_nonnegative(
                0.65 * local_source_for_support
                + 0.35 * local_response_for_support
            )

            norm_displacement = self._normalized_internal_displacement(
                dq_eff=dq_eff,
                q_ref=q_dp,
                coords_flat=coords_flat,
                kinds_flat=kinds_flat,
            )

            support_mask = np.zeros(N, dtype=float)
            if local_coord_score.sum() > eps:
                order_dp = np.argsort(-local_coord_score)
                k_support = max(1, int(np.ceil(support_top_frac * N)))
                support_mask[order_dp[:k_support]] = 1.0

            per_dp_results.append({
                "datapoint": datapoint,
                "dp_weight": float(total_weight_contribution),
                "dq_raw": dq_raw.copy(),
                "dq_eff": dq_eff.copy(),
                "chain": chain.copy(),
                "rho": rho.copy(),
                "response_score": response_score.copy(),
                "source_blame_score": source_blame_score.copy(),
                "helpful_source_score": helpful_source_score.copy(),
                "pair_score": pair_score.copy(),
                "local_coord_score": local_coord_score.copy(),
                "support_mask": support_mask.copy(),
                "norm_displacement": norm_displacement.copy(),
            })

        if not per_dp_results:
            return [], [], [], []

        # ------------------------------------------------------------------
        # Normalize the datapoint weights (every datapoint passed in is used)
        # ------------------------------------------------------------------
        dp_weights = np.array(
            [r["dp_weight"] for r in per_dp_results],
            dtype=float,
        )

        dp_weights = np.maximum(dp_weights, 0.0)
        weight_sum = float(dp_weights.sum())

        if weight_sum < eps:
            dp_weights[:] = 1.0 / max(len(dp_weights), 1)
        else:
            dp_weights /= weight_sum

        # ------------------------------------------------------------------
        # Aggregate global response/source/pair information
        # ------------------------------------------------------------------
        global_response = np.zeros(N, dtype=float)
        global_source_blame = np.zeros(N, dtype=float)
        global_helpful = np.zeros(N, dtype=float)
        global_support = np.zeros(N, dtype=float)
        global_rho = np.zeros(N, dtype=float)
        global_pair = np.zeros((N, N), dtype=float)
        global_displacement = np.zeros(N, dtype=float)

        for idx in range(len(per_dp_results)):
            wk = dp_weights[idx]
            r = per_dp_results[idx]

            global_response += wk * r["response_score"]
            global_source_blame += wk * r["source_blame_score"]
            global_helpful += wk * r["helpful_source_score"]
            global_support += wk * r["support_mask"]
            global_rho += wk * r["rho"]
            global_pair += wk * r["pair_score"]
            global_displacement += wk * r["norm_displacement"]

        # Magnitude-preserving channels have already been aggregated.
        global_response = self._normalize_nonnegative(global_response)
        global_source_blame = self._normalize_nonnegative(global_source_blame)
        global_helpful = self._normalize_nonnegative(global_helpful)

        if global_support.max() > eps:
            global_support = global_support / (global_support.max() + eps)

        if np.max(global_pair) > eps:
            global_pair = global_pair / (np.max(global_pair) + eps)

        # Diagnostic/locality score only.
        regularized_coord_score = (
            w_source * global_source_blame
            + w_response * global_response
            + w_support * global_support
            - helpful_penalty * global_helpful
        )
        regularized_coord_score = self._normalize_nonnegative(
            np.maximum(regularized_coord_score, 0.0)
        )

        # Constraint-selection score.
        global_coord_score = self._normalize_nonnegative(
            0.50 * global_source_blame
            + 0.50 * global_response
        )
        N_coords = N

        d_eff = self._effective_dimension(global_coord_score)

        max_constraints_to_return = int(
            np.clip(
                np.ceil(1.25 * d_eff),
                1,
                min(12, max(3, int(np.ceil(np.sqrt(N_coords)))))
            )
        )

        ranked_blocks = []
        seen_blocks = set()

        components = self._build_faulty_components(
            global_coord_score=global_coord_score,
            global_pair=global_pair,
            coords_flat=coords_flat,
            kinds_flat=kinds_flat,
            constraints_to_exclude=constraints_to_exclude,
            score_z=component_score_z,
            pair_z=component_pair_z,
        )

        if not components:
            best_idx = int(np.argmax(global_coord_score))
            components = [[best_idx]]

        for component in components:

            component = [
                int(i) for i in component
                if tuple(int(x) for x in coords_flat[i]) not in constraints_to_exclude
            ]
            if not component:
                continue

            steering_score = {
                i: float(0.8 * global_source_blame[i] + 0.2 * global_response[i])
                for i in component
            }
            anchor_idx = max(component, key=lambda i: steering_score[i])
            anchor_coord = tuple(int(x) for x in coords_flat[anchor_idx])

            pair_to_anchor = np.asarray(global_pair[anchor_idx], dtype=float)

            component_sorted = sorted(
                component,
                key=lambda i: (
                    0.55 * float(global_coord_score[i])
                    + 0.30 * float(steering_score[i])
                    + 0.15 * float(pair_to_anchor[i])
                ),
                reverse=True,
            )
            component_sorted = [anchor_idx] + [
                int(i) for i in component_sorted if int(i) != int(anchor_idx)
            ]

            total_component_score = float(np.sum(global_coord_score[component_sorted]))
            block_indices = [int(anchor_idx)]
            cumulative_score = float(global_coord_score[anchor_idx])

            for i in component_sorted:
                i = int(i)
                if i == anchor_idx:
                    continue

                if len(block_indices) >= max_component_size:
                    break

                if (
                    total_component_score > eps
                    and cumulative_score >= component_coverage_mass * total_component_score
                    and len(block_indices) >= 1
                ):
                    break

                block_indices.append(i)
                cumulative_score += float(global_coord_score[i])

            key = tuple(sorted(block_indices))
            if key in seen_blocks:
                continue
            seen_blocks.add(key)

            coord_mass = float(np.sum(global_coord_score[block_indices]))
            support_mass = float(np.mean(global_support[block_indices])) if block_indices else 0.0

            pair_sum = 0.0
            n_pairs = 0

            for a in range(len(block_indices)):
                for b in range(a + 1, len(block_indices)):
                    i = block_indices[a]
                    j = block_indices[b]
                    pair_sum += float(global_pair[i, j])
                    n_pairs += 1

            pair_density = pair_sum / max(1, n_pairs)

            acq_score = (
                coord_mass
                + pair_weight * pair_density
                - component_size_penalty * len(block_indices)
            )

            block_coords = [tuple(int(x) for x in coords_flat[i]) for i in block_indices]
            block_types = [kinds_flat[i] for i in block_indices]

            ranked_blocks.append({
                "anchor_idx": int(anchor_idx),
                "anchor_coord": anchor_coord,
                "anchor_type": kinds_flat[anchor_idx],
                "block_indices": [int(i) for i in block_indices],
                "block_coords": block_coords,
                "block_types": block_types,
                "coord_mass": float(coord_mass),
                "pair_mass": float(pair_density),   # now normalized pair density
                "pair_sum": float(pair_sum),
                "support_mass": float(support_mass),
                "acq_score": float(acq_score),
            })

        ranked_blocks = sorted(ranked_blocks, key=lambda x: x["acq_score"], reverse=True)

        if not ranked_blocks:
            return [], [], [], []

        best_block = ranked_blocks[0]

        candidate_constraints, primary_constraint = self._select_constraints_from_block(
            block_indices=best_block["block_indices"],
            global_coord_score=global_coord_score,
            global_source_blame=global_source_blame,
            global_response_score=global_response,
            global_pair=global_pair,
            z_matrix=coords_flat,
            global_displacement=global_displacement,
            max_constraints=max_constraints_to_return,
            singleton_dominance=singleton_dominance,
            singleton_source_frac=singleton_source_frac,
            singleton_response_frac=singleton_response_frac,
            secondary_score_frac=secondary_score_frac,
            min_secondary_pair_frac=min_secondary_pair_frac,
        )

        fallback_constraints = self._select_fallback_locality_constraints(
            block_indices=best_block["block_indices"],
            candidate_constraints=candidate_constraints,
            primary_constraint=primary_constraint,
            coords_flat=coords_flat,
            kinds_flat=kinds_flat,
            constraints_to_exclude=constraints_to_exclude,
            global_coord_score=regularized_coord_score,
            global_pair=global_pair,
            global_displacement=global_displacement,
        )

        self.print_best_repair_block(
            best_block=best_block,
            coords_flat=coords_flat,
            kinds_flat=kinds_flat,
            global_coord_score=global_coord_score,
            global_source_blame=global_source_blame,
            global_response=global_response,
            global_support=global_support,
            global_rho=global_rho,
            global_helpful=global_helpful,
            global_displacement=global_displacement,
            primary_constraint=primary_constraint,
            candidate_constraints=candidate_constraints,
            fallback_constraints=fallback_constraints,
        )

        return primary_constraint, candidate_constraints, best_block["block_coords"], fallback_constraints

    # Helper methods
    def print_best_repair_block(
        self,
        best_block,
        coords_flat,
        kinds_flat,
        global_coord_score,
        global_source_blame,
        global_response,
        global_support,
        global_rho,
        global_helpful,
        global_displacement,
        primary_constraint,
        candidate_constraints,
        fallback_constraints,
    ):
        """
        Prints the block selected by determine_important_internal_coordinates
        with the scores of its coordinates.
        """

        title = 'Selected Internal-Coordinate Repair Block'

        self.ostream.print_blank()
        self.ostream.print_header(title)
        self.ostream.print_header('=' * len(title))
        self.ostream.print_blank()

        self.ostream.print_header(f"{'Anchor':<24}{self.format_coord(best_block['anchor_coord'])}")
        self.ostream.print_header(f"{'Anchor type':<24}{best_block['anchor_type']}")
        self.ostream.print_header(f"{'Block size':<24}{len(best_block['block_indices'])}")
        self.ostream.print_header(f"{'Coordinate mass':<24}{best_block['coord_mass']:.6e}")
        self.ostream.print_header(f"{'Pair mass':<24}{best_block['pair_mass']:.6e}")
        self.ostream.print_header(f"{'Support mass':<24}{best_block['support_mass']:.6e}")
        self.ostream.print_header(f"{'Acquisition score':<24}{best_block['acq_score']:.6e}")
        self.ostream.print_blank()

        header = (
            f"{'Rank':>4} "
            f"{'Idx':>5} "
            f"{'Type':>10} "
            f"{'Coord':>18} "
            f"{'Score':>12} "
            f"{'Source':>12} "
            f"{'Response':>12} "
            f"{'Support':>12} "
            f"{'Disp':>12} "
            f"{'Rho':>12} "
            f"{'Helpful':>12}"
        )

        self.ostream.print_header(header)
        self.ostream.print_header('-' * len(header))

        for rank, idx in enumerate(best_block['block_indices']):
            idx = int(idx)
            coord = tuple(int(x) for x in coords_flat[idx])
            coord_type = self.coord_type(coord, kinds_flat[idx])

            self.ostream.print_header(
                f"{rank:4d} "
                f"{idx:5d} "
                f"{coord_type:>10} "
                f"{self.format_coord(coord):>18} "
                f"{global_coord_score[idx]:12.6e} "
                f"{global_source_blame[idx]:12.6e} "
                f"{global_response[idx]:12.6e} "
                f"{global_support[idx]:12.6e} "
                f"{global_displacement[idx]:12.6e} "
                f"{global_rho[idx]:12.6e} "
                f"{global_helpful[idx]:12.6e}"
            )

        self.ostream.print_blank()
        self.ostream.print_header(f"{'Primary constraint':<24}{primary_constraint}")
        self.ostream.print_header(f"{'Candidate constraints':<24}{candidate_constraints}")

        if fallback_constraints:
            self.ostream.print_header('Fallback locality constraints')
            self.ostream.print_block(str(fallback_constraints))

        self.ostream.print_blank()
        self.ostream.flush()

    def _robust_positive_threshold(
        self,
        values: np.ndarray,
        z: float = 2.0,
        eps: float = 1e-12,
    ) -> float:
        """
        Robust threshold for positive values using median + z * MAD.

        This is more size-consistent than fixed fractions because the threshold
        adapts to the actual score distribution.
        """
        values = np.asarray(values, dtype=float)
        positive = values[values > eps]

        if positive.size == 0:
            return np.inf

        med = np.median(positive)
        mad = np.median(np.abs(positive - med)) + eps

        return float(med + z * 1.4826 * mad)

    def _is_coord_excluded_for_constraint(
        self,
        idx: int,
        coords_flat: List[Tuple[int, ...]],
        kinds_flat: List[str],
        constraints_to_exclude: List[Tuple[int, ...]],
    ) -> bool:
        """
        Centralized exclusion logic for direct constraint coordinates.
        """
        coord = tuple(int(x) for x in coords_flat[idx])

        if coord in constraints_to_exclude:
            return True

        return False

    def _build_faulty_components(
        self,
        *,
        global_coord_score: np.ndarray,
        global_pair: np.ndarray,
        coords_flat: List[Tuple[int, ...]],
        kinds_flat: List[str],
        constraints_to_exclude: List[Tuple[int, ...]],
        score_z: float = 2.0,
        pair_z: float = 1.5,
        min_atom_overlap: float = 0.25,
        bridge_score_frac: float = 0.20,
        bridge_score_z: float = 0.5,
        eps: float = 1e-12,
    ) -> List[List[int]]:
        """
        Build connected faulty coordinate components using core and bridge nodes.

        Core nodes:
            coordinates with clearly high error-source score.

        Bridge nodes:
            coordinates with moderate score that are chemically or blame-coupled
            to a core node.

        This avoids missing coordinates that are necessary to lock the failed
        geometry during constrained optimization but do not have top direct
        gradient-error scores.
        """
        global_coord_score = np.asarray(global_coord_score, dtype=float)
        global_pair = np.asarray(global_pair, dtype=float)

        N = len(global_coord_score)

        if N == 0:
            return []

        score_max = float(np.max(global_coord_score))

        if score_max <= eps:
            return []

        # ------------------------------------------------------------
        # Pair threshold
        # ------------------------------------------------------------
        upper = global_pair[np.triu_indices(N, k=1)]
        pair_threshold = self._robust_positive_threshold(
            upper,
            z=pair_z,
            eps=eps,
        )

        # If the pair matrix is almost empty, avoid pair_threshold = inf
        # blocking all chemically local bridge connections.
        if not np.isfinite(pair_threshold):
            pair_threshold = np.inf

        # ------------------------------------------------------------
        # Core and bridge score thresholds
        # ------------------------------------------------------------
        core_threshold = self._robust_positive_threshold(
            global_coord_score,
            z=score_z,
            eps=eps,
        )

        if not np.isfinite(core_threshold):
            core_threshold = score_max

        bridge_threshold = max(
            bridge_score_frac * score_max,
            self._robust_positive_threshold(
                global_coord_score,
                z=bridge_score_z,
                eps=eps,
            ),
        )

        if not np.isfinite(bridge_threshold):
            bridge_threshold = bridge_score_frac * score_max

        # Make sure bridge threshold is lower than core threshold.
        bridge_threshold = min(bridge_threshold, 0.75 * core_threshold)

        # ------------------------------------------------------------
        # Collect core and bridge-candidate nodes
        # ------------------------------------------------------------
        core_nodes = []
        bridge_candidates = []

        for i in range(N):
            if self._is_coord_excluded_for_constraint(
                i,
                coords_flat,
                kinds_flat,
                constraints_to_exclude,
            ):
                continue

            score_i = float(global_coord_score[i])

            if score_i >= core_threshold:
                core_nodes.append(int(i))

            elif score_i >= bridge_threshold:
                bridge_candidates.append(int(i))

        # Always retain at least the best non-excluded coordinate as a core.
        if not core_nodes:
            order = np.argsort(-global_coord_score)

            for i in order:
                i = int(i)

                if self._is_coord_excluded_for_constraint(
                    i,
                    coords_flat,
                    kinds_flat,
                    constraints_to_exclude,
                ):
                    continue

                core_nodes = [i]
                break

        if not core_nodes:
            return []

        # ------------------------------------------------------------
        # Accept only bridge nodes connected to a core node
        # ------------------------------------------------------------
        accepted_bridges = []

        for b in bridge_candidates:
            atoms_b = set(coords_flat[b])

            for c in core_nodes:
                atoms_c = set(coords_flat[c])

                atom_overlap = len(atoms_b.intersection(atoms_c)) / max(
                    1,
                    len(atoms_b.union(atoms_c)),
                )

                pair_ok = float(global_pair[b, c]) >= pair_threshold

                # If the pair graph is sparse, atom overlap is still allowed
                # to keep chemically local bridge nodes.
                overlap_ok = atom_overlap >= min_atom_overlap

                if pair_ok or overlap_ok:
                    accepted_bridges.append(int(b))
                    break

        nodes = set(core_nodes + accepted_bridges)

        if not nodes:
            return []

        # ------------------------------------------------------------
        # Build graph between all accepted nodes
        # ------------------------------------------------------------
        adjacency = {i: set() for i in nodes}
        node_list = sorted(nodes)

        for a, i in enumerate(node_list):
            atoms_i = set(coords_flat[i])

            for j in node_list[a + 1:]:
                atoms_j = set(coords_flat[j])

                pair_ij = float(global_pair[i, j])

                atom_overlap = len(atoms_i.intersection(atoms_j)) / max(
                    1,
                    len(atoms_i.union(atoms_j)),
                )

                pair_ok = pair_ij >= pair_threshold

                overlap_ok = atom_overlap >= min_atom_overlap

                # Additional score compatibility condition:
                # avoid connecting very weak bridge nodes to everything.
                score_compat_ok = (
                    min(global_coord_score[i], global_coord_score[j])
                    >= bridge_threshold
                )

                if pair_ok or (overlap_ok and score_compat_ok):
                    adjacency[i].add(j)
                    adjacency[j].add(i)

        # ------------------------------------------------------------
        # Connected components
        # ------------------------------------------------------------
        components = []
        visited = set()

        for start in node_list:
            if start in visited:
                continue

            stack = [start]
            visited.add(start)
            component = []

            while stack:
                u = stack.pop()
                component.append(int(u))

                for v in adjacency[u]:
                    if v not in visited:
                        visited.add(v)
                        stack.append(v)

            components.append(sorted(component))

        # ------------------------------------------------------------
        # Rank components by score mass
        # ------------------------------------------------------------
        components = sorted(
            components,
            key=lambda comp: float(np.sum(global_coord_score[comp])),
            reverse=True,
        )

        return components

    def _flatten_internal_coordinates(
        self,
        z_matrix: Dict[str, List[Tuple[int, ...]]]
    ) -> Tuple[List[Tuple[int, ...]], List[str]]:
        """
        Flatten z_matrix into a single coordinate list with explicit kinds.

        Order is fixed and must remain consistent everywhere.
        """
        coords = []
        kinds = []

        for coord in z_matrix.get("bonds", []):
            coords.append(tuple(int(x) for x in coord))
            kinds.append("bond")

        for coord in z_matrix.get("angles", []):
            coords.append(tuple(int(x) for x in coord))
            kinds.append("angle")

        for coord in z_matrix.get("dihedrals", []):
            coords.append(tuple(int(x) for x in coord))
            kinds.append("dihedral")

        for coord in z_matrix.get("impropers", []):
            coords.append(tuple(int(x) for x in coord))
            kinds.append("improper")

        return coords, kinds

    def _normalize_nonnegative(self, x: np.ndarray, eps: float = 1e-12) -> np.ndarray:
        """
        Clips x at zero and normalises it to unit sum (zeros if the sum is
        below eps).
        """
        x = np.asarray(x, dtype=float)
        x = np.maximum(x, 0.0)
        s = x.sum()
        if s < eps:
            return np.zeros_like(x)
        return x / s

    def _principal_torsion_delta(self, delta: float) -> float:
        """Wrap angular difference to [-pi, pi)."""
        return (delta + np.pi) % (2.0 * np.pi) - np.pi

    def _improper_dihedral_displacement(self, delta):
        """
        Improper dihedral displacement used in the local Taylor model:
            2 * tan(delta / 2)
        """
        d = self._principal_torsion_delta(delta)
        return 2.0 * np.tan(0.5 * d)

    def _improper_dihedral_chain(self, delta):
        """
        Chain factor for d[2*tan(delta/2)]/d(delta):
            sec^2(delta / 2)
        """
        d = self._principal_torsion_delta(delta)
        cos_half = np.cos(0.5 * d)
        return 1.0 / np.maximum(cos_half * cos_half, 1.0e-12)

    def coord_type(self, coord: Tuple[int, ...], kind: Optional[str] = None) -> str:
        """
        Prefer explicit kind if available. Fall back to tuple-length classification.
        """
        if kind is not None:
            return kind

        n = len(coord)
        if n == 2:
            return "bond"
        elif n == 3:
            return "angle"
        elif n == 4:
            return "dihedral"
        return "unknown"

    def format_coord(self, coord: Tuple[int, ...]) -> str:
        """
        Returns the coordinate as a string of a tuple of ints.
        """
        return str(tuple(int(x) for x in coord))

    def _compute_internal_differences(
        self,
        q_current: np.ndarray,
        q_ref: np.ndarray,
        z_matrix: Dict[str, List[Tuple[int, ...]]],
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Returns (dq_raw, dq_eff, chain) of q_current - q_ref: the raw
        differences (torsions wrapped to [-pi, pi)), the Taylor-model
        displacements (2 sin(Delta/2) for dihedrals, 2 tan(Delta/2) for
        impropers) and their chain factors.
        """

        q_current = np.asarray(q_current, dtype=float)
        q_ref = np.asarray(q_ref, dtype=float)

        coords_flat, kinds_flat = self._flatten_internal_coordinates(z_matrix)
        N = len(coords_flat)

        dq_raw = np.zeros(N, dtype=float)
        dq_eff = np.zeros(N, dtype=float)
        chain = np.ones(N, dtype=float)

        for idx, (_, kind) in enumerate(zip(coords_flat, kinds_flat)):
            if kind == "dihedral":
                d = self._principal_torsion_delta(q_current[idx] - q_ref[idx])
                dq_raw[idx] = d
                dq_eff[idx] = 2.0 * np.sin(0.5 * d)
                chain[idx] = np.cos(0.5 * d)

            elif kind == "improper":
                d = self._principal_torsion_delta(q_current[idx] - q_ref[idx])
                dq_raw[idx] = d
                dq_eff[idx] = self._improper_dihedral_displacement(d)
                chain[idx] = self._improper_dihedral_chain(d)

            else:
                d = q_current[idx] - q_ref[idx]
                dq_raw[idx] = d
                dq_eff[idx] = d
                chain[idx] = 1.0

        return dq_raw, dq_eff, chain

    def _predict_internal_gradient_runtime(self, datapoint, molecule, org_int_coords):
        """
        Returns the internal gradient predicted by the Taylor model of a
        datapoint: the Cartesian model gradient of compute_potential at
        org_int_coords, projected back with
        transform_gradient_to_internal_coordinates (the same B-matrix and
        metric as the projected QM gradient). bond_rmsd, angle_rmsd and
        dihedral_rmsd are left unchanged.

        :param datapoint:
            The InterpolationDatapoint.
        :param molecule:
            The molecule at the current geometry.
        :param org_int_coords:
            The internal coordinates at which the model is evaluated.
        """
        org_int_coords = np.asarray(org_int_coords, dtype=np.float64)

        # Avoid polluting outer diagnostics buffers.
        old_bond = self.bond_rmsd
        old_angle = self.angle_rmsd
        old_dihedral = self.dihedral_rmsd
        self.bond_rmsd = []
        self.angle_rmsd = []
        self.dihedral_rmsd = []

        try:
            _, pred_g_mw_cart, _ = self.compute_potential(datapoint, org_int_coords)
        finally:
            self.bond_rmsd = old_bond
            self.angle_rmsd = old_angle
            self.dihedral_rmsd = old_dihedral

        pred_g_int = self.transform_gradient_to_internal_coordinates(
            molecule,
            np.asarray(pred_g_mw_cart, dtype=np.float64),
            self.impes_coordinate.b_matrix,
        )
        return np.asarray(pred_g_int, dtype=np.float64).reshape(-1)

    def compute_internal_gradient_diagnostics_runtime_model(
        self,
        *,
        dq_raw,
        dq_eff,
        chain,
        datapoint,
        molecule,
        z_matrix,
        q_current,
        q_ref,
        g_qm_int,
        max_sources=None,
        eps=1.0e-12,
    ):
        """
        Source-blame diagnostics that mirror the actual interpolation model.
        For each source j, remove its displacement and re-evaluate the full model.

        :param dq_raw:
            The raw differences of _compute_internal_differences (unused).
        :param dq_eff:
            The Taylor-model displacements (used to rank the sources).
        :param chain:
            The chain factors (unused).
        :param datapoint:
            The InterpolationDatapoint.
        :param molecule:
            The molecule at the current geometry.
        :param z_matrix:
            The z-matrix dictionary.
        :param q_current:
            The internal coordinates of the current geometry.
        :param q_ref:
            The internal coordinates of the datapoint.
        :param g_qm_int:
            The projected QM gradient in internal coordinates.
        :param max_sources:
            Only the max_sources largest displacements are tested (None:
            all).

        :return:
            A dictionary with the per-coordinate scores response_score,
            source_blame_score, helpful_source_score, rho, the pair_score
            matrix, the blame matrices and the model gradient g_model.
        """
        coords_flat, kinds_flat = self._flatten_internal_coordinates(z_matrix)
        N = len(coords_flat)

        q_current = np.asarray(q_current, dtype=np.float64).reshape(-1)
        q_ref = np.asarray(q_ref, dtype=np.float64).reshape(-1)
        g_qm = np.asarray(g_qm_int, dtype=np.float64).reshape(-1)

        if q_current.size != N or q_ref.size != N or g_qm.size != N:
            raise ValueError("runtime diagnostics: size mismatch")

        g_model = self._predict_internal_gradient_runtime(
            datapoint,
            molecule,
            q_current,
        )

        delta_g_err = g_qm - g_model
        abs_err = np.abs(delta_g_err)

        source_indices = np.arange(N, dtype=np.int64)
        if isinstance(max_sources, int) and 0 < max_sources < N:
            source_indices = np.argsort(-np.abs(dq_eff))[:max_sources]

        blame_matrix = np.zeros((N, N), dtype=np.float64)

        for j in source_indices:
            q_wo_j = q_current.copy()

            if kinds_flat[j] in ("dihedral", "improper"):
                d = self._principal_torsion_delta(q_current[j] - q_ref[j])
                q_wo_j[j] = q_current[j] - d
            else:
                q_wo_j[j] = q_ref[j]

            g_model_wo_j = self._predict_internal_gradient_runtime(
                datapoint,
                molecule,
                q_wo_j,
            )

            err_wo_j = np.abs(g_qm - g_model_wo_j)
            scaled_err_wo_j = err_wo_j

            blame_matrix[:, j] = abs_err - scaled_err_wo_j

        harmful_matrix = np.maximum(blame_matrix, 0.0)
        helpful_matrix = np.maximum(-blame_matrix, 0.0)

        # Raw, magnitude-preserving channels.
        response_score = abs_err.copy()
        source_blame_score = np.sum(harmful_matrix, axis=0)
        helpful_source_score = np.sum(helpful_matrix, axis=0)

        pair_score = harmful_matrix + harmful_matrix.T
        np.fill_diagonal(pair_score, 0.0)

        rho = abs_err / (np.abs(g_model) + eps)

        return {
            "dq_eff": dq_eff,
            "chain": chain,
            "Hdq": np.zeros(N, dtype=np.float64),
            "diag_resp": np.zeros(N, dtype=np.float64),
            "offdiag_resp": np.zeros(N, dtype=np.float64),
            "g_model": g_model,
            "delta_g_err": delta_g_err,
            "rho": rho,
            "pE_total": np.zeros(N, dtype=np.float64),
            "source_contrib": blame_matrix,
            "blame_matrix": blame_matrix,
            "harmful_matrix": harmful_matrix,
            "helpful_matrix": helpful_matrix,
            "response_score": response_score,
            "source_blame_score": source_blame_score,
            "helpful_source_score": helpful_source_score,
            "pair_score": pair_score,
        }

    def _select_constraints_from_block(
        self,
        block_indices: List[int],
        global_coord_score: np.ndarray,
        global_source_blame: np.ndarray,
        global_response_score: np.ndarray,
        global_pair: np.ndarray,
        z_matrix: List[Tuple[int, ...]],
        global_displacement=None,
        max_constraints: int = 3,
        singleton_dominance: float = 1.35,
        singleton_source_frac: float = 0.85,
        singleton_response_frac: float = 0.85,
        secondary_score_frac: float = 0.60,
        min_secondary_pair_frac: float = 0.35,
        min_secondary_displacement: float = 0.5,
        min_secondary_relative_displacement: float = 0.10,
    ):
        """
        Choose final constraints from the winning block.

        Critical behavior:
        ------------------
        If one coordinate clearly dominates as both
        - a source of error and
        - a response channel,
        then return only that coordinate.

        This means:
        if a dihedral is truly the problematic piece, it is returned alone.
        No artificial non-dihedral is added in that case.
        """
        block_indices = list(block_indices)
        if not block_indices:
            return [], []

        eps = 1e-12

        block_scores = {i: float(global_coord_score[i]) for i in block_indices}
        steering_score = {
            i: float(0.5 * global_source_blame[i] + 0.5 * global_response_score[i])
            for i in block_indices
        }

        primary_idx = max(block_indices, key=lambda i: steering_score[i])
        primary_coord = tuple(int(x) for x in z_matrix[primary_idx])
        primary_constraint = [primary_coord]

        selected = [primary_coord]
        remaining = [i for i in block_indices if i != primary_idx]

        second_best = max([block_scores[i] for i in remaining], default=0.0)
        block_max_source = max([float(global_source_blame[i]) for i in block_indices], default=0.0)
        block_max_response = max([float(global_response_score[i]) for i in block_indices], default=0.0)

        # singleton rule:
        # if one coordinate clearly dominates, return only that coordinate
        if (
            block_scores[primary_idx] >= singleton_dominance * max(second_best, eps)
            and float(global_source_blame[primary_idx]) >= singleton_source_frac * max(block_max_source, eps)
            and float(global_response_score[primary_idx]) >= singleton_response_frac * max(block_max_response, eps)
        ):
            return selected, primary_constraint

        pair_row = np.asarray(global_pair[primary_idx], dtype=float)
        best_pair = max([float(pair_row[i]) for i in remaining], default=0.0)

        def secondary_rank_key(i):
            return 0.70 * float(pair_row[i]) + 0.30 * float(global_coord_score[i])

        remaining_sorted = sorted(remaining, key=secondary_rank_key, reverse=True)

        for i in remaining_sorted:
            if len(selected) >= max_constraints:
                break

            if block_scores[i] < secondary_score_frac * max(block_scores[primary_idx], eps):
                continue

            if best_pair > eps:
                pair_ok = float(pair_row[i]) >= min_secondary_pair_frac * best_pair
                if not pair_ok and block_scores[i] < 0.60 * max(block_scores[primary_idx], eps):
                    continue
            if global_displacement is not None:
                primary_disp = float(global_displacement[primary_idx])
                disp_i = float(global_displacement[i])

                required_disp = min_secondary_displacement

                if primary_disp > eps:
                    required_disp = max(
                        required_disp,
                        min_secondary_relative_displacement * primary_disp,
                    )

                if disp_i < required_disp:
                    continue

            coord_i = tuple(int(x) for x in z_matrix[i])
            if coord_i not in selected:
                selected.append(coord_i)

        return selected, primary_constraint

    def _section_from_kind(self, kind):
        """
        Maps 'bond', 'angle', 'dihedral', 'improper' to the z-matrix keys.
        """
        return {
            "bond": "bonds",
            "angle": "angles",
            "dihedral": "dihedrals",
            "improper": "impropers",
        }.get(kind)

    def _normalized_internal_displacement(self, dq_eff, q_ref, coords_flat, kinds_flat):
        """
        Returns |dq_eff| of every coordinate in units of its tc tolerance
        (_imp_coordinate_sigma).
        """
        disp = np.zeros(len(coords_flat), dtype=float)

        for idx, kind in enumerate(kinds_flat):
            section = self._section_from_kind(kind)
            if section is None:
                continue

            sigma = self._imp_coordinate_sigma(section, idx, q_ref)
            disp[idx] = abs(float(dq_eff[idx])) / max(float(sigma), 1.0e-12)

        return disp

    def _select_fallback_locality_constraints(
        self,
        *,
        block_indices,
        candidate_constraints,
        primary_constraint,
        coords_flat,
        kinds_flat,
        constraints_to_exclude,
        global_coord_score,
        global_pair,
        global_displacement,
        min_displacement_sigma=1.0,
        min_score_frac=0.10,
        min_pair_frac=0.10,
    ):
        """
        Returns further coordinates near the selected constraints (same
        block, shared atoms or coupled by the pair score) with a large
        displacement or score, ranked by displacement, score and coupling.
        """
        eps = 1.0e-12
        selected = {tuple(int(x) for x in c)
                    for c in list(candidate_constraints) + list(primary_constraint)}

        selected_atoms = {a for coord in selected for a in coord}
        block_set = set(int(i) for i in block_indices)

        primary_idx = None
        if primary_constraint:
            primary_coord = tuple(int(x) for x in primary_constraint[0])
            for idx, coord in enumerate(coords_flat):
                if tuple(int(x) for x in coord) == primary_coord:
                    primary_idx = idx
                    break

        pair_row = np.zeros(len(coords_flat), dtype=float)
        best_pair = 0.0
        if primary_idx is not None:
            pair_row = np.asarray(global_pair[primary_idx], dtype=float)
            best_pair = max([float(pair_row[i]) for i in block_set if i != primary_idx], default=0.0)

        disp_score = self._normalize_nonnegative(global_displacement)
        block_score_ref = max([float(global_coord_score[i]) for i in block_set], default=float(np.max(global_coord_score)))

        ranked = []
        for idx, coord_raw in enumerate(coords_flat):
            coord = tuple(int(x) for x in coord_raw)
            if coord in selected:
                continue
            if self._is_coord_excluded_for_constraint(idx, coords_flat, kinds_flat, constraints_to_exclude):
                continue

            atoms = set(coord)
            atom_overlap = bool(selected_atoms.intersection(atoms)) if selected_atoms else False
            pair_ok = best_pair > eps and float(pair_row[idx]) >= min_pair_frac * best_pair
            local_ok = idx in block_set or atom_overlap or pair_ok

            disp_ok = float(global_displacement[idx]) >= min_displacement_sigma
            score_ok = float(global_coord_score[idx]) >= min_score_frac * max(block_score_ref, eps)

            if not local_ok or not (disp_ok or score_ok):
                continue

            rank_value = (
                0.50 * float(disp_score[idx])
                + 0.30 * float(global_coord_score[idx])
                + 0.20 * float(pair_row[idx])
            )
            ranked.append((idx, rank_value))

        ranked = sorted(ranked, key=lambda x: x[1], reverse=True)
        return [tuple(int(x) for x in coords_flat[i]) for i, _ in ranked[:]]

    def transform_gradient_to_internal_coordinates(self, molecule, gradient, b_matrix, tol=1e-6):
        """
        Projects a gradient onto the internal coordinates,
        g_q = G^+ B g with G = B B^T. G^+ keeps the 3N-6 (3N-5 for two
        atoms) largest singular values among those above tol (the largest
        ones if fewer are above tol). The gradient must be in the metric of
        b_matrix (mass-weighted if b_matrix is).

        :param molecule:
            Unused; kept for the call signature.
        :param gradient:
            The gradient, (natoms, 3) or flat.
        :param b_matrix:
            The B-matrix.
        :param tol:
            The singular value threshold.
        """
        grad = np.asarray(gradient, dtype=np.float64).reshape(-1)
        B = np.asarray(b_matrix, dtype=np.float64)

        natm = grad.size // 3
        target_rank = 3 * natm - 6
        if natm == 2:
            target_rank += 1
        target_rank = max(1, int(target_rank))

        G = B @ B.T
        U, s, Vt = np.linalg.svd(G, full_matrices=False)

        s_inv = np.zeros_like(s)
        nz = np.where(s > tol)[0]

        if nz.size >= target_rank:
            keep = nz[np.argsort(s[nz])[-target_rank:]]
        else:
            keep = np.argsort(s)[-target_rank:]

        s_inv[keep] = 1.0 / s[keep]
        G_pinv = (U * s_inv) @ Vt
        return G_pinv @ (B @ grad)

    def get_energy(self):
        """ Returns the potential energy obtained by interpolation.
        """
        return self.impes_coordinate.energy

    def get_gradient(self):
        """ Returns the gradient obtained by interpolation.
        """
        return self.impes_coordinate.gradient
