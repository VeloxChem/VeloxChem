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
"""
The QM phase of the metal site force field builder: the three expensive
steps -- the constrained optimization, the Hessian and the RESP charges --
and the fit of the metal terms that reads them.
"""

from contextlib import contextmanager
from pathlib import Path
import numpy as np

from ..veloxchemlib import mpi_master
from ..molecularbasis import MolecularBasis
from ..respchargesdriver import RespChargesDriver
from ..scfrestdriver import ScfRestrictedDriver
from ..scfunrestdriver import ScfUnrestrictedDriver
from ..scfhessiandriver import ScfHessianDriver
from ..optimizationdriver import OptimizationDriver
from ..errorhandler import assert_msg_critical
from ..molecule import Molecule
from . import core
from .util import (Shell, on_master, collective, param, get_metal_keys,
                   extract_pairs, d4_charges, _folder_file, DONOR_ELEMENTS,
                   PARTIAL_HESSIAN_CUTOFF, WEAK_BRIDGE_TOLERANCE,
                   FITTED_COMMENT, GEOMETRY_FILE, HESSIAN_FILE, CHARGES_FILE)


@contextmanager
def _muted(driver, mute_scf=True):
    """
    Mutes the output stream of a driver for the duration of a run.

    OutputStream.mute is reference counted, so a mute that is not balanced
    by its unmute swallows every later warning of that stream outright.
    Doing it here means the balance holds even when the driver raises.

    :param driver:
        The driver whose stream to mute.
    :param mute_scf:
        Whether to mute at all.
    """

    if mute_scf:
        driver.ostream.mute()
    try:
        yield
    finally:
        if mute_scf:
            driver.ostream.unmute()


class QmParameterizer(Shell):
    """
    Pays for the QM behind a metal site force field and fits the metal
    terms from it.

    The three QM steps are collective: every rank enters them with the
    same molecule, and the drivers parallelize inside. Everything else --
    the pairs a partial Hessian is restricted to, the fit, the pruning,
    the resolution of what an earlier run left in a folder -- runs on the
    master rank and is broadcast.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    """

    # ------------------------------------------------------------------
    # the QM steps: collective, every rank enters them
    # ------------------------------------------------------------------

    @collective
    def _get_scf_driver(self,
                        molecule,
                        basis_set_label='def2-svp',
                        scf_drv=None,
                        xcfun='PBE0'):
        """
        Returns the SCF driver and basis set for the active site

        The driver assigned to scf_drv is used exactly as given, so it carries
        every QM setting beyond the functional and the basis set.

        :param molecule:
            The active site molecule.

        :return:
            The tuple of the SCF driver and the basis set.
        """

        if scf_drv is None:
            if molecule.get_multiplicity() != 1:
                scf_drv = ScfUnrestrictedDriver(self.comm, self.ostream)
            else:
                scf_drv = ScfRestrictedDriver(self.comm, self.ostream)
            if xcfun is not None:
                scf_drv.xcfun = xcfun

        basis = MolecularBasis.read(molecule, basis_set_label)

        return scf_drv, basis

    @collective
    def _run_scf(self, scf_drv, molecule, basis, mute_scf=True):
        """
        Runs the SCF for a geometry, whatever the driver already holds.

        Both the gradient and the Hessian driver reuse whatever is already in
        scf_drv.scf_results and only fall back to running an SCF when it is
        empty. They do not check that those results belong to the geometry they
        were handed, so the SCF is run explicitly whenever the geometry may have
        moved since the driver was last used.

        :param scf_drv:
            The SCF driver.
        :param molecule:
            The active site molecule.
        :param basis:
            The basis set.
        :param mute_scf:
            Whether to mute the driver while it runs.
        """

        with _muted(scf_drv, mute_scf):
            scf_drv.compute(molecule, basis)

    @collective
    def optimize_active_site(self,
                             active_site,
                             frozen_indices=None,
                             mute_scf=True,
                             basis_set_label='def2-svp',
                             scf_drv=None,
                             xcfun='PBE0',
                             constrain_capping_hydrogens=False):
        """
        Optimizes the active site with the beta carbons frozen.

        Freezing them keeps the spatial arrangement imposed by the protein
        backbone. Without it the site relaxes to a gas-phase geometry, and
        since the Seminario method takes the equilibrium values straight from
        the geometry, every fitted bond length and angle would then describe
        the wrong structure.

        :param frozen_indices:
            The zero-based active site indices to freeze. Defaults to
            constrained_indices().

        :return:
            The tuple of the optimized molecule and the results of the
            optimization driver. The caller decides whether to put the molecule
            back into the active site.
        """

        if frozen_indices is None:
            frozen_indices = core.constrained_indices(
                active_site,
                constrain_capping_hydrogens=constrain_capping_hydrogens)

        molecule = active_site['molecule']

        constraint = core.freeze_constraints(active_site,
                                             frozen_indices=frozen_indices)
        constraints = None if constraint is None else [constraint]

        self._print_muted_notice(
            f'the constrained optimization with {len(frozen_indices)} '
            'atom(s) frozen',
            mute_scf=mute_scf)

        scf_drv, basis = self._get_scf_driver(molecule,
                                              basis_set_label=basis_set_label,
                                              scf_drv=scf_drv,
                                              xcfun=xcfun)

        opt_drv = OptimizationDriver(scf_drv)
        opt_drv.constraints = constraints

        with _muted(opt_drv, mute_scf):
            opt_results = opt_drv.compute(molecule, basis)

        optimized = Molecule.read_xyz_string(opt_results['final_geometry'])
        optimized.set_charge(molecule.get_charge())
        optimized.set_multiplicity(molecule.get_multiplicity())

        return optimized, opt_results

    @collective
    def compute_hessian(self,
                        active_site,
                        atom_pairs=None,
                        mute_scf=True,
                        basis_set_label='def2-svp',
                        scf_drv=None,
                        xcfun='PBE0'):
        """
        Computes the nuclear Hessian of the active site.

        With atom_pairs the analytical Hessian is restricted to those blocks
        and everything else is left at zero, which is all the Seminario method
        needs when only the metal terms are being fitted. The diagonal blocks
        are added by the Hessian driver itself, so only the off-diagonal pairs
        need to be given.

        compute() passes the pairs of hessian_pairs unless
        calculate_partial_hessian is off, in which case it asks for the whole
        Hessian instead.

        :param atom_pairs:
            The list of zero-based (i, j) tuples, typically from
            hessian_pairs. None computes the full Hessian.

        :return:
            The Hessian as a (3N, 3N) numpy array in Hartree per Bohr squared.
        """

        molecule = active_site['molecule']

        scf_drv, basis = self._get_scf_driver(molecule,
                                              basis_set_label=basis_set_label,
                                              scf_drv=scf_drv,
                                              xcfun=xcfun)

        assert_msg_critical(
            scf_drv.solvation_model is None, 'compute_hessian: ScfHessianDriver '
            'does not support a solvation model')

        # The Hessian driver reuses scf_drv.scf_results without checking
        # which geometry they belong to, so the SCF is run here for the
        # current one. This costs nothing: the driver would otherwise run the
        # same SCF itself.
        self._print_muted_notice('the SCF for the Hessian', mute_scf=mute_scf)
        self._run_scf(scf_drv, molecule, basis, mute_scf=mute_scf)

        self._print_muted_notice('the Hessian', mute_scf=mute_scf)

        hessian_drv = ScfHessianDriver(scf_drv)
        # the numerical path ignores atom_pairs entirely and would silently
        # compute the full Hessian instead
        hessian_drv.numerical = False
        if atom_pairs is None:
            hessian_drv.atom_pairs = None
        else:
            hessian_drv.atom_pairs = [tuple(pair) for pair in atom_pairs]

        with _muted(hessian_drv, mute_scf):
            hessian_drv.compute(molecule, basis)

        hessian = np.copy(hessian_drv.hessian)

        return hessian

    @collective
    def compute_resp_charges(self, active_site, mute_scf=True):
        """
        Computes RESP charges for the active site.

        :return:
            The partial charges as an (N,) numpy array.
        """

        molecule = active_site['molecule']

        self._print_muted_notice(
            'the RESP charge fit at Hartree-Fock/6-31G*', mute_scf=mute_scf)

        resp_drv = RespChargesDriver(self.comm, self.ostream)

        # Neither a basis nor SCF results are passed: the driver then defaults
        # to Hartree-Fock with 6-31G*, which is what RESP charges are meant to
        # be fitted to, and runs its own SCF. Handing it the active site's own
        # functional and basis would silently fit the charges at a level the
        # RESP parameters were never derived for.
        with _muted(resp_drv, mute_scf):
            charges = resp_drv.compute(molecule)

        charges = self.comm.bcast(charges, root=mpi_master())
        charges = np.array(charges)

        return charges

    # ------------------------------------------------------------------
    # what is fitted from them, and what a run reads back
    # ------------------------------------------------------------------

    @on_master
    def hessian_pairs(self,
                      active_site,
                      bond_count=2,
                      partial_hessian_cutoff=PARTIAL_HESSIAN_CUTOFF):
        """
        Finds the atom pairs a partial Hessian has to hold blocks for.

        The pairs extract_pairs walks out of the connectivity cover exactly the
        metal terms the fit makes on the geometry as it stands. That is one
        geometry's worth of perception, and a coordination the QM optimization
        opened up past metal_bond_cutoff is precisely what add_metal_bond is
        reached for afterwards -- at which point the Hessian holds an all-zero
        block for the new bond and Seminario gives it no force constant, with
        the whole Hessian to pay for again to repair it. The walk therefore
        starts from a connectivity that additionally treats every donor atom
        within partial_hessian_cutoff of a metal as bonded to it, so a bond
        added later is already covered.

        Only what is computed is widened. Nothing reads this connectivity back:
        the active site keeps the bonding it was perceived with, and the fit
        keeps fitting that.

        :param active_site:
            The active site whose Hessian is being computed. Not modified.
        :param bond_count:
            The number of bonds to walk, as extract_pairs takes it.
        :param partial_hessian_cutoff:
            The distance in Angstrom within which an unbonded donor atom is
            covered along with the metal anyway. None walks the connectivity as
            it stands.

        :return:
            The tuple of the sorted pair list and the sorted atom list.
        """

        matrix = np.array(active_site['connectivity_matrix'], dtype=bool)
        metals = list(active_site['metal_indices'])

        if partial_hessian_cutoff is not None:
            molecule = active_site['molecule']
            coordinates = molecule.get_coordinates_in_angstrom()
            labels = molecule.get_labels()

            donors = [
                index for index, label in enumerate(labels)
                if label in DONOR_ELEMENTS
            ]

            for metal in metals:
                for donor in donors:
                    if matrix[metal, donor]:
                        continue
                    distance = np.linalg.norm(coordinates[donor] -
                                              coordinates[metal])
                    if distance <= partial_hessian_cutoff:
                        matrix[metal, donor] = True
                        matrix[donor, metal] = True

        return extract_pairs(matrix, metals, bond_count=bond_count)

    @on_master
    def d4_charges(self, active_site):
        """
        Returns D4 partial charges for the active site, and says so.

        The fallback of a run whose charges were neither fitted nor handed
        in; see util.d4_charges for what they are worth.

        :param active_site:
            The active site.

        :return:
            The charges as an (N,) numpy array, capping hydrogens included.
        """

        charges = d4_charges(active_site)

        self.ostream.print_info(
            f'Using D4 partial charges: {charges.size} atoms summing to '
            f'{charges.sum():+.3f} e. No RESP charges were supplied.')
        self.ostream.flush()

        return charges

    @on_master
    def fit_forcefield(self,
                       active_site,
                       forcefield,
                       hessian,
                       protected_bonds=None,
                       average_metal_terms=False,
                       metal_hessian_fitting_method='seminario',
                       prune_weak_bridge_bonds=True,
                       reparameterize_metal_angles=True,
                       weak_bridge_tolerance=WEAK_BRIDGE_TOLERANCE):
        """
        Fits the metal terms of a force field against a Hessian.

        Takes the force field core.build_forcefield seeded and replaces every
        metal bond and angle in it with the Seminario fit, then prunes the
        weak arm of any bridging residue and reports the force constants the
        fit could not determine. Everything else in the force field -- the
        typing, the charges, the annotations, the planarity impropers -- is
        the seeded build's and is left as it is; an improper across a pruned
        bond goes with the bond.

        :param active_site:
            The active site the force field was built for.
        :param forcefield:
            The seeded force field generator. Fitted in place and returned.
        :param hessian:
            The Hessian as a (3N, 3N) numpy array.
        :param protected_bonds:
            The metal bond keys the weak bridge pruning must leave alone,
            which is what a bond added by hand needs: manual_bond_keys turns
            the records of the binding modes into them.
        :param average_metal_terms:
            Whether equivalent metal terms are averaged, as
            MMForceFieldGenerator.reparameterize takes it.
        :param metal_hessian_fitting_method:
            The method MMForceFieldGenerator.reparameterize fits the metal
            bonds and angles with: 'seminario' (default), 'improved-seminario'
            or 'phf'/'phf(k)'.
        :param prune_weak_bridge_bonds:
            Whether the weak arm of a bridging residue is dropped after the
            fit; see _prune_weak_bridges.
        :param reparameterize_metal_angles:
            Whether the metal angles are fitted as well as the bonds. Off,
            they stay at whatever the generator guessed.
        :param weak_bridge_tolerance:
            The distance in Angstrom a bridge arm has to be longer than its
            residue's shortest metal bond by before it can be dropped.

        :return:
            The force field generator, fitted.
        """

        assert_msg_critical(
            not isinstance(hessian, str), 'fit_forcefield: the Hessian must be '
            'a numpy array.')

        n_atoms = active_site['molecule'].number_of_atoms()
        hessian = np.asarray(hessian)
        assert_msg_critical(
            hessian.shape == (3 * n_atoms, 3 * n_atoms),
            'fit_forcefield: Hessian shape '
            f'{hessian.shape} does not match {(3 * n_atoms, 3 * n_atoms)}')

        bonds, angles = get_metal_keys(forcefield, active_site)
        if not reparameterize_metal_angles:
            angles = []

        forcefield.reparameterize(hessian,
                                  reparameterize_keys=bonds + angles,
                                  average_metal_terms=average_metal_terms,
                                  method=metal_hessian_fitting_method)

        # the generator appends to the comment the seeding wrote
        for key in bonds:
            forcefield.bonds[key]['comment'] = FITTED_COMMENT
        for key in angles:
            forcefield.angles[key]['comment'] = FITTED_COMMENT

        if prune_weak_bridge_bonds:
            bonds, angles = self._prune_weak_bridges(
                forcefield,
                active_site,
                bonds,
                angles,
                weak_bridge_tolerance=weak_bridge_tolerance,
                protected=protected_bonds)

        self._check_force_constants(forcefield, active_site, bonds, angles,
                                    hessian)

        return forcefield

    @on_master
    def _prune_weak_bridges(self,
                            forcefield,
                            active_site,
                            bonds,
                            angles,
                            weak_bridge_tolerance=WEAK_BRIDGE_TOLERANCE,
                            protected=None):
        """
        Drops the long arm of a bridging residue that the fit gave no force
        constant.

        A residue that reaches two metals can hold one of them far more
        weakly than the other, whether it does the bridging through one atom
        that binds both (mu-1,1) or through two atoms of the same group that
        take one metal each (a mu-1,3 carboxylate). Both are the same
        situation seen from the residue, so the residue is what is looked at:
        the metal bonds of everything the metals leave connected together.

        Two independent things have to agree before an arm is dropped: the
        Hessian, by giving it no force constant at all, and the geometry, by
        holding it at least weak_bridge_tolerance further out than the
        shortest metal bond of that same residue. A zero on its own says
        nothing here - it can equally mean a geometry that is not stationary,
        or a Hessian that never covered the pair - which is why the distance
        has to agree.

        A bond with no stiffness is not the same thing as no bond: it leaves
        the pair at an equilibrium the dynamics never restores while its
        angles, torsions and exclusions all still act as though the two were
        bonded. So the bond goes, and with it every angle, torsion and
        improper that crosses it, the entry in the generator's own
        connectivity matrix and the 1-4 pairs derived from it. The generator
        is left describing one topology rather than two: a stale matrix would
        put the bond straight back the next time anything called
        create_topology on it, and a stale pair list keeps a scaled 1-4
        interaction for atoms that are no longer 1-4 at all.

        The active site's connectivity is not this function's to rewrite --
        it must not mutate its arguments -- so a caller that keeps a site
        beside the force field lifts it across with
        connectivity_from_forcefield.

        Only the longer arms are ever dropped, so a bridging residue can
        never lose the contact it is held by.

        :param forcefield:
            The force field being built, whose terms are removed in place.
        :param active_site:
            The active site, for the metal indices and the geometry the fit
            was made on.
        :param bonds:
            The metal bond keys.
        :param angles:
            The metal angle keys.

        :return:
            The metal bond and angle keys that are left.
        """

        # a bond asked for by hand is a decision, and the two things this
        # reads - a zero force constant and a long distance - are exactly what
        # a hand-added bond looks like when the Hessian does not cover it, so
        # leaving it in reach of the heuristic would take it straight back out
        protected = {tuple(key) for key in (protected or ())}

        metals = set(active_site['metal_indices'])
        coordinates = active_site['molecule'].get_coordinates_in_angstrom()
        labels = active_site['molecule'].get_labels()
        fragments = self._ligand_fragments(forcefield, metals)

        # the metal bonds of each residue, keyed by the fragment it is
        reached = {}
        for key in bonds:
            ligands = [index for index in key if index not in metals]
            if len(ligands) != 1:
                # a metal-metal bond belongs to no residue
                continue
            reached.setdefault(fragments[ligands[0]], []).append(key)

        removed = []
        for keys in reached.values():
            touched = {index for key in keys for index in key if index in metals}
            if len(touched) < 2:
                # one metal is a grip, however many atoms it is made with
                continue

            def length(key):
                return float(
                    np.linalg.norm(coordinates[key[0]] - coordinates[key[1]]))

            lengths = {key: length(key) for key in keys}
            shortest = min(lengths.values())

            for key in keys:
                if key in protected:
                    continue
                if forcefield.bonds[key]['force_constant'] != 0.0:
                    continue
                if lengths[key] - shortest < weak_bridge_tolerance:
                    continue
                removed.append((key, lengths[key], shortest))

        if not removed:
            return bonds, angles

        for key, length, shortest in removed:
            pair = set(key)
            names = '-'.join(labels[index] for index in key)

            crossing = {
                'angle': self._terms_crossing(forcefield.angles, pair),
                'torsion': self._terms_crossing(forcefield.dihedrals, pair),
                'improper': self._terms_crossing(forcefield.impropers, pair, path=False),
            }

            del forcefield.bonds[key]
            for table, found in ((forcefield.angles, crossing['angle']),
                                 (forcefield.dihedrals, crossing['torsion']),
                                 (forcefield.impropers, crossing['improper'])):
                for crossed in found:
                    del table[crossed]

            forcefield.connectivity_matrix[key[0], key[1]] = 0
            forcefield.connectivity_matrix[key[1], key[0]] = 0

            counts = ', '.join(f'{len(found)} {name}(s)'
                               for name, found in crossing.items())
            self.ostream.print_info(
                f'Removed the bridging bond {key} {names}: it was fitted to '
                f'no force constant and is {length:.2f} A long against '
                f'{shortest:.2f} A for the shortest metal bond of the same '
                f'residue. {counts} crossing it were removed with it.')

        # The 1-4 pairs are a function of the connectivity and are not a table
        # a term can be deleted out of: a pair reached by a second dihedral path
        # is still 1-4 once this bond is gone, so filtering by key would drop
        # pairs that must stay. Re-deriving them off the pruned matrix gives
        # exactly what create_topology would have produced had the bond never
        # been there, which is the point.
        *_, forcefield.pairs = forcefield.generate_topology_indices(len(labels))

        self.ostream.flush()

        # filtered rather than re-derived, so that an angle list the caller
        # emptied on purpose stays empty
        return ([key for key in bonds if key in forcefield.bonds],
                [key for key in angles if key in forcefield.angles])

    @on_master
    def _ligand_fragments(self, forcefield, metals):
        """
        Numbers the residues a force field holds, as what the metals leave
        connected together.

        The bonds of the force field are walked rather than the connectivity
        matrix, so that what counts as one residue is what this force field
        is actually wired as. A metal is not part of any of them, which is
        the whole point: it is what would otherwise join two residues into
        one.

        :param forcefield:
            The force field.
        :param metals:
            The indices of the metal centers.

        :return:
            A dictionary from atom index to fragment number, holding every
            atom that is not a metal.
        """

        neighbors = {}
        for first, second in forcefield.bonds:
            if first in metals or second in metals:
                continue
            neighbors.setdefault(first, set()).add(second)
            neighbors.setdefault(second, set()).add(first)

        fragments = {}
        count = 0

        for atom in forcefield.atoms:
            if atom in metals or atom in fragments:
                continue

            stack = [atom]
            while stack:
                current = stack.pop()
                if current in fragments:
                    continue
                fragments[current] = count
                stack.extend(neighbors.get(current, set()) - fragments.keys())

            count += 1

        return fragments

    @on_master
    def _terms_crossing(self, table, pair, path=True):
        """
        Finds the terms of one table that act across a bond.

        A bond, an angle and a torsion are paths, so a term crosses the bond
        when the two atoms are neighbours in its key. An improper is not a
        path - its central atom is written first, with the substituents after
        it - so there it is enough that both atoms are in the key at all.

        :param table:
            The bonds, angles, dihedrals or impropers of a force field.
        :param pair:
            The two atoms of the bond, as a set.
        :param path:
            Whether the keys of the table are paths.

        :return:
            The keys that cross the bond.
        """

        if not path:
            return [key for key in table if pair <= set(key)]

        return [
            key for key in table if any({key[index], key[index + 1]} == pair
                                        for index in range(len(key) - 1))
        ]

    @on_master
    def _check_force_constants(self,
                               forcefield,
                               active_site,
                               bonds,
                               angles,
                               hessian=None):
        """
        Warns about metal terms whose force constant was fitted to zero, and
        says which of the two reasons put it there.

        A term reads its own blocks of the Hessian and nothing else: a bond
        (i, j) reads the (i, j) block and an angle (i, j, k) reads the (i, j)
        and (j, k) blocks. So a zero comes from one of two places, and the fix
        is not the same:

        - **the Hessian does not cover the term.** compute_hessian is
          restricted to the pairs hessian_pairs walks out of the connectivity,
          and everything else is left at zero. A bond that was not there when
          those pairs were taken, and further out than
          partial_hessian_cutoff was forgiving of, has nothing to project,
          which is what a Hessian reused from a folder or supplied by hand
          runs into when the coordination has moved on since. Only recomputing
          it on this coordination fixes that.
        - **the projection was negative and Seminario clamped it**, which
          means the geometry is not stationary along that coordinate. On an
          unrelaxed structure this typically wipes out the long, strained
          metal-ligand bonds, which are exactly the ones that matter for a
          bridged binuclear site.

        Without a Hessian the two cannot be told apart and everything is
        reported as the second.

        :param bonds:
            The metal bond keys.
        :param angles:
            The metal angle keys.
        :param hessian:
            The Hessian the terms were fitted to, for telling an uncovered
            term from a clamped one.
        """

        labels = active_site['molecule'].get_labels()

        zero = [(key, 'bond') for key in bonds
                if forcefield.bonds[key]['force_constant'] == 0.0]
        zero += [(key, 'angle') for key in angles
                 if forcefield.angles[key]['force_constant'] == 0.0]

        if not zero:
            return

        n_bonds = sum(1 for _, kind in zero if kind == 'bond')
        n_angles = len(zero) - n_bonds

        uncovered = [(key, kind) for key, kind in zero
                     if not self._hessian_covers(hessian, key)]
        clamped = [pair for pair in zero if pair not in uncovered]

        def report(terms):
            for key, kind in terms:
                names = '-'.join(labels[index] for index in key)
                self.ostream.print_warning(
                    f'  zero force constant: {kind} {key} {names}')

        self.ostream.print_warning(
            f'{n_bonds} of {len(bonds)} metal bond(s) and {n_angles} of '
            f'{len(angles)} metal angle(s) got a zero force constant.')

        if uncovered:
            self.ostream.print_warning(
                f'{len(uncovered)} of them are not covered by the Hessian at '
                'all: it was restricted to the atom pairs of a different '
                'coordination, so it holds no data for these terms. Recompute '
                'it on this one rather than reusing it.')
            report(uncovered)

        if clamped:
            self.ostream.print_warning(
                f'{len(clamped)} of them were clamped from a negative '
                'projection. The Hessian is most likely not evaluated at a '
                'stationary point; run optimize_active_site first.')
            report(clamped)

        self.ostream.flush()

    @on_master
    def _hessian_covers(self, hessian, key):
        """
        Says whether a Hessian holds anything for one term.

        The blocks a term reads are the bonded pairs of its key, which is what
        extract_pairs walks out of the connectivity: (i, j) for a bond and
        (i, j) together with (j, k) for an angle. A block left at exactly zero
        was never computed, since compute_hessian fills only the pairs it is
        given.

        :param hessian:
            The Hessian, or None when there is none to look at.
        :param key:
            The bond or angle key.

        :return:
            True when every block the term reads holds something, and when
            there is no Hessian to say otherwise.
        """

        if hessian is None:
            return True

        hessian = np.asarray(hessian)
        pairs = [(key[index], key[index + 1]) for index in range(len(key) - 1)]

        for first, second in pairs:
            # the driver fills the pairs it is given, and a key is not stored
            # in any particular order, so both orientations are looked at
            block = hessian[3 * first:3 * first + 3, 3 * second:3 * second + 3]
            mirror = hessian[3 * second:3 * second + 3, 3 * first:3 * first + 3]
            if not np.any(block) and not np.any(mirror):
                return False

        return True

    @on_master
    def _resolve_source(self, supplied, filename, label, folder=None):
        """
        Applies the precedence the three resolvers share: what the caller
        handed in beats what an earlier run left in the working folder, and
        neither means there is nothing to use.

        Stated once here rather than three times, so the rule cannot come to
        mean different things for a geometry, a Hessian and a set of charges.

        Where a value came from is announced here for the same reason. This is
        the only place that knows the answer, and a caller that guessed said
        'given to build_forcefield' about a file a resumed run had left in the
        folder -- misattributing the single most useful fact when a reused
        folder produces a wrong fit.

        :param supplied:
            What the caller passed, or None.
        :param filename:
            The name the step writes its result under, one of the file name
            constants.
        :param label:
            What the value is, for the announcement.
        :param folder:
            The working folder to fall back on.

        :return:
            The tuple of the source and the flag that is True when it came
            from the folder rather than from the caller. The source is None
            when there is nothing to use.
        """

        if supplied is not None:
            self.ostream.print_info(f'Using the {label} supplied by the caller.')
            self.ostream.flush()
            return supplied, False

        source = _folder_file(filename, folder=folder)

        if source is None:
            return None, False

        self.ostream.print_info(f'Reusing {filename} from {folder}.')
        self.ostream.flush()

        return source, True

    @on_master
    def _resolve_optimized_geometry(self,
                                    active_site,
                                    folder=None,
                                    optimized_geometry=None):
        """
        Validates a geometry supplied through optimized_geometry, or left
        behind in the working folder by an earlier run. The element sequence is checked against the extracted active site
        :param active_site:
            The active site, to validate against.

        :return:
            The molecule, or None if there is nothing to use.
        """

        source, _ = self._resolve_source(optimized_geometry,
                                         GEOMETRY_FILE,
                                         'geometry',
                                         folder=folder)

        if source is None:
            return None

        if isinstance(source, Molecule):
            molecule = Molecule(source)
        else:
            path = Path(source)
            assert_msg_critical(
                path.is_file(), '_resolve_optimized_geometry: the geometry file '
                f'{path} not found')
            molecule = Molecule.read_xyz_file(str(path))

        active_site = active_site['molecule']

        assert_msg_critical(
            molecule.number_of_atoms() == active_site.number_of_atoms(),
            '_resolve_optimized_geometry: the geometry has '
            f'{molecule.number_of_atoms()} atoms but the extracted active site '
            f'has {active_site.number_of_atoms()}')

        assert_msg_critical(
            list(molecule.get_labels()) == list(active_site.get_labels()),
            '_resolve_optimized_geometry: the elements of '
            'optimized_geometry do not match the extracted active site, so it '
            'describes a different structure')

        molecule.set_charge(active_site.get_charge())
        molecule.set_multiplicity(active_site.get_multiplicity())

        return molecule

    @on_master
    def _resolve_hessian(self, active_site, folder=None, hessian=None):
        """
        Validates a Hessian supplied through the hessian setting, or left
        behind in the working folder by an earlier run.

        :param active_site:
            The active site, to validate the shape against.

        :return:
            The Hessian, or None if there is nothing to use.
        """

        source, reused = self._resolve_source(hessian,
                                              HESSIAN_FILE,
                                              'Hessian',
                                              folder=folder)

        if source is None:
            return None

        assert_msg_critical(
            not isinstance(source, (str, Path)) or Path(source).is_file(),
            f'_resolve_hessian: hessian file {source} not found')

        if isinstance(source, (str, Path)):
            hessian = np.loadtxt(source)
        else:
            hessian = np.asarray(source)

        n_atoms = active_site['molecule'].number_of_atoms()
        expected = (3 * n_atoms, 3 * n_atoms)

        assert_msg_critical(
            hessian.shape == expected,
            f'_resolve_hessian: hessian has shape {hessian.shape} '
            f'but the extracted active site needs {expected}')

        # The shape and the elements match any site of the same composition, so
        # a file left behind by a run whose coordination differed passes both
        # and then fits zeros. A block the metal terms read that was never
        # filled is what says so, and is worth recomputing rather than warning
        # about: an explicitly supplied Hessian is an instruction, a file lying
        # in the folder is a guess.
        if reused and not self._hessian_covers_site(hessian, active_site):
            self.ostream.print_warning(
                f'{HESSIAN_FILE} in {folder} holds nothing for some of the '
                'metal terms of this active site, so it was computed for a '
                'different coordination. Ignoring it.')
            self.ostream.flush()
            return None

        return hessian

    @on_master
    def _hessian_covers_site(self, hessian, active_site):
        """
        Says whether a Hessian holds data for every metal term of a site.

        The blocks the metal terms read are the pairs extract_pairs walks out
        of the connectivity as it stands. Deliberately not hessian_pairs: what
        the terms need is what has to be there, while the donors
        partial_hessian_cutoff was forgiving of are a surplus a file computed
        under another setting is not worth rejecting over.

        :param hessian:
            The Hessian.
        :param active_site:
            The active site whose metal terms have to be covered.

        :return:
            True when every pair the terms read holds something.
        """

        pairs, _ = extract_pairs(active_site['connectivity_matrix'],
                                 active_site['metal_indices'],
                                 bond_count=2)

        return all(self._hessian_covers(hessian, pair) for pair in pairs)

    @on_master
    def _resolve_partial_charges(self,
                                 active_site,
                                 folder=None,
                                 partial_charges=None):
        """
        Validates charges supplied through the partial_charges setting, or
        left behind in the working folder by an earlier run.

        :param active_site:
            The active site, to validate the count against.

        :return:
            The partial charges, or None if there is nothing to use.
        """

        source, _ = self._resolve_source(partial_charges,
                                         CHARGES_FILE,
                                         'partial charges',
                                         folder=folder)

        if source is None:
            return None

        if isinstance(source, (str, Path)):
            assert_msg_critical(
                Path(source).is_file(),
                '_resolve_partial_charges: the charges file '
                f'{source} not found')
            charges = np.loadtxt(source)
        else:
            charges = np.asarray(source)

        n_atoms = active_site['molecule'].number_of_atoms()

        assert_msg_critical(
            charges.shape == (n_atoms, ),
            f'_resolve_partial_charges: the charges have shape '
            f'{charges.shape} but the extracted active site has {n_atoms} atoms')

        total = float(np.sum(charges))
        expected = int(active_site['molecule'].get_charge())

        if abs(total - expected) > 1.0e-3:
            self.ostream.print_warning(
                f'The supplied partial charges sum to {total:+.3f}, but the '
                f'active site charge is {expected:+d}')
            self.ostream.flush()

        return charges

    @on_master
    def print_partial_charges(self,
                              topology,
                              active_site,
                              partial_charges,
                              corrected_charges,
                              residue_labels):
        """
        Prints the fitted charges and what the capping correction did to them.

        :param topology:
            The protonated topology, for the residue each active site atom belongs
            to.
        :param active_site:
            The active site.
        :param partial_charges:
            The charges as fitted, capping hydrogens included.
        :param corrected_charges:
            The same charges after redistribute_cap_charges.
        :param residue_labels:
            The ASP130-style label of every topology residue, by residue index.
        """

        charges = np.asarray(partial_charges)
        caps = sorted(active_site['cap_indices'])
        metals = active_site['metal_indices']
        labels = active_site['molecule'].get_labels()
        n_atoms = len(charges)
        rest = [index for index in range(n_atoms) if index not in caps]
        cap_charge = float(sum(charges[index] for index in caps))

        self.ostream.print_blank()
        self.ostream.print_header('Partial charges')
        self.ostream.print_header(15 * '-')
        self.ostream.print_header(
            param('active site charge',
                  f'{int(active_site["molecule"].get_charge()):+d}'))
        self.ostream.print_header(param('fitted total', f'{charges.sum():+.4f} e'))
        self.ostream.print_header(param('on capping hydrogens', f'{cap_charge:+.4f} e'))
        self.ostream.print_header(
            param(f'spread over {len(rest)} atoms',
                  f'{cap_charge / len(rest):+.4f} e each'))
        self.ostream.print_blank()

        # group what is left by the residue each atom came from
        atoms = list(topology.atoms())
        by_residue = {}
        for index in rest:
            residue = atoms[active_site['atom_map'][index]].residue
            by_residue.setdefault(residue, []).append(index)

        corrected = np.asarray(corrected_charges)

        valstr = '{:>16} {:>7} | {:>12}'.format('fragment', 'atoms', 'charge')
        self.ostream.print_header(valstr)
        self.ostream.print_header(45 * '-')

        for residue, indices in by_residue.items():
            total = sum(corrected[index] for index in indices)
            if len(indices) == 1 and indices[0] in metals:
                name = f'{labels[indices[0]]} (metal)'
            else:
                name = residue_labels[residue.index]
            valstr = '{:>16} {:>7} | {:>12.4f}'.format(name, len(indices), total)
            self.ostream.print_header(valstr)

        self.ostream.print_blank()
        self.ostream.flush()

    @on_master
    def print_metal_parameters(self, active_site, forcefield, metal_keys):
        """
        Prints the fitted metal bonds and angles.

        :param active_site:
            The active site.
        :param forcefield:
            The fitted force field.
        :param metal_keys:
            The (bonds, angles) keys of the metal terms, from get_metal_keys.
        """

        labels = active_site['molecule'].get_labels()
        metals = set(active_site['metal_indices'])
        coords = active_site['molecule'].get_coordinates_in_angstrom()
        bonds, angles = metal_keys

        self.ostream.print_blank()
        self.ostream.print_header('Metal bonds')
        self.ostream.print_header(11 * '-')
        valstr = '{:>12} {:>7} | {:>9} | {:>21}'.format('atoms', 'elements',
                                                        'r0 (A)',
                                                        'k (kcal/mol/A^2)')
        self.ostream.print_header(valstr)
        self.ostream.print_header(60 * '-')

        for key in bonds:
            params = forcefield.bonds[key]
            names = '-'.join(labels[index] for index in key)
            # kJ/mol/nm^2 to kcal/mol/A^2
            force_constant = params['force_constant'] / 100.0 / 4.184
            valstr = '{:>12} {:>7} | {:>9.3f} | {:>21.1f}'.format(
                str(key), names, params['equilibrium'] * 10.0, force_constant)
            self.ostream.print_header(valstr)

        self.ostream.print_blank()
        self.ostream.print_header('Metal angles')
        self.ostream.print_header(12 * '-')
        valstr = '{:>14} {:>9} | {:>12} | {:>19}'.format('atoms', 'elements',
                                                         'theta0 (deg)',
                                                         'k (kJ/mol/rad^2)')
        self.ostream.print_header(valstr)
        self.ostream.print_header(60 * '-')

        for key in angles:
            params = forcefield.angles[key]
            names = '-'.join(labels[index] for index in key)
            bridging = key[0] in metals and key[2] in metals
            # the marker gets a fixed-width field of its own, otherwise the
            # centering of print_header would shift the marked line
            valstr = '{:>14} {:>9} | {:>12.1f} | {:>19.1f} {:<9}'.format(
                str(key), names, params['equilibrium'], params['force_constant'],
                'bridging' if bridging else '')
            self.ostream.print_header(valstr)

        self.ostream.print_blank()

        for metal_a in sorted(metals):
            for metal_b in sorted(metals):
                if metal_a >= metal_b:
                    continue
                distance = np.linalg.norm(coords[metal_a] - coords[metal_b])
                self.ostream.print_header(
                    param(f'{labels[metal_a]}-{labels[metal_b]} distance',
                          f'{distance:.3f} A'))

        self.ostream.print_blank()
        self.ostream.flush()

    @on_master
    def _print_muted_notice(self, step, mute_scf=True):
        """
        Announces a long calculation whose output is being suppressed.

        :param step:
            A description of the step about to run.
        """

        if mute_scf:
            self.ostream.print_info(
                f'Running {step} with muted QM output. Set mute_scf to False '
                'to follow it.')
        else:
            self.ostream.print_info(f'Running {step}.')
        self.ostream.flush()
