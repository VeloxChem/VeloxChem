# This class is based and originated by the turtlemap programm by Lukas Lampe and 
# he is being the main contributer to that class which have been implemented based on
# his work!!

import numpy as np
import math
from copy import copy
from scipy.optimize import linear_sum_assignment

# from .utils import fix_groups
# from .utils import linear_assignment_solver
# from .utils import permute
# from .utils import output


class AtomMapper:

    """ An object that maps atoms between two molecule objects. """
    def __init__(self, reactant, product, tau=0, max_maps=1000, oop_thresh=100.0, sphere_guess=False, n_spheres=4, max_iter=100, print_output=False):
        """ Initializes the optimizer.

        Args:
            reactant:   The molecule object of the reactant.
            product:    The molecule object of the product.
            tau:        The upper limit for the difference between chemical distances and minimal chemical distance to be
                        tolerated for the reaction center detection, prescreening of atom maps and final evaluation of atom
                        maps.
            max_maps:   The maximum number of atom maps to be processed.
            oop_threshold:
                        Threshold in square degrees for an element in the dihedral cost matrix that has to be exceeded for
                        using out-of-plane distortions instead of dihedral angles. This avoids using out-of-plane distortions
                        for planar functional groups.
            sphere_guess:
                        The elements in coordination spheres are counted to generate an initial guess if True.
            n_spheres:  The number of coordination spheres used for the sphere guess.
            max_iter:   The maximum number of iterations for the steepest descent procedure.
            print_output:
                        Output is shown if True.
        """
        if len(reactant.get_labels()) != len(product.get_labels()):
            exit("Molecules need to contain the same number atoms.")
        if (np.sort(reactant.get_labels()) != np.sort(product.get_labels())).any():
            exit("Molecules need to contain the same number of elements.")
        self.n_atoms = len(reactant.get_labels())
        self.reactant = reactant; self.product = product
        self.labels_r = np.array(reactant.get_labels(), dtype='<U1'); self.labels_p = np.array(product.get_labels(), dtype='<U1')
        self.con_r = reactant.get_connectivity_matrix(); self.con_p = product.get_connectivity_matrix()
        self.tau = tau
        self.max_maps = max_maps
        self.oop_thresh = oop_thresh
        self.sphere_guess = sphere_guess
        self.n_spheres = n_spheres
        self.max_iter = max_iter
        self.print_output = print_output
        self.mono_el = np.intersect1d(reactant.identify_monovalent_elements(), product.identify_monovalent_elements())
    
    def determine_symmetry_group(self):

        atom_map, symmetry_groups, rc_groups = self.detect_reaction_center()
        _, refined_symmetry_groups = self.refine_environment(atom_map, symmetry_groups, rc_groups)
    
        return refined_symmetry_groups
    
    def perform(self, write_xyz=False):
        """ Performs the atom mapping. """
        # Construct initial atom map and detects reaction center.
        atom_map, env_groups, rc_groups = self.detect_reaction_center()
        # Refine environment groups.
        atom_map, env_groups = self.refine_environment(atom_map, env_groups, rc_groups)
        # Permute reaction center groups and optimize mappings for all atom maps.
        atom_maps = self.run_optimization(atom_map, env_groups, rc_groups)
        # Evaluate and print final atom maps.
        atom_maps = self.evaluate_results(atom_maps, write_xyz)
        return atom_maps

    def detect_reaction_center(self):
        """ Constructs initial atom map and detects reaction center.

        Returns:
            atom_map:   The atom map.
            env_groups: The environment groups.
            rc_groups:  The reaction center groups.
        """
        atom_map = np.arange(self.n_atoms)
        env_groups = []
        rc_groups = []
        rc_atoms = []
        for el, n_el in zip(*np.unique(self.labels_r, return_counts=True)):
            if el in self.mono_el:
                for i, j in zip(np.where(self.labels_r == el)[0], np.where(self.labels_p == el)[0]):
                    atom_map[i] = j
            else:
                list_r = np.where(self.labels_r == el)[0]
                list_p = np.where(self.labels_p == el)[0]
                cost = np.zeros((n_el, n_el), dtype=int)
                for i, a in enumerate(list_r):
                    for j, b in enumerate(list_p):
                        cost[i, j] += self.initial_guess_cost(a, b)
                el_map = self.linear_assignment_solver(cost)
                cost = cost[:,el_map]
                for i, j in enumerate(list_r):
                    atom_map[j] = list_p[el_map[i]]
                groups = []
                for i in range(n_el):
                    for j in range(i + 1, n_el):
                        add_to_groups = False
                        if cost[i, i] == 0 and cost[j, j] == 0:
                            if cost[i, j] == 0 and cost[j, i] == 0:
                                add_to_groups = True
                        elif (cost[i, j] + cost[j, i] - cost[i, i] - cost[j, j]) <= self.tau:
                            add_to_groups = True
                        if add_to_groups:
                            in_group = False
                            for g, group in enumerate(groups):
                                if list_r[i] in group and list_r[j] not in group:
                                    groups[g].append(list_r[j])
                                    in_group = True
                                elif list_r[i] not in group and list_r[j] in group:
                                    groups[g].append(list_r[i])
                                    in_group = True
                                elif list_r[i] in group and list_r[j] in group:
                                    in_group = True
                            if not in_group:
                                groups.append([list_r[i], list_r[j]])
                groups = self.fix_groups(groups)
                for group in groups:
                    cost_tot = 0
                    for a in group:
                        i = np.where(list_r == a)[0][0]
                        cost_tot += cost[i, i]
                    if cost_tot == 0:
                        env_groups.append(group)
                    else:
                        rc_groups.append(group)
                for i in range(n_el):
                    if cost[i, i] == 0:
                        continue
                    in_group = False
                    for group in rc_groups:
                        if list_r[i] in group:
                            in_group = True
                            break
                    if not in_group:
                        rc_atoms.append(list_r[i])
        self.output("            Reaction Center Analysis  \n"
             + "          ----------------------------\n", self.print_output)
        if len(self.mono_el) > 0:
            self.output("     Element(s) treated as monovalent: " + " ".join(self.mono_el), self.print_output)
            self.output("     (These elements are excluded from the\n"
                 + "      reaction center analysis.)          \n", self.print_output)
        if len(rc_groups) > 0 or len(rc_atoms) > 0:
            self.output("   Atoms associated with the reaction center:    \n\n"
                 + "  Group | Element | No. reactant | No. product     \n"
                 + "  --------------------------------------------       ", self.print_output)
            for g, group in enumerate(rc_groups):
                g_str = str(g + 1)
                for i in group:
                    self.output("{:>5}{:>9}{:>14}{:>14}".format(g_str, self.labels_r[i], i + 1, atom_map[i] + 1), self.print_output)
                    g_str = ""
            for i in rc_atoms:
                self.output("{:>5}{:>9}{:>14}{:>14}".format("-", self.labels_r[i], i + 1, atom_map[i] + 1), self.print_output)
        else:
            self.output("          No reaction center was found.", self.print_output)
        for el in self.mono_el:
            group = []
            for i in range(self.n_atoms):
                if self.labels_r[i] == el:
                    group.append(i)
            env_groups.append(group)
        env_groups = sorted(env_groups, key=len)
        return atom_map, env_groups, rc_groups
    
    def initial_guess_cost(self, i, j):
        """ Calculates a lower limit for the number of formed and broken bonds for
            mapping atom i of the reactant structure onto atom j of the product structure.

        Args:
            i:          The atom in the reactant structure.
            j:          The atom in the product structure.

        Returns:
            final_cost: The final cost.
        """
        final_cost = 0
        if self.sphere_guess:
            env_r = []
            exclude_r = [i]
            env_p = []
            exclude_p = [j]
            n_sphere = 0
            while n_sphere < self.n_spheres:
                n_sphere += 1
                sphere = []
                for k in copy(exclude_r):
                    for l, bonded in enumerate(self.con_r[k]):
                        if l not in exclude_r and bonded:
                            sphere.append(self.labels_r[l])
                            exclude_r.append(l)
                env_r.append(sphere)
                sphere = []
                for k in copy(exclude_p):
                    for l, bonded in enumerate(self.con_p[k]):
                        if l not in exclude_p and bonded:
                            sphere.append(self.labels_p[l])
                            exclude_p.append(l)
                env_p.append(sphere)
            for (sphere_r, sphere_p) in zip(env_r, env_p):
                for e in np.unique(self.labels_r):
                    final_cost += abs(sphere_r.count(e) - sphere_p.count(e))
        else:
            env_r = np.where(self.con_r[i])[0]
            env_p = np.where(self.con_p[j])[0]
            len_env = max([len(env_r), len(env_p)])
            cost = np.ones((len_env, len_env), dtype=int)
            for k, a in enumerate(env_r):
                for l, b in enumerate(env_p):
                    if self.labels_r[a] != self.labels_p[b]:
                        cost[k, l] = 2
                    else:
                        el_a = []
                        el_b = []
                        for m in range(self.n_atoms):
                            if self.con_r[a, m]:
                                el_a.append(self.labels_r[m])
                            if self.con_p[b, m]:
                                el_b.append(self.labels_p[m])
                        cost[k, l] = 0
                        for e in np.unique(self.labels_r):
                            cost[k, l] += abs(el_a.count(e) - el_b.count(e))
                    if cost[k, l] > 2:
                        cost[k, l] = 2
            env_map = self.linear_assignment_solver(cost)
            final_cost = cost[np.arange(len_env), env_map].sum()
        return final_cost
        
    def refine_environment(self, atom_map, env_groups, rc_groups):
        """ Refines environment groups.

        Args:
            atom_map:   The atom map.
            env_groups: The environment groups.
            rc_groups:  The reaction center groups.

        Returns:
            atom_map:   Updated atom map.
            env_groups: Updated environment groups.
        """
        update = True
        while update:
            update = False
            con_p = np.copy(self.con_p)
            con_p[:, :] = con_p[atom_map, :]
            con_p[:, :] = con_p[:, atom_map]
            group_atoms = []
            for group in (env_groups + rc_groups):
                group_atoms.extend(group)
            non_group_atoms = []
            for i in range(self.n_atoms):
                if i not in group_atoms and self.labels_r[i] not in self.mono_el:
                    non_group_atoms.append(i)
            for group in env_groups:
                cost = np.zeros((len(group), len(group)), dtype=int)
                for i, a in enumerate(group):
                    for j, b in enumerate(group):
                        n_group = np.zeros((len(env_groups + rc_groups)), dtype=int)
                        for c in range(self.n_atoms):
                            if not self.con_r[a, c] and not con_p[b, c] and c in [a, b]:
                                continue
                            elif c in group_atoms:
                                for g, other_group in enumerate(env_groups + rc_groups):
                                    if c in other_group:
                                        n_group[g] += int(self.con_r[a, c]) - int(con_p[b, c])
                                        break
                            else:
                                cost[i, j] += abs(int(self.con_r[a, c]) - int(con_p[b, c]))
                        cost[i, j] += abs(n_group).sum()
                if cost.sum() == 0:
                    continue
                sub_map = self.linear_assignment_solver(cost)
                tmp = np.copy(atom_map)
                for i in range(len(group)):
                    atom_map[group[i]] = tmp[group[sub_map[i]]]
                cost[:, :] = cost[:, sub_map]
                sub_groups = []
                for i in range(len(group)):
                    for j in range(i + 1, len(group)):
                        if (cost[i, j] + cost[j, i] - cost[i, i] - cost[j, j]) == 0:
                            in_group = False
                            for g, sub_group in enumerate(sub_groups):
                                if group[i] in sub_group and group[j] not in sub_group:
                                    sub_groups[g].append(group[j])
                                    in_group = True
                                elif group[j] in sub_group and group[i] not in sub_group:
                                    sub_groups[g].append(group[i])
                                    in_group = True
                                elif group[i] in sub_group and group[j] in sub_group:
                                    in_group = True
                            if not in_group:
                                sub_groups.append([group[i], group[j]])
                sub_groups = self.fix_groups(sub_groups)
                if len(sub_groups) == 1 and group == sub_groups[0]:
                    continue
                update = True
                env_groups.remove(group)
                env_groups.extend(sub_groups)
                break
        env_groups = sorted(env_groups, key=len)
        return atom_map, env_groups
    
    def run_optimization(self, atom_map, env_groups, rc_groups):
        """ Permutes reaction center groups and optimizes mappings for all atom maps.

        Args:
            atom_map:   The atom map.
            env_groups: The environment groups.
            rc_groups:   The reaction center groups.

        Returns:
            atom_maps:  Constructed atom maps.
        """
        self.output("\n"
             + "            Atom Mapping Optimization  \n"
             + "          -----------------------------\n", self.print_output)
        atom_maps = []
        if len(rc_groups) > 0:
            self.output("         Searching for atom maps through\n"
                 + "         permutation of groups.\n", self.print_output)
            for g, rc_group in enumerate(rc_groups):
                if math.factorial(len(rc_group)) > self.max_maps:
                    self.output("         Maximum number of atom maps is\n"
                         + "         exceeded for group {}.".format(g + 1), self.print_output)
                    if self.tau > 0:
                        self.output("         You could decrease --tau.\n", self.print_output)
                    else:
                        self.output("         You could increase --max_maps.\n", self.print_output)
                else:
                    groups = rc_groups + env_groups
                    groups.remove(rc_group)
                    tmp_atom_maps = self.permute_reaction_center(atom_map, groups, rc_group)
                    for i in range(len(tmp_atom_maps)):
                        tmp_atom_maps[i] = self.optimize_mapping(tmp_atom_maps[i], groups)
                    atom_maps.extend(tmp_atom_maps)
        if len(atom_maps) == 0:
            if len(rc_groups) > 0:
                self.output("     Continue without permutation of groups!\n", self.print_output)
            groups = rc_groups + env_groups
            atom_maps.append(self.optimize_mapping(atom_map, groups))
        return atom_maps
        
    def permute_reaction_center(self, atom_map, env_groups, rc_group):
        """ Permutes reaction center group with highest cost.

        Args:
            atom_map:   The atom map.
            env_groups: The environment groups.
            rc_group:   The reaction center group to be permuted.

        Returns:
            atom_maps:  Constructed atom maps.
        """
        atom_maps = []
        pers = self.permute(rc_group)
        tmp_atom_maps = []
        for p in pers:
            new_atom_map = np.copy(atom_map)
            for (i, j) in zip(rc_group, p):
                new_atom_map[i] = atom_map[j]
            tmp_atom_maps.append(new_atom_map)
        group_atoms = []
        for group in env_groups:
            group_atoms.extend(group)
        costs = []
        for atom_map in tmp_atom_maps:
            cost = 0
            con_p = np.copy(self.con_p)
            con_p[:, :] = con_p[atom_map, :]
            con_p[:, :] = con_p[:, atom_map]
            for i in rc_group:
                n_group = np.zeros((len(env_groups)), dtype=int)
                for j in range(self.n_atoms):
                    if (not self.con_r[i, j] and not con_p[i, j]) or i == j:
                        continue
                    elif j in group_atoms:
                        for g, other_group in enumerate(env_groups):
                            if j in other_group:
                                n_group[g] += int(self.con_r[i, j]) - int(con_p[i, j])
                                break
                    else:
                        cost += abs(int(self.con_r[i, j]) - int(con_p[i, j]))
                cost += abs(n_group).sum()
            costs.append(cost)
        for i in range(len(tmp_atom_maps)):
            if (costs[i] - min(costs)) <= self.tau:
                atom_maps.append(tmp_atom_maps[i])
        return atom_maps

    def optimize_mapping(self, atom_map, env_groups):
        """ Optimizes atom mapping.

        Args:
            atom_map:   The atom map.
            env_groups: The environment groups.

        Returns:
            atom_map:   Updated atom map.
        """
        groups = copy(env_groups)
        while len(groups) > 0:
            update = False
            con_p = np.copy(self.con_p)
            con_p[:, :] = con_p[atom_map, :]
            con_p[:, :] = con_p[:, atom_map]
            group_atoms = []
            for group in groups:
                group_atoms.extend(group)
            non_group_atoms = []
            for i in range(self.n_atoms):
                if i not in group_atoms:
                    non_group_atoms.append(i)
            for group in groups:
                cost = np.zeros((len(group), len(group)), dtype=int)
                for i, a in enumerate(group):
                    for j, b in enumerate(group):
                        n_group = np.zeros((len(groups)), dtype=int)
                        for c in range(self.n_atoms):
                            if not self.con_r[a, c] and not con_p[b, c] and c in [a, b]:
                                continue
                            elif c in group_atoms:
                                for g, other_group in enumerate(groups):
                                    if c in other_group:
                                        n_group[g] += int(self.con_r[a, c]) - int(con_p[b, c])
                                        break
                            else:
                                cost[i, j] += abs(int(self.con_r[a, c]) - int(con_p[b, c]))
                        cost[i, j] += abs(n_group).sum()
                sub_map = self.linear_assignment_solver(cost)
                tmp = np.copy(atom_map)
                for i in range(len(group)):
                    atom_map[group[i]] = tmp[group[sub_map[i]]]
                cost[:, :] = cost[:, sub_map]
                sub_groups = []
                for i in range(len(group)):
                    for j in range(i + 1, len(group)):
                        if (cost[i, j] + cost[j, i] - cost[i, i] - cost[j, j]) == 0:
                            in_group = False
                            for g, sub_group in enumerate(sub_groups):
                                if group[i] in sub_group and group[j] not in sub_group:
                                    sub_groups[g].append(group[j])
                                    in_group = True
                                elif group[j] in sub_group and group[i] not in sub_group:
                                    sub_groups[g].append(group[i])
                                    in_group = True
                                elif group[i] in sub_group and group[j] in sub_group:
                                    in_group = True
                            if not in_group:
                                sub_groups.append([group[i], group[j]])
                sub_groups = self.fix_groups(sub_groups)
                amb_dihs = False
                amb_sub_group = []
                for sub_group in sub_groups:
                    cost = self.get_dihedral_cost(atom_map, sub_group, non_group_atoms)
                    if np.sum(cost == 0.0) > len(sub_group):
                        if group == sub_group:
                            amb_dihs = True
                        elif sub_group not in amb_sub_group:
                            amb_sub_group.append(sub_group)
                    else:
                        sub_map = self.linear_assignment_solver(cost)
                        old_atom_map = np.copy(atom_map)
                        for i in range(len(sub_group)):
                            atom_map[sub_group[i]] = old_atom_map[sub_group[sub_map[i]]]
                if amb_dihs:
                    continue
                groups.remove(group)
                for sub_group in amb_sub_group:
                    groups.append(sub_group)
                update = True
                break
            if not update:
                groups = sorted(groups, key=len)
                per_group = groups[0]
                if math.factorial(len(per_group)) > self.max_maps:
                    groups.remove(per_group)
                    atom_map = self.steepest_descent(atom_map, per_group, groups)
                    continue
                else:
                    pers = self.permute(per_group)
                    groups.remove(per_group)
                    atom_maps = []
                    for p in pers:
                        new_atom_map = np.copy(atom_map)
                        for (i, j) in zip(per_group, p):
                            new_atom_map[i] = atom_map[j]
                        atom_maps.append(new_atom_map)
                    min_rmsd = 0.0
                    min_cd = 0
                    for i in range(len(atom_maps)):
                        atom_maps[i] = self.optimize_mapping(atom_maps[i], groups)
                        con_p = np.copy(self.con_p)
                        con_p[:, :] = con_p[atom_maps[i], :]
                        con_p[:, :] = con_p[:, atom_maps[i]]
                        cd = abs(1*con_p - 1*self.con_r).sum()
                        dih_rmsd = self.get_dihedral_rmsd(atom_maps[i], per_group)
                        if i == 0:
                            min_rmsd = copy(dih_rmsd)
                            min_cd = copy(cd)
                            atom_map = atom_maps[i]
                        elif cd < min_cd:
                            min_rmsd = copy(dih_rmsd)
                            min_cd = copy(cd)
                            atom_map = atom_maps[i]
                        elif dih_rmsd < min_rmsd and cd == min_cd:
                            min_rmsd = copy(dih_rmsd)
                            atom_map = atom_maps[i]
                    break
        return atom_map

    def get_dihedral_cost(self, atom_map, sub_group, allowed_atoms):
        """ Calculates the cost matrix to minimize the squared error of dihedrals.

        Args:
            atom_map:   The atom map.
            sub_group:  The group which has to be mapped according to the dihedrals.
            allowed_atoms:
                        Atoms that are already mapped and hence allowed to be used for dihedrals.

        Returns:
            cost:       The cost matrix with squared errors of dihedrals.
        """
        con_p = np.copy(self.con_p)
        con_p[:, :] = con_p[atom_map, :]
        con_p[:, :] = con_p[:, atom_map]
        cost = np.zeros((len(sub_group), len(sub_group)), dtype=float)
        cost_oops = np.zeros((len(sub_group), len(sub_group)), dtype=float)
        for i in range(len(sub_group)):
            for j in range(len(sub_group)):
                for k in allowed_atoms:
                    if not self.con_r[sub_group[i], k] or not con_p[sub_group[j], k]:
                        continue
                    for l in allowed_atoms:
                        if not self.con_r[k, l] or not con_p[k, l]:
                            continue
                        for m in allowed_atoms:
                            if k != m and self.con_r[l, m] and con_p[l, m]:
                                iklm = self.reactant.get_dihedral_in_degrees([sub_group[i] + 1, k + 1, l + 1, m + 1])
                                jklm = self.product.get_dihedral_in_degrees([atom_map[sub_group[j]] + 1, atom_map[k] + 1, atom_map[l] + 1, atom_map[m] + 1])
                                if iklm and jklm:
                                    cost[i, j] += ((abs(iklm - jklm) > 180.0) *  360.0 - abs(iklm - jklm))**2
                            elif len(sub_group) == 2 and l != m and self.con_r[k, m] and con_p[k, m]:
                                ilmk = self.reactant.get_dihedral_in_degrees([sub_group[i] + 1, l + 1, m + 1, k + 1])
                                jlmk = self.product.get_dihedral_in_degrees([atom_map[sub_group[j]] + 1, atom_map[l] + 1, atom_map[m] + 1, atom_map[k] + 1])
                                if ilmk and jlmk:
                                    cost_oops[i, j] += ((abs(ilmk - jlmk) > 180.0) *  360.0 - abs(ilmk - jlmk))**2
        if np.sum(cost_oops > self.oop_thresh) > 1:
            cost = np.copy(cost_oops)
        return cost

    def steepest_descent(self, atom_map, opt_group, env_groups):
        """ Optimizes atom mapping.

        Args:
            atom_map:   The atom map.
            opt_group:  The group to be optimized.
            env_groups: The environment groups.

        Returns:
            atom_map:   Updated atom map.
        """
        group_atoms = []
        for group in env_groups:
            group_atoms.extend(group)
        it = 0
        while it < self.max_iter:
            it += 1
            con_p = np.copy(self.con_p)
            con_p[:, :] = con_p[atom_map, :]
            con_p[:, :] = con_p[:, atom_map]
            cost = np.zeros((len(opt_group), len(opt_group)), dtype=int)
            for i, a in enumerate(opt_group):
                for j, b in enumerate(opt_group):
                    n_group = np.zeros((len(env_groups)), dtype=int)
                    for c in range(self.n_atoms):
                        if not self.con_r[a, c] and not con_p[b, c] and c in [a, b]:
                            continue
                        elif c in group_atoms:
                            for g, other_group in enumerate(env_groups):
                                if c in other_group:
                                    n_group[g] += int(self.con_r[a, c]) - int(con_p[b, c])
                                    break
                        else:
                            cost[i, j] += abs(int(self.con_r[a, c]) - int(con_p[b, c]))
                    cost[i, j] += abs(n_group).sum()
            per = []
            max_grad = 0
            for i, a in enumerate(opt_group):
                for j, b in enumerate(opt_group):
                    if i < j:
                        grad = cost[i, i] + cost[j, j] - cost[i, j] - cost[j, i]
                        if grad > max_grad:
                            max_grad = copy(grad)
                            per = [a, b]
            if len(per) == 2:
                atom_map[per[0]], atom_map[per[1]] = atom_map[per[1]], atom_map[per[0]]
            else:
                break
        if it == self.max_iter:
            exit("Maximum number of iterations reached. Steepest descent not converged!")
        return atom_map
    
    def evaluate_results(self, atom_maps, write_xyz):
        """ Evaluates all remaining atom maps.

        Args:
            atom_maps:  The remaining atom maps.
            write_xyz:  If XYZ files are written.

        Returns:
            final_atom_maps:
                        The final atom maps.
        """
        cds = []
        deltas = []
        dih_rmsds = []
        inv_centers = []
        for atom_map in atom_maps:
            con_p = np.copy(self.con_p)
            con_p[:, :] = con_p[atom_map, :]
            con_p[:, :] = con_p[:, atom_map]
            delta = 1*con_p - 1*self.con_r
            cds.append(abs(delta).sum() / 2)
            deltas.append(delta)
            dih_rmsds.append(self.get_dihedral_rmsd(atom_map, np.arange(self.n_atoms)))
            inv_centers.append(self.check_stereocenters(atom_map))
        final_atom_maps = []
        final_deltas = []
        final_dih_rmsds = []
        final_inv_centers = []
        min_cd = min(cds)
        for i in range(len(atom_maps)):
            skip = False
            for j in range(len(atom_maps)):
                if np.all(atom_maps[i] == atom_maps[j]) and i < j:
                    skip = True
                    break
                elif np.all(deltas[i] == deltas[j]) and len(inv_centers[i]) > len(inv_centers[j]):
                    skip = True
                    break
                elif np.all(deltas[i] == deltas[j]) and len(inv_centers[i]) == len(inv_centers[j]) and dih_rmsds[i] > dih_rmsds[j]:
                    skip = True
                    break
            if (cds[i] - min_cd) <= self.tau and not skip:
                final_atom_maps.append(atom_maps[i])
                final_deltas.append(deltas[i])
                final_dih_rmsds.append(dih_rmsds[i])
                final_inv_centers.append(inv_centers[i])
        for atom_map, delta, dih_rmsd, inv_centers in zip(final_atom_maps, final_deltas, final_dih_rmsds, final_inv_centers):
            if (write_xyz):
                filename = self.product.filename.split(".xyz")[0] + "_" + str(self.product.file_index) + ".xyz"
                self.product.write_xyz(filename, atom_map)
                self.product.file_index += 1
                self.output("  => {} has been created.\n\n".format(filename)
                     + "      Properties related to this atom map:", self.print_output)
            if abs(delta).sum() > 0:
                self.output("\n"
                     + "  Change | Element | No. reactant | No. product\n"
                     + "  ---------------------------------------------  ", self.print_output)
                for i in range(self.n_atoms):
                    for j in range(i + 1, self.n_atoms):
                        if delta[i, j] == 1:
                            self.output("  Formed {:>5}-{:<5} {:>6}-{:<5} {:>8}-{:<5}".format(
                                self.labels_r[i], self.labels_r[j], i + 1, j + 1, atom_map[i] + 1, atom_map[j] + 1), self.print_output)
                        elif delta[i, j] == -1:
                            self.output("  Broke  {:>5}-{:<5} {:>6}-{:<5} {:>8}-{:<5}".format(
                                self.labels_r[i], self.labels_r[j], i + 1, j + 1, atom_map[i] + 1, atom_map[j] + 1), self.print_output)
            self.output("\n      RMSD of dihedrals in degrees: {:5.2f}\n".format(dih_rmsd), self.print_output)
            if len(inv_centers) > 0:
                self.output("    Inverted stereocenters have been traced!\n\n"
                     + "      Element | No. reactant | No. product\n"
                     + "      ------------------------------------", self.print_output)
                for i in inv_centers:
                    self.output("{:>11}{:>12}{:>15}".format(self.labels_r[i], i + 1, atom_map[i] + 1), self.print_output)
                self.output("", self.print_output)
        return final_atom_maps
    
    def get_dihedral_rmsd(self, atom_map, allowed_atoms):
        """ Calculates the RMSD of dihedrals for a given atom map.

        Args:
            atom_map:   The atom map.
            allowed_atoms:
                        Atoms that are already mapped and hence allowed to be used for dihedrals.

        Returns:
            dih_rmsd:   The RMSD of dihedrals.
        """
        n_dihs = 0
        dih_error = 0.0
        con_p = np.copy(self.con_p)
        con_p[:, :] = con_p[atom_map, :]
        con_p[:, :] = con_p[:, atom_map]
        for i in allowed_atoms:
            for j in range(self.n_atoms):
                if not self.con_r[i, j] or not con_p[i, j]:
                    continue
                for k in range(self.n_atoms):
                    if i == k or not self.con_r[j, k] or not con_p[j, k]:
                        continue
                    for l in range(self.n_atoms):
                        if l == j or not self.con_r[k, l] or not con_p[k, l]:
                            continue
                        ijkl_r = self.reactant.get_dihedral_in_degrees([i + 1, j + 1, k + 1, l + 1])
                        ijkl_p = self.product.get_dihedral_in_degrees([atom_map[i] + 1, atom_map[j] + 1, atom_map[k] + 1, atom_map[l] + 1])
                        if ijkl_r and ijkl_p:
                            dih_error += (ijkl_r - ijkl_p)**2
                        n_dihs += 1
        if n_dihs == 0:
            return float('nan')
        dih_rmsd = np.sqrt(dih_error / n_dihs)
        return dih_rmsd

    def check_stereocenters(self, atom_map):
        """ Checks if stereocenters have been inverted.

        Args:
            atom_map:   The atom map.

        Returns:
            inverted_centers:
                        The atoms with inverted stereocenters.
        """
        inverted_centers = []
        con_p = np.copy(self.con_p)
        con_p[:, :] = con_p[atom_map, :]
        con_p[:, :] = con_p[:, atom_map]
        for i in range(self.n_atoms):
            intersect = np.intersect1d(np.where(self.con_r[i])[0], np.where(con_p[i])[0])
            if len(intersect) == 4:
                [j, k, l, m] = intersect.tolist()
                jlmi_r = self.reactant.get_dihedral_in_degrees([j + 1, l + 1, m + 1, i + 1])
                klmi_r = self.reactant.get_dihedral_in_degrees([k + 1, l + 1, m + 1, i + 1])
                jlmi_p = self.product.get_dihedral_in_degrees([atom_map[j] + 1, atom_map[l] + 1, atom_map[m] + 1, atom_map[i] + 1])
                klmi_p = self.product.get_dihedral_in_degrees([atom_map[k] + 1, atom_map[l] + 1, atom_map[m] + 1, atom_map[i] + 1])
                if not jlmi_r or not klmi_r or not jlmi_p or not klmi_p:
                    continue
                elif jlmi_r < klmi_r and jlmi_p < klmi_p:
                    continue
                elif jlmi_r > klmi_r and jlmi_p > klmi_p:
                    continue
                else:
                    inverted_centers.append(i)
        return inverted_centers

    def linear_assignment_solver(self, cost):
        """ This tool is used to solve linear assignment problems.

        Args:
            cost:       A cost matrix.

        Returns:
            col:        The solution vector to the assignment problem.
        """
        row, col = linear_sum_assignment(cost)
        if not (row == np.arange(cost.shape[0])).all():
            exit("Error in linear assignment solver.")
        while True:
            lap_check = True
            new_cost = np.copy(cost)[:, col]
            for i in range(cost.shape[0]):
                for j in range(cost.shape[0]):
                    if (new_cost[i, j] + new_cost[j, i]) < (new_cost[i, i] + new_cost[j, j]):
                        col[i], col[j] = col[j], col[i]
                        lap_check = False
                        break
                else:
                    continue
                break
            if lap_check:
                break
        return col
        
    def fix_groups(self, groups):
        """ This tool is used to prevent incorrect group definitions.

        Args:
            groups:     Groups of atoms that are treated indistinguishable.

        Returns:
            groups:     Groups of atoms that are treated indistinguishable.
        """
        fixed = False
        while not fixed:
            fixed = True
            for g in range(len(groups)):
                for a in groups[g]:
                    for h in range(g + 1, len(groups)):
                        if g != h and a in groups[h]:
                            groups[g].extend(groups[h])
                            groups[g] = sorted(list(set(groups[g])))
                            groups.remove(groups[h])
                            fixed = False
                            break
                    else:
                        continue
                    break
                else:
                    continue
                break
        for g in range(len(groups)):
            groups[g] = sorted(groups[g])
        groups = sorted(groups, key=len)
        return groups

    def permute(self, group):
        """ This tool is used to construct all possible permutations of a group.

        Args:
            group:      A group of atoms that are treated indistinguishable.

        Returns:
            pers:       All possible permutations of this group.
        """
        if len(group) == 0:
            return []
        if len(group) == 1:
            return [group]
        pers = []
        for i in range(len(group)):
            j = group[i]
            reduced_group = group[:i] + group[i+1:]
            for p in self.permute(reduced_group):
                pers.append([j] + p)
        return pers

    def output(self, string, level):
        """ This tool is used to control the print level. """
        if level:
            print(string)