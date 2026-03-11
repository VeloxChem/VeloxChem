# Four edits to mmforcefieldgenerator.py to add option-B force field selection
# Apply them in order.  Each block shows the OLD lines and the NEW lines.
# ============================================================================

# EDIT 1 — __init__: declare the attribute
# ----------------------------------------
# Find this line (last line of __init__):
        self.fitting_summary = None

# Replace with:
        self.fitting_summary = None

        # Force field selection: 'gaff' (default) or 'opls'
        # Set via   mmgen.force_field = 'opls'   before calling create_topology.
        self.force_field = 'gaff'


# ============================================================================

# EDIT 2 — update_settings: register the keyword so it can be set from input
# ---------------------------------------------------------------------------
# Find this block inside the ffg_keywords dict:
            'keep_files': 'bool',
        }

# Replace with:
            'keep_files': 'bool',
            'force_field': 'str',
        }


# ============================================================================

# EDIT 3 — populate_atoms: add OPLS branch before the GAFF block
# ---------------------------------------------------------------
# Find (first few lines of the per-atom loop):
        for i, atom_type in enumerate(self.atom_types_dict.values()):
            atom_type_found = False

            if 'gaff' in atom_type:
                # Note: need strip() for converting e.g. 'c ' to 'c'
                gafftype = atom_type['gaff'].strip()

# Replace with:
        for i, atom_type in enumerate(self.atom_types_dict.values()):
            atom_type_found = False

            if self.force_field == 'opls' and 'opls' in atom_type:
                oplstype = atom_type.get('opls')
                if oplstype is None:
                    # No OPLS-AA equivalent exists for this chemical environment.
                    # The 'opls' key was set to None in atomtypeidentifier.py to
                    # flag this gap explicitly.  Fall back to GAFF parameters and
                    # warn the user.
                    gafftype = atom_type['gaff'].strip() if 'gaff' in atom_type else '?'
                    warnmsg = (
                        f'MMForceFieldGenerator: no OPLS-AA type for GAFF type '
                        f'{gafftype!r} (atom index {i + 1}). '
                        'Falling back to GAFF sigma/epsilon parameters.')
                    self.ostream.print_warning(warnmsg)
                # Whether oplstype is a valid string or None, we fall through to
                # the GAFF block below for the actual sigma/epsilon lookup.
                # TODO: Once an OPLS-AA parameter file is bundled with VeloxChem,
                # add a sigma/epsilon lookup keyed on oplstype here and set
                # atom_type_found = True so the GAFF block is skipped.

            if 'gaff' in atom_type:
                # Note: need strip() for converting e.g. 'c ' to 'c'
                gafftype = atom_type['gaff'].strip()


# ============================================================================

# EDIT 4 — write_top: use comb_rule=3 (geometric mixing) for OPLS-AA
# -------------------------------------------------------------------
# OPLS-AA uses geometric-mean mixing (comb_rule 3 in GROMACS), whereas GAFF/
# AMBER uses Lorentz-Berthelot (comb_rule 2).  Without this change, an OPLS
# topology would use the wrong combining rule.
#
# Find:
                f_top.write('{}{:16}{:>18}{:21.6f}{:10.6f}\n'.format(
                    self.nbfunc, self.comb_rule, gen_pairs, self.fudgeLJ,
                    self.fudgeQQ))

# Replace with:
                comb_rule = 3 if self.force_field == 'opls' else self.comb_rule
                f_top.write('{}{:16}{:>18}{:21.6f}{:10.6f}\n'.format(
                    self.nbfunc, comb_rule, gen_pairs, self.fudgeLJ,
                    self.fudgeQQ))
