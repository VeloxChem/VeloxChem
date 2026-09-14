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
Formatting for the metal site force field code.

A printer computes nothing. Everything it prints is handed to it, so the
table always describes what the run actually did rather than the settings
that asked for it, and so nothing here has to be kept in step with the
caller's state.

This is the leaf of the package: it imports nothing from core, the builder
or the manager, which is what lets all three import it.
"""

from pathlib import Path
import numpy as np
import math

from ..outputstream import OutputStream

# ----------------------------------------------------------------------
# the shared line shapes
# ----------------------------------------------------------------------


def stream(ostream):
    """
    Returns the stream to report through.

    A function reports only when it is given somewhere to report to, so a
    caller that wants the numbers and not the commentary simply leaves the
    argument out.

    :param ostream:
        The output stream, or None.

    :return:
        The given stream, or a silent one.
    """

    return OutputStream(None) if ostream is None else ostream


def param(label, value, label_width=26, value_width=20):
    """
    Formats one parameter line with fixed label and value widths.

    print_header centers text, so all lines need the same total length to
    appear left-aligned relative to each other.

    :param label:
        The parameter name.
    :param value:
        The parameter value.
    :param label_width:
        The width of the label field.
    :param value_width:
        The width of the value field.

    :return:
        The formatted line.
    """

    return f'{label:<{label_width}} : {str(value):>{value_width}}'


def print_param_list(label, items, value_width=20, ostream=None):
    """
    Prints a list of values as parameter lines of uniform width.

    print_header centers each line, so a value that overflows the field
    would make its line start further left than the others. Long lists are
    therefore wrapped over several lines, with the label only on the
    first.

    :param label:
        The parameter name.
    :param items:
        The values to list.
    :param value_width:
        The width of the value field.
    :param ostream:
        The output stream, or None to print nothing.
    """

    ostream = stream(ostream)

    chunks = []
    current = ''

    for item in items:
        candidate = item if not current else f'{current}, {item}'
        if len(candidate) > value_width and current:
            chunks.append(current + ',')
            current = item
        else:
            current = candidate

    if current:
        chunks.append(current)

    for i, chunk in enumerate(chunks):
        ostream.print_header(param(label if i == 0 else '', chunk))


def print_section(title, ostream=None):
    """
    Prints a title underlined to its own length.

    Hand-counting that length is how an underline ends up one character
    short of the title it sits under.

    :param title:
        The title.
    :param ostream:
        The output stream, or None to print nothing.
    """

    ostream = stream(ostream)

    ostream.print_header(title)
    ostream.print_header(len(title) * '-')


# ----------------------------------------------------------------------
# templates
# ----------------------------------------------------------------------


def print_template(template, bonds, angles, ostream=None):
    """
    Prints what one template holds.

    :param template:
        The template.
    :param bonds:
        Its metal bond keys.
    :param angles:
        Its metal angle keys.
    :param ostream:
        The output stream, or None to print nothing.
    """

    ostream = stream(ostream)

    labels = template['molecule'].get_labels()
    metals = ', '.join(labels[index] for index in template['metal_indices'])

    ostream.print_blank()
    print_section(f'Template {template["name"]}', ostream)
    ostream.print_header(param('geometry', template['geometry_kind']))
    ostream.print_header(param('atoms', template['molecule'].number_of_atoms()))
    ostream.print_header(param('metal centers', metals))
    ostream.print_header(
        param('capping hydrogens', len(template['cap_indices'])))
    ostream.print_header(param('metal bonds', len(bonds)))
    ostream.print_header(param('metal angles', len(angles)))
    ostream.print_header(
        param('total charge', f'{float(np.sum(template["charges"])):+.3f}'))
    ostream.print_blank()
    ostream.print_info(f'Loaded from {template["folder"]}')
    ostream.flush()


def print_templates(templates, ostream=None):
    """
    Prints every template that is loaded.

    :param templates:
        The templates, by name.
    :param ostream:
        The output stream, or None to print nothing.
    """

    ostream = stream(ostream)

    ostream.print_header(f'Loaded templates ({len(templates)})')
    ostream.print_header(60 * '-')
    valstr = '{:>24} | {:>7} | {:>7} | {:>13}'.format('name', 'atoms', 'metals',
                                                      'geometry')
    ostream.print_header(valstr)
    ostream.print_header(60 * '-')

    for name, template in templates.items():
        valstr = '{:>24} | {:>7} | {:>7} | {:>13}'.format(
            name[:24], template['molecule'].number_of_atoms(),
            len(template['metal_indices']), template['geometry_kind'])
        ostream.print_header(valstr)

    ostream.print_blank()
    ostream.flush()


# ----------------------------------------------------------------------
# shoehorning
# ----------------------------------------------------------------------


def print_shoehorn_header(template_name,
                          n_template_residues,
                          site_residues,
                          max_include_radius,
                          ostream=None):
    """
    Says what a shoehorning is starting from.

    :param template_name:
        The name of the template being edited onto.
    :param n_template_residues:
        How many residues the template coordinates.
    :param site_residues:
        The labels of the residues the site holds.
    :param max_include_radius:
        How far out from a metal center a residue may be picked up from.
    :param ostream:
        The output stream, or None to print nothing.
    """

    ostream = stream(ostream)

    ostream.print_blank()
    print_section(f'Shoehorning the site into {template_name}', ostream)
    ostream.print_blank()
    ostream.print_header(param('template residues', n_template_residues))
    ostream.print_header(param('site residues', len(site_residues)))
    ostream.print_header(param('search radius', f'{max_include_radius:.1f} A'))
    ostream.print_blank()
    ostream.print_info(f'Site holds: {", ".join(site_residues)}')
    ostream.print_blank()
    ostream.flush()


def print_shoehorn_summary(template_name, modes, variants, ostream=None):
    """
    Says what the site was made into.

    The coordination is what a shoehorning is for, so it is printed the way
    the builder prints it rather than left to be asked for: which residue
    ended up on which metal center, how far out, and what it is protonated
    as.

    :param template_name:
        The name of the template it was edited onto.
    :param modes:
        The binding modes, bound once by the caller.
    :param variants:
        The protonation of the site's own residues, as (label, variant)
        pairs in the order they should read.
    :param ostream:
        The output stream, or None to print nothing.
    """

    ostream = stream(ostream)

    ostream.print_blank()
    ostream.print_info(
        f'The site is now built the way {template_name} is built.')
    ostream.print_blank()

    for metal in modes['metals']:
        bonds = []
        for ligand in modes['ligands']:
            for index, distance in zip(ligand['metals'], ligand['distances']):
                if index != metal['index']:
                    continue
                bonds.append(f'{ligand["residue"]} {ligand["atom"]} '
                             f'{distance:.2f} A')
        ostream.print_info(f'  {metal["element"]} {metal["index"]}: ' +
                           (', '.join(sorted(bonds)) if bonds else 'nothing'))

    protonation = ', '.join(f'{label} {variant}' for label, variant in variants)

    ostream.print_info(f'  protonation: {protonation}')
    ostream.print_blank()
    ostream.print_info('Call compare_active_site() to measure it, then '
                       'build_ff_from_template.')
    ostream.print_blank()
    ostream.flush()


# ----------------------------------------------------------------------
# comparison
# ----------------------------------------------------------------------


def print_spec(title, described, bridging, ostream=None):
    """
    Prints what a site is made of: which residues it holds and which of them
    coordinate which metal.

    The residues are named by their formula, which is for reading; two sites
    are compared on the keys behind them.

    :param title:
        What the block is describing.
    :param described:
        A described active site or a template.
    :param bridging:
        Its residue nodes that coordinate more than one metal center, from
        matching.bridging_nodes -- the same answer the spec two sites are
        compared by is built from.
    :param ostream:
        The output stream, or None to print nothing.
    """

    ostream = stream(ostream)

    labels = described['molecule'].get_labels()
    coarse = described['coarse_topology']

    def named(nodes):
        return ', '.join(
            sorted(f'{coarse.nodes[node]["formula"]}/'
                   f'{coarse.nodes[node]["key"][:6]}' for node in nodes))

    ostream.print_blank()
    ostream.print_info(f'Site spec, {title}:')

    for metal in described['metal_indices']:
        node = ('metal', metal)
        ostream.print_info(
            f'  {labels[metal]}{metal}: {named(coarse.neighbors(node))}')

    if bridging:
        ostream.print_info(f'  bridging: {named(bridging)}')

    ostream.flush()


def ic_cell(ic_rmsd, name):
    """
    Formats one internal coordinate type for a table cell.

    :param ic_rmsd:
        The deviations, as get_ic_rmsd reports them.
    :param name:
        The type to format.

    :return:
        The cell.
    """

    if ic_rmsd is None:
        return ''

    found = ic_rmsd.get(name)

    if found is None:
        return ''

    return f'{found["rms"]:.2f} / {found["max"]:.2f}'


def print_template_comparison(name, entry, spec=None, ostream=None):
    """
    Prints the numbers and the verdicts of one template.

    :param name:
        The name of the template.
    :param entry:
        What compare_active_site measured for it.
    :param spec:
        The template and its residue nodes, for the case where the site
        coordinates a different set of residues and the two specs are worth
        reading side by side. None when nothing was measured for another
        reason.
    :param ostream:
        The output stream, or None to print nothing.
    """

    ostream = stream(ostream)

    ostream.print_blank()

    if entry['status'] == 'composition':
        ostream.print_info(
            f'{name}: holds different atoms, so nothing was measured.')
        ostream.flush()
        return

    if entry['status'] == 'spec':
        ostream.print_info(f'{name}: coordinates a different set of residues, '
                           'so nothing was measured.')
        if spec is not None:
            print_spec(f'{name} holds', *spec, ostream=ostream)
        ostream.flush()
        return

    bonds = entry['metal_bonds']
    summary = (f'{bonds["shared"]} metal bond(s) shared, within '
               f'{bonds["deviation"]:.3f} A')
    if bonds['template_only']:
        summary += f', {bonds["template_only"]} only in the template'
    if bonds['query_only']:
        summary += f', {bonds["query_only"]} only in the structure'

    ostream.print_info(f'{name}: {entry["n_mappings"]} atom mapping(s) from '
                       f'{entry["n_coarse_mappings"]} coarse mapping(s), '
                       f'{summary}')
    ostream.print_blank()

    row = '{:>19} | {:>5} | {:>7} | {:>7} | {:>14} | {:>14} | {:>14}'
    ostream.print_header(
        row.format('region', 'atoms', 'RMSD', 'heavy', 'bonds rms/max',
                   'angles rms/max', 'dihed rms/max'))
    ostream.print_header(98 * '-')

    for region, found in entry['regions'].items():
        ostream.print_header(
            row.format(region, found['atoms'], f'{found["rmsd"]:.3f}',
                       f'{found["rmsd_heavy"]:.3f}',
                       ic_cell(found['ic_rmsd'], 'bonds'),
                       ic_cell(found['ic_rmsd'], 'angles'),
                       ic_cell(found['ic_rmsd'], 'dihedrals')))

    ostream.print_blank()
    ostream.flush()


def print_comparison_summary(results, ranked_on, scores, ostream=None):
    """
    Ranks the templates on what a selection is decided by.

    :param results:
        The last comparison, from compare_active_site.
    :param ranked_on:
        The (region, ic type, measure) a selection is ranked on.
    :param scores:
        That measure per template, computed by the caller.
    :param ostream:
        The output stream, or None to print nothing.
    """

    ostream = stream(ostream)

    region, ic_type, measure = ranked_on
    heavy = not results['include_hydrogens']

    def rmsd(entry):
        found = entry['regions'].get(region)
        if found is None:
            return None
        return found['rmsd_heavy'] if heavy else found['rmsd']

    order = sorted(results['templates'].items(),
                   key=lambda item: scores[item[0]])

    ostream.print_blank()
    print_section(f'Ranked on the {region} {ic_type} {measure}', ostream)

    row = '{:>24} | {:>14} | {:>11} | {:>9} | {:>10}'
    ostream.print_header(
        row.format('template', 'status', f'{ic_type} {measure}', 'RMSD',
                   'metal bond'))
    ostream.print_header(78 * '-')

    for name, entry in order:
        found = rmsd(entry)
        if found is None:
            ostream.print_header(
                row.format(name[:24], entry['status'], '', '', ''))
            continue
        ostream.print_header(
            row.format(name[:24], entry['status'], f'{scores[name]:.3f}',
                       f'{found:.3f}',
                       f'{entry["metal_bonds"]["deviation"]:.3f}'))

    ostream.print_blank()
    ostream.flush()


def print_no_selection(decision, ostream=None):
    """
    Says which template came closest when none of them was good enough.

    :param decision:
        The decision, as _select_template makes it.
    """
    ostream = stream(ostream)
    closest = min(decision['scores'],
                  key=lambda name: decision['scores'][name],
                  default=None)

    if closest is not None and math.isfinite(decision['scores'][closest]):
        ostream.print_info(
            f'No template is within the {decision["criteria_name"]} '
            f'criteria. The closest is {closest}: '
            f'{decision["verdicts"][closest]}.')
    else:
        ostream.print_info(
            'No template describes this site: none of them maps onto all '
            'of its atoms.')

    ostream.print_info(
        "Set selection_criteria to 'loose' to widen what counts as a "
        'match, or build this site with MetalSiteForceFieldBuilder.')
    ostream.flush()


def print_comparison(results, specs, ranked_on, scores, ostream=None):
    """
    Prints everything compare_active_site measured.

    One table of numbers and one of verdicts per template that could be
    measured, and a closing summary ranking the templates by the region that
    is configured, so the closest one is visible without reading every table.

    :param results:
        The last comparison, from compare_active_site.
    :param specs:
        The (described, residue nodes) pair per site whose spec is printed:
        the structure under the key None, and a template under its name.
    :param ranked_on:
        The (region, ic type, measure) a selection is ranked on.
    :param scores:
        That measure per template, computed by the caller.
    :param ostream:
        The output stream, or None to print nothing.
    """

    ostream = stream(ostream)

    active_site = results['active_site']
    labels = active_site['molecule'].get_labels()
    metals = ', '.join(labels[index] for index in active_site['metal_indices'])

    ostream.print_blank()
    print_section('Comparison against every template', ostream)
    ostream.print_header(param('source', Path(results['source']).name))
    ostream.print_header(
        param('active site atoms', active_site['molecule'].number_of_atoms()))
    ostream.print_header(param('metal centers', metals))
    ostream.print_header(param('geometry', results['geometry']))
    ostream.print_header(
        param('measured over',
              'all atoms' if results['include_hydrogens'] else 'heavy atoms'))
    ostream.print_header(param('templates', len(results['templates'])))
    ostream.print_blank()
    ostream.print_info(f'Residues: {", ".join(active_site["residues"])}')

    print_spec('the structure holds', *specs[None], ostream=ostream)

    for name, entry in results['templates'].items():
        print_template_comparison(name,
                                  entry,
                                  spec=specs.get(name),
                                  ostream=ostream)

    print_comparison_summary(results, ranked_on, scores, ostream=ostream)


def print_selection(comparison,
                    decision,
                    rmsd_regions,
                    ic_types,
                    ranked_on,
                    ostream=None):
    """
    Prints how every template stands against the criteria, and which one was
    taken.

    The whole field is printed rather than the winner alone: whether the
    others are near misses or a long way off is what says how much the chosen
    one is worth.

    :param comparison:
        The last comparison, from compare_active_site.
    :param decision:
        The decision, as _select_template makes it.
    :param rmsd_regions:
        Every region there can be a column for.
    :param ic_types:
        The internal coordinate types, with their units.
    :param ranked_on:
        The (region, ic type, measure) a selection is ranked on.
    :param ostream:
        The output stream, or None to print nothing.
    """

    ostream = stream(ostream)

    regions = [
        region for region in rmsd_regions if decision['criteria'].get(region)
    ]

    ostream.print_blank()
    print_section(
        f'Choosing a template on the {decision["criteria_name"]} criteria',
        ostream)
    ostream.print_blank()

    for region in regions:
        thresholds = decision['criteria'][region]
        # only the measures the set actually holds, since either of them may
        # be left out of one
        measures = {
            name:
            ' / '.join(f'{measure} {limit:.2f}'
                       for measure, limit in given.items() if limit is not None)
            for name, given in thresholds.items() if given
        }
        limits = '; '.join(f'{name} {shown} {ic_types[name]}'
                           for name, shown in measures.items())
        ostream.print_header(param(region, limits, value_width=44))

    ostream.print_blank()

    # one column per region the criteria name, so a custom set of them prints
    # as readably as the two that come with the class
    row = ' | '.join(['{:>22}'] + ['{:>13}'] * len(regions) +
                     ['{:>26}', '{:>5}'])
    header = row.format('template', *[region[:13] for region in regions],
                        'verdict', 'taken')
    ostream.print_header(header)
    ostream.print_header(len(header) * '-')

    order = sorted(comparison['templates'],
                   key=lambda name: (decision['verdicts'][name] is not None,
                                     decision['scores'][name], name))

    for name in order:
        entry = comparison['templates'][name]
        cells = []
        for region in regions:
            found = entry['regions'].get(region)
            cells.append('' if found is
                         None else ic_cell(found['ic_rmsd'], 'bonds'))

        verdict = decision['verdicts'][name] or 'within the criteria'
        ostream.print_header(
            row.format(name[:22], *cells, verdict[:26],
                       'yes' if name == decision['name'] else ''))

    ostream.print_blank()

    if decision['name'] is None:
        ostream.print_info('No template was taken.')
    elif decision['forced']:
        ostream.print_info(
            f'{decision["name"]} was named rather than chosen, so the '
            'criteria were measured but did not decide.')
    else:
        ranked = ' '.join(ranked_on)
        ostream.print_info(
            f'{len(decision["candidates"])} of '
            f'{len(comparison["templates"])} template(s) are within the '
            f'criteria. Taking {decision["name"]}, whose {ranked} of '
            f'{decision["score"]:.3f} is the lowest of them.')

    ostream.print_blank()
    ostream.flush()


# ----------------------------------------------------------------------
# reporting from the core
# ----------------------------------------------------------------------


def print_binding_modes(binding_modes, ostream=None):
    """
    Prints the detected coordination sphere.
    """

    ostream = stream(ostream)

    ostream.print_blank()
    ostream.print_header('Coordination sphere')
    ostream.print_header(19 * '-')

    for metal in binding_modes['metals']:
        ostream.print_header(
            param(f'metal {metal["element"]} (index {metal["index"]})',
                  f'charge {metal["formal_charge"]:+d}'))

    ostream.print_blank()
    valstr = '{:>10} {:>9} | {:>18} | {:>16}'.format('residue', 'atoms',
                                                     'distances (A)', 'mode')
    ostream.print_header(valstr)
    ostream.print_header(60 * '-')

    by_residue = {}
    for ligand in binding_modes['ligands']:
        by_residue.setdefault(ligand['res_index'], []).append(ligand)

    for group in by_residue.values():
        # a residue binding through several atoms is one ligand, so it gets
        # one row listing them side by side, the same way an atom bridging
        # two metals lists both of its distances. Merging is only
        # unambiguous while every atom contributes exactly one distance;
        # otherwise the distances could not be read back onto their atoms,
        # so that group stays one row per atom
        if any(len(ligand['distances']) != 1 for ligand in group):
            rows = [[ligand] for ligand in group]
        else:
            rows = [group]

        for row in rows:
            atoms = ', '.join(ligand['atom'] for ligand in row)
            distances = ', '.join(f'{d:.2f}' for ligand in row
                                  for d in ligand['distances'])
            modes = '/'.join(dict.fromkeys(ligand['mode'] for ligand in row))
            valstr = '{:>10} {:>9} | {:>18} | {:>16}'.format(
                row[0]['residue'], atoms, distances, modes)
            ostream.print_header(valstr)

    for note in binding_modes['notes']:
        ostream.print_warning(note)

    ostream.print_blank()
    ostream.flush()


def print_binding_mode_update(changes,
                              largest_shift,
                              dropped_residues=None,
                              ostream=None):
    """
    Prints what re-detecting the coordination on a new geometry did.

    :param changes:
        The list of (kind, atom, detail) tuples describing the contacts
        that were gained, lost or reclassified. Empty when the coordination
        is unchanged.
    :param largest_shift:
        The largest change in Angstrom that the recorded metal-ligand
        bonds underwent.
    :param dropped_residues:
        The residues that no longer coordinate at all.
    """

    ostream = stream(ostream)

    ostream.print_blank()
    ostream.print_header('Coordination update')
    ostream.print_header(19 * '-')
    ostream.print_header(param('largest bond change', f'{largest_shift:.2f} A'))

    if not changes:
        ostream.print_header(param('coordination', 'unchanged'))
        ostream.print_blank()
        ostream.print_info(
            'The new geometry gives the same coordination sphere; the '
            'binding modes and the connectivity matrix are kept as they '
            'are.')
        ostream.print_blank()
        ostream.flush()
        return

    ostream.print_header(param('contacts changed', len(changes)))
    ostream.print_blank()

    valstr = '{:>10} {:>12} | {:>46}'.format('change', 'atom', 'detail')
    ostream.print_header(valstr)
    ostream.print_header(72 * '-')

    for kind, atom, detail in changes:
        ostream.print_header('{:>10} {:>12} | {:>46}'.format(
            kind, atom, detail))

    ostream.print_blank()
    ostream.print_info(
        'The binding modes and the connectivity matrix were updated to '
        'the new geometry. Overwrite the ones you hold, or the fit will '
        'use a coordination the geometry no longer has.')

    for residue in dropped_residues or []:
        ostream.print_warning(
            f'{residue} no longer coordinates a metal, but it is still '
            'part of the truncated active site; extract the active site '
            'again to leave it out')

    ostream.print_blank()
    ostream.flush()


def print_active_site(active_site, binding_modes, ostream=None):
    """
    Prints the composition of the truncated active site.
    """

    ostream = stream(ostream)

    molecule = active_site['molecule']

    ostream.print_blank()
    ostream.print_header('Truncated active site')
    ostream.print_header(21 * '-')
    ostream.print_header(param('atoms', molecule.number_of_atoms()))
    ostream.print_header(param('charge', f'{int(molecule.get_charge()):+d}'))
    ostream.print_header(param('multiplicity',
                               int(molecule.get_multiplicity())))
    ostream.print_header(
        param('capping hydrogens', len(active_site['cap_indices'])))
    ostream.print_header(
        param('bonds', int(active_site['connectivity_matrix'].sum() // 2)))
    print_param_list('residues', active_site['residues'], ostream=ostream)

    variants = sorted(binding_modes['variants'].values())
    print_param_list('protonation', variants, ostream=ostream)

    ostream.print_blank()
    ostream.flush()


def print_partial_charges(topology,
                          active_site,
                          partial_charges,
                          corrected_charges,
                          residue_labels,
                          ostream=None):
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
    :param ostream:
        The output stream, or None to print nothing.
    """

    ostream = stream(ostream)

    charges = np.asarray(partial_charges)
    caps = sorted(active_site['cap_indices'])
    metals = active_site['metal_indices']
    labels = active_site['molecule'].get_labels()
    n_atoms = len(charges)
    rest = [index for index in range(n_atoms) if index not in caps]
    cap_charge = float(sum(charges[index] for index in caps))

    ostream.print_blank()
    ostream.print_header('Partial charges')
    ostream.print_header(15 * '-')
    ostream.print_header(
        param('active site charge',
              f'{int(active_site["molecule"].get_charge()):+d}'))
    ostream.print_header(param('fitted total', f'{charges.sum():+.4f} e'))
    ostream.print_header(param('on capping hydrogens', f'{cap_charge:+.4f} e'))
    ostream.print_header(
        param(f'spread over {len(rest)} atoms',
              f'{cap_charge / len(rest):+.4f} e each'))
    ostream.print_blank()

    # group what is left by the residue each atom came from
    atoms = list(topology.atoms())
    by_residue = {}
    for index in rest:
        residue = atoms[active_site['atom_map'][index]].residue
        by_residue.setdefault(residue, []).append(index)

    corrected = np.asarray(corrected_charges)

    valstr = '{:>16} {:>7} | {:>12}'.format('fragment', 'atoms', 'charge')
    ostream.print_header(valstr)
    ostream.print_header(45 * '-')

    for residue, indices in by_residue.items():
        total = sum(corrected[index] for index in indices)
        if len(indices) == 1 and indices[0] in metals:
            name = f'{labels[indices[0]]} (metal)'
        else:
            name = residue_labels[residue.index]
        valstr = '{:>16} {:>7} | {:>12.4f}'.format(name, len(indices), total)
        ostream.print_header(valstr)

    ostream.print_blank()
    ostream.flush()


def print_mm_optimization(active_site,
                          forcefield,
                          relaxed,
                          frozen_indices,
                          metal_keys,
                          equilibrium_labels,
                          bond_change_warning=0.25,
                          ostream=None):
    """
    Prints what the crude MM relaxation did to the coordination sphere.

    The metal-ligand distances before and against after are the point of
    the table: the pass is there to clean up contacts and hydrogens, and
    a coordination sphere that moved more than a few hundredths of an
    Angstrom is the sign that it did something else instead.

    What the seeding did is read off the force field rather than off the
    settings that asked for it, so the table describes the terms the
    relaxation actually ran with.

    :param active_site:
        The active site, holding the geometry the pass started from.
    :param forcefield:
        The seeded force field the relaxation ran on.
    :param relaxed:
        The relaxed molecule.
    :param frozen_indices:
        The indices that were held fixed.
    :param metal_keys:
        The (bonds, angles) keys of the metal terms, from get_metal_keys.
    :param equilibrium_labels:
        The table from a seeded term's comment to the label its equilibrium
        source is printed as (SEEDED_EQUILIBRIUM_LABELS).
    :param bond_change_warning:
        How far a metal-ligand bond may move before it is reported.
    :param ostream:
        The output stream, or None to print nothing.
    """

    def _seeded_constant(table, keys):
        """
        Returns the force constant the seeding put on a set of terms.

        :param table:
            The bond or angle table of the force field.
        :param keys:
            The keys of the metal terms.

        :return:
            The constant as a string, or 'varies' when they are not all the
            same, which the crude pass never makes them.
        """

        constants = {round(table[key]['force_constant'], 6) for key in keys}

        if len(constants) != 1:
            return 'varies'

        return f'{constants.pop():.0f}'

    def _seeded_equilibria(table, keys):
        """
        Returns where the seeding took its equilibrium values from.

        _seed_metal_terms writes that into the comment of every term it
        touches, so the force field says it without being asked again.

        :param table:
            The bond or angle table of the force field.
        :param keys:
            The keys of the metal terms.

        :return:
            'requested', 'given', 'measured', or 'mixed'.
        """

        sources = {table[key].get('comment') for key in keys}

        if len(sources) != 1:
            return 'mixed'

        return equilibrium_labels.get(sources.pop(), 'mixed')

    ostream = stream(ostream)

    labels = active_site['molecule'].get_labels()
    molecule = active_site['molecule']
    metals = sorted(active_site['metal_indices'])
    bonds, angles = metal_keys

    before = molecule.get_coordinates_in_angstrom()
    after = relaxed.get_coordinates_in_angstrom()
    shift = np.linalg.norm(after - before, axis=1)

    ostream.print_blank()
    ostream.print_header('Crude MM relaxation')
    ostream.print_header(19 * '-')

    ostream.print_header(param('frozen atoms', len(frozen_indices)))
    ostream.print_header(
        param('metal centers',
              'frozen' if set(metals) <= set(frozen_indices) else 'free'))
    ostream.print_header(
        param('metal bonds',
              f'{len(bonds)}, k = {_seeded_constant(forcefield.bonds, bonds)}'))
    ostream.print_header(
        param(
            'metal angles', f'{len(angles)}, k = '
            f'{_seeded_constant(forcefield.angles, angles)}'
            if angles else 'left untouched'))
    ostream.print_header(
        param('bond equilibria', _seeded_equilibria(forcefield.bonds, bonds)))
    if angles:
        ostream.print_header(
            param('angle equilibria',
                  _seeded_equilibria(forcefield.angles, angles)))
    ostream.print_blank()

    valstr = '{:>12} {:>8} | {:>10} | {:>9} | {:>8}'.format(
        'atoms', 'elements', 'before (A)', 'after (A)', 'change')
    ostream.print_header(valstr)
    ostream.print_header(60 * '-')

    worst_bond = None
    for key in bonds:
        one_based = [index + 1 for index in key]
        was = molecule.get_distance_in_angstroms(one_based)
        now = relaxed.get_distance_in_angstroms(one_based)
        names = '-'.join(labels[index] for index in key)
        valstr = '{:>12} {:>8} | {:>10.2f} | {:>9.2f} | {:>+8.2f}'.format(
            str(key), names, was, now, now - was)
        ostream.print_header(valstr)
        if worst_bond is None or abs(now - was) > abs(worst_bond[1]):
            worst_bond = (names, now - was)

    for first in range(len(metals)):
        for second in range(first + 1, len(metals)):
            pair = (metals[first], metals[second])
            one_based = [index + 1 for index in pair]
            was = molecule.get_distance_in_angstroms(one_based)
            now = relaxed.get_distance_in_angstroms(one_based)
            names = '-'.join(labels[index] for index in pair)
            valstr = ('{:>12} {:>8} | {:>10.2f} | {:>9.2f} | '
                      '{:>+8.2f}').format(str(pair), names, was, now, now - was)
            ostream.print_header(valstr)

    ostream.print_blank()

    largest = int(np.argmax(shift))
    ostream.print_header(
        param('largest shift', f'{shift[largest]:.2f} A on '
              f'{labels[largest]} {largest}'))
    ostream.print_header(param('mean shift', f'{shift.mean():.2f} A'))

    # A ligand swinging around its metal moves a long way in Cartesian
    # terms while the coordination sphere itself is untouched, so what
    # the pass has to be held to is the bond lengths, not the shifts.
    if worst_bond is not None:
        ostream.print_header(
            param('largest bond change',
                  f'{worst_bond[1]:+.2f} A on {worst_bond[0]}'))
        if abs(worst_bond[1]) > bond_change_warning:
            ostream.print_warning(
                f'The crude relaxation changed a {worst_bond[0]} bond by '
                f'{worst_bond[1]:+.2f} A. Check the metal terms it was '
                'given before trusting the geometry it produced.')

    ostream.print_blank()
    ostream.flush()


def print_metal_parameters(active_site, forcefield, metal_keys, ostream=None):
    """
    Prints the fitted metal bonds and angles.

    :param active_site:
        The active site.
    :param forcefield:
        The fitted force field.
    :param metal_keys:
        The (bonds, angles) keys of the metal terms, from get_metal_keys.
    :param ostream:
        The output stream, or None to print nothing.
    """

    ostream = stream(ostream)

    labels = active_site['molecule'].get_labels()
    metals = set(active_site['metal_indices'])
    coords = active_site['molecule'].get_coordinates_in_angstrom()
    bonds, angles = metal_keys

    ostream.print_blank()
    ostream.print_header('Metal bonds')
    ostream.print_header(11 * '-')
    valstr = '{:>12} {:>7} | {:>9} | {:>21}'.format('atoms', 'elements',
                                                    'r0 (A)',
                                                    'k (kcal/mol/A^2)')
    ostream.print_header(valstr)
    ostream.print_header(60 * '-')

    for key in bonds:
        params = forcefield.bonds[key]
        names = '-'.join(labels[index] for index in key)
        # kJ/mol/nm^2 to kcal/mol/A^2
        force_constant = params['force_constant'] / 100.0 / 4.184
        valstr = '{:>12} {:>7} | {:>9.3f} | {:>21.1f}'.format(
            str(key), names, params['equilibrium'] * 10.0, force_constant)
        ostream.print_header(valstr)

    ostream.print_blank()
    ostream.print_header('Metal angles')
    ostream.print_header(12 * '-')
    valstr = '{:>14} {:>9} | {:>12} | {:>19}'.format('atoms', 'elements',
                                                     'theta0 (deg)',
                                                     'k (kJ/mol/rad^2)')
    ostream.print_header(valstr)
    ostream.print_header(60 * '-')

    for key in angles:
        params = forcefield.angles[key]
        names = '-'.join(labels[index] for index in key)
        bridging = key[0] in metals and key[2] in metals
        # the marker gets a fixed-width field of its own, otherwise the
        # centering of print_header would shift the marked line
        valstr = '{:>14} {:>9} | {:>12.1f} | {:>19.1f} {:<9}'.format(
            str(key), names, params['equilibrium'], params['force_constant'],
            'bridging' if bridging else '')
        ostream.print_header(valstr)

    ostream.print_blank()

    for metal_a in sorted(metals):
        for metal_b in sorted(metals):
            if metal_a >= metal_b:
                continue
            distance = np.linalg.norm(coords[metal_a] - coords[metal_b])
            ostream.print_header(
                param(f'{labels[metal_a]}-{labels[metal_b]} distance',
                      f'{distance:.3f} A'))

    ostream.print_blank()
    ostream.flush()


def print_muted_notice(step, mute_scf=True, ostream=None):
    """
    Announces a long calculation whose output is being suppressed.

    :param step:
        A description of the step about to run.
    """

    ostream = stream(ostream)

    if mute_scf:
        ostream.print_info(
            f'Running {step} with muted QM output. Set mute_scf to False '
            'to follow it.')
    else:
        ostream.print_info(f'Running {step}.')
    ostream.flush()
