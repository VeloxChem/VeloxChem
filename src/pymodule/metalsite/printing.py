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
    ostream.print_header(
        param('atoms', template['molecule'].number_of_atoms()))
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
    valstr = '{:>24} | {:>7} | {:>7} | {:>13}'.format('name', 'atoms',
                                                      'metals', 'geometry')
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


def print_spec(title, described, residue_nodes, ostream=None):
    """
    Prints what a site is made of: which residues it holds and which of them
    coordinate which metal.

    The residues are named by their formula, which is for reading; two sites
    are compared on the keys behind them.

    :param title:
        What the block is describing.
    :param described:
        A described active site or a template.
    :param residue_nodes:
        Its residue nodes, metals left out.
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

    bridging = [node for node in residue_nodes if coarse.degree(node) > 1]
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


def print_selection(comparison, decision, rmsd_regions, ic_types, ranked_on,
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
            name: ' / '.join(f'{measure} {limit:.2f}'
                             for measure, limit in given.items()
                             if limit is not None)
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
            cells.append('' if found is None else ic_cell(
                found['ic_rmsd'], 'bonds'))

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
