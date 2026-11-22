"""Reads the json files run.py writes.

    python benchmarks/compare.py results/RUN.json
    python benchmarks/compare.py results/BASE.json results/NEW.json

With one file: each case's RI row against its four-center row, from the same
run on the same machine. The time ratio is the speedup of RI; the differences
of the results are the error of the fitting and are reported, not judged.

With two files: each row of NEW against the row of BASE with the same case
and variant. Results are judged against the tolerances below; wall time and
peak memory are reported as ratios and marked where they grew by more than
--time-tol and --mem-tol. Timings are only comparable between runs with the
same machine, threads, ranks and launcher; the differences are printed first.

Exits with 1 when a result is outside its tolerance or a row of the (last)
run has a problem (a failed process, missing results, an RI row that fell
back to the four-center integrals).
"""

import argparse
import json
import sys

import numpy as np

# NOTE: (absolute, relative) tolerances for a result of the same variant in
# two runs: a difference is accepted when it is at most atol + rtol * |base|.
# Energies and gradients are in Hartree and Hartree/Bohr; excitation and
# photon energies in Hartree; the strengths in the units the output prints.
TOLERANCES = {
    'scf_energy': (1.0e-6, 0.0),
    'opt_final_energy': (1.0e-6, 0.0),
    'gradient': (1.0e-5, 0.0),
    'excitation_energies': (1.0e-5, 0.0),
    'tpa_photon_energies': (1.0e-5, 0.0),
    'oscillator_strengths': (1.0e-4, 1.0e-3),
    'tpa_strengths_linear': (1.0e-2, 1.0e-3),
    'tpa_strengths_circular': (1.0e-2, 1.0e-3),
    'tpa_cross_sections': (0.0, 1.0e-3),
    '3pa_strengths_linear': (1.0e-2, 1.0e-3),
    '3pa_strengths_circular': (1.0e-2, 1.0e-3),
}

# NOTE: counts that may change between runs without being wrong; a change is
# printed and not judged.
NOTED = ('opt_steps',)


def load(path):
    with open(path) as fh:
        document = json.load(fh)
    rows = {(r['suite'], r['case'], r['variant']): r for r in document['rows']}
    return document, rows


def difference(a, b):
    """
    The largest absolute difference, and the same divided by the largest
    |b|. The divisor is the largest value of the array rather than each
    element's own, because a strength can be zero by symmetry and print as
    1e-30 in one run and 1e-31 in the other.
    """

    a = np.atleast_1d(np.asarray(a, dtype=float))
    b = np.atleast_1d(np.asarray(b, dtype=float))
    if a.shape != b.shape:
        return None, None
    diff = float(np.abs(a - b).max(initial=0.0))
    scale = float(np.abs(b).max(initial=0.0))
    return diff, (diff / scale if scale > 0 else 0.0)


def within(key, a, b):
    atol, rtol = TOLERANCES[key]
    a = np.atleast_1d(np.asarray(a, dtype=float))
    b = np.atleast_1d(np.asarray(b, dtype=float))
    if a.shape != b.shape:
        return False
    return bool(np.all(np.abs(a - b) <= atol + rtol * np.abs(b)))


def gib(n):
    return n / 1024**3


def scf_iterations(row):
    stage = row.get('stages', {}).get('SCF')
    return stage['iterations'] if stage else None


def print_provenance_diff(base, new):
    """Prints what differs between two runs and changes what a time means."""

    def view(doc):
        p = doc['provenance']
        return {
            'host': p['machine'].get('host'),
            'cpu': p['machine'].get('cpu'),
            'launcher': p.get('launcher'),
            'OMP_NUM_THREADS': p.get('env', {}).get('OMP_NUM_THREADS'),
            'veloxchem': p.get('veloxchem'),
            'commit': p.get('commit'),
            'dirty': p.get('dirty'),
        }

    a, b = view(base), view(new)
    for key in a:
        mark = '' if a[key] == b[key] else '   <-- differs'
        print(f'  {key:16s} {str(a[key]):30s} {str(b[key]):30s}{mark}')
    for key in ('host', 'cpu', 'launcher', 'OMP_NUM_THREADS'):
        if a[key] != b[key]:
            print('  The machine or its configuration differs: the time and '
                  'memory ratios do not measure the code alone.')
            break
    print()


def summarize(document, rows):
    """One run: RI against four-center, case by case."""

    print(f'{document["run_id"]}  ({document["provenance"].get("commit")})')
    print()
    header = (f'{"case":36s} {"4c wall":>9s} {"RI wall":>9s} {"speedup":>8s} '
              f'{"4c GiB":>7s} {"RI GiB":>7s}  RI - 4c')
    print(header)
    print('-' * len(header))

    failed = False
    cases = sorted({(s, c) for s, c, _ in rows})
    for suite, case in cases:
        ri = rows.get((suite, case, 'ri'))
        fc = rows.get((suite, case, 'four_center'))

        def wall(r):
            return f'{r["wall"]:9.1f}' if r else f'{"-":>9s}'

        def mem(r):
            return f'{gib(r["peak_rss"][0]):7.2f}' if r else f'{"-":>7s}'

        speedup = (f'{fc["wall"] / ri["wall"]:8.1f}'
                   if ri and fc and ri['wall'] > 0 else f'{"-":>8s}')

        diffs = []
        if ri and fc:
            for key in TOLERANCES:
                if key in ri['results'] and key in fc['results']:
                    absdiff, rel = difference(ri['results'][key],
                                              fc['results'][key])
                    if absdiff is None:
                        diffs.append(f'{key} shape differs')
                    elif TOLERANCES[key][1] > 0:
                        diffs.append(f'{key} {rel:.1e} rel')
                    else:
                        diffs.append(f'{key} {absdiff:.1e}')

        print(f'{suite + "/" + case:36s} {wall(fc)} {wall(ri)} {speedup} '
              f'{mem(fc)} {mem(ri)}  {"; ".join(diffs)}')

        for r in (fc, ri):
            if r and r['problems']:
                failed = True
                print(f'    {r["variant"]}: {"; ".join(r["problems"])}')

    return failed


def regress(base_doc, base, new_doc, new, time_tol, mem_tol):
    """Two runs: each row of new against the same row of base."""

    print(f'base {base_doc["run_id"]}')
    print(f'new  {new_doc["run_id"]}')
    print()
    print_provenance_diff(base_doc, new_doc)

    header = (f'{"case":36s} {"variant":12s} {"wall":>7s} {"memory":>7s} '
              f'{"SCF it":>9s}  results')
    print(header)
    print('-' * len(header))

    failed = False
    for key in sorted(new):
        row = new[key]
        ref = base.get(key)
        name = f'{key[0]}/{key[1]}'
        if ref is None:
            print(f'{name:36s} {key[2]:12s}  not in base')
            continue

        wall = row['wall'] / ref['wall'] if ref['wall'] > 0 else float('nan')
        memory = row['peak_rss'][0] / ref['peak_rss'][0]
        wall_mark = '!' if wall > 1.0 + time_tol else ' '
        mem_mark = '!' if memory > 1.0 + mem_tol else ' '

        it_new, it_ref = scf_iterations(row), scf_iterations(ref)
        iterations = (f'{it_ref}->{it_new}'
                      if it_new != it_ref else f'{it_new}')

        notes = []
        for result in TOLERANCES:
            in_new = result in row['results']
            in_ref = result in ref['results']
            if in_new != in_ref:
                failed = True
                notes.append(f'{result} {"missing" if in_ref else "new"}')
            elif in_new and not within(result, row['results'][result],
                                       ref['results'][result]):
                failed = True
                absdiff, rel = difference(row['results'][result],
                                          ref['results'][result])
                if absdiff is None:
                    notes.append(f'{result} FAIL (shape differs)')
                else:
                    notes.append(f'{result} FAIL ({absdiff:.1e}, '
                                 f'{rel:.1e} rel)')
        for result in NOTED:
            a = ref['results'].get(result)
            b = row['results'].get(result)
            if a != b:
                notes.append(f'{result} {a}->{b}')

        if row['problems']:
            failed = True
            notes.extend(row['problems'])

        print(f'{name:36s} {key[2]:12s} {wall:6.2f}{wall_mark} '
              f'{memory:6.2f}{mem_mark} {iterations:>9s}  '
              f'{"; ".join(notes) or "ok"}')

    missing = sorted(set(base) - set(new))
    if missing:
        print(f'{len(missing)} rows of base are not in new: ' +
              ', '.join(f'{s}/{c} ({v})' for s, c, v in missing))

    print()
    print(f'wall and memory are new / base; ! marks growth above '
          f'{time_tol:.0%} (wall) and {mem_tol:.0%} (memory).')
    return failed


def main():
    parser = argparse.ArgumentParser(
        description='Summarize one benchmark run or compare two.')
    parser.add_argument('runs', nargs='+', help='one or two run json files')
    parser.add_argument('--time-tol',
                        type=float,
                        default=0.10,
                        help='relative growth of wall time that is marked '
                        '(default: 0.10)')
    parser.add_argument('--mem-tol',
                        type=float,
                        default=0.10,
                        help='relative growth of peak memory that is marked '
                        '(default: 0.10)')
    args = parser.parse_args()

    if len(args.runs) == 1:
        failed = summarize(*load(args.runs[0]))
    elif len(args.runs) == 2:
        failed = regress(*load(args.runs[0]), *load(args.runs[1]),
                         args.time_tol, args.mem_tol)
    else:
        parser.error('one or two run files')

    sys.exit(1 if failed else 0)


if __name__ == '__main__':
    main()
