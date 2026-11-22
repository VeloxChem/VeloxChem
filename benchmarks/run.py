"""Runs benchmark cases through the VeloxChem command line and records them.

    python benchmarks/run.py quick                     # one suite, ri: yes and no
    python benchmarks/run.py quick standard --ri yes   # RI only
    python benchmarks/run.py cases/quick/scf_hf.inp    # one case
    python benchmarks/run.py quick --launcher 'mpirun -n 2'
    python benchmarks/run.py --list

Each case is an input file under cases/<suite>/. It is run once with
`ri: yes` and once with `ri: no` added to its @method settings group, as its
own process, in a fresh directory, so that no checkpoint file of an earlier
run is read as a restart and a killed case costs one row and not the run.

One json file is written per run, to results/<run id>.json, and rewritten
after every case, so a run stopped early leaves a valid file. Read it with
compare.py.
"""

import argparse
import json
import os
import platform
import re
import shlex
import shutil
import socket
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
CASES = HERE / 'cases'
GEOMETRIES = HERE / 'geometries'
RESULTS = HERE / 'results'

SCHEMA = 1

# NOTE: the h5 datasets read from the final results file (<stem>.h5). Only
# quantities that do not depend on the phase of a response vector are listed;
# transition dipoles and transition moments change sign with the phase.
H5_SCALARS = {
    'scf_energy': 'scf/scf_energy',
    'c6': 'rsp/c6',
}
H5_FLAGS = {
    'ri': 'scf/ri',
    'ri_jk': 'scf/ri_jk',
    'ri_coulomb': 'scf/ri_coulomb',
    'ri_fallback': 'scf/ri_fallback',
    'force_ri': 'scf/force_ri',
}
H5_ARRAYS = {
    'excitation_energies': 'rsp/eigenvalues',
    'oscillator_strengths': 'rsp/oscillator_strengths',
    'tpa_photon_energies': 'rsp/photon_energies',
    'tpa_strengths_linear': 'rsp/tpa_strengths/linear',
    'tpa_strengths_circular': 'rsp/tpa_strengths/circular',
    'tpa_cross_sections': 'rsp/cross_sections',
    'opt_energies': 'opt/opt_energies',
    'vib_frequencies': 'vib/vib_frequencies',
    'ir_intensities': 'vib/ir_intensities',
}
# NOTE: a dictionary of linear response functions <<A;B>>(w), keyed by the
# two operator components and the frequency; read into a dictionary of real
# and imaginary parts.
H5_RESPONSE_FUNCTIONS = 'rsp/response_functions'

RE_BUILD = re.compile(r'\(build commit (\w+)\)')
RE_CONVERGED = re.compile(
    r'\*\*\* (.+?) converged in (\d+) iterations\. Time: ([\d.]+) sec')
RE_TIME_SPENT = re.compile(r'\*\*\* Time spent in (.+?): ([\d.]+) sec')
RE_TOTAL = re.compile(r'Total execution time is ([\d.]+) sec')
RE_WARNING = re.compile(r'^\s*(?:\* Warning \*|\*\*\* Warning:)\s*(.*?)\s*$')
RE_POLGRAD_ATOM = re.compile(r'\*\* Atom #(\d+): (\w+) \*\*')

# NOTE: the result tables read from the output, as (title, result key, row
# pattern, values of a row). A table starts at its title and ends at the
# next title, a "Time spent" line or the end of the output. A title seen
# again starts the table again, so the last one is kept (an optimization
# prints a gradient at every step). The polarizability gradient is printed
# per atom and only for components above 1e-6 in magnitude, so it is read
# into a dictionary keyed by atom, component and frequency.
OUTPUT_TABLES = (
    ('Analytical Gradient (Hartree/Bohr)', 'gradient',
     re.compile(r'^\s*[A-Z][a-z]?\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)'
                r'\s+(-?\d+\.\d+)\s*$'),
     lambda m: [float(g) for g in m.groups()]),
    ('3PA Strength (Linear Polarization)', '3pa_strengths_linear',
     re.compile(r'^\s*\d+\s+[\d.]+ eV\s+(-?[\d.]+) a\.u\.\s*$'),
     lambda m: [float(m.group(1))]),
    ('3PA Strength (Circular Polarization)', '3pa_strengths_circular',
     re.compile(r'^\s*\d+\s+[\d.]+ eV\s+(-?[\d.]+) a\.u\.\s*$'),
     lambda m: [float(m.group(1))]),
    ('SHG hyperpolarizability', 'shg_beta',
     re.compile(r'^\s*beta\s+[\d.]+\s+(-?[\d.]+)\s+(-?[\d.]+)j\s*$'),
     lambda m: [float(m.group(1)), float(m.group(2))]),
    ('THG Summary', 'thg_gamma',
     re.compile(r'^\s*[\d.]+ eV\s+(-?[\d.]+)\s+(-?[\d.]+)\s+[\d.]+\s*$'),
     lambda m: [float(m.group(1)), float(m.group(2))]),
    ('Analytical Polarizability Gradients', 'polarizability_gradient',
     re.compile(r'^\s*(d<<\w;\w>>/d\w)\s+([\d.]+)\s+(-?[\d.]+)\s*$'), None),
)
RE_RI_KEY = re.compile(r'^\s*ri\s*:', re.IGNORECASE | re.MULTILINE)
RE_XYZFILE = re.compile(r'^(\s*xyzfile\s*:\s*)(\S+)\s*$',
                        re.IGNORECASE | re.MULTILINE)
RE_VALUE = re.compile(r'^\s*[\w ]+?\s*:\s*(\S+)\s*$', re.MULTILINE)
RE_KILLED = re.compile(r'signal 9|\bkilled\b|out of memory|oom', re.IGNORECASE)


def find_cases(names):
    """
    Resolves suite names, directories and input files to a list of cases.

    :param names:
        Suite names (directories under cases/), paths to directories of
        .inp files, or paths to .inp files.

    :return:
        A list of (suite, path) pairs, in the order given and sorted within
        a suite. The suite of a case is the name of its directory.
    """

    found = []
    for name in names:
        path = Path(name)
        if path.suffix == '.inp' and path.is_file():
            found.append((path.resolve().parent.name, path.resolve()))
            continue
        suite = path if path.is_dir() else CASES / name
        if not suite.is_dir():
            sys.exit(f'run.py: no suite, directory or input file named '
                     f'{name}; suites: {", ".join(list_suites())}')
        suite = suite.resolve()
        found.extend((suite.name, p) for p in sorted(suite.glob('*.inp')))
    return found


def list_suites():
    return sorted(p.name for p in CASES.iterdir() if p.is_dir())


def geometry_folders(case_path):
    """
    The folders a geometry named by a case is looked up in, in order: the
    folder of the case, a geometries/ folder beside it or beside one of its
    two parent folders, and benchmarks/geometries.
    """

    folders = [case_path.parent]
    for parent in list(case_path.parents)[:3]:
        folders.append(parent / 'geometries')
    folders.append(GEOMETRIES)
    return folders


def describe(path):
    """The first comment line of a case file."""

    for line in path.read_text().splitlines():
        if line.startswith('#'):
            return line.lstrip('#').strip()
    return ''


def prepare_input(case_path, ri, workdir):
    """
    Writes the input of one variant of a case into its working directory.

    Adds `ri: yes` or `ri: no` to the @method settings group and copies each
    geometry named by `xyzfile:` beside the input, from the first of
    geometry_folders that has it. Any other value that names a file in the
    folder of the case (point charges, parameter files) is copied as well.

    :return:
        The path of the written input.
    """

    text = case_path.read_text()

    if RE_RI_KEY.search(text):
        sys.exit(f'run.py: {case_path} sets ri itself; the runner adds it')

    def copy_geometry(match):
        name = match.group(2)
        for folder in geometry_folders(case_path):
            source = folder / name
            if source.is_file():
                shutil.copy(source, workdir / Path(name).name)
                return match.group(1) + Path(name).name
        sys.exit(f'run.py: geometry {name} of {case_path} not found')

    text = RE_XYZFILE.sub(copy_geometry, text)

    for match in RE_VALUE.finditer(text):
        name = match.group(1)
        source = case_path.parent / name
        if '/' not in name and source.is_file() and source != case_path:
            shutil.copy(source, workdir / name)

    ri_line = f'ri: {"yes" if ri else "no"}'
    lines = text.splitlines()
    for i, line in enumerate(lines):
        # NOTE: the input parser names a group by its words joined with '_',
        # so '@method settings' and '@method_settings' are one group.
        if '_'.join(line.lower().split()) == '@method_settings':
            lines.insert(i + 1, ri_line)
            break
    else:
        lines = ['@method settings', ri_line, '@end', ''] + lines

    path = workdir / f'{case_path.stem}.inp'
    path.write_text('\n'.join(lines) + '\n')
    return path


def run_process(command, cwd, env, log_path):
    """
    Runs one process and returns its exit code, wall time and peak memory.

    The peak resident set size comes from wait4: the largest process among
    the child and the descendants it waited for, on this host. With mpirun
    that is the largest rank; ranks on other hosts are not seen.

    :return:
        (exit code, wall time in seconds, peak RSS in bytes)
    """

    with open(log_path, 'w') as log:
        t0 = time.time()
        proc = subprocess.Popen(command,
                                cwd=cwd,
                                env=env,
                                stdout=log,
                                stderr=subprocess.STDOUT)
        _, status, usage = os.wait4(proc.pid, 0)
        wall = time.time() - t0
    # NOTE: the process was reaped by wait4 above, so Popen is told its exit
    # code rather than left to wait for it again.
    proc.returncode = os.waitstatus_to_exitcode(status)
    scale = 1 if platform.system() == 'Darwin' else 1024
    return proc.returncode, wall, usage.ru_maxrss * scale


def read_h5(path):
    """The results of one calculation that the final h5 file holds."""

    record = {'results': {}, 'flags': {}}
    if not path.is_file():
        return record

    with h5py.File(path, 'r') as h5:
        for key, name in H5_SCALARS.items():
            if name in h5:
                record['results'][key] = float(h5[name][()])
        for key, name in H5_FLAGS.items():
            if name in h5:
                record['flags'][key] = bool(h5[name][()])
        for key, name in H5_ARRAYS.items():
            if name in h5:
                record['results'][key] = np.ravel(h5[name][()]).tolist()
        if H5_RESPONSE_FUNCTIONS in h5:
            record['results']['response_functions'] = read_response_functions(
                h5[H5_RESPONSE_FUNCTIONS])

    energies = record['results'].pop('opt_energies', None)
    if energies:
        record['results']['opt_final_energy'] = energies[-1]
        record['results']['opt_steps'] = len(energies)

    return record


def read_response_functions(group):
    """
    The response functions of an h5 group as {'x,y,0.3000 re': value, ...}.

    Each entry holds a key, (operator component, operator component,
    frequency), and a complex value. The frequency is printed to four
    decimals, so that 0.36 and 0.36000000000000004 give one key.
    """

    values = {}
    for entry in group.values():
        parts = [entry['__key__'][k][()] for k in sorted(entry['__key__'])]
        parts = [
            p.decode() if isinstance(p, bytes) else f'{float(p):.4f}'
            for p in parts
        ]
        value = complex(entry['__value__'][()])
        label = ','.join(parts)
        values[f'{label} re'] = value.real
        values[f'{label} im'] = value.imag
    return dict(sorted(values.items()))


def read_output(path):
    """
    What the output text holds that the h5 file does not.

    The stages are the solver summaries ("... converged in N iterations.
    Time: T sec") and the "Time spent in ..." lines, summed by label. The
    results are the tables of OUTPUT_TABLES.
    """

    record = {'stages': {}, 'warnings': {}, 'results': {}}
    if not path.is_file():
        return record

    lines = path.read_text(errors='replace').splitlines()

    def add_stage(label, seconds, iterations=None):
        stage = record['stages'].setdefault(label, {
            'count': 0,
            'seconds': 0.0,
            'iterations': 0
        })
        stage['count'] += 1
        stage['seconds'] = round(stage['seconds'] + seconds, 3)
        if iterations is not None:
            stage['iterations'] += iterations

    tables = {}
    table = None
    atom = None
    for line in lines:

        match = RE_BUILD.search(line)
        if match and 'build_commit' not in record:
            record['build_commit'] = match.group(1)

        match = RE_CONVERGED.search(line)
        if match:
            add_stage(match.group(1), float(match.group(3)),
                      int(match.group(2)))

        match = RE_TIME_SPENT.search(line)
        if match:
            add_stage(f'Time spent in {match.group(1)}', float(match.group(2)))
            table = None

        match = RE_TOTAL.search(line)
        if match:
            record['total_execution_time'] = float(match.group(1))

        match = RE_WARNING.match(line)
        if match and match.group(1):
            text = match.group(1)
            record['warnings'][text] = record['warnings'].get(text, 0) + 1

        title = next((t for t in OUTPUT_TABLES if t[0] in line), None)
        if title is not None:
            table = title
            tables[table[1]] = {} if table[3] is None else []
            continue
        if table is None:
            continue

        _, key, pattern, values = table
        if values is None:
            match = RE_POLGRAD_ATOM.search(line)
            if match:
                atom = f'{match.group(1)} {match.group(2)}'
                continue
        match = pattern.match(line)
        if match and values is None:
            tables[key][f'{atom} {match.group(1)} {match.group(2)}'] = float(
                match.group(3))
        elif match:
            tables[key].extend(values(match))

    record['results'].update({k: v for k, v in tables.items() if v})

    return record


def clean(workdir, stem):
    """Removes checkpoint and scratch files, keeping input and output."""

    keep = {f'{stem}.inp', f'{stem}.out', f'{stem}.log'}
    for path in workdir.iterdir():
        if path.name not in keep and not path.name.endswith('.xyz'):
            if path.is_dir():
                shutil.rmtree(path)
            else:
                path.unlink()


def memory_limit_prefix(limit):
    """
    The command that runs a process under a memory limit, without swap: a
    transient systemd scope of the user (Linux with systemd only). The
    program reads the limit through the control group, as it reads one set
    by a batch system.
    """

    if shutil.which('systemd-run') is None:
        sys.exit('run.py: --memory-limit needs systemd-run (Linux)')
    return [
        'systemd-run', '--user', '--scope', '--quiet', '-p',
        f'MemoryMax={limit}', '-p', 'MemorySwapMax=0'
    ]


def was_killed(rc, log_path):
    """
    Whether a process was killed rather than ended: exit code 137 or signal
    9, or a launcher's message about a killed task or the memory running out
    in the log of standard output and error.
    """

    if rc in (137, -9):
        return True
    try:
        return bool(RE_KILLED.search(log_path.read_text(errors='replace')))
    except OSError:
        return False


def run_case(suite, case_path, ri, rundir, launcher, repeat, keep_files,
             memory_limit):
    """Runs one variant of one case, repeat times, and returns its record."""

    variant = 'ri' if ri else 'four_center'
    record = {
        'suite': suite,
        'case': case_path.stem,
        'description': describe(case_path),
        'variant': variant,
        'memory_limit': memory_limit,
        'walls': [],
        'peak_rss': [],
        'exit_codes': [],
        'killed': False,
    }
    prefix = memory_limit_prefix(memory_limit) if memory_limit else []

    for rep in range(repeat):
        workdir = rundir / suite / case_path.stem / variant
        if repeat > 1:
            workdir = workdir / f'rep{rep + 1}'
        workdir.mkdir(parents=True)

        inp = prepare_input(case_path, ri, workdir)
        out = workdir / f'{case_path.stem}.out'
        command = prefix + launcher + [
            sys.executable, '-m', 'veloxchem', inp.name, out.name
        ]

        log = workdir / f'{case_path.stem}.log'
        rc, wall, rss = run_process(command, workdir, os.environ.copy(), log)
        record['walls'].append(round(wall, 3))
        record['peak_rss'].append(rss)
        record['exit_codes'].append(rc)
        if rc != 0 and was_killed(rc, log):
            record['killed'] = True

        # NOTE: the results, stages and warnings are read from the first
        # repetition; later ones only add a wall time and a peak memory.
        if rep == 0:
            h5 = read_h5(workdir / f'{case_path.stem}.h5')
            text = read_output(out)
            record['results'] = {**h5['results'], **text['results']}
            record['flags'] = h5['flags']
            record['stages'] = text['stages']
            record['warnings'] = text['warnings']
            record['build_commit'] = text.get('build_commit')
            record['total_execution_time'] = text.get('total_execution_time')
            record['workdir'] = str(workdir.relative_to(rundir))

        if not keep_files:
            clean(workdir, case_path.stem)

        if rc != 0:
            break

    record['wall'] = min(record['walls'])
    record['problems'], record['notes'] = problems(record, ri, memory_limit)
    return record


def problems(record, ri, memory_limit):
    """
    What makes a record unusable as a measurement of its variant, and what
    is worth noting about it.

    An ri row that fell back to the four-center integrals measured the
    four-center path; the SCF reports that in the h5 file, and the gradient
    and response drivers in a warning that names the four-center integrals.
    Under a memory limit a fallback is what the limit is there to bring
    about: it is a note, and a killed process is the problem.

    :return:
        (problems, notes), two lists of strings.
    """

    found = []
    if record['killed']:
        found.append(f'killed (exit code {record["exit_codes"][-1]})')
    elif any(rc != 0 for rc in record['exit_codes']):
        found.append(f'exit code {record["exit_codes"]}')
    if not record.get('results'):
        found.append('no results read')

    fallbacks = []
    flags = record.get('flags', {})
    if ri:
        if flags and not (flags.get('ri_jk') or flags.get('ri_coulomb')):
            found.append('ri asked for, the SCF used neither RI-J nor RI-JK')
        if flags.get('ri_fallback'):
            fallbacks.append('the SCF fell back to the four-center integrals')
        if any('four-center' in w for w in record.get('warnings', {})):
            fallbacks.append('a step fell back to the four-center integrals')
    elif flags.get('ri'):
        found.append('ri: no asked for, the SCF has ri set')

    if memory_limit:
        return found, fallbacks
    return found + fallbacks, []


def git(*args):
    try:
        return subprocess.run(['git', *args],
                              cwd=HERE,
                              capture_output=True,
                              text=True,
                              check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def cpu_name():
    try:
        if platform.system() == 'Darwin':
            return subprocess.run(['sysctl', '-n', 'machdep.cpu.brand_string'],
                                  capture_output=True,
                                  text=True,
                                  check=True).stdout.strip()
        for line in Path('/proc/cpuinfo').read_text().splitlines():
            if line.startswith('model name'):
                return line.split(':', 1)[1].strip()
    except (OSError, subprocess.CalledProcessError):
        pass
    return platform.processor() or platform.machine()


def runtime_info():
    """
    The veloxchem the runner's Python imports, and the BLAS libraries and
    thread pools of a process that has imported it (threadpoolctl, when
    installed).
    """

    probe = ('import json, veloxchem\n'
             'info = {"veloxchem": veloxchem.__version__,\n'
             '        "veloxchem_path": veloxchem.__file__}\n'
             'try:\n'
             '    from threadpoolctl import threadpool_info\n'
             '    info["blas"] = [{"lib": p.get("filepath"),\n'
             '                     "api": p.get("internal_api"),\n'
             '                     "threads": p.get("num_threads")}\n'
             '                    for p in threadpool_info()]\n'
             'except ImportError:\n'
             '    info["blas"] = None\n'
             'print(json.dumps(info))\n')
    try:
        out = subprocess.run([sys.executable, '-c', probe],
                             capture_output=True,
                             text=True,
                             check=True)
        return json.loads(out.stdout.strip().splitlines()[-1])
    except (OSError, subprocess.CalledProcessError, ValueError, IndexError):
        return {}


def provenance(machine, launcher):
    env_keys = ('OMP_', 'OPENBLAS_', 'MKL_', 'BLIS_', 'VLX_', 'SLURM_NTASKS',
                'SLURM_CPUS_PER_TASK')
    now = datetime.now(timezone.utc).astimezone()
    return {
        'date': now.isoformat(timespec='seconds'),
        'commit': git('rev-parse', '--short', 'HEAD'),
        'dirty': bool(git('status', '--porcelain', '--untracked-files=no')),
        'machine': {
            'name': machine,
            'host': socket.gethostname(),
            'cpu': cpu_name(),
            'cores': os.cpu_count(),
            'os': f'{platform.system()} {platform.release()}',
        },
        'launcher': launcher,
        'python': sys.executable,
        'env': {
            k: v
            for k, v in sorted(os.environ.items())
            if k.startswith(env_keys)
        },
        **runtime_info(),
    }


def check_provenance(info):
    """Prints what makes the timings of a run hard to read."""

    if 'OMP_NUM_THREADS' not in info['env']:
        print('run.py: OMP_NUM_THREADS is not set; the thread count is '
              "whatever the OpenMP runtime picks")
    threads = info['env'].get('OMP_NUM_THREADS', '')
    if threads.isdigit() and info.get('blas'):
        for lib in info['blas']:
            if lib['threads'] is not None and lib['threads'] < int(threads):
                print(f'run.py: {Path(lib["lib"]).name} serves '
                      f'{lib["threads"]} threads, fewer than '
                      f'OMP_NUM_THREADS={threads}')


def main():
    parser = argparse.ArgumentParser(
        description='Run benchmark cases with and without RI.')
    parser.add_argument('targets',
                        nargs='*',
                        help='suite names (directories under cases/), '
                        'directories of .inp files, or .inp files')
    parser.add_argument('--memory-limit',
                        default=None,
                        help='run each process under this memory limit, e.g. '
                        '8G (systemd-run); an RI fallback is then a note and '
                        'a killed process the problem')
    parser.add_argument('--ri',
                        choices=('both', 'yes', 'no'),
                        default='both',
                        help='which variants to run (default: both)')
    parser.add_argument('--machine',
                        default=socket.gethostname().split('.')[0],
                        help='a label for the machine (default: host name)')
    parser.add_argument('--launcher',
                        default='',
                        help="command put before python, e.g. 'mpirun -n 2'")
    parser.add_argument('--repeat',
                        type=int,
                        default=1,
                        help='runs of each variant; the smallest wall time is '
                        'compared (default: 1)')
    parser.add_argument('--out',
                        default=str(RESULTS),
                        help='directory for results (default: results/)')
    parser.add_argument('--name',
                        default=None,
                        help='run id (default: <date>_<machine>_<targets>)')
    parser.add_argument('--keep-files',
                        action='store_true',
                        help='keep checkpoint and h5 files of each run')
    parser.add_argument('--list',
                        action='store_true',
                        help='list suites and cases and exit')
    args = parser.parse_args()

    if args.list or not args.targets:
        listed = find_cases(args.targets or list_suites())
        for i, (suite, path) in enumerate(listed):
            if i == 0 or listed[i - 1][0] != suite:
                print(suite)
            print(f'  {path.stem:32s} {describe(path)}')
        return

    cases = find_cases(args.targets)
    if args.memory_limit:
        memory_limit_prefix(args.memory_limit)
    launcher = shlex.split(args.launcher)
    variants = {'both': (True, False), 'yes': (True,), 'no': (False,)}[args.ri]

    label = '+'.join(Path(t).stem for t in args.targets)
    run_id = args.name or (
        f'{datetime.now():%Y-%m-%d_%H%M}_{args.machine}_{label}')
    out_dir = Path(args.out)
    rundir = out_dir / run_id
    json_path = out_dir / f'{run_id}.json'
    if rundir.exists() or json_path.exists():
        sys.exit(f'run.py: {rundir} exists; a run is never appended to')
    rundir.mkdir(parents=True)

    info = provenance(args.machine, args.launcher)
    info['memory_limit'] = args.memory_limit
    check_provenance(info)
    print(f'run.py: {len(cases)} cases x {len(variants)} variants at '
          f'{info["commit"]}{" (dirty)" if info["dirty"] else ""}, '
          f'writing {json_path}')

    document = {
        'schema': SCHEMA,
        'run_id': run_id,
        'targets': args.targets,
        'repeat': args.repeat,
        'provenance': info,
        'rows': [],
    }

    build_checked = False
    for suite, case_path in cases:
        for ri in variants:
            record = run_case(suite, case_path, ri, rundir, launcher,
                              args.repeat, args.keep_files, args.memory_limit)
            document['rows'].append(record)
            json_path.write_text(json.dumps(document, indent=1) + '\n')

            status = '; '.join(record['problems']) or 'ok'
            if record['notes']:
                status += f' ({"; ".join(record["notes"])})'
            print(
                f'{suite}/{case_path.stem:28s} {record["variant"]:12s} '
                f'{record["wall"]:9.1f} s '
                f'{record["peak_rss"][0] / 1024**3:6.2f} GiB  {status}',
                flush=True)

            # NOTE: the output names its build commit only for a build
            # configured with ENABLE_COMMIT_HASH.
            build = record.get('build_commit')
            commit = info.get('commit')
            if build and commit and not build_checked:
                build_checked = True
                short = min(len(build), len(commit))
                if build[:short] != commit[:short]:
                    print(f'run.py: the build commit {build} is not the '
                          f'checked-out commit {commit}')

    print(f'run.py: wrote {json_path}')


if __name__ == '__main__':
    main()
