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
QM levels of the grouped MSI construction and record-level data validity
(ADR-009; PROMPT-M3B Part A).

A grouped QM level supplies the energy and the Cartesian gradient at x, the
numerical Hessian at x by two-point central differences of fresh single-point
gradients, and its level-of-theory record. Two implementations:

- `GxtbLevel`: the xtb binary of the g-xTB development release (`--gxtb`;
  `--gfn 2` is also available), one subprocess per single point in a fresh
  directory with `--norestart`, Turbomole `coord` input in bohr, the energy
  from the `:: total energy` summary line and the gradient from the
  `gradient` file. The binary's own `--hess` is never used (it restarts every
  displaced SCC from the reference wavefunction and its g-xTB Hessian is
  wrong by 2.8e-2 per element at its default SCC accuracy; M3B report).
- `XtbPythonLevel`: GFN2-xTB through xtb-python with `set_accuracy`, a fresh
  `Calculator` per single point (the legacy `XtbDriver` does not expose the
  SCC accuracy, and at its default accuracy a few percent of the single
  points stop with gradient errors of 1e-4; M3B report section 3).

Every single point of a level runs at its SCC accuracy (0.01 by decision E2).
The displaced single points of a Hessian run in parallel worker processes
(subprocesses for the binary, persistent worker interpreters for xtb-python);
the construction keeps one writer (ADR-006).

The validity battery (`hessian_battery`, decision E4) compares a Hessian
with QM values at the record's own geometry and level that were not used to
compute it, with check single points at a tighter SCC accuracy (1e-4):

  V0  identity: the record's energy and gradient against a recomputation;
  V1  symmetry of the unsymmetrised Hessian;
  V2  central gradient differences along two random rigid-body-free
      directions (step 0.003 bohr, not the Hessian's own) against H d;
  V3  the curvature from energy second differences along the two softest
      projected modes of H and, at relaxed nodes, the softest constrained
      mode of the SOSC check (whose sign is judged on the constrained
      curvatures);
  V4  the order of the Taylor remainder |g(x + a u) - g(x) - a H u| at
      a = 0.04, 0.02, 0.01 bohr (ratio 4 for the right H, 2 otherwise);
  V5  for numerical Hessians, the pair check r_k = max |g(x + h e_k) -
      2 g(x) + g(x - h e_k)| over the displaced pairs against the Hessian's
      own median; the pairs above the limit are re-checked, their two
      single points recomputed at the check accuracy (a large r_k of a stiff
      coordinate reproduces, a failed single point does not).

Internal consistency (a model reproducing its node data) is not data
validity; only comparisons with independent QM values are.
"""

import hashlib
import os
import re
import shutil
import subprocess
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from .errorhandler import assert_msg_critical
from .groupedcoordinates import rigid_body_projector
from .groupeddatabase import level_of_theory_record
from .groupeddefinitions import hessian_validity_parameters
from .groupedrelaxation import rigid_body_modes

# the binary of the g-xTB development release on the development machine (PROMPT-M3B section 1 item 4); any other
# location is given by VLX_GXTB_BINARY or the `binary` argument
DEFAULT_GXTB_BINARY = '/home/vlind06/Programs/g-xtb/binaries/xtb-6.7.1/bin/xtb'

# the battery's protocol (ADR-009): part of the level record, because its single points are part of how the data
# are validated; the tolerances are validation parameters (groupeddefinitions.HESSIAN_VALIDITY_PARAMETERS)
BATTERY_PROTOCOL = {
    'check_scc_accuracy': 1.0e-4,
    'directional_step_bohr': 3.0e-3,
    'directional_directions': 2,
    'curvature_step_bohr': 3.0e-3,
    'curvature_soft_modes': 2,
    'taylor_amplitudes_bohr': [0.04, 0.02, 0.01],
}
HESSIAN_STEP_BOHR = 5.0e-3
# E2: SCC accuracy 0.01 for every single point; for g-xTB the user chose 0.001 on 2026-10-13, because at 0.01 its gradient
# differs from the derivative of its energy by 2.1e-6 hartree/bohr at the melatonin anchor, above the precision contract
# (0.1 x relaxation_gradient_threshold), and at 0.001 by 1.2e-7 at the same cost
SCC_ACCURACY = 1.0e-2
GXTB_SCC_ACCURACY = 1.0e-3

_ENERGY = re.compile(r'::\s+total energy\s+(-?\d+\.\d+)\s+Eh')
_VERSION = re.compile(r'(xtb version[^\n]*)')


# ------------------------------------------------------------------ helpers

def displaced_geometries(x, step):
    """The 6N geometries of two-point central differences: +h then -h for every Cartesian coordinate in turn."""
    x = np.asarray(x, dtype=np.float64).reshape(-1)
    out = []
    for k in range(x.size):
        for sign in (1.0, -1.0):
            y = x.copy()
            y[k] += sign * step
            out.append(y.reshape(-1, 3))
    return out


def central_difference_hessian(results, step, gradient_at_x=None):
    """
    The unsymmetrised Hessian with rows k = (g(x + h e_k) - g(x - h e_k)) /
    2h from the single-point results (energy, gradient) of
    `displaced_geometries`, and the pair check of every displaced pair
    r_k = max |g(x + h e_k) - 2 g(x) + g(x - h e_k)| when the gradient at x
    (same accuracy) is given: on a smooth surface h^2 times a third
    derivative, while one failed single point shows its full error.
    """

    n = len(results) // 2
    hessian = np.zeros((n, n))
    pair = None
    if gradient_at_x is not None:
        pair = np.zeros(n)
        g0 = np.asarray(gradient_at_x, dtype=np.float64).reshape(-1)
    for k in range(n):
        g_plus = np.asarray(results[2 * k][1], dtype=np.float64).reshape(-1)
        g_minus = np.asarray(results[2 * k + 1][1], dtype=np.float64).reshape(-1)
        hessian[k] = (g_plus - g_minus) / (2.0 * step)
        if pair is not None:
            pair[k] = float(np.abs(g_plus - 2.0 * g0 + g_minus).max())
    return hessian, pair


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def _atomic_numbers(labels):
    from .molecule import Molecule
    mol = Molecule(list(labels), np.zeros((len(labels), 3)), 'bohr')
    return np.array(mol.get_identifiers(), dtype=int)


# ------------------------------------------------------------------ levels

class GroupedQMLevel:
    """
    The interface a grouped construction calls (ADR-009): energy and
    gradient at x, the numerical Hessian at x, the level-of-theory record.

    :param scc_accuracy: the SCC accuracy of every single point (E2)
    :param hessian_step: the central-difference step in bohr (E3)
    :param workers: parallel workers for batches of single points (the
        displaced geometries of a Hessian, the battery)
    :param threads: OpenMP threads of a single sequential single point (a
        relaxation step); the workers of a batch run on one thread each
    """

    driver = None

    def __init__(self, scc_accuracy=SCC_ACCURACY, hessian_step=HESSIAN_STEP_BOHR, workers=16, threads=None):
        assert_msg_critical(float(scc_accuracy) > 0.0, 'grouped MSI QM level: the SCC accuracy must be positive.')
        assert_msg_critical(float(hessian_step) > 0.0, 'grouped MSI QM level: the Hessian step must be positive.')
        self.scc_accuracy = float(scc_accuracy)
        self.check_accuracy = float(BATTERY_PROTOCOL['check_scc_accuracy'])
        self.hessian_step = float(hessian_step)
        self.workers = max(1, int(workers))
        self.threads = int(threads) if threads is not None else int(os.environ.get('OMP_NUM_THREADS', '1') or 1)
        self.counts = {}
        self.wall = {}
        self._record = None
        self._record_key = None
        # one batch at a time per level (review of 2026-10-13: two threads sharing the workers mixed their results)
        self._batch_lock = threading.Lock()

    # the subclasses implement _single_point(labels, x, charge, multiplicity, accuracy, threads) and _settings()

    def _tick(self, accuracy, n, wall):
        key = f'single_point_acc_{accuracy:g}'
        self.counts[key] = self.counts.get(key, 0) + int(n)
        self.wall[key] = self.wall.get(key, 0.0) + float(wall)

    def level_record(self):
        """The restricted level-of-theory record (ADR-006 amendment 1, ADR-009); its hash separates the levels."""
        key = (self.scc_accuracy, self.hessian_step)
        if self._record is None or self._record_key != key:
            self._record_key = key
            settings = dict(self._settings())
            settings.update({'scc_accuracy': self.scc_accuracy, 'hessian': 'numerical',
                             'hessian_protocol': {'scheme': 'central_two_point', 'step_bohr': self.hessian_step,
                                                  'fresh_single_points': True},
                             'battery': dict(BATTERY_PROTOCOL)})
            self._record = level_of_theory_record(self.driver, **settings)
        return self._record

    def energy_gradient(self, labels, x, charge=0, multiplicity=1, accuracy=None):
        """(energy in hartree, gradient (N, 3) in hartree/bohr) at x, one single point at the level's SCC accuracy."""
        acc = self.scc_accuracy if accuracy is None else float(accuracy)
        t0 = time.time()
        out = self._single_point(list(labels), np.asarray(x, dtype=np.float64).reshape(-1, 3), int(charge), int(multiplicity),
                                 acc, self.threads)
        self._tick(acc, 1, time.time() - t0)
        return out

    def energy_gradient_many(self, labels, xs, charge=0, multiplicity=1, accuracy=None):
        """The single points of a batch in parallel (one thread each); the order of `xs` is kept."""
        acc = self.scc_accuracy if accuracy is None else float(accuracy)
        xs = [np.asarray(y, dtype=np.float64).reshape(-1, 3) for y in xs]
        t0 = time.time()
        with self._batch_lock:
            out = self._many(list(labels), xs, int(charge), int(multiplicity), acc)
        self._tick(acc, len(xs), time.time() - t0)
        return out

    def _many(self, labels, xs, charge, multiplicity, accuracy):
        with ThreadPoolExecutor(max_workers=self.workers) as pool:
            return list(pool.map(lambda y: self._single_point(labels, y, charge, multiplicity, accuracy, 1), xs))

    def hessian(self, labels, x, charge=0, multiplicity=1, accuracy=None):
        """
        The unsymmetrised numerical Hessian at x (E3: two-point central
        differences of fresh single-point gradients at `hessian_step`, every
        single point at `accuracy`, default the level's SCC accuracy) and its
        protocol record: step, accuracy, the number of single points, wall
        time, and the pair check (V5) with the gradient at x.
        """

        acc = self.scc_accuracy if accuracy is None else float(accuracy)
        x = np.asarray(x, dtype=np.float64).reshape(-1, 3)
        t0 = time.time()
        results = self.energy_gradient_many(labels, [x] + displaced_geometries(x, self.hessian_step), charge, multiplicity,
                                            accuracy=acc)
        hessian, pair = central_difference_hessian(results[1:], self.hessian_step, gradient_at_x=results[0][1])
        protocol = {'scheme': 'central_two_point', 'step_bohr': self.hessian_step, 'scc_accuracy': acc,
                    'single_points': len(results), 'wall_time': time.time() - t0,
                    'pair_check': pair_check_record(pair, step=self.hessian_step,
                                                    displaced=[(r[1], s[1]) for r, s in zip(results[1::2], results[2::2])])}
        return hessian, protocol

    def precision_contract(self, labels, x, charge=0, multiplicity=1, steps=(2.0e-3, 1.0e-3)):
        """
        ADR-009: the analytic gradient against central differences of the
        energy at the level's SCC accuracy, Richardson-extrapolated over the
        two steps (the h^2 truncation removed), over every component. Returns
        the record (max and median |g - dE/dx|, the steps, the single points).
        """

        x = np.asarray(x, dtype=np.float64).reshape(-1, 3)
        t0 = time.time()
        e0, g0 = self.energy_gradient(labels, x, charge, multiplicity)
        diffs = []
        for h in steps:
            res = self.energy_gradient_many(labels, displaced_geometries(x, h), charge, multiplicity)
            diffs.append(np.array([(res[2 * k][0] - res[2 * k + 1][0]) / (2.0 * h) for k in range(x.size)]))
        ratio = (steps[0] / steps[1]) ** 2
        extrapolated = (ratio * diffs[1] - diffs[0]) / (ratio - 1.0)
        err = np.abs(extrapolated - g0.reshape(-1))
        plain = np.abs(diffs[1] - g0.reshape(-1))
        return {'max': float(err.max()), 'median': float(np.median(err)), 'argmax': int(np.argmax(err)),
                'plain_central_difference_max': float(plain.max()), 'steps_bohr': list(steps), 'scc_accuracy': self.scc_accuracy,
                'single_points': 1 + 2 * x.size * len(steps), 'wall_time': time.time() - t0}

    def close(self):
        """Releases worker resources (subclasses)."""


class GxtbLevel(GroupedQMLevel):
    """
    The xtb binary of the g-xTB development release (ADR-009).

    :param method: 'gxtb' (`--gxtb`, default) or 'gfn2' (`--gfn 2`)
    :param scc_accuracy: default 0.001 for g-xTB (the user's decision of
        2026-10-13, the precision contract), 0.01 for GFN2
    :param binary: the xtb binary; default VLX_GXTB_BINARY, else
        DEFAULT_GXTB_BINARY
    :param work_directory: where the per-call directories are made (default:
        a new temporary directory, removed by `close`)
    """

    driver = 'GxtbLevel'

    def __init__(self, method='gxtb', binary=None, scc_accuracy=None, hessian_step=HESSIAN_STEP_BOHR, workers=16, threads=None,
                 work_directory=None):
        if scc_accuracy is None:
            scc_accuracy = GXTB_SCC_ACCURACY if method == 'gxtb' else SCC_ACCURACY
        super().__init__(scc_accuracy=scc_accuracy, hessian_step=hessian_step, workers=workers, threads=threads)
        assert_msg_critical(method in ('gxtb', 'gfn2'), f'grouped MSI QM level: unsupported xtb binary method {method!r}.')
        self.method = method
        self.binary = os.path.abspath(binary or os.environ.get('VLX_GXTB_BINARY') or DEFAULT_GXTB_BINARY)
        assert_msg_critical(os.path.isfile(self.binary) and os.access(self.binary, os.X_OK),
                            f'grouped MSI QM level: the xtb binary {self.binary} does not exist or is not executable.')
        self._own_work = work_directory is None
        self.work_directory = os.path.abspath(work_directory) if work_directory else tempfile.mkdtemp(prefix='vlx_grouped_xtb_')
        os.makedirs(self.work_directory, exist_ok=True)
        self.version_line = self._version()
        self.sha256 = _sha256(self.binary)

    def _version(self):
        p = subprocess.run([self.binary, '--version'], capture_output=True, text=True)
        found = _VERSION.search(p.stdout + p.stderr)
        assert_msg_critical(p.returncode == 0 and found is not None,
                            f'grouped MSI QM level: `{self.binary} --version` failed: {(p.stdout + p.stderr)[-300:]}')
        return found.group(1).strip()

    def _settings(self):
        return {'program': 'xtb', 'version': self.version_line, 'sha256': self.sha256,
                'method': 'g-xTB' if self.method == 'gxtb' else 'GFN2-xTB', 'flags': self._flags()}

    def _flags(self):
        return ['--gxtb'] if self.method == 'gxtb' else ['--gfn', '2']

    def _single_point(self, labels, x, charge, multiplicity, accuracy, threads):
        directory = tempfile.mkdtemp(dir=self.work_directory)
        try:
            with open(os.path.join(directory, 'coord'), 'w') as f:
                f.write('$coord\n')
                for label, r in zip(labels, x):
                    f.write(f'{r[0]:26.18f}{r[1]:26.18f}{r[2]:26.18f}  {label.lower()}\n')
                f.write('$end\n')
            env = dict(os.environ, OMP_NUM_THREADS=str(threads), MKL_NUM_THREADS=str(threads), OPENBLAS_NUM_THREADS=str(threads),
                       OMP_STACKSIZE='1G')
            cmd = [self.binary, 'coord', *self._flags(), '--grad', '--norestart', '--acc', repr(float(accuracy)),
                   '--chrg', str(int(charge)), '--uhf', str(int(multiplicity) - 1)]
            p = subprocess.run(cmd, cwd=directory, env=env, capture_output=True, text=True)
            ok = p.returncode == 0 and 'normal termination' in p.stderr
            assert_msg_critical(ok, f'grouped MSI QM level: the xtb binary did not terminate normally (exit {p.returncode}, '
                                f'{self.method}, SCC accuracy {accuracy:g}): {p.stderr[-400:]}')
            found = _ENERGY.search(p.stdout)
            assert_msg_critical(found is not None, 'grouped MSI QM level: no total energy in the xtb output.')
            energy = float(found.group(1))
            n = len(labels)
            with open(os.path.join(directory, 'gradient')) as f:
                lines = f.read().split('\n')
            gradient = np.array([[float(t.replace('D', 'E')) for t in ln.split()] for ln in lines[2 + n:2 + 2 * n]], dtype=np.float64)
            assert_msg_critical(gradient.shape == (n, 3) and np.all(np.isfinite(gradient)),
                                'grouped MSI QM level: the xtb gradient file could not be read.')
            return energy, gradient
        finally:
            shutil.rmtree(directory, ignore_errors=True)

    def close(self):
        if self._own_work:
            shutil.rmtree(self.work_directory, ignore_errors=True)


def _xtb_python_single_point(numbers, x, charge, uhf, accuracy):
    """One GFN2-xTB single point through xtb-python with a fresh Calculator."""
    from xtb.interface import Calculator, Param
    from xtb.libxtb import VERBOSITY_MUTED
    calc = Calculator(Param.GFN2xTB, np.asarray(numbers), np.asarray(x, dtype=np.float64).reshape(-1, 3), charge=float(charge),
                      uhf=int(uhf))
    calc.set_verbosity(VERBOSITY_MUTED)
    calc.set_accuracy(float(accuracy))
    res = calc.singlepoint()
    return float(res.get_energy()), np.asarray(res.get_gradient(), dtype=np.float64).reshape(-1, 3)


# The worker of XtbPythonLevel: a separate interpreter that imports only numpy and xtb-python (not veloxchem, and not
# the caller's main module, which a multiprocessing spawn would re-execute), reads (numbers, x, charge, uhf, accuracy)
# requests from stdin and writes its replies to the file descriptor given as its argument, so that anything the
# Fortran library prints cannot corrupt them.
_XTB_PYTHON_WORKER = r"""
import os, pickle, sys
import numpy as np
from xtb.interface import Calculator, Param
from xtb.libxtb import VERBOSITY_MUTED
requests = sys.stdin.buffer
replies = os.fdopen(int(sys.argv[1]), 'wb')
while True:
    try:
        numbers, x, charge, uhf, accuracy = pickle.load(requests)
    except EOFError:
        break
    try:
        calc = Calculator(Param.GFN2xTB, np.asarray(numbers), np.asarray(x, dtype=np.float64).reshape(-1, 3),
                          charge=float(charge), uhf=int(uhf))
        calc.set_verbosity(VERBOSITY_MUTED)
        calc.set_accuracy(float(accuracy))
        res = calc.singlepoint()
        out = ('ok', float(res.get_energy()), np.asarray(res.get_gradient(), dtype=np.float64).reshape(-1, 3))
    except Exception as exc:
        out = ('error', repr(exc))
    pickle.dump(out, replies)
    replies.flush()
"""


class _XtbPythonWorker:
    """One persistent worker interpreter with one OpenMP thread."""

    def __init__(self):
        import pickle
        import sys
        self._pickle = pickle
        read_fd, write_fd = os.pipe()
        env = dict(os.environ, OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', OMP_STACKSIZE='1G')
        self.process = subprocess.Popen([sys.executable, '-c', _XTB_PYTHON_WORKER, str(write_fd)], stdin=subprocess.PIPE,
                                        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, pass_fds=(write_fd,), env=env)
        os.close(write_fd)
        self.replies = os.fdopen(read_fd, 'rb')

    def call(self, request):
        try:
            self._pickle.dump(request, self.process.stdin)
            self.process.stdin.flush()
            out = self._pickle.load(self.replies)
        except (EOFError, BrokenPipeError, OSError, self._pickle.UnpicklingError) as exc:
            assert_msg_critical(False, f'grouped MSI QM level: an xtb-python worker died (exit {self.process.poll()}): {exc!r}')
        assert_msg_critical(out[0] == 'ok', f'grouped MSI QM level: the xtb-python single point failed in a worker: {out[-1]}')
        assert_msg_critical(np.isfinite(out[1]) and np.all(np.isfinite(out[2])),
                            'grouped MSI QM level: an xtb-python single point returned a non-finite energy or gradient.')
        return out[1], out[2]

    def close(self):
        try:
            self.process.stdin.close()
            self.process.wait(timeout=10)
        except Exception:
            self.process.kill()
        self.replies.close()


class XtbPythonLevel(GroupedQMLevel):
    """
    GFN2-xTB through xtb-python with SCC accuracy control (ADR-009). A
    sequential single point runs in this process; the batches run in a pool
    of persistent worker interpreters with one OpenMP thread each (started on
    first use, released by `close`).
    """

    driver = 'XtbPythonLevel'

    def __init__(self, scc_accuracy=SCC_ACCURACY, hessian_step=HESSIAN_STEP_BOHR, workers=16, threads=None):
        super().__init__(scc_accuracy=scc_accuracy, hessian_step=hessian_step, workers=workers, threads=threads)
        import xtb
        self.version_line = f'xtb-python {getattr(xtb, "__version__", "unknown")}'
        self._numbers = {}
        self._workers = None

    def _settings(self):
        return {'program': 'xtb-python', 'version': self.version_line, 'method': 'GFN2-xTB'}

    def _numbers_of(self, labels):
        key = tuple(labels)
        if key not in self._numbers:
            self._numbers[key] = _atomic_numbers(labels)
        return self._numbers[key]

    def _single_point(self, labels, x, charge, multiplicity, accuracy, threads):
        try:
            energy, gradient = _xtb_python_single_point(self._numbers_of(labels), x, charge, multiplicity - 1, accuracy)
        except Exception as exc:   # an SCC failure of xtb-python surfaces as a QM failure
            assert_msg_critical(False, f'grouped MSI QM level: the xtb-python single point failed (SCC accuracy {accuracy:g}): {exc}')
        assert_msg_critical(np.isfinite(energy) and np.all(np.isfinite(gradient)),
                            'grouped MSI QM level: an xtb-python single point returned a non-finite energy or gradient.')
        return energy, gradient

    def _many(self, labels, xs, charge, multiplicity, accuracy):
        import queue
        if self._workers is None:
            self._workers = [_XtbPythonWorker() for _ in range(self.workers)]
        idle = queue.Queue()
        for w in self._workers:
            idle.put(w)
        numbers = self._numbers_of(labels)

        def run(y):
            w = idle.get()
            try:
                return w.call((numbers, np.asarray(y, dtype=np.float64), int(charge), int(multiplicity) - 1, float(accuracy)))
            finally:
                idle.put(w)

        with ThreadPoolExecutor(max_workers=len(self._workers)) as pool:
            return list(pool.map(run, xs))

    def close(self):
        if self._workers is not None:
            for w in self._workers:
                w.close()
            self._workers = None


# ------------------------------------------------------------------ the validity battery

def pair_check_record(pair, step=None, displaced=None):
    """
    The V5 summary of the pair check values r_k: max, median, the coordinate
    of the max, max/median, the ten largest (coordinate, r_k) and the step.
    With `displaced` (the (g(x + h e_k), g(x - h e_k)) of every k) the record
    also keeps them under the transient key '_displaced' for the re-check of
    flagged pairs; it is dropped before the result is stored.
    """

    if pair is None:
        return None
    pair = np.asarray(pair, dtype=np.float64)
    med = float(np.median(pair))
    order = np.argsort(pair)[::-1][:10]
    out = {'max': float(pair.max()), 'median': med, 'argmax': int(np.argmax(pair)),
           'ratio': float(pair.max() / med) if med > 0.0 else None,
           'largest': [[int(k), float(pair[k])] for k in order], 'step_bohr': None if step is None else float(step)}
    if displaced is not None:
        out['_displaced'] = displaced
    return out


def _rigid_free_unit(rng, x):
    p = rigid_body_projector(np.asarray(x, dtype=np.float64).reshape(-1, 3))
    v = p @ rng.normal(size=p.shape[0])
    return v / np.linalg.norm(v)


def constrained_soft_mode(coordinates, x, gradient, hessian, frozen_rows=()):
    """
    The softest mode of the SOSC check at a phase-constrained stationary point
    (the algebra of groupedrelaxation.well_posedness): with A the constraint
    Jacobian (phases, frozen rows, rigid modes), lambda the multipliers of
    grad U = A^T lambda and K = H - sum lambda grad2 phi - sum lambda B2 the
    Lagrangian Hessian, the lowest eigenpair of Z^T K Z. Returns lambda_min,
    the unit Cartesian direction v = Z e_0, v^T H v and the multiplier part
    v^T (K - H) v (from the gradient only), so that the constrained curvature
    from energies is (the energy curvature along v) + that part.
    """

    model = coordinates
    x = np.asarray(x, dtype=np.float64).reshape(model.n_atoms, 3)
    q, b_mat, b2 = model.ic.values_b_b2(x)
    _, grad_phi, hess_phi, _ = model.phases_cartesian(x, order=2, qbb2=(q, b_mat, b2))
    frozen = np.array(sorted(int(r) for r in frozen_rows), dtype=int)
    jac = [grad_phi] + ([b_mat[frozen]] if frozen.size else []) + [rigid_body_modes(x)]
    jac = np.vstack(jac)
    lam = np.linalg.lstsq(jac.T, np.asarray(gradient, dtype=np.float64).reshape(-1), rcond=None)[0]
    n_g = model.n_groups
    h_sym = 0.5 * (np.asarray(hessian) + np.asarray(hessian).T)
    shift = -np.einsum('g,gxy->xy', lam[:n_g], hess_phi)
    if frozen.size:
        shift = shift - np.einsum('k,kxy->xy', lam[n_g:n_g + frozen.size], b2[frozen])
    k_mat = h_sym + shift
    _, s, vt = np.linalg.svd(jac, full_matrices=True)
    rank = int((s > 1.0e-10 * s[0]).sum())
    z = vt[rank:].T
    reduced = z.T @ k_mat @ z
    evals, evecs = np.linalg.eigh(0.5 * (reduced + reduced.T))
    v = z @ evecs[:, 0]
    v = v / np.linalg.norm(v)
    return {'lambda_min': float(evals[0]), 'direction': v, 'hessian_curvature': float(v @ h_sym @ v),
            'multiplier_curvature': float(v @ shift @ v)}


def _finite_or_none(obj):
    """The result with every non-finite float replaced by None (canonical JSON refuses NaN and inf)."""
    if isinstance(obj, dict):
        return {k: _finite_or_none(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_finite_or_none(v) for v in obj]
    if isinstance(obj, (float, np.floating)):
        return float(obj) if np.isfinite(obj) else None
    return obj


def battery_verdict(result, parameters):
    """
    The verdict of a battery result from its stored values and the given
    tolerances, without QM (also used to re-judge a stored result under the
    current tolerances). Sets 'verdict', 'failed', 'valid' and 'parameters';
    returns the result.
    """

    tol = dict(parameters)
    verdict = {}
    if result.get('finite') is False:
        verdict['finite'] = False
    if 'V0' in result:
        v0 = result['V0']
        verdict['V0'] = (v0['energy'] is not None and v0['gradient'] is not None and v0['energy'] <= tol['hessian_validity_energy_max']
                         and v0['gradient'] <= tol['hessian_validity_gradient_max'])
    if 'V1' in result:
        verdict['V1'] = result['V1'] is not None and result['V1'] <= tol['hessian_validity_asymmetry_max']
    if 'V3' in result:
        ok3 = True
        floor = tol['hessian_validity_curvature_sign_floor']
        for name, e in result['V3'].items():
            c_h, c_e = e['hessian'], e['energy']
            if c_h is None or c_e is None:
                e['pass'] = False
                ok3 = False
                continue
            a, b = (e['lambda_min'], e['lambda_energy']) if 'lambda_min' in e else (c_h, c_e)
            close = abs(c_e - c_h) <= max(tol['hessian_validity_curvature_abs'], tol['hessian_validity_curvature_rel'] * abs(c_h))
            same_sign = (b is not None) and ((np.sign(a) == np.sign(b)) or (abs(a) <= floor and abs(b) <= floor))
            e['pass'] = bool(close and same_sign)
            ok3 = ok3 and e['pass']
        verdict['V3'] = ok3
    if 'V2' in result:
        verdict['V2'] = all(v is not None for v in result['V2']) and max(result['V2']) <= tol['hessian_validity_directional_max']
    if 'V4' in result:
        r = result['V4']['residuals']
        ratios = result['V4']['ratios']
        verdict['V4'] = (r[-1] is not None and r[-1] <= tol['hessian_validity_taylor_floor']) or \
            all(q is not None and q >= tol['hessian_validity_taylor_ratio_min'] for q in ratios)
    if isinstance(result.get('V5'), dict):
        v5 = result['V5']
        limit = max(tol['hessian_validity_pair_factor'] * v5['median'], tol['hessian_validity_pair_floor'])
        v5['limit'] = float(limit)
        flagged = [int(k) for k, r_k in v5.get('largest', []) if r_k > limit]
        if v5['max'] > limit and not v5.get('largest'):
            flagged = [int(v5['argmax'])]
        rechecked = {int(e['coordinate']): e for e in v5.get('recheck', [])}
        # a flagged pair passes when its two single points are reproduced at the check accuracy to the V0 gradient
        # tolerance: then its large r_k is a property of the surface (a stiff coordinate), not a failed single point
        verdict['V5'] = all(k in rechecked and rechecked[k]['deviation'] is not None
                            and rechecked[k]['deviation'] <= tol['hessian_validity_gradient_max'] for k in flagged)
        v5['flagged'] = flagged
    result['verdict'] = {k: bool(v) for k, v in verdict.items()}
    result['failed'] = [k for k, v in verdict.items() if not v]
    result['valid'] = not result['failed']
    result['parameters'] = tol
    return result


def hessian_battery(hessian, x, check_many, record=None, constrained=None, pair_check=None, parameters=None, seed=0):
    """
    The validity battery of one Hessian at x (ADR-009, decision E4).

    :param hessian: the unsymmetrised Hessian (3N, 3N)
    :param x: the record's geometry (N, 3), bohr
    :param check_many: callable list of geometries -> list of (energy,
        gradient) at the record's level with the check SCC accuracy
    :param record: the record's (energy, gradient) for V0, or None
    :param constrained: the result of `constrained_soft_mode` at a relaxed
        node or anchor, or None
    :param pair_check: the V5 record of a numerical Hessian (with its
        displaced gradients for the re-check of flagged pairs), or None
    :param parameters: the tolerances (hessian_validity_parameters)
    :param seed: the seed of the random directions
    :returns: a JSON-ready dict with every value, the verdict per item,
        'valid', 'failed' and the tolerances used ('parameters')
    """

    tol = hessian_validity_parameters() if parameters is None else dict(parameters)
    t0 = time.time()
    x = np.asarray(x, dtype=np.float64).reshape(-1)
    h = np.asarray(hessian, dtype=np.float64)
    if not np.all(np.isfinite(h)):
        out = {'finite': False, 'single_points': 0, 'wall_time': time.time() - t0}
        return battery_verdict(out, tol)
    hs = 0.5 * (h + h.T)
    proj = rigid_body_projector(x.reshape(-1, 3))
    w, basis = np.linalg.eigh(proj)
    basis = basis[:, w > 0.5]
    lam, vec = np.linalg.eigh(basis.T @ hs @ basis)
    directions = {}
    for k in range(int(BATTERY_PROTOCOL['curvature_soft_modes'])):
        directions[f'soft {k + 1}'] = basis @ vec[:, k]
    if constrained is not None:
        directions['constrained'] = np.asarray(constrained['direction'], dtype=np.float64).reshape(-1)
    rng = np.random.default_rng(int(seed))
    rand = [_rigid_free_unit(rng, x) for _ in range(int(BATTERY_PROTOCOL['directional_directions']))]
    u = _rigid_free_unit(rng, x)
    s3 = float(BATTERY_PROTOCOL['curvature_step_bohr'])
    s2 = float(BATTERY_PROTOCOL['directional_step_bohr'])
    amps = [float(a) for a in BATTERY_PROTOCOL['taylor_amplitudes_bohr']]
    geoms = [x]
    for d in directions.values():
        geoms += [x + s3 * d, x - s3 * d]
    for d in rand:
        geoms += [x + s2 * d, x - s2 * d]
    geoms += [x + a * u for a in amps]
    # V5: the flagged pairs' own displaced geometries, recomputed at the check accuracy
    v5 = None
    recheck = []
    if pair_check is not None:
        v5 = {k: v for k, v in pair_check.items() if not k.startswith('_')}
        limit = max(tol['hessian_validity_pair_factor'] * v5['median'], tol['hessian_validity_pair_floor'])
        displaced = pair_check.get('_displaced')
        step = v5.get('step_bohr')
        if displaced is not None and step is not None:
            for k, r_k in v5.get('largest', [])[:6]:
                if r_k > limit:
                    e_k = np.zeros(x.size)
                    e_k[int(k)] = float(step)
                    recheck.append(int(k))
                    geoms += [x + e_k, x - e_k]
    res = check_many([y.reshape(-1, 3) for y in geoms])
    finite = all(np.isfinite(r[0]) and np.all(np.isfinite(r[1])) for r in res)
    e0 = float(res[0][0])
    g0 = np.asarray(res[0][1], dtype=np.float64).reshape(-1)
    out = {'lowest_projected_eigenvalue': float(lam[0]), 'single_points': len(geoms), 'finite': bool(finite)}
    if record is not None:
        out['V0'] = {'energy': abs(e0 - float(record[0])),
                     'gradient': float(np.abs(g0 - np.asarray(record[1], dtype=np.float64).reshape(-1)).max())}
    out['V1'] = float(np.abs(h - h.T).max())
    i = 1
    v3 = {}
    for name, d in directions.items():
        (ep, gp), (em, gm) = res[i], res[i + 1]
        i += 2
        c_h = float(d @ hs @ d)
        c_e = (float(ep) - 2.0 * e0 + float(em)) / s3 ** 2
        c_g = float(d @ (np.asarray(gp) - np.asarray(gm)).reshape(-1)) / (2.0 * s3)
        entry = {'hessian': c_h, 'energy': c_e, 'gradient': c_g}
        if name == 'constrained':
            # the SOSC check is judged on the constrained curvatures: lambda_min from the Hessian, and the energy curvature
            # plus the multiplier part (gradient data only); their difference is the same as that of c_h and c_e
            entry.update({'lambda_min': float(constrained['lambda_min']), 'lambda_energy': c_e + float(constrained['multiplier_curvature'])})
        v3[name] = entry
    out['V3'] = v3
    v2 = []
    for d in rand:
        (_, gp), (_, gm) = res[i], res[i + 1]
        i += 2
        fd = (np.asarray(gp) - np.asarray(gm)).reshape(-1) / (2.0 * s2)
        hd = hs @ d
        v2.append(float(np.linalg.norm(fd - hd) / max(np.linalg.norm(hd), 1.0e-300)))
    out['V2'] = v2
    r = [float(np.linalg.norm(np.asarray(res[i + k][1]).reshape(-1) - g0 - a * (hs @ u))) for k, a in enumerate(amps)]
    i += len(amps)
    ratios = [r[k] / r[k + 1] if r[k + 1] > 0.0 else None for k in range(len(r) - 1)]
    out['V4'] = {'residuals': r, 'ratios': ratios}
    if v5 is not None:
        entries = []
        for k in recheck:
            g_plus, g_minus = pair_check['_displaced'][k]
            (_, cp), (_, cm) = res[i], res[i + 1]
            i += 2
            dev = max(float(np.abs(np.asarray(g_plus).reshape(-1) - np.asarray(cp).reshape(-1)).max()),
                      float(np.abs(np.asarray(g_minus).reshape(-1) - np.asarray(cm).reshape(-1)).max()))
            entries.append({'coordinate': int(k), 'r': float(dict((int(a), b) for a, b in v5['largest'])[k]), 'deviation': dev})
        v5['recheck'] = entries
        out['V5'] = v5
    out['wall_time'] = time.time() - t0
    return battery_verdict(_finite_or_none(out), tol)


def reevaluate_battery(result, parameters):
    """A stored battery result judged again under the given tolerances (no QM); a copy is returned."""
    import copy
    return battery_verdict(copy.deepcopy(result), parameters)


def sosc_verdict(result, curvature_min):
    """
    The SOSC classification of decision E4 from a battery result with the
    constrained mode: 'invalid' when the battery fails (a Hessian whose
    constrained curvature the energies contradict is invalid data, not a
    bifurcation); otherwise 'positive' (lambda_min >= curvature_min),
    'bifurcation' (lambda_min < -curvature_min, confirmed in sign by the
    energies) or 'below the margin'. None without a constrained mode.
    """

    c = result.get('V3', {}).get('constrained')
    if c is None:
        return None
    if not result['valid']:
        return 'invalid'
    if c['lambda_min'] >= curvature_min:
        return 'positive'
    if c['lambda_min'] < -curvature_min and c.get('lambda_energy') is not None and c['lambda_energy'] < 0.0:
        return 'bifurcation'
    return 'below the margin'


def battery_summary(result):
    """The compact form of a battery result for a ledger row and the validation record."""
    if result is None:
        return None
    v3 = result.get('V3', {})
    worst3 = max((abs(e['energy'] - e['hessian']) for e in v3.values() if e['energy'] is not None and e['hessian'] is not None),
                 default=None)
    ratios = [q for q in result.get('V4', {}).get('ratios', []) if q is not None]
    v2 = [v for v in result.get('V2', []) if v is not None]
    out = {'valid': bool(result['valid']), 'failed': list(result['failed']), 'V1': result.get('V1'),
           'V2': max(v2) if v2 else None, 'V3': worst3, 'V4': min(ratios) if ratios else None,
           'single_points': result.get('single_points')}
    if 'V0' in result:
        out['V0_energy'] = result['V0']['energy']
        out['V0_gradient'] = result['V0']['gradient']
    if isinstance(result.get('V5'), dict):
        out['V5_max'] = result['V5']['max']
        out['V5_ratio'] = result['V5']['ratio']
        if result['V5'].get('recheck'):
            out['V5_rechecked'] = [[e['coordinate'], e['deviation']] for e in result['V5']['recheck']]
    if 'constrained' in v3:
        out['lambda_min'] = v3['constrained']['lambda_min']
        out['lambda_energy'] = v3['constrained']['lambda_energy']
    return out
