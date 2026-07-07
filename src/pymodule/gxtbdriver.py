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

from dataclasses import dataclass
from pathlib import Path
from mpi4py import MPI
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np

from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .gradientdriver import GradientDriver
from .hessiandriver import HessianDriver


_FLOAT_RE = re.compile(
    r"[-+]?(?:(?:\d+\.\d*)|(?:\.\d+)|(?:\d+))(?:[EeDd][-+]?\d+)?")


@dataclass
class GxtbResult:
    """
    Result container returned by :class:`GxtbDriver`.

    Energies are in Hartree, gradients in Hartree/Bohr, and Hessians in
    Hartree/Bohr**2.
    """

    energy: float
    gradient: np.ndarray | None
    hessian: np.ndarray | None
    stdout: str
    stderr: str
    workdir: str | None
    metadata: dict


class GxtbDriver:
    """
    Subprocess-backed driver for the g-xTB-modified ``xtb`` executable.

    g-xTB is distributed as a modified ``xtb`` binary and activated with
    ``--gxtb``. This driver intentionally does not use xtb-python.
    """

    def __init__(self,
                 comm=None,
                 ostream=None,
                 xtb_binary=None,
                 charge=None,
                 uhf=None,
                 nthreads=1,
                 timeout=300,
                 keep_workdir=False,
                 verify_binary=True,
                 hessian_acc=1.0,
                 solvent_model=None,
                 solvent=None):
        """
        Initializes the g-xTB driver.

        Binary discovery order:
        explicit ``xtb_binary``, ``VELOXCHEM_GXTB_BINARY``, ``GXTB_BINARY``,
        then ``shutil.which("xtb")``.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        if solvent_model is not None or solvent is not None:
            raise NotImplementedError(
                'GxtbDriver: solvation is disabled. g-xTB solvation support '
                'is limited; explicit --gbe/--cosmo handling is not '
                'implemented in this wrapper.')

        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()
        self.ostream = ostream

        self.xtb_binary = self._discover_binary(xtb_binary)
        self.charge = charge
        self.uhf = uhf
        self.nthreads = int(nthreads)
        self.timeout = timeout
        self.keep_workdir = bool(keep_workdir)
        self.hessian_acc = float(hessian_acc)

        self._last_result = None
        self.energy_metadata = None

        if verify_binary:
            self._verify_binary()

    @staticmethod
    def _resolve_binary(candidate):
        if candidate is None:
            return None

        candidate = str(candidate).strip()
        if not candidate:
            return None

        path_candidate = Path(candidate).expanduser()
        if path_candidate.is_absolute() or os.sep in candidate:
            return str(path_candidate)

        resolved = shutil.which(candidate)
        return resolved if resolved is not None else candidate

    @classmethod
    def _discover_binary(cls, explicit_binary=None):
        for candidate in (
                explicit_binary,
                os.environ.get('VELOXCHEM_GXTB_BINARY'),
                os.environ.get('GXTB_BINARY'),
                shutil.which('xtb'),
        ):
            resolved = cls._resolve_binary(candidate)
            if resolved:
                return resolved

        raise RuntimeError(
            'GxtbDriver: could not find an xtb executable. Set '
            'VELOXCHEM_GXTB_BINARY or GXTB_BINARY to the g-xTB-modified xtb '
            'binary.')

    @classmethod
    def _binary_accepts_gxtb(cls, binary, timeout=20):
        try:
            help_proc = subprocess.run(
                [binary, '--help'],
                capture_output=True,
                text=True,
                timeout=timeout,
                check=False,
            )
        except (OSError, subprocess.TimeoutExpired):
            return False

        help_text = (help_proc.stdout or '') + (help_proc.stderr or '')
        if '--gxtb' in help_text.lower():
            return True

        with tempfile.TemporaryDirectory(prefix='veloxchem-gxtb-probe-') as tmpdir:
            xyz_path = Path(tmpdir) / 'probe.xyz'
            xyz_path.write_text(
                '2\nprobe\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n',
                encoding='utf-8',
            )
            env = cls._make_subprocess_env(1)
            try:
                proc = subprocess.run(
                    [
                        binary, 'probe.xyz', '--gxtb', '--chrg', '0',
                        '--uhf', '0', '--silent'
                    ],
                    cwd=tmpdir,
                    capture_output=True,
                    text=True,
                    timeout=timeout,
                    check=False,
                    env=env,
                )
            except (OSError, subprocess.TimeoutExpired):
                return False

        text = ((proc.stdout or '') + (proc.stderr or '')).lower()
        unknown_gxtb = (
            'gxtb' in text and
            ('unknown option' in text or 'unrecognized option' in text or
             'invalid option' in text)
        )

        return proc.returncode == 0 and not unknown_gxtb

    @classmethod
    def is_available(cls, xtb_binary=None):
        """
        Returns whether a g-xTB-capable ``xtb`` executable is available.
        """

        try:
            binary = cls._discover_binary(xtb_binary)
        except RuntimeError:
            return False

        return cls._binary_accepts_gxtb(binary)

    def _verify_binary(self):
        if not self._binary_accepts_gxtb(self.xtb_binary):
            raise RuntimeError(
                'GxtbDriver: the selected xtb executable does not support '
                "'--gxtb'. A g-xTB-modified xtb binary from grimme-lab/g-xtb "
                'is required. Set VELOXCHEM_GXTB_BINARY or GXTB_BINARY to '
                'that binary.')

    @staticmethod
    def _make_subprocess_env(nthreads):
        env = os.environ.copy()
        env['OMP_NUM_THREADS'] = str(int(nthreads))
        env['MKL_NUM_THREADS'] = str(int(nthreads))
        env['OMP_MAX_ACTIVE_LEVELS'] = '1'
        env['OMP_STACKSIZE'] = '4G'
        return env

    @staticmethod
    def _float_tokens(text):
        return [
            float(token.replace('D', 'E').replace('d', 'e'))
            for token in _FLOAT_RE.findall(text)
        ]

    @staticmethod
    def _format_command(command):
        return ' '.join(str(item) for item in command)

    @staticmethod
    def _tail(text, max_chars=4000):
        text = text or ''
        if len(text) <= max_chars:
            return text
        return text[-max_chars:]

    def _charge_and_uhf(self, molecule, charge=None, uhf=None):
        calc_charge = self.charge if charge is None else charge
        calc_uhf = self.uhf if uhf is None else uhf

        if calc_charge is None:
            calc_charge = molecule.get_charge()
        if calc_uhf is None:
            calc_uhf = max(int(molecule.get_multiplicity()) - 1, 0)

        return int(calc_charge), int(calc_uhf)

    def _build_command(self, xyz_name, gradient, hessian, charge, uhf):
        command = [
            self.xtb_binary,
            xyz_name,
            '--gxtb',
            '--chrg',
            str(charge),
            '--uhf',
            str(uhf),
            '--silent',
        ]

        if gradient:
            command.append('--grad')
        if hessian:
            command.extend(['--hess', '--acc', f'{self.hessian_acc:g}'])

        return command

    def _parse_energy_from_engrad(self, workdir):
        engrad_files = sorted(Path(workdir).glob('*.engrad'))
        for path in engrad_files:
            values = self._float_tokens(path.read_text(encoding='utf-8'))
            if len(values) >= 2:
                return float(values[1])
        return None

    def _parse_energy_file(self, workdir):
        energy_path = Path(workdir) / 'energy'
        if not energy_path.is_file():
            return None

        in_block = False
        for line in energy_path.read_text(encoding='utf-8').splitlines():
            stripped = line.strip()
            if stripped.lower().startswith('$energy'):
                in_block = True
                continue
            if stripped.lower().startswith('$end'):
                break
            if not in_block:
                continue
            values = self._float_tokens(line)
            if len(values) >= 2:
                return float(values[1])
            if len(values) == 1:
                return float(values[0])

        return None

    def _parse_energy_from_stdout(self, stdout):
        float_pat = _FLOAT_RE.pattern
        patterns = [
            rf'TOTAL\s+ENERGY\s+({float_pat})\s+Eh',
            rf'total\s+energy\s+({float_pat})\s+Eh',
        ]

        for pattern in patterns:
            matches = re.findall(pattern, stdout or '', flags=re.IGNORECASE)
            if matches:
                return float(matches[-1].replace('D', 'E').replace('d', 'e'))

        return None

    def _parse_energy(self, workdir, stdout):
        for parser in (
                lambda: self._parse_energy_from_engrad(workdir),
                lambda: self._parse_energy_from_stdout(stdout),
                lambda: self._parse_energy_file(workdir),
        ):
            energy = parser()
            if energy is not None:
                return energy

        raise RuntimeError(
            'GxtbDriver: failed to parse total energy from g-xTB output.')

    def _parse_engrad_gradient(self, workdir, natoms):
        engrad_files = sorted(Path(workdir).glob('*.engrad'))
        for path in engrad_files:
            values = self._float_tokens(path.read_text(encoding='utf-8'))
            if len(values) < 2 + 3 * natoms:
                continue

            file_natoms = int(round(values[0]))
            if file_natoms != natoms:
                raise RuntimeError(
                    'GxtbDriver: gradient atom count mismatch in '
                    f'{path.name}: expected {natoms}, got {file_natoms}.')

            grad_values = values[2:2 + 3 * natoms]
            return np.array(grad_values, dtype=np.float64).reshape(natoms, 3)

        return None

    def _parse_legacy_gradient(self, workdir, natoms):
        grad_path = Path(workdir) / 'gradient'
        if not grad_path.is_file():
            return None

        lines = grad_path.read_text(encoding='utf-8').splitlines()
        block = []
        in_block = False

        for line in lines:
            stripped = line.strip()
            if stripped.lower().startswith('$grad'):
                in_block = True
                continue
            if stripped.lower().startswith('$end'):
                break
            if in_block and stripped:
                block.append(line)

        if block and 'cycle' in block[0].lower():
            block = block[1:]

        if len(block) < 2 * natoms:
            return None

        gradient_lines = block[natoms:natoms + natoms]
        values = []
        for line in gradient_lines:
            row = self._float_tokens(line)
            if len(row) != 3:
                raise RuntimeError(
                    'GxtbDriver: failed to parse a gradient row in the '
                    'legacy gradient file.')
            values.extend(row)

        return np.array(values, dtype=np.float64).reshape(natoms, 3)

    def _parse_gradient(self, workdir, natoms):
        gradient = self._parse_engrad_gradient(workdir, natoms)
        if gradient is not None:
            return gradient

        gradient = self._parse_legacy_gradient(workdir, natoms)
        if gradient is not None:
            return gradient

        raise RuntimeError(
            'GxtbDriver: failed to parse analytic gradient from g-xTB files.')

    def _parse_hessian(self, workdir, natoms):
        hessian_path = Path(workdir) / 'hessian'
        if not hessian_path.is_file():
            raise RuntimeError(
                'GxtbDriver: g-xTB Hessian calculation did not write a '
                'hessian file.')

        values = self._float_tokens(hessian_path.read_text(encoding='utf-8'))
        ncart = 3 * natoms
        expected = ncart * ncart

        if len(values) != expected:
            raise RuntimeError(
                'GxtbDriver: invalid Hessian size in hessian file: expected '
                f'{expected} values for a ({ncart}, {ncart}) matrix, got '
                f'{len(values)}.')

        return np.array(values, dtype=np.float64).reshape(ncart, ncart)

    def _run_one(self,
                 molecule,
                 gradient=False,
                 hessian=False,
                 charge=None,
                 uhf=None):
        workdir_obj = None
        if self.keep_workdir:
            workdir = tempfile.mkdtemp(prefix='veloxchem-gxtb-')
        else:
            workdir_obj = tempfile.TemporaryDirectory(
                prefix='veloxchem-gxtb-')
            workdir = workdir_obj.name

        try:
            xyz_path = Path(workdir) / 'input.xyz'
            xyz_path.write_text(
                molecule.get_xyz_string(comment='VeloxChem g-XTB'),
                encoding='utf-8',
            )

            calc_charge, calc_uhf = self._charge_and_uhf(molecule, charge, uhf)
            command = self._build_command(
                xyz_path.name,
                gradient=gradient,
                hessian=hessian,
                charge=calc_charge,
                uhf=calc_uhf,
            )

            proc = subprocess.run(
                command,
                cwd=workdir,
                capture_output=True,
                text=True,
                timeout=self.timeout,
                check=False,
                env=self._make_subprocess_env(self.nthreads),
            )

            text = ((proc.stdout or '') + (proc.stderr or '')).lower()
            unknown_gxtb = (
                'gxtb' in text and
                ('unknown option' in text or 'unrecognized option' in text or
                 'invalid option' in text)
            )

            if proc.returncode != 0 or unknown_gxtb:
                raise RuntimeError(
                    'GxtbDriver: g-xTB subprocess failed.\n'
                    f'Command: {self._format_command(command)}\n'
                    f'Return code: {proc.returncode}\n'
                    f'Stdout tail:\n{self._tail(proc.stdout)}\n'
                    f'Stderr tail:\n{self._tail(proc.stderr)}')

            natoms = molecule.number_of_atoms()
            energy = self._parse_energy(workdir, proc.stdout)
            parsed_gradient = (
                self._parse_gradient(workdir, natoms) if gradient else None
            )
            parsed_hessian = (
                self._parse_hessian(workdir, natoms) if hessian else None
            )

            retained_workdir = workdir if self.keep_workdir else None
            metadata = {
                'backend': 'gxtb',
                'binary_path': self.xtb_binary,
                'command': command,
                'charge': calc_charge,
                'uhf': calc_uhf,
                'nthreads': self.nthreads,
                'timeout': self.timeout,
                'hessian_acc': self.hessian_acc,
                'workdir': retained_workdir,
            }

            return GxtbResult(
                energy=energy,
                gradient=parsed_gradient,
                hessian=parsed_hessian,
                stdout=proc.stdout or '',
                stderr=proc.stderr or '',
                workdir=retained_workdir,
                metadata=metadata,
            )

        except subprocess.TimeoutExpired as exc:
            raise RuntimeError(
                'GxtbDriver: g-xTB subprocess timed out after '
                f'{self.timeout} seconds.') from exc
        finally:
            if workdir_obj is not None:
                workdir_obj.cleanup()

    def _compute_local(self,
                       molecule,
                       gradient=False,
                       hessian=False,
                       charge=None,
                       uhf=None):
        if hessian and gradient:
            grad_result = self._run_one(
                molecule,
                gradient=True,
                hessian=False,
                charge=charge,
                uhf=uhf,
            )
            hess_result = self._run_one(
                molecule,
                gradient=False,
                hessian=True,
                charge=charge,
                uhf=uhf,
            )

            if abs(grad_result.energy - hess_result.energy) > 1.0e-5:
                raise RuntimeError(
                    'GxtbDriver: inconsistent energies from gradient and '
                    'Hessian g-xTB runs.')

            metadata = dict(grad_result.metadata)
            metadata['commands'] = [
                grad_result.metadata['command'],
                hess_result.metadata['command'],
            ]
            metadata['workdirs'] = [
                grad_result.workdir,
                hess_result.workdir,
            ]

            return GxtbResult(
                energy=grad_result.energy,
                gradient=grad_result.gradient,
                hessian=hess_result.hessian,
                stdout=(
                    grad_result.stdout + '\n'
                    '--- g-XTB Hessian run ---\n' + hess_result.stdout
                ),
                stderr=(
                    grad_result.stderr + '\n'
                    '--- g-XTB Hessian run ---\n' + hess_result.stderr
                ),
                workdir=grad_result.workdir,
                metadata=metadata,
            )

        return self._run_one(
            molecule,
            gradient=gradient,
            hessian=hessian,
            charge=charge,
            uhf=uhf,
        )

    def compute(self,
                molecule,
                gradient=False,
                hessian=False,
                point_charges=None,
                charge=None,
                uhf=None):
        """
        Runs a g-xTB calculation for a VeloxChem molecule.
        """

        if point_charges is not None:
            raise NotImplementedError(
                'GxtbDriver: point-charge embedding through $pcem is not '
                'supported by g-xTB.')

        payload = (True, None)
        if self.rank == mpi_master():
            try:
                if not self.ostream.is_muted:
                    self.print_title()
                result = self._compute_local(
                    molecule,
                    gradient=gradient,
                    hessian=hessian,
                    charge=charge,
                    uhf=uhf,
                )
                payload = (True, result)
            except Exception as exc:
                payload = (False, f'{type(exc).__name__}: {exc}')

        ok, result_or_error = self.comm.bcast(payload, root=mpi_master())
        if not ok:
            raise RuntimeError(result_or_error)

        self._last_result = result_or_error
        return self._last_result

    def get_energy(self, molecule=None, **kwargs):
        """
        Returns the last g-xTB energy or computes one for ``molecule``.
        """

        if molecule is not None:
            return self.compute(molecule, **kwargs).energy

        return None if self._last_result is None else self._last_result.energy

    def get_gradient(self, molecule=None, **kwargs):
        """
        Returns the last g-xTB gradient or computes one for ``molecule``.
        """

        if molecule is not None:
            return self.get_energy_and_gradient(molecule, **kwargs).gradient

        return None if self._last_result is None else self._last_result.gradient

    def get_hessian(self, molecule=None, **kwargs):
        """
        Returns the last g-xTB Hessian or computes one for ``molecule``.
        """

        if molecule is not None:
            return self.compute(molecule, hessian=True, **kwargs).hessian

        return None if self._last_result is None else self._last_result.hessian

    def get_energy_and_gradient(self, molecule, **kwargs):
        """
        Computes and returns energy and analytic gradient.
        """

        return self.compute(molecule, gradient=True, **kwargs)

    def get_energy_gradient_hessian(self, molecule, **kwargs):
        """
        Computes and returns energy, analytic gradient, and Hessian.
        """

        return self.compute(molecule, gradient=True, hessian=True, **kwargs)

    def get_metadata(self):
        """
        Returns metadata for the most recent g-xTB calculation.
        """

        if self._last_result is None:
            return None
        return self._last_result.metadata

    def print_title(self):
        """
        Prints title for g-XTB calculation.
        """

        self.ostream.print_blank()
        self.ostream.print_header('g-XTB Driver')
        self.ostream.print_header(13 * '=')
        self.ostream.print_blank()
        self.ostream.print_reference('Reference:')
        self.ostream.print_reference(self.get_reference())
        self.ostream.flush()

    def get_reference(self):
        """
        Gets reference string for g-XTB.
        """

        return 'S. Grimme and co-workers, grimme-lab/g-xtb'


class GxtbGradientDriver(GradientDriver):
    """
    VeloxChem gradient-driver adapter for :class:`GxtbDriver`.
    """

    def __init__(self, gxtb_drv):
        """
        Initializes g-XTB gradient driver.
        """

        super().__init__(gxtb_drv.comm, gxtb_drv.ostream)
        self.gxtb_driver = gxtb_drv
        self.energy = None
        self.metadata = None
        self.flag = 'g-XTB Gradient Driver'

    def compute(self, molecule):
        """
        Performs calculation of the g-XTB analytic gradient.
        """

        self.print_header()
        self.ostream.mute()
        result = self.gxtb_driver.get_energy_and_gradient(molecule)
        self.ostream.unmute()

        self.energy = result.energy
        self.gradient = result.gradient
        self.metadata = result.metadata

        self.print_geometry(molecule)
        self.print_gradient(molecule)
        self.ostream.print_blank()
        self.ostream.flush()

    def compute_energy(self, molecule):
        """
        Performs a g-XTB energy calculation.
        """

        self.ostream.mute()
        energy = self.gxtb_driver.get_energy(molecule)
        self.ostream.unmute()
        return energy

    def get_metadata(self):
        """
        Returns metadata for the most recent g-XTB gradient calculation.
        """

        return self.metadata


class GxtbHessianDriver(HessianDriver):
    """
    VeloxChem Hessian-driver adapter for :class:`GxtbDriver`.
    """

    def __init__(self, gxtb_drv):
        """
        Initializes g-XTB Hessian driver.
        """

        super().__init__(gxtb_drv.comm, gxtb_drv.ostream)
        self.gxtb_driver = gxtb_drv
        self.flag = 'g-XTB Hessian Driver'
        self.numerical = True
        self.elec_energy = None
        self.metadata = None

    def compute(self, molecule, *args):
        """
        Computes the g-XTB Hessian from the CLI backend.
        """

        self.print_header()
        self.ostream.mute()
        result = self.gxtb_driver.compute(molecule, hessian=True)
        self.ostream.unmute()

        self.elec_energy = result.energy
        self.hessian = result.hessian
        self.metadata = result.metadata

        if self.do_print_hessian:
            self.print_geometry(molecule)
            self.ostream.print_blank()
            self.print_hessian(molecule)

        self.ostream.print_blank()
        self.ostream.flush()

    def get_metadata(self):
        """
        Returns metadata for the most recent g-XTB Hessian calculation.
        """

        return self.metadata
