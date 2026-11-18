"""Electrode subpackage (veloxchem.electrode): graphite builder and electrode-system checks.

The checks live in tests/electrode_checks/ as scripts (they share state between steps) and
run here in a subprocess. FEP-driver checks must all execute: --require-fep turns a skipped
FEP check into a failure (the branch ships the 3-stage production driver).
"""
import subprocess
import sys
from pathlib import Path

import pytest

CHECKS = Path(__file__).resolve().parent / 'electrode_checks'


def _run(script, *args, timeout=3600):
    r = subprocess.run([sys.executable, str(CHECKS / script), *args], cwd=str(CHECKS),
                       capture_output=True, text=True, timeout=timeout)
    assert r.returncode == 0, r.stdout[-4000:] + r.stderr[-4000:]
    return r.stdout


def test_electrode_graphite():
    assert 'all graphite checks passed' in _run('graphite_checks.py')


def test_electrode_system():
    out = _run('system_checks.py', '--require-fep')
    assert 'all system checks passed (no skips)' in out


@pytest.mark.timeconsuming
def test_electrode_system_tiny_fep():
    out = _run('system_checks.py', '--fep', '--require-fep')
    assert 'all system checks passed (no skips)' in out
    assert 'FEP result records the driver file and SHA-256' in out
