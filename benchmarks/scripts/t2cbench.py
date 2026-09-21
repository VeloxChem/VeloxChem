"""One two-center Coulomb matrix, one record: the plain driver against the one which
forms the range separated matrix beside it.

    python t2cbench.py <geometry> <fitting set> <omega> [--check]

Run as its own process and prints one json record. The laptop suite spawns it that
way rather than looping in one process, for the reason the existing table was taken
that way: the largest case holds tens of gigabytes, and a process which is killed by
the system for asking too much should cost one row and not the run.
"""
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import veloxchem as vlx
from veloxchem.veloxchemlib import (SimdTwoCenterElectronRepulsionDriver,
                                    SimdTwoCenterElectronRepulsionRsDriver)


def best_of(call, repeats=3):
    """The best of several calls after one which is not counted.

    The first call warms the pages the matrix is written into and the tables the
    Boys function reads, neither of which a second call pays for.
    """
    call()

    walls = []
    for _ in range(repeats):
        t0 = time.time()
        call()
        walls.append(time.time() - t0)

    return min(walls), walls


def main():
    geometry, fitting, omega = sys.argv[1], sys.argv[2], float(sys.argv[3])
    check = "--check" in sys.argv

    molecule = vlx.Molecule.read_xyz_file(geometry)
    basis = vlx.MolecularBasis.read(molecule, fitting.upper(), ostream=None)

    nao = basis.get_dimensions_of_basis()

    plain_drv = SimdTwoCenterElectronRepulsionDriver()
    rs_drv = SimdTwoCenterElectronRepulsionRsDriver()

    plain_wall, plain_walls = best_of(lambda: plain_drv.compute(molecule, basis))
    rs_wall, rs_walls = best_of(lambda: rs_drv.compute(molecule, basis, omega))

    record = {
        "molecule": Path(geometry).stem,
        "atoms": molecule.number_of_atoms(),
        "basis": fitting,
        "nao": nao,
        "lmax": basis.max_angular_momentum(),
        "omega": omega,
        "threads": int(os.environ.get("OMP_NUM_THREADS", os.cpu_count())),
        # NOTE: the lower triangle of one matrix. The range separated driver holds
        # two of these, which is what the memory of the largest case is about.
        "packed_gb": round(nao * (nao + 1) / 2 * 8 / 1024 ** 3, 4),
        "plain_wall": round(plain_wall, 6),
        "rs_wall": round(rs_wall, 6),
        "plain_walls": [round(w, 6) for w in plain_walls],
        "rs_walls": [round(w, 6) for w in rs_walls],
    }

    if check:
        # NOTE: the plain block of the range separated driver must be exactly what
        # the plain driver gives -- same kernels, same path -- so anything but zero
        # is a mistake in the new driver and not an accuracy question.
        one = plain_drv.compute(molecule, basis).to_numpy()
        two, att = rs_drv.compute(molecule, basis, omega)
        record["plain_block_vs_plain_driver"] = float(np.abs(two.to_numpy() - one).max())
        record["attenuated_largest"] = float(np.abs(att.to_numpy()).max())

    print(json.dumps(record))


if __name__ == "__main__":
    main()
