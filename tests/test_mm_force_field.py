from pathlib import Path
import numpy as np
import pytest
import sys

from veloxchem.mpitask import MpiTask
from veloxchem.mmforcefieldgenerator import MMForceFieldGenerator

try:
    import scipy
except ImportError:
    pass


@pytest.mark.filterwarnings(
    'ignore:.*tostring.*tobytes:DeprecationWarning:geometric')
class TestForceField:

    @pytest.mark.skipif('scipy' not in sys.modules,
                        reason='scipy not available')
    def test_force_field(self):

        # vlxtag: RKS, MM_Force_Field_Generation

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'butane.inp')

        task = MpiTask([inpfile, None])
        partial_charges = np.array([
            -0.10519, 0.023018, 0.023018, 0.023018, 0.022809, 0.006663,
            0.006663, 0.022809, 0.006663, 0.006663, -0.105188, 0.023018,
            0.023018, 0.023018
        ])

        ff_gen = MMForceFieldGenerator(task.mpi_comm, task.ostream)
        ff_gen.partial_charges = partial_charges
        ff_gen.create_topology(task.molecule)

        scan_file = str(here / 'data' / '1-5-8-11.xyz')

        ff_gen.reparametrize_dihedrals([5, 8],
                                       scan_file=scan_file,
                                       fit_extremes=False)

        fitted_barriers = []
        for (i, j, k, l), dih in ff_gen.dihedrals.items():
            if (j + 1, k + 1) == (5, 8) or (k + 1, j + 1) == (5, 8):
                if dih['multiple']:
                    fitted_barriers += list(dih['barrier'])
                else:
                    fitted_barriers.append(dih['barrier'])
        fitted_barriers = np.array(fitted_barriers)

        ref_barriers = np.array([
            2.33528300e-01, 2.33528300e-01, 5.99127082e-01, 8.62776008e-15,
            1.85250447e+00, 6.88993804e-01, 6.88993804e-01, 2.33528300e-01,
            6.88993804e-01, 6.88993804e-01, 2.33528300e-01
        ])

        assert np.max(np.abs(fitted_barriers - ref_barriers)) < 1.0e-6
