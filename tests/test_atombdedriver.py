from mpi4py import MPI
from pathlib import Path

from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import Molecule
from veloxchem.atombdedriver import AtomBdeDriver


class TestAtomBdeDriver:

    def run_atom_bde(self, molecules):

        bde = AtomBdeDriver()
        bde.ostream.mute()
        bde.mute_output = True

        bde.compute(molecules)

        # check all ranks
        assert abs(bde.mols_bdes_list[0]["unique_BDEs_hartree"][0] -
                   0.17657796026960781) < 1.0e-4

    def cleanup_h5files(self):

        files = [
            'bde_mol_1_1_final_scf.h5',
            'bde_mol_1_1_final.h5',
            'bde_mol_1_1_opt_scf.h5',
            'bde_mol_1_1_opt.h5',
            'bde_mol_1_final_scf.h5',
            'bde_mol_1_final.h5',
            'bde_mol_1_opt_scf.h5',
            'bde_mol_1_opt.h5',
            'bde_mol_1_target_H_final_scf.h5',
            'bde_mol_1_target_H_final.h5',
        ]

        # unlink files on master rank
        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            for f in files:
                p = Path(f)
                if p.is_file():
                    p.unlink()
        MPI.COMM_WORLD.barrier()

    def test_atom_bde_driver(self):

        self.cleanup_h5files()

        # use a methane molecule as test case
        xyz_str = """5

        C              1.529985000000        -0.102981000000        -0.995481000000
        H              2.583412000000        -0.313470000000        -0.718431000000
        H              1.140736000000         0.724598000000        -0.367493000000
        H              0.915727000000        -1.012123000000        -0.831416000000
        H              1.480064000000         0.189072000000        -2.064585000000
        """
        self.run_atom_bde(Molecule.read_xyz_string(xyz_str))

        self.cleanup_h5files()
