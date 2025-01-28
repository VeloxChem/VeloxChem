import pickle

from mpi4py import MPI

from veloxchem import MolecularBasis
from veloxchem import Molecule
from veloxchem import BlockedGtoPairBlock
from veloxchem import T4CScreener
from veloxchem import make_gto_blocks
from veloxchem import make_gto_pair_blocks


class TestT4CScreener:

    def get_data_h2o(self):

        h2ostr = """
            O    0.000000000000        0.000000000000        0.000000000000
            H    0.000000000000        0.740848095288        0.582094932012
            H    0.000000000000       -0.740848095288        0.582094932012
        """

        mol = Molecule.read_str(h2ostr, 'angstrom')
        bas = MolecularBasis.read(mol, 'sto-3g')

        return mol, bas

    def test_pickle(self):

        mol_h2o, bas_sto3g = self.get_data_h2o()

        t4c_drv_a = T4CScreener()
        t4c_drv_a.partition(bas_sto3g, mol_h2o, 'eri')
        bobj = pickle.dumps(t4c_drv_a)
        t4c_drv_b = pickle.loads(bobj)
        assert t4c_drv_a == t4c_drv_b

    def test_partition(self):

        mol_h2o, bas_sto3g = self.get_data_h2o()

        t4c_drv = T4CScreener()
        t4c_drv.partition(bas_sto3g, mol_h2o, 'eri')

        # compare GTO pair blocks
        agp_pair_blocks = t4c_drv.gto_pair_blocks()
        b_pair_blocks = make_gto_pair_blocks(bas_sto3g, mol_h2o)

        assert len(agp_pair_blocks) == len(b_pair_blocks)
        tints_ss = [
            4.78506575181571, 0.136873367356134, 0.817206295835749,
            0.00785752036126304, 0.168705263230473, 0.774605944211488,
            0.00785752036126304, 0.168705263230473, 0.0395957018967669,
            0.774605944211488
        ]
        assert agp_pair_blocks[0] == BlockedGtoPairBlock(
            b_pair_blocks[0], tints_ss)
        tints_sp = [
            0.0244774107811672, 0.180518387631544, 0.112067640667896,
            0.112067640667896
        ]
        assert agp_pair_blocks[1] == BlockedGtoPairBlock(
            b_pair_blocks[1], tints_sp)
        tints_pp = [
            0.880159089647115,
        ]
        assert agp_pair_blocks[2] == BlockedGtoPairBlock(
            b_pair_blocks[2], tints_pp)

    def test_mpi_bcast(self):

        comm = MPI.COMM_WORLD

        mol_h2o, bas_sto3g = self.get_data_h2o()

        t4c_drv_a = None
        if comm.Get_rank() == 0:
            t4c_drv_a = T4CScreener()
            t4c_drv_a.partition(bas_sto3g, mol_h2o, 'eri')
        t4c_drv_a = comm.bcast(t4c_drv_a, 0)
        t4c_drv_b = T4CScreener()
        t4c_drv_b.partition(bas_sto3g, mol_h2o, 'eri')
        assert t4c_drv_a == t4c_drv_b
