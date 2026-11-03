from mpi4py import MPI
import pytest

from veloxchem.veloxchemlib import XCIntegrator
from veloxchem.veloxchemlib import parse_xc_func
from veloxchem.griddriver import GridDriver
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis


@pytest.mark.solvers
class TestXCBoxAOCounts:

    @staticmethod
    def water_chain(nwaters, spacing):

        lines = []
        for i in range(nwaters):
            x = i * spacing
            lines.append(f'O {x} 0.000 0.000')
            lines.append(f'H {x} 0.757 0.587')
            lines.append(f'H {x} -0.757 0.587')

        return Molecule.read_str('\n'.join(lines))

    @staticmethod
    def grid(molecule, comm):

        grid_drv = GridDriver(comm)
        grid_drv.set_level(4)

        return grid_drv.generate(molecule)

    @staticmethod
    def every_box(molecule, basis, grid, xcfun):

        xc_drv = XCIntegrator()

        return xc_drv.get_largest_box_ao_counts_for_response_fock_xc(
            molecule, basis, grid, parse_xc_func(xcfun), 10**9)

    @pytest.mark.parametrize('xcfun', ['slater', 'blyp', 'tpssh'])
    def test_a_box_of_a_small_molecule_keeps_every_function(self, xcfun):

        molecule = self.water_chain(1, 0.0)
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)
        grid = self.grid(molecule, MPI.COMM_WORLD)

        counts = self.every_box(molecule, basis, grid, xcfun)

        largest = MPI.COMM_WORLD.allreduce(max(counts, default=0), op=MPI.MAX)

        assert largest == basis.get_dimension_of_basis()

    @pytest.mark.parametrize('xcfun', ['slater', 'blyp'])
    def test_no_box_keeps_the_functions_of_a_distant_molecule(self, xcfun):

        one = self.water_chain(1, 0.0)
        nao_one = MolecularBasis.read(one, 'def2-svp',
                                      ostream=None).get_dimension_of_basis()

        two = self.water_chain(2, 30.0)
        basis = MolecularBasis.read(two, 'def2-svp', ostream=None)
        grid = self.grid(two, MPI.COMM_WORLD)

        counts = self.every_box(two, basis, grid, xcfun)

        largest = MPI.COMM_WORLD.allreduce(max(counts, default=0), op=MPI.MAX)

        assert basis.get_dimension_of_basis() == 2 * nao_one
        assert largest == nao_one

    def test_the_counts_are_largest_first_and_cut_to_the_count(self):

        molecule = self.water_chain(4, 4.0)
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)
        grid = self.grid(molecule, MPI.COMM_WORLD)
        xcfun = parse_xc_func('blyp')

        xc_drv = XCIntegrator()

        counts = xc_drv.get_largest_box_ao_counts_for_response_fock_xc(
            molecule, basis, grid, xcfun, 10**9)

        assert counts == sorted(counts, reverse=True)
        assert all(0 <= c <= basis.get_dimension_of_basis() for c in counts)

        for count in [0, 1, 5]:
            assert xc_drv.get_largest_box_ao_counts_for_response_fock_xc(
                molecule, basis, grid, xcfun, count) == counts[:count]

    def test_lda_keeps_no_more_functions_than_gga(self):

        # NOTE: the response Fock XC integrations screen LDA by the values of
        # the functions and GGA by their values and first derivatives, so in
        # each box GGA keeps every function LDA keeps.
        molecule = self.water_chain(4, 4.0)
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)
        grid = self.grid(molecule, MPI.COMM_WORLD)

        lda = self.every_box(molecule, basis, grid, 'slater')
        gga = self.every_box(molecule, basis, grid, 'blyp')

        assert len(lda) == len(gga)
        assert all(a <= b for a, b in zip(lda, gga))
        assert lda[0] < basis.get_dimension_of_basis()

    def test_the_ranks_count_their_own_boxes(self):

        # NOTE: the boxes of all ranks are the boxes of the grid made on one
        # rank, so the counts gathered over the ranks are the counts of one.
        molecule = self.water_chain(4, 4.0)
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)

        mine = self.every_box(molecule, basis,
                              self.grid(molecule, MPI.COMM_WORLD), 'blyp')
        whole = self.every_box(molecule, basis,
                               self.grid(molecule, MPI.COMM_SELF), 'blyp')

        gathered = sorted(
            [c for counts in MPI.COMM_WORLD.allgather(mine) for c in counts],
            reverse=True)

        assert gathered == whole
