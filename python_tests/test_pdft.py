import numpy as np
import pytest

from veloxchem.veloxchemlib import (GridDriver, MolecularGrid, XCNewIntegrator)
from veloxchem.veloxchemlib import parse_xc_func, is_single_node
from veloxchem.veloxchemlib import denmat, mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestopendriver import ScfRestrictedOpenDriver
from veloxchem.aodensitymatrix import AODensityMatrix


class TestPDFT:

    def run_RODFT(self, func):

        O2_xyz = """
            O 0.0 0.0 -0.6
            O 0.0 0.0  0.6
        """
        molecule = Molecule.read_str(O2_xyz)
        molecule.set_multiplicity(3)
        basis = MolecularBasis.read(molecule, 'cc-pvdz', ostream=None)

        # Optimize ROHF wavefunction
        scfdrv = ScfRestrictedOpenDriver()
        scfdrv.ostream.state = False
        scfdrv.compute(molecule, basis)

        # Compute SLDA correction
        grid_drv = GridDriver()
        molgrid = grid_drv.generate(molecule)

        xcfun = parse_xc_func(func)
        xc_drv = XCNewIntegrator()
        molgrid_new = MolecularGrid(molgrid)
        molgrid_new.partition_grid_points()
        vxc_mat = xc_drv.integrate_vxc_fock(molecule, basis, scfdrv.density,
                                                molgrid_new,
                                                xcfun.get_func_label())

        comm = scfdrv.comm
        rank = scfdrv.rank

        # Compute total and on-top pair densities
        if rank == mpi_master():
            total_density = (scfdrv.scf_tensors['D_alpha'] +
                             scfdrv.scf_tensors['D_beta'])
            den_mat = AODensityMatrix([total_density], denmat.rest)
        else:
            den_mat = AODensityMatrix()
        den_mat.broadcast(rank, comm)

        # Only in the 2 singly occupied orbitals
        Dact = np.identity(2)
        D2act = -0.5 * np.einsum('mt,np->mntp', Dact,
                                 Dact)  # True density minus "closed shell"

        if rank == mpi_master():
            mo_act = scfdrv.mol_orbs.alpha_to_numpy()[:, [7, 8]]
        else:
            mo_act = None
        mo_act = comm.bcast(mo_act, root=mpi_master())

        pdft_ene = xc_drv.integrate_vxc_pdft(den_mat, D2act, mo_act.T.copy(),
                                         molecule, basis, molgrid,
                                         xcfun.get_func_label())

        return vxc_mat.get_energy(), pdft_ene

    @pytest.mark.skipif(not is_single_node(), reason='single node only')
    def test_O2_ROLDA(self):
        ksdft, pdft = self.run_RODFT("SLDA")
        assert abs(ksdft - pdft) < 1.0e-6

    @pytest.mark.skipif(not is_single_node(), reason='single node only')
    def test_O2_ROGGA(self):
        ksdft, pdft = self.run_RODFT("BLYP")
        # do not match in GGA case when using Li-Manni's formulation
        assert abs(-17.006969151998145 - pdft) < 1.0e-6
