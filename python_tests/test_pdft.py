import numpy as np
from mpi4py import MPI

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestopendriver import ScfRestrictedOpenDriver
from veloxchem.veloxchemlib import GridDriver
from veloxchem.veloxchemlib import parse_xc_func
from veloxchem.veloxchemlib import XCIntegrator
from veloxchem.aodensitymatrix import AODensityMatrix
from veloxchem.veloxchemlib import denmat


class TestPDFT:

    def run_RODFT(self,func):

        O2_xyz= """
O 0.0 0.0 -0.6
O 0.0 0.0  0.6
"""
        molecule=Molecule.read_str(O2_xyz)
        molecule.set_multiplicity(3)
        basis = MolecularBasis.read(molecule,"cc-pvdz")

        #Optimize ROHF wavefunction
        scfdrv = ScfRestrictedOpenDriver()
        scfdrv.compute(molecule, basis)

        #Compute SLDA correction
        grid_drv = GridDriver()
        molgrid = grid_drv.generate(molecule)

        xcfun = parse_xc_func(func)
        xc_drv = XCIntegrator()
        vxc_mat = xc_drv.integrate(scfdrv.density, molecule,basis, molgrid, xcfun.get_func_label())

        comm=MPI.COMM_WORLD
        den_mat=AODensityMatrix()
        #Compute total and on-top pair densities
        if comm.rank==0:
            total_density = scfdrv.density.alpha_to_numpy(0) + scfdrv.density.beta_to_numpy(0)
            den_mat=AODensityMatrix([total_density],denmat.rest)
        den_mat.broadcast(comm.rank, comm)

        #Only in the 2 singly occupied orbitals
        Dact=np.identity(2)
        D2act=-0.5*np.einsum('mt,np->mntp', Dact, Dact) #True density minus "closed shell"

        mo_act=None
        if comm.rank==0:
            mo_act=scfdrv.mol_orbs.alpha_to_numpy()[:,[7,8]]
        mo_act = comm.bcast(mo_act, root=0)

        pdft_ene = xc_drv.integrate_pdft( den_mat, np.array(D2act), np.array(mo_act.transpose()), molecule, basis, molgrid, xcfun.get_func_label())
        return vxc_mat.get_energy(), pdft_ene

    def test_O2_ROLDA(self):
        ksdft, pdft = self.run_RODFT("SLDA")
        assert abs(ksdft- pdft) < 1.0e-6

    def test_O2_ROGGA(self):
        ksdft, pdft = self.run_RODFT("BLYP")
        # do not match in GGA case when using Li-Manni's formulation
        assert abs(-17.006969151998145 - pdft) < 1.0e-6
