import numpy as np
import time as tm

from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import XCIntegrator
from veloxchem.griddriver import GridDriver
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.dftutils import get_default_grid_level


class TestXCMolGrad:

    def run_xc_new(self, molecule, basis, xcfun):

        scf_drv = ScfRestrictedDriver()

        scf_drv.xcfun = xcfun
        scf_drv.conv_thresh = 1.0e-6
        scf_drv.ri_coulomb = True

        scf_drv.compute(molecule, basis)
        density = scf_drv.density
        density = density.broadcast(scf_drv.comm, root=mpi_master())

        grid_drv = GridDriver(scf_drv.comm)
        grid_level = (get_default_grid_level(scf_drv.xcfun)
                      if scf_drv.grid_level is None else scf_drv.grid_level)
        grid_drv.set_level(grid_level)
        mol_grid = grid_drv.generate(molecule)
        
        xc_drv = XCIntegrator()
        
        start = tm.time()
        kx_old = xc_drv.integrate_vxc_fock(molecule, basis,
                                          [density.alpha_to_numpy(0)],
                                          mol_grid, 'slda')
        end = tm.time()
        print("* OLD VXC INTEGRATION: ", end - start, " sec.  Energy: ", kx_old.get_energy(), " a.u. Electrons : ", kx_old.get_electrons())
        
        start = tm.time()
        kx_new = xc_drv.new_integrate_vxc_fock(molecule, basis,
                                               [density.alpha_to_numpy(0)],
                                               mol_grid, 'slda')
        end = tm.time()
        print("* NEW VXC INTEGRATION: ", end - start, " sec.  Energy: ", kx_new.get_energy(), " a.u. Electrons : ", kx_new.get_electrons())
        
        print("Max. diff. : ", np.max(np.abs(kx_old.alpha_to_numpy() - kx_new.alpha_to_numpy())))

    def test_xc_new_lda(self):

        molecule = Molecule.read_xyz_file('caffeine.xyz')
        basis = MolecularBasis.read(molecule, 'def2-tzvpd', ostream=None)
        
        self.run_xc_new(molecule, basis, 'slda')
        
        assert False
