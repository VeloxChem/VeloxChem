from copy import deepcopy

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.optimizationengine import OptimizationEngine
from veloxchem.optimizationdriver import OptimizationDriver


class TestOptimizeMiscellaneous:

    @staticmethod
    def get_ch3_molecule_and_basis():

        xyz_string = """
        4
        ch3
        C   -1.85334300   -0.63945100    1.29623300
        H   -2.40884500   -1.56570200    1.04276400
        H   -2.24160900   -0.22442700    2.24900500
        H   -1.98830700    0.10613700    0.48589200
        """
        molecule = Molecule.read_xyz_string(xyz_string)
        molecule.set_multiplicity(2)
        basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)
        return molecule, basis

    @staticmethod
    def run_unrestricted_opt(molecule, basis, xcfun='hf'):

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun

        opt_drv = OptimizationDriver(scf_drv)
        opt_results = opt_drv.compute(molecule, basis)

        return opt_drv, opt_results

    def test_optimization_engine_and_driver_deepcopy(self):

        molecule, basis = self.get_ch3_molecule_and_basis()

        opt_drv, opt_results = self.run_unrestricted_opt(
            molecule, basis, 'b3lyp')

        opt_engine = OptimizationEngine(opt_drv.grad_drv, molecule, basis)

        opt_engine_copy = deepcopy(opt_engine)

        assert opt_engine_copy.molecule == opt_engine.molecule
        assert opt_engine_copy.opt_unparsed_input == opt_engine.opt_unparsed_input
        assert opt_engine_copy.grad_drv.xcfun == opt_engine.grad_drv.xcfun
        assert opt_engine_copy.args[0] == opt_engine.args[0]

        opt_drv_copy = deepcopy(opt_drv)

        assert opt_drv_copy.grad_drv.xcfun == opt_drv.grad_drv.xcfun
        assert opt_drv_copy.grad_drv.ostream is opt_drv.grad_drv.ostream
