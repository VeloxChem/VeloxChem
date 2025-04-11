from mpi4py import MPI
from pathlib import Path
import pytest
import sys

from veloxchem.molecule import Molecule
from veloxchem.mmforcefieldgenerator import MMForceFieldGenerator

try:
    from openmm import app
    from openmm import LangevinIntegrator, Platform
    from openmm.app import NoCutoff, Simulation, PDBFile, GromacsTopFile
    from openmm.unit import nanometer, md_unit_system, kelvin, picoseconds, picosecond
except ImportError:
    pass


class TestEnergyMinimization:

    def run_energy_minimization(self, molecule, partial_charges, filename, tol):

        # generate force field files

        ffgen = MMForceFieldGenerator()
        ffgen.ostream.mute()

        ffgen.partial_charges = partial_charges
        ffgen.create_topology(molecule)
        ffgen.write_openmm_files(filename, 'MOL')
        ffgen.write_gromacs_files(filename, 'MOL')

        em_tolerance = 1.0
        coords_nm = molecule.get_coordinates_in_angstrom() * 0.1
        fpath = Path(filename)

        # run EM using gromacs files

        top = GromacsTopFile(str(fpath.with_suffix('.top')))
        system = top.createSystem(NoCutoff)

        platform = Platform.getPlatformByName("CPU")
        platform.setPropertyDefaultValue("Threads", "1")
        integrator = LangevinIntegrator(300 * kelvin, 1 / picosecond,
                                        0.001 * picoseconds)

        simulation = Simulation(top.topology, system, integrator, platform)
        simulation.context.setPositions(coords_nm * nanometer)
        simulation.minimizeEnergy(tolerance=em_tolerance)

        state = simulation.context.getState(getPositions=True,
                                            getEnergy=True,
                                            getForces=True)

        energy_gmx = state.getPotentialEnergy().value_in_unit_system(
            md_unit_system)

        # run EM using openmm files

        pdb = PDBFile(str(fpath.with_suffix('.pdb')))
        forcefield = app.ForceField(str(fpath.with_suffix('.xml')))
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff)

        platform = Platform.getPlatformByName("CPU")
        platform.setPropertyDefaultValue("Threads", "1")
        integrator = LangevinIntegrator(300 * kelvin, 1 / picosecond,
                                        0.001 * picoseconds)

        simulation = Simulation(pdb.topology, system, integrator, platform)
        simulation.context.setPositions(coords_nm * nanometer)
        simulation.minimizeEnergy(tolerance=em_tolerance)

        state = simulation.context.getState(getPositions=True,
                                            getEnergy=True,
                                            getForces=True)

        energy_omm = state.getPotentialEnergy().value_in_unit_system(
            md_unit_system)

        for suffix in ['.top', '.itp', '.gro', '.xml', '.pdb']:
            fpath = Path(filename).with_suffix(suffix)
            if fpath.is_file():
                fpath.unlink()

        assert abs(energy_gmx - energy_omm) < tol

    @pytest.mark.skipif("openmm" not in sys.modules,
                        reason="openmm not available")
    def test_energy_minimization(self):

        mol_xyz_1 = """21
        xyz
        C             -1.18850         1.75976        -0.61127
        C             -0.48079         0.90103         0.45056
        C             -1.46584         0.32453         1.43433
        O             -1.69674         0.93020         2.51520
        O             -2.15886        -0.84717         1.14101
        N              0.31797        -0.15681        -0.18029
        C              1.58586         0.10684        -0.79871
        O              2.09488         1.25926        -0.74655
        C              2.34028        -1.00310        -1.47194
        N              3.10306        -1.78190        -0.50000
        S             -2.23066         0.76855        -1.73578
        H             -0.42624         2.28472        -1.22450
        H             -1.81327         2.53689        -0.12028
        H              0.20397         1.56497         1.02454
        H             -2.83532        -1.22789         1.79174
        H             -0.05555        -1.13200        -0.22710
        H              1.62817        -1.67100        -2.00245
        H              3.03476        -0.58181        -2.23077
        H              2.44257        -2.17802         0.20823
        H              3.76085        -1.14191         0.00255
        H             -3.29610         0.69310        -0.81898
        """

        mol_xyz_2 = """21
        xyz
        O             -2.94267        -0.66700        -0.99123
        C             -2.59052        -1.80018        -0.56501
        C             -3.47696        -2.97712        -0.79084
        O             -1.43526        -1.95304         0.21232
        C             -0.39619        -1.01428         0.31735
        C             -0.05669        -0.18991        -0.77366
        C              0.98591         0.73108        -0.67116
        C              1.71057         0.83588         0.51363
        C              1.40144         0.01415         1.59842
        C              0.35227        -0.92344         1.51719
        C              0.05199        -1.78232         2.69357
        O              0.75352        -1.69446         3.73790
        O             -1.00186        -2.69177         2.66739
        H             -4.24873        -3.01688         0.00560
        H             -2.87794        -3.91141        -0.76511
        H             -3.97406        -2.89731        -1.78066
        H             -0.57720        -0.27758        -1.71754
        H              1.23842         1.35866        -1.51618
        H              2.52012         1.55038         0.58990
        H              1.98723         0.11438         2.50365
        H             -1.20631        -3.27593         3.46941
        """

        mol_xyz_3 = """41
        xyz
        C             -2.73759        -2.08419        -0.21779
        C             -1.55997        -2.11414        -1.20369
        C             -1.73617        -1.01917        -2.30582
        N             -0.41126        -0.48695        -2.50834
        C              0.72617        -0.92337        -1.75863
        S              0.00180        -1.72083        -0.28214
        C              1.27608         0.40980        -1.64356
        C              0.14535         0.88529        -2.50176
        O             -0.22044         1.99774        -2.96293
        N              1.21211         0.97578        -0.29400
        C              2.24581         0.73654         0.66976
        O              3.26457         0.06158         0.35663
        C              2.14022         1.32041         2.04513
        C              2.49021         2.78264         2.01430
        C              1.48044         3.75226         1.88759
        C              1.81147         5.10932         1.81988
        C              3.15054         5.50694         1.87144
        C              4.16136         4.54806         1.98632
        C              3.83495         3.18965         2.05346
        C             -2.33179        -1.49740        -3.60168
        O             -1.59123        -1.90016        -4.53938
        O             -3.71440        -1.53328        -3.74738
        C             -1.40780        -3.53726        -1.77835
        H             -3.69590        -2.25852        -0.75366
        H             -2.61749        -2.86973         0.55981
        H             -2.79176        -1.09707         0.29261
        H             -2.37432        -0.19088        -1.92416
        H              1.35769        -1.61500        -2.35932
        H              2.26325         0.55377        -2.13612
        H              0.39235         1.56370        -0.01646
        H              2.82230         0.78141         2.73933
        H              1.10649         1.17428         2.42823
        H              0.43943         3.45765         1.83899
        H              1.02935         5.85267         1.72474
        H              3.40505         6.55775         1.81696
        H              5.19847         4.85696         2.01908
        H              4.62758         2.45547         2.13282
        H             -4.13818        -1.85702        -4.60904
        H             -2.34476        -3.86599        -2.27681
        H             -0.57650        -3.57421        -2.51422
        H             -1.18111        -4.25990        -0.96564
        """

        xyz_list = [mol_xyz_1, mol_xyz_2, mol_xyz_3]
        tol = 0.05

        here = Path(__file__).parent
        filename = str(here / 'data' / 'vlx_em_mol_test')

        if MPI.COMM_WORLD.Get_rank() == 0:
            for xyzstr in xyz_list:
                molecule = Molecule.read_xyz_string(xyzstr)
                # quick and likely inaccurate charges just for testing purposes
                partial_charges = molecule.get_partial_charges(
                    molecule.get_charge())
                self.run_energy_minimization(molecule, partial_charges,
                                             filename, tol)
