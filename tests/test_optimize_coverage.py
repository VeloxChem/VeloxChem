from copy import deepcopy

import h5py
import numpy as np
from mpi4py import MPI

from veloxchem.veloxchemlib import mpi_master
from veloxchem.inputparser import read_unparsed_input_from_hdf5
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.optimizationdriver import OptimizationDriver
from veloxchem.optimizationengine import OptimizationEngine
from veloxchem.outputstream import OutputStream
from veloxchem.scfgradientdriver import ScfGradientDriver
from veloxchem.scfrestdriver import ScfRestrictedDriver


class TestOptimizationDriverCoverage:

    @staticmethod
    def make_molecule(xyz):
        return Molecule.read_xyz_string(xyz)

    @staticmethod
    def make_basis(molecule, basis_label="sto-3g"):
        return MolecularBasis.read(molecule, basis_label.upper(), ostream=None)

    @staticmethod
    def make_scf_opt_driver(molecule, basis_label="sto-3g", filename=None):
        ostream = OutputStream(None)
        scf_drv = ScfRestrictedDriver(MPI.COMM_WORLD, ostream)
        scf_drv.filename = filename
        grad_drv = ScfGradientDriver(scf_drv)
        opt_drv = OptimizationDriver(grad_drv)
        basis = MolecularBasis.read(molecule, basis_label.upper(), ostream=None)
        return opt_drv, basis

    def test_optimizationdriver_helper_methods(self, tmp_path):
        molecule = self.make_molecule("""3
            water
            O  0.000000  0.000000  0.000000
            H  0.000000  0.000000  0.950000
            H  0.000000  0.750000 -0.250000
        """)
        opt_drv, basis_not_used = self.make_scf_opt_driver(molecule)

        assert opt_drv.is_scf is True

        opt_drv.transition = True
        opt_drv.update_settings({"conv_energy": 1.0e-6, "filename": "opt_job"})
        assert opt_drv.filename == "opt_job"
        assert opt_drv.hessian == "first"

        opt_drv.conv_grms = 2.0e-5
        opt_drv.conv_gmax = 3.0e-5
        opt_drv.conv_drms = 4.0e-5
        opt_drv.conv_dmax = 5.0e-5
        opt_drv.conv_maxiter = True
        assert opt_drv.conv_flags() == [
            "energy", "1e-06", "grms", "2e-05", "gmax", "3e-05", "drms",
            "4e-05", "dmax", "5e-05", "maxiter"
        ]

        opt_drv.print_keywords()
        opt_drv.print_header()
        assert "Wang" in opt_drv.get_geometric_reference()

        missing_xyz = tmp_path / "missing.xyz"
        invalid_msg = opt_drv.get_ic_rmsd(molecule, str(missing_xyz))
        assert invalid_msg == "*** Note: invalid reference xyz file!"

        other = self.make_molecule("""2
            other
            H 0.0 0.0 0.0
            H 0.0 0.0 0.7
        """)
        invalid_labels = opt_drv.get_ic_rmsd(molecule, other)
        assert invalid_labels == "*** Note: invalid reference molecule!"

        opt_drv.print_ic_rmsd(molecule, str(missing_xyz))

        vdata_file = tmp_path / "opt.vdata_last"
        vdata_file.write_text("# mode 1\nnot printed\n# mode 2\n")
        opt_drv.print_vib_analysis(str(tmp_path / "opt"), "vdata_last")

        removable = tmp_path / "remove.me"
        removable.write_text("x")
        opt_drv.clean_up_file(removable)
        assert not removable.exists()

        opt_copy = deepcopy(opt_drv)
        assert opt_copy is not opt_drv
        assert opt_copy.grad_drv is not opt_drv.grad_drv
        assert opt_copy.filename == opt_drv.filename

    def test_optimizationdriver_write_final_hdf5_scan_and_opt(self, tmp_path):
        molecule = self.make_molecule("""2
            H2
            H 0.0 0.0 0.0
            H 0.0 0.0 0.7
        """)
        opt_drv, basis = self.make_scf_opt_driver(molecule)

        opt_h5 = tmp_path / "opt_results.h5"
        with h5py.File(opt_h5, "w"):
            pass

        xyz = molecule.get_xyz_string()
        opt_drv._write_final_hdf5(
            str(opt_h5),
            molecule,
            {
                "opt_energies": [-1.0, -1.1],
                "opt_geometries": [xyz, xyz],
            },
            basis=basis,
        )

        scan_h5 = tmp_path / "scan_results.h5"
        with h5py.File(scan_h5, "w"):
            pass

        opt_drv._write_final_hdf5(
            str(scan_h5),
            molecule,
            {
                "scan_energies": [-1.0, -0.9],
                "scan_geometries": [xyz, xyz],
            },
            basis=basis,
        )

        with h5py.File(opt_h5) as h5f:
            assert "opt/opt_energies" in h5f
            assert "opt/opt_coordinates_au" in h5f
            assert "opt/nuclear_repulsion_energies" in h5f

        with h5py.File(scan_h5) as h5f:
            assert "opt/scan_energies" in h5f
            assert "opt/scan_coordinates_au" in h5f
            assert "opt/nuclear_repulsion_energies" in h5f

    def test_optimizationdriver_restart_uses_checkpoint_geometry(
            self, tmp_path):
        checkpoint_molecule = self.make_molecule("""2
            H2
            H 0.0 0.0 0.0
            H 0.0 0.0 0.7
        """)
        input_molecule = self.make_molecule("""2
            H2
            H 0.0 0.0 0.0
            H 0.0 0.0 1.2
        """)
        basis = self.make_basis(checkpoint_molecule)

        filename = str(tmp_path / "opt_restart_real")
        scf_drv = ScfRestrictedDriver(MPI.COMM_WORLD, OutputStream(None))
        scf_drv.filename = filename
        scf_results = scf_drv.compute(checkpoint_molecule, basis)

        opt_drv = OptimizationDriver(scf_drv)
        opt_drv.ostream.mute()
        opt_drv.max_iter = 1
        opt_drv.conv_maxiter = True

        opt_results = opt_drv.compute(input_molecule, basis, scf_results)
        opt_mol = Molecule.read_xyz_string(opt_results["final_geometry"])

        input_dist = np.linalg.norm(
            input_molecule.get_coordinates_in_bohr()[0] -
            input_molecule.get_coordinates_in_bohr()[1])
        checkpoint_dist = np.linalg.norm(
            checkpoint_molecule.get_coordinates_in_bohr()[0] -
            checkpoint_molecule.get_coordinates_in_bohr()[1])
        final_dist = np.linalg.norm(opt_mol.get_coordinates_in_bohr()[0] -
                                    opt_mol.get_coordinates_in_bohr()[1])

        if opt_drv.rank == mpi_master():
            assert opt_drv.filename == scf_results["filename"]
        assert final_dist != input_dist
        assert abs(final_dist - checkpoint_dist) < abs(final_dist - input_dist)

    def test_optimizationdriver_real_scan(self, tmp_path):
        molecule = self.make_molecule("""2
            H2
            H 0.0 0.0 0.0
            H 0.0 0.0 0.8
        """)
        filename = str(tmp_path / "opt_scan_real")

        opt_drv, basis = self.make_scf_opt_driver(molecule, filename=filename)
        opt_drv.ostream.mute()
        opt_drv.filename = filename
        opt_drv.constraints = ["scan distance 1 2 0.6 1.0 3"]
        opt_drv.conv_maxiter = True

        opt_results = opt_drv.compute(molecule, basis)

        if opt_drv.rank == mpi_master():
            assert "scan_energies" in opt_results
            assert "scan_geometries" in opt_results
            assert len(opt_results["scan_energies"]) == 3
            assert len(opt_results["scan_geometries"]) == 3

            with h5py.File(filename + ".h5") as h5f:
                assert "opt/scan_energies" in h5f
                assert "opt/scan_coordinates_au" in h5f

    def test_optimizationengine_calc_new_and_deepcopy(self, tmp_path):
        molecule = self.make_molecule("""2
            H2
            H 0.0 0.0 0.0
            H 0.0 0.0 0.7
        """)
        basis = self.make_basis(molecule)
        filename = str(tmp_path / "opt_engine_real")

        scf_drv = ScfRestrictedDriver(MPI.COMM_WORLD, OutputStream(None))
        scf_drv.filename = filename
        scf_results = scf_drv.compute(molecule, basis)

        grad_drv = ScfGradientDriver(scf_drv)
        engine = OptimizationEngine(grad_drv, molecule, basis, scf_results)
        engine.opt_unparsed_input = {"coordsys": "tric"}

        coords = molecule.get_coordinates_in_bohr().reshape(-1)
        result = engine.calc_new(coords, ".")

        checkpoint_file = scf_drv.get_checkpoint_file()
        if grad_drv.rank == mpi_master():
            opt_settings = read_unparsed_input_from_hdf5(
                checkpoint_file, group_name="opt_settings")
            assert opt_settings == {"coordsys": "tric"}

        assert isinstance(result["energy"], float)
        assert result["gradient"].shape == (3 * molecule.number_of_atoms(),)
        assert engine.opt_current_step == 1
        assert engine.lower() == "custom engine"
        assert engine.copy_scratch("src", "dest") is None

        engine_copy = deepcopy(engine)
        assert engine_copy is not engine
        assert engine_copy.molecule == engine.molecule
        assert len(engine_copy.args) == len(engine.args)
        assert engine_copy.args[0] == engine.args[0]
        if grad_drv.rank == mpi_master():
            assert engine_copy.args[1]["filename"] == engine.args[1]["filename"]
