from mpi4py import MPI
from pathlib import Path
import pytest
import sys

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.mmforcefieldgenerator import MMForceFieldGenerator

try:
    import rdkit
except ImportError:
    pass


class TestMMForceFieldGenerator:

    @pytest.mark.skipif('rdkit' not in sys.modules,
                        reason='rdkit not available')
    def test_ffgen_water(self):

        ff_gen = MMForceFieldGenerator()
        ff_gen.ostream.mute()
        ff_gen.create_water(water_model='spce')

        fname = '_vlx_test_water_'
        ff_gen.write_openmm_files(fname)

        here = Path(__file__).parent
        xmlfile_ref = str(here / 'data' / 'mmffgen_water.xml')

        assert self.match_xml(f'{fname}.xml', xmlfile_ref, 1e-10)

        self.clean_up_xml_and_pdb(fname)

    def test_ffgen_transition_metal(self):

        xyzstr = """25
        xyz
        Cl   -2.1166   -0.1503   -4.1637
        Ni   -1.8895    0.2693   -2.0458
        N    -1.6050    0.6134   -0.1842
        C    -2.4664    1.0748    0.7434
        C    -2.0815    1.2718    2.0886
        C    -0.7782    0.9909    2.5055
        C     0.1541    0.5027    1.5504
        C    -0.3211    0.3348    0.2191
        C     0.5368   -0.1488   -0.8081
        C     1.8942   -0.4781   -0.5331
        C     2.6710   -0.9531   -1.6243
        C     2.0665   -1.0663   -2.8785
        C     0.7098   -0.7158   -3.0608
        N    -0.0508   -0.2627   -2.0451
        Cl   -3.9666    0.8936   -1.9495
        H    -3.4816    1.2810    0.3624
        H    -2.8315    1.6511    2.7988
        H     2.6351   -1.4301   -3.7475
        H     0.1854   -0.7864   -4.0295
        C     1.5316    0.1656    1.8165
        C     2.3641   -0.3035    0.8199
        H     3.4123   -0.5541    1.0477
        H     1.9135    0.2905    2.8421
        H    -0.4707    1.1419    3.5522
        H     3.7280   -1.2245   -1.4752
        """

        tm_complex = Molecule.read_xyz_string(xyzstr)

        ff_gen = MMForceFieldGenerator()
        ff_gen.ostream.mute()
        ff_gen.create_topology(tm_complex, resp=False)

        fname = '_vlx_test_tm_complex_'
        ff_gen.write_openmm_files(fname)

        here = Path(__file__).parent
        xmlfile_ref = str(here / 'data' / 'mmffgen_tm_noresp.xml')

        assert self.match_xml(f'{fname}.xml', xmlfile_ref, 1e-10)

        self.clean_up_xml_and_pdb(fname)

    def clean_up_xml_and_pdb(self, fname):

        if MPI.COMM_WORLD.Get_rank() == mpi_master():

            ffgen_xml_file = Path(f'{fname}.xml')
            if ffgen_xml_file.is_file():
                ffgen_xml_file.unlink()

            ffgen_pdb_file = Path(f'{fname}.pdb')
            if ffgen_pdb_file.is_file():
                ffgen_pdb_file.unlink()

    def match_xml(self, fname_data, fname_ref, tolerance):

        data = self.read_xml(fname_data)
        data_ref = self.read_xml(fname_ref)

        if sorted(list(data.keys())) != sorted(list(data_ref.keys())):
            return False

        for key in data:
            assert len(data[key]) == len(data_ref[key])
            for val, val_ref in zip(data[key], data_ref[key]):
                if isinstance(val, str):
                    if val != val_ref:
                        return False
                else:
                    if abs(val - val_ref) >= tolerance:
                        return False

        return True

    @staticmethod
    def read_xml(fname):

        # This is a simple function that reads an xml file and collects the
        # double-quoted text.

        with open(fname, "r") as fh:
            xml_lines = fh.readlines()

        data = {}

        for line in xml_lines:
            if 'xml version' in line:
                continue

            if '="' not in line:
                continue

            key = line.split()[0].strip().replace('<', '')
            assert len(key) > 0

            if key not in data:
                data[key] = []

            content = line.split('"')
            assert len(content) > 0

            for idx in range((len(content) - 1) // 2):
                valstr = content[idx * 2 + 1]
                try:
                    val = float(valstr)
                    data[key].append(val)
                except ValueError:
                    data[key].append(valstr)

        return data
