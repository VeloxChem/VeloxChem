from veloxchem.evbdriver import EvbDriver
from veloxchem.evbffbuilder import EvbForceFieldBuilder
from veloxchem.evbsystembuilder import EvbSystemBuilder
from veloxchem.evbdataprocessing import EvbDataProcessing
from pathlib import Path
import veloxchem as vlx
import openmm as mm
import numpy as np


class TestEvb:

    def test_forcefield_builder(self):
        # build reactant and product forcefields from unordered xyz inputs and compare outputs with reference
        ffbuilder = EvbForceFieldBuilder()

        ethanol_xyz = """
        9

        C         -2.15527        1.01475        0.02411
        C         -3.06542       -0.21462        0.05176
        H         -1.13231        0.72209       -0.26791
        H         -2.09671        1.44552        1.03309
        H         -2.66108       -0.95762        0.77380
        H         -3.10001       -0.63761       -0.95125
        O         -4.35015        0.17497        0.44422
        H         -2.54504        1.76624       -0.66695
        H         -4.87808       -0.68100        0.42754
        """

        ethene_xyz = """
        6

        C         -1.45521       -1.62344        0.00000
        C         -0.12198       -1.66092        0.00000
        H         -1.98694       -0.94219       -0.65611
        H         -2.02437       -2.27373        0.65611
        H          0.44719       -1.01062       -0.65611
        H          0.40976       -2.34216        0.65611
        """

        water_xyz = """
        3

        H         -0.59707        0.49671        0.00000
        H         -1.89040        1.13075       -0.65905
        O         -1.56707        0.49671        0.00000
        """

        reactant_input = {
            "molecule": vlx.Molecule.from_xyz_string(ethanol_xyz),
            "optimise": False,
            "forcefield": None,
            "hessian": None,
            "charges": None,
        }
        product_input = [
            {
                "molecule": vlx.Molecule.from_xyz_string(ethene_xyz),
                "optimise": False,
                "forcefield": None,
                "hessian": None,
                "charges": None,
            }
            ,
            {
                "molecule": vlx.Molecule.from_xyz_string(water_xyz),
                "optimise": False,
                "forcefield": None,
                "hessian": None,
                "charges": None,
            }
        ]

        reactant, product = ffbuilder.build_forcefields(
            reactant_input,
            product_input,
            reactant_charge=0,
            product_charge=[0,0],
            reactant_multiplicity=1,
            product_multiplicity=[1,1]
        )

        here = Path(__file__).parent
        reapath = str(here / 'data' / 'evb_ethanol_ff_data.json')
        propath = str(here / 'data' / 'evb_ethene_H2O_ff_data.json')
        reactant_ref = EvbDriver.load_forcefield_from_json(reapath)
        product_ref = EvbDriver.load_forcefield_from_json(propath)

        assert TestEvb._compare_dict(reactant.bonds, reactant_ref.bonds)
        assert TestEvb._compare_dict(reactant.angles, reactant_ref.angles)
        assert TestEvb._compare_dict(reactant.dihedrals, reactant_ref.dihedrals)
        assert TestEvb._compare_dict(reactant.impropers, reactant_ref.impropers)

        assert TestEvb._compare_dict(product.bonds, product_ref.bonds)
        assert TestEvb._compare_dict(product.angles, product_ref.angles)
        assert TestEvb._compare_dict(product.dihedrals, product_ref.dihedrals)
        assert TestEvb._compare_dict(product.impropers, product_ref.impropers)
    
    """
    Compare two dictionaries recursively while rounding floats to 1e-6. Keys need to match exactly
    """
    @staticmethod
    def _compare_dict(dict1, dict2, float_tol=1e-2):
        if len(dict1.keys()) != len(dict2.keys()):
            print(f"Key count mismatch: {len(dict1.keys())} != {len(dict2.keys())}")
            return False

        for (key, val1) in dict1.items():
            if key not in dict2:
                print(f"Could not find key {key} in second dictionary")
                return False

            type1 = type(val1)
            val2 = dict2[key]

            if type1 != type(val2):
                try:
                    val2 = type1(val2)
                except (ValueError, TypeError):
                    print(f"Type mismatch: {type1} != {type(val2)} for key {key}")
                    return False
            if key == 'comment':
                pass
            elif type1 == dict:
                TestEvb._compare_dict(val1, val2)
            elif type1 == float:
                if abs(val1 - val2) > float_tol:
                    print(f"Value mismatch: {val1} != {val2} for key {key}")
                    return False
            elif type1 == list or type1 == np.ndarray:
                if len(val1) != len(val2):
                    print(f"Length mismatch: {len(val1)} != {len(val2)} for key {key}")
                    return False
                if not np.allclose(val1, val2, atol=float_tol):
                    print(f"Array mismatch: {val1} != {val2} for key {key}")
                    return False
            else:
                if val1 != val2:
                    print(f"Value mismatch: {val1} != {val2} for key {key}")
                    return False
        return True

    def test_system_builder(self):
        data_path = Path(__file__).parent / 'data'
        # load forcefields
        reapath = str(data_path / 'evb_ethanol_ff_data.json')
        propath = str(data_path / 'evb_ethene_H2O_ff_data.json')

        # build systems in water and vacuum
        ethanol_xyz = """
        9

        C         -2.15527        1.01475        0.02411
        C         -3.06542       -0.21462        0.05176
        H         -1.13231        0.72209       -0.26791
        H         -2.09671        1.44552        1.03309
        H         -2.66108       -0.95762        0.77380
        H         -3.10001       -0.63761       -0.95125
        O         -4.35015        0.17497        0.44422
        H         -2.54504        1.76624       -0.66695
        H         -4.87808       -0.68100        0.42754
        """
        reactant_mol = vlx.Molecule.from_xyz_string(ethanol_xyz)
        reactant = EvbDriver.load_forcefield_from_json(reapath)
        reactant.molecule = reactant_mol
        product = EvbDriver.load_forcefield_from_json(propath)

        EVB = EvbDriver()
        vac_conf = EVB.default_system_configurations("vacuum")
        wat_conf = EVB.default_system_configurations("water")

        # 0.4 is chosen instead of 0.5 because for lambda=0.4, 1-lambda=/=lambda
        Lambda = [0,0.4,1]
        system_builder = EvbSystemBuilder()
        vac_systems, vac_topology, vac_positions = system_builder.build_systems(reactant, product, Lambda, vac_conf)
        wat_systems, wat_topology, wat_positions = system_builder.build_systems(reactant, product, Lambda, wat_conf)

        # compare outputs with reference
        EVB.Lambda = Lambda
        assert TestEvb._compare_systems(vac_systems[0], str(data_path / 'evb_ethanol_vac_0.000_sys.xml'))
        assert TestEvb._compare_systems(vac_systems[0.4], str(data_path / 'evb_ethanol_vac_0.400_sys.xml'))
        assert TestEvb._compare_systems(vac_systems[1], str(data_path / 'evb_ethanol_vac_1.000_sys.xml'))
        assert TestEvb._compare_systems(vac_systems['reactant'], str(data_path / 'evb_ethanol_vac_recalc_reactant_sys.xml'))
        assert TestEvb._compare_systems(vac_systems['product'], str(data_path / 'evb_ethanol_vac_recalc_product_sys.xml'))

        assert TestEvb._compare_systems(wat_systems[0], str(data_path / 'evb_ethanol_solv_0.000_sys.xml'))
        assert TestEvb._compare_systems(wat_systems[0.4], str(data_path / 'evb_ethanol_solv_0.400_sys.xml'))
        assert TestEvb._compare_systems(wat_systems[1], str(data_path / 'evb_ethanol_solv_1.000_sys.xml'))
        assert TestEvb._compare_systems(wat_systems['reactant'], str(data_path / 'evb_ethanol_solv_recalc_reactant_sys.xml'))
        assert TestEvb._compare_systems(wat_systems['product'], str(data_path / 'evb_ethanol_solv_recalc_product_sys.xml'))

    @staticmethod
    def _compare_systems(system: mm.System, path: str):
        # Compare strings of serialised systems instead of the systems themselves because the systems are swig proxy's
        sys_string = mm.XmlSerializer.serialize(system)
        with open(path, 'r') as input:
            ref_string = input.read()
            sys_lines = sys_string.splitlines()
            ref_lines = ref_string.splitlines()
            result = True
            if len(sys_lines) != len(ref_lines):
                print(f"The amount of lines in the test and reference system mismatch: {len(sys_lines)} != {len(ref_lines)}")
                result = False
            min_len = min(len(sys_lines), len(ref_lines))
            for i, (sys_line, ref_line) in enumerate(zip(sys_lines[:min_len], ref_lines[:min_len])):
                if sys_line != ref_line:
                    print(f"Line mismatch on line {i}: {sys_line} != {ref_line}")
                    result = False
            return result
        return False

    def test_data_processing(self):
        # Load simulation data
        input_results = {}
        EVB = EvbDriver()
        folder = Path(__file__).parent / 'data'

        specific_results = {}

        E_file = folder / 'evb_RuVbdaO_Br_2_vacuum_Energies.dat'
        data_file = folder / 'evb_RuVbdaO_Br_2_vacuum_Data_combined.dat'
        options_file = folder / 'evb_options.json'
        specific, common = EVB._load_output_files(E_file, data_file, options_file)
        specific_results.update({'vacuum': specific})

        E_file = folder / 'evb_RuVbdaO_Br_2_water_Energies.dat'
        data_file = folder / 'evb_RuVbdaO_Br_2_water_Data_combined.dat'
        options_file = folder / 'evb_options.json'
        specific, common = EVB._load_output_files(E_file, data_file, options_file)
        specific_results.update({'water': specific})

        input_results.update(common)
        input_results.update({"configuration_results": specific_results})


        # EVB.load_initialisation(str(vac_folder), 'vacuum', skip_systems=True, skip_pdb=True)
        # EVB.load_initialisation(str(water_folder), 'water', skip_systems=True, skip_pdb=True)
        # do data processing
        dp = EvbDataProcessing()
        
        comp_results = dp.compute(input_results, 5, 10)

        # compare with final results
        reference_results = EVB._load_dict_from_h5(folder / "evb_reference_results.h5")
        assert TestEvb._compare_dict(comp_results, reference_results)
        
