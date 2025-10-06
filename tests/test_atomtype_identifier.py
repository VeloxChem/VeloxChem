from veloxchem.atomtypeidentifier import AtomTypeIdentifier
from veloxchem.molecule import Molecule


class TestAtomTypeIdentifier:

    def run_atomtypeidentifier(self, xyz_string, expected_atomtypes,
                               expected_equal_charges_list,
                               expected_equiv_atoms):

        atomtypeidentifier = AtomTypeIdentifier()
        atomtypeidentifier.ostream.mute()

        molecule = Molecule.read_xyz_string(xyz_string)
        gaff_atomtypes = atomtypeidentifier.generate_gaff_atomtypes(molecule)
        atom_types_dict = atomtypeidentifier.atom_types_dict

        atomtypeidentifier.identify_equivalences()
        equal_charges_list = []
        equal_charges = atomtypeidentifier.equivalent_charges
        if equal_charges:
            for eq_types in equal_charges.split(','):
                equal_charges_list.append(
                    sorted([int(eq_at) for eq_at in eq_types.split('=')]))
            equal_charges_list.sort()

        assert equal_charges_list == expected_equal_charges_list

        assert atomtypeidentifier.equivalent_atoms == expected_equiv_atoms

        assert atom_types_dict == expected_atomtypes

    def test_atomtypeidentifier_1(self):

        xyz_string = """11
        (1E)-1-bromo-3-isocyanato-1-propene
        Br -2.267871 -0.291954 -0.040277
        O 3.663540 -1.567304 0.108174
        N 2.809816 0.618208 -0.205898
        C 1.638478 1.423302 0.082027
        C 0.360464 0.765580 -0.357880
        C -0.627875 0.510467 0.482479
        C 3.186038 -0.505612 -0.014415
        H 1.767571 2.367862 -0.447493
        H 1.606371 1.652917 1.148662
        H 0.279107 0.508641 -1.406687
        H -0.597210 0.737527 1.537844
        """
        expected_atomtypes = {
            'Br1': {
                'opls': 'opls_XXX',
                'gaff': 'br',
                'uff': 'Br'
            },
            'O2': {
                'opls': 'opls_XXX',
                'gaff': 'o',
                'uff': 'O'
            },
            'N3': {
                'opls': 'opls_XXX',
                'gaff': 'n2',
                'uff': 'N'
            },
            'C4': {
                'opls': 'opls_135',
                'gaff': 'c3',
                'uff': 'C'
            },
            'C5': {
                'opls': 'opls_141',
                'gaff': 'c2',
                'uff': 'C'
            },
            'C6': {
                'opls': 'opls_141',
                'gaff': 'c2',
                'uff': 'C'
            },
            'C7': {
                'opls': 'opls_235',
                'gaff': 'c1',
                'uff': 'C'
            },
            'H8': {
                'opls': 'opls_140',
                'gaff': 'h1',
                'uff': 'H'
            },
            'H9': {
                'opls': 'opls_140',
                'gaff': 'h1',
                'uff': 'H'
            },
            'H10': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'H11': {
                'opls': 'opls_146',
                'gaff': 'h4',
                'uff': 'H'
            }
        }
        expected_equal_charges_list = [[8, 9]]
        expected_equiv_atoms = [
            'br_00', 'o_00', 'n2_00', 'c3_00', 'c2_00', 'c2_01', 'c1_00',
            'h1_00', 'h1_00', 'ha_00', 'h4_00'
        ]

        self.run_atomtypeidentifier(xyz_string, expected_atomtypes,
                                    expected_equal_charges_list,
                                    expected_equiv_atoms)

    def test_atomtypeidentifier_2(self):

        xyz_string = """16
        bicyclopropyl
        C 0.000000 1.843020 0.753782
        C 0.427743 0.611952 0.000000
        C 0.000000 1.843020 -0.753782
        C -0.427743 -0.611952 0.000000
        C 0.000000 -1.843020 0.753782
        C 0.000000 -1.843020 -0.753782
        H 0.754138 2.423213 1.266843
        H -0.956630 1.813429 1.257842
        H 1.492677 0.407118 0.000000
        H -0.956630 1.813429 -1.257842
        H 0.754138 2.423213 -1.266843
        H -1.492677 -0.407118 0.000000
        H -0.754138 -2.423213 1.266843
        H 0.956630 -1.813429 1.257842
        H 0.956630 -1.813429 -1.257842
        H -0.754138 -2.423213 -1.266843
        """
        expected_atomtypes = {
            'C1': {
                'opls': 'opls_CX',
                'gaff': 'cx',
                'uff': 'C'
            },
            'C2': {
                'opls': 'opls_CX',
                'gaff': 'cx',
                'uff': 'C'
            },
            'C3': {
                'opls': 'opls_CX',
                'gaff': 'cx',
                'uff': 'C'
            },
            'C4': {
                'opls': 'opls_CX',
                'gaff': 'cx',
                'uff': 'C'
            },
            'C5': {
                'opls': 'opls_CX',
                'gaff': 'cx',
                'uff': 'C'
            },
            'C6': {
                'opls': 'opls_CX',
                'gaff': 'cx',
                'uff': 'C'
            },
            'H7': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'H8': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'H9': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'H10': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'H11': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'H12': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'H13': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'H14': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'H15': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'H16': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            }
        }
        expected_equal_charges_list = [[1, 3, 5, 6], [2, 4],
                                       [7, 8, 10, 11, 13, 14, 15, 16], [9, 12]]
        expected_equiv_atoms = [
            'cx_00', 'cx_01', 'cx_00', 'cx_01', 'cx_00', 'cx_00', 'hc_00',
            'hc_00', 'hc_01', 'hc_00', 'hc_00', 'hc_01', 'hc_00', 'hc_00',
            'hc_00', 'hc_00'
        ]

        self.run_atomtypeidentifier(xyz_string, expected_atomtypes,
                                    expected_equal_charges_list,
                                    expected_equiv_atoms)

    def test_atomtypeidentifier_3(self):

        xyz_string = """8
        carbonochloridodithioic acid
        Cl -0.142353 1.682318 0.000006
        S 0.908059 -1.112252 -0.000003
        S -1.983399 -0.651574 -0.000005
        C -0.475009 -0.060877 0.000000
        C 2.395894 -0.059942 0.000003
        H 2.437172 0.555719 0.893120
        H 3.225785 -0.764724 0.000001
        H 2.437174 0.555726 -0.893110
        """
        expected_atomtypes = {
            'Cl1': {
                'opls': 'opls_XXX',
                'gaff': 'cl',
                'uff': 'Cl'
            },
            'S2': {
                'opls': 'opls_SS',
                'gaff': 'ss',
                'uff': 'S'
            },
            'S3': {
                'opls': 'opls_920S',
                'gaff': 's',
                'uff': 'S'
            },
            'C4': {
                'opls': 'opls_cs',
                'gaff': 'cs',
                'uff': 'C'
            },
            'C5': {
                'opls': 'opls_135',
                'gaff': 'c3',
                'uff': 'C'
            },
            'H6': {
                'opls': 'opls_140',
                'gaff': 'h1',
                'uff': 'H'
            },
            'H7': {
                'opls': 'opls_140',
                'gaff': 'h1',
                'uff': 'H'
            },
            'H8': {
                'opls': 'opls_140',
                'gaff': 'h1',
                'uff': 'H'
            }
        }
        expected_equal_charges_list = [[6, 7, 8]]
        expected_equiv_atoms = [
            'cl_00', 'ss_00', 's_00', 'cs_00', 'c3_00', 'h1_00', 'h1_00',
            'h1_00'
        ]

        self.run_atomtypeidentifier(xyz_string, expected_atomtypes,
                                    expected_equal_charges_list,
                                    expected_equiv_atoms)

    def test_atomtypeidentifier_4(self):

        xyz_string = """8
        hexachloroethane
        C 0.000000 0.000000 0.792268
        C 0.000000 0.000000 -0.792268
        Cl 0.000000 1.679244 1.403393
        Cl -1.454268 -0.839622 1.403393
        Cl 1.454268 -0.839622 1.403393
        Cl 0.000000 -1.679244 -1.403393
        Cl -1.454268 0.839622 -1.403393
        Cl 1.454268 0.839622 -1.403393
        """
        expected_atomtypes = {
            'C1': {
                'opls': 'opls_135',
                'gaff': 'c3',
                'uff': 'C'
            },
            'C2': {
                'opls': 'opls_135',
                'gaff': 'c3',
                'uff': 'C'
            },
            'Cl3': {
                'opls': 'opls_XXX',
                'gaff': 'cl',
                'uff': 'Cl'
            },
            'Cl4': {
                'opls': 'opls_XXX',
                'gaff': 'cl',
                'uff': 'Cl'
            },
            'Cl5': {
                'opls': 'opls_XXX',
                'gaff': 'cl',
                'uff': 'Cl'
            },
            'Cl6': {
                'opls': 'opls_XXX',
                'gaff': 'cl',
                'uff': 'Cl'
            },
            'Cl7': {
                'opls': 'opls_XXX',
                'gaff': 'cl',
                'uff': 'Cl'
            },
            'Cl8': {
                'opls': 'opls_XXX',
                'gaff': 'cl',
                'uff': 'Cl'
            }
        }
        expected_equal_charges_list = [[1, 2], [3, 4, 5, 6, 7, 8]]
        expected_equiv_atoms = [
            'c3_00', 'c3_00', 'cl_00', 'cl_00', 'cl_00', 'cl_00', 'cl_00',
            'cl_00'
        ]

    def test_atomtypeidentifier_5(self):

        xyz_string = """9
        propanedial
        C 1.070710 -0.368234 0.146305
        C -1.282848 0.323971 -0.320122
        C 0.020532 0.716170 0.346563
        O 2.179848 -0.153436 -0.260150
        O -1.960910 -0.597247 0.054811
        H 0.734564 -1.388064 0.407927
        H -1.572330 0.919367 -1.206311
        H 0.403972 1.666523 -0.018102
        H -0.168069 0.776197 1.422723
        """
        expected_atomtypes = {
            'C1': {
                'opls': 'opls_235',
                'gaff': 'c',
                'uff': 'C'
            },
            'C2': {
                'opls': 'opls_235',
                'gaff': 'c',
                'uff': 'C'
            },
            'C3': {
                'opls': 'opls_135',
                'gaff': 'c3',
                'uff': 'C'
            },
            'O4': {
                'opls': 'opls_XXX',
                'gaff': 'o',
                'uff': 'O'
            },
            'O5': {
                'opls': 'opls_XXX',
                'gaff': 'o',
                'uff': 'O'
            },
            'H6': {
                'opls': 'opls_146',
                'gaff': 'h4',
                'uff': 'H'
            },
            'H7': {
                'opls': 'opls_146',
                'gaff': 'h4',
                'uff': 'H'
            },
            'H8': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'H9': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            }
        }
        expected_equal_charges_list = [[1, 2], [4, 5], [6, 7], [8, 9]]
        expected_equiv_atoms = [
            'c_00', 'c_00', 'c3_00', 'o_00', 'o_00', 'h4_00', 'h4_00', 'hc_00',
            'hc_00'
        ]

        self.run_atomtypeidentifier(xyz_string, expected_atomtypes,
                                    expected_equal_charges_list,
                                    expected_equiv_atoms)

    def test_atomtypeidentifier_6(self):

        xyz_string = """10
        2-bromoethylisocyanate
        Br -1.348575 -0.496729 0.033420
        O 2.744713 -1.356928 -0.005847
        N 1.761073 0.788135 -0.179037
        C 0.676014 1.566183 0.358859
        C -0.650037 1.308808 -0.332914
        C 2.204873 -0.321207 -0.046476
        H 0.918673 2.620383 0.205614
        H 0.579485 1.404172 1.433564
        H -0.560282 1.378536 -1.411868
        H -1.408068 1.998180 0.026210
        """
        expected_atomtypes = {
            'Br1': {
                'opls': 'opls_XXX',
                'gaff': 'br',
                'uff': 'Br'
            },
            'O2': {
                'opls': 'opls_XXX',
                'gaff': 'o',
                'uff': 'O'
            },
            'N3': {
                'opls': 'opls_XXX',
                'gaff': 'n2',
                'uff': 'N'
            },
            'C4': {
                'opls': 'opls_135',
                'gaff': 'c3',
                'uff': 'C'
            },
            'C5': {
                'opls': 'opls_135',
                'gaff': 'c3',
                'uff': 'C'
            },
            'C6': {
                'opls': 'opls_235',
                'gaff': 'c1',
                'uff': 'C'
            },
            'H7': {
                'opls': 'opls_140',
                'gaff': 'h1',
                'uff': 'H'
            },
            'H8': {
                'opls': 'opls_140',
                'gaff': 'h1',
                'uff': 'H'
            },
            'H9': {
                'opls': 'opls_140',
                'gaff': 'h1',
                'uff': 'H'
            },
            'H10': {
                'opls': 'opls_140',
                'gaff': 'h1',
                'uff': 'H'
            }
        }
        expected_equal_charges_list = [[7, 8], [9, 10]]
        expected_equiv_atoms = [
            'br_00', 'o_00', 'n2_00', 'c3_00', 'c3_01', 'c1_00', 'h1_00',
            'h1_00', 'h1_01', 'h1_01'
        ]

        self.run_atomtypeidentifier(xyz_string, expected_atomtypes,
                                    expected_equal_charges_list,
                                    expected_equiv_atoms)

    def test_atomtypeidentifier_7(self):

        xyz_string = """12
        2-pyridinol
        N -1.176140 0.295001 0.000000
        C 0.000000 0.900809 0.000000
        C 1.223975 0.226351 0.000000
        C 1.190077 -1.154667 0.000000
        C -0.039417 -1.813282 0.000000
        C -1.189267 -1.042689 0.000000
        O 0.005822 2.254025 0.000000
        H 2.149101 0.783058 0.000000
        H 2.113772 -1.717700 0.000000
        H -0.101257 -2.891524 0.000000
        H -2.168354 -1.506394 0.000000
        H -0.919060 2.536214 0.000000
        """
        expected_atomtypes = {
            'N1': {
                'opls': 'opls_520',
                'gaff': 'nb',
                'uff': 'N'
            },
            'C2': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'C3': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'C4': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'C5': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'C6': {
                'opls': 'opls_521',
                'gaff': 'ca',
                'uff': 'C'
            },
            'O7': {
                'opls': 'opls_154',
                'gaff': 'oh',
                'uff': 'O'
            },
            'H8': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'H9': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'H10': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'H11': {
                'opls': 'opls_146',
                'gaff': 'h4',
                'uff': 'H'
            },
            'H12': {
                'opls': 'opls_155',
                'gaff': 'ho',
                'uff': 'H'
            }
        }
        expected_equal_charges_list = []
        expected_equiv_atoms = [
            'nb_00', 'ca_00', 'ca_01', 'ca_02', 'ca_03', 'ca_04', 'oh_00',
            'ha_00', 'ha_01', 'ha_02', 'h4_00', 'ho_00'
        ]

        self.run_atomtypeidentifier(xyz_string, expected_atomtypes,
                                    expected_equal_charges_list,
                                    expected_equiv_atoms)

    def test_atomtypeidentifier_8(self):

        xyz_string = """19
        3-methyl-1H-indole
        C 2.400164 0.371989 0.000000
        C 2.378036 -1.032831 0.000000
        C 1.225975 1.104452 0.000000
        C 1.182367 -1.731751 0.000000
        C -2.128810 -0.282841 0.000000
        C 0.000000 0.426694 0.000000
        C -1.373591 0.855619 0.000000
        C 0.002265 -0.990913 0.000000
        C -1.867833 2.265725 0.000000
        N -1.311240 -1.397591 0.000000
        H 3.352447 0.884837 0.000000
        H 3.311976 -1.578318 0.000000
        H 1.256605 2.186234 0.000000
        H 1.169380 -2.813921 0.000000
        H -3.200520 -0.392058 0.000000
        H -1.518424 2.813683 0.878543
        H -1.518424 2.813683 -0.878543
        H -2.957161 2.300824 0.000000
        H -1.628625 -2.348693 0.000000
        """
        expected_atomtypes = {
            'C1': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'C2': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'C3': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'C4': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'C5': {
                'opls': 'opls_508',
                'gaff': 'cc',
                'uff': 'C'
            },
            'C6': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'C7': {
                'opls': 'opls_508',
                'gaff': 'cd',
                'uff': 'C'
            },
            'C8': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'C9': {
                'opls': 'opls_135',
                'gaff': 'c3',
                'uff': 'C'
            },
            'N10': {
                'opls': 'opls_na',
                'gaff': 'na',
                'uff': 'N'
            },
            'H11': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'H12': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'H13': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'H14': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'H15': {
                'opls': 'opls_146',
                'gaff': 'h4',
                'uff': 'H'
            },
            'H16': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'H17': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'H18': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'H19': {
                'opls': 'opls_240',
                'gaff': 'hn',
                'uff': 'H'
            }
        }
        expected_equal_charges_list = [[16, 17, 18]]
        expected_equiv_atoms = [
            'ca_00', 'ca_01', 'ca_02', 'ca_03', 'cc_00', 'ca_04', 'cc_01',
            'ca_05', 'c3_00', 'na_00', 'ha_00', 'ha_01', 'ha_02', 'ha_03',
            'h4_00', 'hc_00', 'hc_00', 'hc_00', 'hn_00'
        ]

        self.run_atomtypeidentifier(xyz_string, expected_atomtypes,
                                    expected_equal_charges_list,
                                    expected_equiv_atoms)

    def test_atomtypeidentifier_9(self):

        xyz_string = """22
        4-propylphenol
        C -0.358008 -1.192497 -0.280254
        C -0.343198 1.192831 -0.282150
        C -1.721958 -1.189747 -0.010477
        C -1.704067 1.215219 -0.013027
        C 0.359997 -0.006025 -0.420220
        C -2.398815 0.017812 0.125458
        C 4.192562 -0.004905 0.338317
        C 1.847021 -0.015398 -0.674427
        C 2.689971 0.002765 0.610921
        O -3.741651 0.088900 0.388088
        H 0.152278 -2.141809 -0.388206
        H 0.183030 2.133209 -0.392492
        H -2.257935 -2.126629 0.088883
        H -2.239210 2.149244 0.086996
        H 4.489024 -0.895456 -0.219839
        H 4.492677 0.865907 -0.248259
        H 4.763947 0.008939 1.267358
        H 2.112002 -0.900851 -1.258622
        H 2.117837 0.849206 -1.286581
        H 2.424357 0.886004 1.197601
        H 2.420111 -0.860563 1.224708
        H -4.105938 -0.798727 0.458901
        """
        expected_atomtypes = {
            'C1': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'C2': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'C3': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'C4': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'C5': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'C6': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'C7': {
                'opls': 'opls_135',
                'gaff': 'c3',
                'uff': 'C'
            },
            'C8': {
                'opls': 'opls_135',
                'gaff': 'c3',
                'uff': 'C'
            },
            'C9': {
                'opls': 'opls_135',
                'gaff': 'c3',
                'uff': 'C'
            },
            'O10': {
                'opls': 'opls_154',
                'gaff': 'oh',
                'uff': 'O'
            },
            'H11': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'H12': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'H13': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'H14': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'H15': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'H16': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'H17': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'H18': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'H19': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'H20': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'H21': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'H22': {
                'opls': 'opls_155',
                'gaff': 'ho',
                'uff': 'H'
            }
        }
        expected_equal_charges_list = [[1, 2], [3, 4], [11, 12], [13, 14],
                                       [15, 16, 17], [18, 19], [20, 21]]
        expected_equiv_atoms = [
            'ca_00', 'ca_00', 'ca_01', 'ca_01', 'ca_02', 'ca_03', 'c3_00',
            'c3_01', 'c3_02', 'oh_00', 'ha_00', 'ha_00', 'ha_01', 'ha_01',
            'hc_00', 'hc_00', 'hc_00', 'hc_01', 'hc_01', 'hc_02', 'hc_02',
            'ho_00'
        ]

        self.run_atomtypeidentifier(xyz_string, expected_atomtypes,
                                    expected_equal_charges_list,
                                    expected_equiv_atoms)

    def test_atomtypeidentifier_10(self):

        xyz_string = """35
        hs276
        C  4.373000  0.283000 -0.025000
        C  4.415000  1.595000 -0.512000
        H  3.512000  2.197000 -0.625000
        N  5.570000  2.158000 -0.868000
        C  6.709000  1.437000 -0.751000
        H  7.623000  1.952000 -1.052000
        C  6.768000  0.125000 -0.292000
        H  7.720000 -0.402000 -0.234000
        C  5.560000 -0.488000  0.083000
        C  5.181000 -1.787000  0.570000
        H  5.834000 -2.627000  0.778000
        C  3.821000 -1.767000  0.740000
        H  3.152000 -2.538000  1.110000
        N  3.307000 -0.521000  0.391000
        C  1.957000 -0.159000  0.437000
        S  0.730000 -1.238000 -0.213000
        C  1.403000  1.005000  0.911000
        H  1.989000  1.800000  1.365000
        C -0.012000  1.036000  0.774000
        H -0.623000  1.870000  1.113000
        C -0.550000 -0.103000  0.208000
        C -1.939000 -0.434000 -0.042000
        C -2.470000 -1.612000 -0.545000
        H -1.859000 -2.469000 -0.822000
        C -3.881000 -1.591000 -0.658000
        H -4.487000 -2.414000 -1.028000
        S -3.211000  0.726000  0.293000
        C -4.440000 -0.400000 -0.246000
        C -5.874000 -0.081000 -0.247000
        O -6.746000 -0.844000 -0.616000
        O -6.116000  1.171000  0.215000
        C -7.501000  1.564000  0.245000
        H -8.073000  0.894000  0.899000
        H -7.504000  2.584000  0.637000
        H -7.928000  1.534000 -0.765000
        """
        expected_atomtypes = {
            'C1': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'C2': {
                'opls': 'opls_521',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H3': {
                'opls': 'opls_146',
                'gaff': 'h4',
                'uff': 'H'
            },
            'N4': {
                'opls': 'opls_520',
                'gaff': 'nb',
                'uff': 'N'
            },
            'C5': {
                'opls': 'opls_521',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H6': {
                'opls': 'opls_146',
                'gaff': 'h4',
                'uff': 'H'
            },
            'C7': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H8': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C9': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'C10': {
                'opls': 'opls_508',
                'gaff': 'cc',
                'uff': 'C'
            },
            'H11': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C12': {
                'opls': 'opls_508',
                'gaff': 'cd',
                'uff': 'C'
            },
            'H13': {
                'opls': 'opls_146',
                'gaff': 'h4',
                'uff': 'H'
            },
            'N14': {
                'opls': 'opls_na',
                'gaff': 'na',
                'uff': 'N'
            },
            'C15': {
                'opls': 'opls_508',
                'gaff': 'cc',
                'uff': 'C'
            },
            'S16': {
                'opls': 'opls_SS',
                'gaff': 'ss',
                'uff': 'S'
            },
            'C17': {
                'opls': 'opls_508',
                'gaff': 'cd',
                'uff': 'C'
            },
            'H18': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C19': {
                'opls': 'opls_508',
                'gaff': 'cd',
                'uff': 'C'
            },
            'H20': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C21': {
                'opls': 'opls_508',
                'gaff': 'cc',
                'uff': 'C'
            },
            'C22': {
                'opls': 'opls_508',
                'gaff': 'cc',
                'uff': 'C'
            },
            'C23': {
                'opls': 'opls_508',
                'gaff': 'cd',
                'uff': 'C'
            },
            'H24': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C25': {
                'opls': 'opls_508',
                'gaff': 'cd',
                'uff': 'C'
            },
            'H26': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'S27': {
                'opls': 'opls_SS',
                'gaff': 'ss',
                'uff': 'S'
            },
            'C28': {
                'opls': 'opls_508',
                'gaff': 'cc',
                'uff': 'C'
            },
            'C29': {
                'opls': 'opls_235',
                'gaff': 'c',
                'uff': 'C'
            },
            'O30': {
                'opls': 'opls_XXX',
                'gaff': 'o',
                'uff': 'O'
            },
            'O31': {
                'opls': 'opls_XXX',
                'gaff': 'os',
                'uff': 'O'
            },
            'C32': {
                'opls': 'opls_135',
                'gaff': 'c3',
                'uff': 'C'
            },
            'H33': {
                'opls': 'opls_140',
                'gaff': 'h1',
                'uff': 'H'
            },
            'H34': {
                'opls': 'opls_140',
                'gaff': 'h1',
                'uff': 'H'
            },
            'H35': {
                'opls': 'opls_140',
                'gaff': 'h1',
                'uff': 'H'
            }
        }
        expected_equal_charges_list = [[33, 34, 35]]
        expected_equiv_atoms = [
            'ca_00', 'ca_01', 'h4_00', 'nb_00', 'ca_02', 'h4_01', 'ca_03',
            'ha_00', 'ca_04', 'cc_00', 'ha_01', 'cc_01', 'h4_02', 'na_00',
            'cc_02', 'ss_00', 'cc_03', 'ha_02', 'cc_04', 'ha_03', 'cc_05',
            'cc_06', 'cc_07', 'ha_04', 'cc_08', 'ha_05', 'ss_01', 'cc_09',
            'c_00', 'o_00', 'os_00', 'c3_00', 'h1_00', 'h1_00', 'h1_00'
        ]

        self.run_atomtypeidentifier(xyz_string, expected_atomtypes,
                                    expected_equal_charges_list,
                                    expected_equiv_atoms)

    def test_atomtypeidentifier_hexatriene(self):

        xyz_string = """14
        hexatriene
        C     -4.25966    -1.52840     0.04900
        C     -2.93241    -1.66109    -0.01724
        H     -4.74560    -1.00165     0.87382
        H     -4.90399    -1.95396    -0.72224
        C     -2.03936    -1.12221     0.98343
        H     -2.49052    -2.19802    -0.86169
        C     -0.70820    -1.27308     0.92617
        H     -2.48902    -0.57469     1.81552
        H     -0.25517    -1.82532     0.10111
        C      0.18073    -0.73586     1.93157
        C      1.50990    -0.89760     1.88986
        H     -0.26288    -0.18121     2.75110
        H      2.00545    -1.44399     1.09891
        H      2.14185    -0.48639     2.66959
        """
        expected_atomtypes = {
            'C1': {
                'opls': 'opls_141',
                'gaff': 'c2',
                'uff': 'C'
            },
            'C2': {
                'opls': 'opls_XXX',
                'gaff': 'ce',
                'uff': 'C'
            },
            'H3': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'H4': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C5': {
                'opls': 'opls_XXX',
                'gaff': 'ce',
                'uff': 'C'
            },
            'H6': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C7': {
                'opls': 'opls_XXX',
                'gaff': 'cf',
                'uff': 'C'
            },
            'H8': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'H9': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C10': {
                'opls': 'opls_XXX',
                'gaff': 'cf',
                'uff': 'C'
            },
            'C11': {
                'opls': 'opls_141',
                'gaff': 'c2',
                'uff': 'C'
            },
            'H12': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'H13': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'H14': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            }
        }
        expected_equal_charges_list = [[1, 11], [2, 10], [3, 4, 13, 14], [5, 7],
                                       [6, 12], [8, 9]]
        expected_equiv_atoms = [
            'c2_00', 'ce_00', 'ha_00', 'ha_00', 'ce_01', 'ha_01', 'ce_01',
            'ha_02', 'ha_02', 'ce_00', 'c2_00', 'ha_01', 'ha_00', 'ha_00'
        ]

        self.run_atomtypeidentifier(xyz_string, expected_atomtypes,
                                    expected_equal_charges_list,
                                    expected_equiv_atoms)

    def test_atomtypeidentifier_water(self):

        xyz_string = """3
        water
        O      1.035010   -0.055790    0.042900
        H      2.002900   -0.067860    0.087740
        H      0.756160   -0.297060    0.939000
        """
        expected_atomtypes = {
            'O1': {
                'opls': 'opls_111',
                'gaff': 'ow',
                'uff': 'O'
            },
            'H2': {
                'opls': 'opls_112',
                'gaff': 'hw',
                'uff': 'H'
            },
            'H3': {
                'opls': 'opls_112',
                'gaff': 'hw',
                'uff': 'H'
            }
        }
        expected_equal_charges_list = [[2, 3]]
        expected_equiv_atoms = ['ow_00', 'hw_00', 'hw_00']

        self.run_atomtypeidentifier(xyz_string, expected_atomtypes,
                                    expected_equal_charges_list,
                                    expected_equiv_atoms)

    def test_atomtypeidentifier_butane(self):

        xyz_string = """14
        butane
        C     -1.962049    0.121016    0.000050
        H     -2.109257    0.754214    0.888026
        H     -2.109343    0.754296   -0.887853
        H     -2.751897   -0.643609    0.000052
        C     -0.568377   -0.514088   -0.000056
        H     -0.463006   -1.168632    0.880682
        H     -0.463086   -1.168544   -0.880871
        C      0.568377    0.514088   -0.000056
        H      0.463086    1.168542   -0.880872
        H      0.463006    1.168634    0.880681
        C      1.962049   -0.121016    0.000051
        H      2.109336   -0.754313   -0.887842
        H      2.751897    0.643609    0.000031
        H      2.109264   -0.754197    0.888037
        """
        expected_atomtypes = {
            'C1': {
                'opls': 'opls_135',
                'gaff': 'c3',
                'uff': 'C'
            },
            'H2': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'H3': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'H4': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'C5': {
                'opls': 'opls_135',
                'gaff': 'c3',
                'uff': 'C'
            },
            'H6': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'H7': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'C8': {
                'opls': 'opls_135',
                'gaff': 'c3',
                'uff': 'C'
            },
            'H9': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'H10': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'C11': {
                'opls': 'opls_135',
                'gaff': 'c3',
                'uff': 'C'
            },
            'H12': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'H13': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'H14': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            }
        }
        expected_equal_charges_list = [[1, 11], [2, 3, 4, 12, 13, 14], [5, 8],
                                       [6, 7, 9, 10]]
        expected_equiv_atoms = [
            'c3_00', 'hc_00', 'hc_00', 'hc_00', 'c3_01', 'hc_01', 'hc_01',
            'c3_01', 'hc_01', 'hc_01', 'c3_00', 'hc_00', 'hc_00', 'hc_00'
        ]

        self.run_atomtypeidentifier(xyz_string, expected_atomtypes,
                                    expected_equal_charges_list,
                                    expected_equiv_atoms)

    def test_atomtypeidentifier_noradrenaline(self):

        xyz_string = """23
        noradrenaline
        H       -1.8193999674   -2.5097458050   -0.5507950455
        C       -1.4355452194   -1.5284769930   -0.2951571832
        C       -2.3089191215   -0.4593458985   -0.1886635299
        O       -3.6690464443   -0.5422210733   -0.3841138405
        C       -1.8255516449    0.8097372093    0.1469046451
        C       -0.4687604354    0.9845502431    0.3694347370
        H       -0.1169650009    1.9731850398    0.6359161766
        C        0.4203521635   -0.0850019011    0.2563051158
        C        1.9048341672    0.1432961237    0.4528861030
        C       -0.0738036281   -1.3439593880   -0.0746469636
        H        0.6049466138   -2.1813113202   -0.1482676236
        O       -2.6753261092    1.8678164911    0.2657837532
        H       -3.9133744651   -1.4429007059   -0.6142762974
        H       -3.5694739457    1.5581900821    0.0765607930
        O        2.5513450005   -0.9759063254    1.0721777983
        H        2.0561900870    1.0385736809    1.0650086348
        C        2.6349951121    0.3620412906   -0.8693646417
        H        2.1445626765    1.1787081547   -1.4017612725
        N        4.0273077512    0.7299013711   -0.6192540306
        H        2.5065073570   -0.5476233229   -1.4727537911
        H        4.4869804996   -0.0238192548   -0.1223947057
        H        4.5220904955    0.8575433726   -1.4930596178
        H        2.1653934540   -1.1046663754    1.9436295799
        """
        expected_atomtypes = {
            'H1': {
                'uff': 'H',
                'opls': 'opls_146',
                'gaff': 'ha'
            },
            'C2': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'C3': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'O4': {
                'opls': 'opls_154',
                'gaff': 'oh',
                'uff': 'O'
            },
            'C5': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'C6': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H7': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C8': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'C9': {
                'opls': 'opls_135',
                'gaff': 'c3',
                'uff': 'C'
            },
            'C10': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H11': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'O12': {
                'opls': 'opls_154',
                'gaff': 'oh',
                'uff': 'O'
            },
            'H13': {
                'opls': 'opls_155',
                'gaff': 'ho',
                'uff': 'H'
            },
            'H14': {
                'opls': 'opls_155',
                'gaff': 'ho',
                'uff': 'H'
            },
            'O15': {
                'opls': 'opls_154',
                'gaff': 'oh',
                'uff': 'O'
            },
            'H16': {
                'opls': 'opls_140',
                'gaff': 'h1',
                'uff': 'H'
            },
            'C17': {
                'opls': 'opls_135',
                'gaff': 'c3',
                'uff': 'C'
            },
            'H18': {
                'opls': 'opls_140',
                'gaff': 'h1',
                'uff': 'H'
            },
            'N19': {
                'opls': 'opls_300',
                'gaff': 'n8',
                'uff': 'N'
            },
            'H20': {
                'opls': 'opls_140',
                'gaff': 'h1',
                'uff': 'H'
            },
            'H21': {
                'opls': 'opls_240',
                'gaff': 'hn',
                'uff': 'H'
            },
            'H22': {
                'opls': 'opls_240',
                'gaff': 'hn',
                'uff': 'H'
            },
            'H23': {
                'opls': 'opls_155',
                'gaff': 'ho',
                'uff': 'H'
            }
        }
        expected_equal_charges_list = [[18, 20], [21, 22]]
        expected_equiv_atoms = [
            'ha_00', 'ca_00', 'ca_01', 'oh_00', 'ca_02', 'ca_03', 'ha_01',
            'ca_04', 'c3_00', 'ca_05', 'ha_02', 'oh_01', 'ho_00', 'ho_01',
            'oh_02', 'h1_00', 'c3_01', 'h1_01', 'n8_00', 'h1_01', 'hn_00',
            'hn_00', 'ho_02'
        ]

        self.run_atomtypeidentifier(xyz_string, expected_atomtypes,
                                    expected_equal_charges_list,
                                    expected_equiv_atoms)

    def test_atomtypeidentifier_cp_cq(self):

        xyz_string = """32
        xyz
        H        1.1673441590   -0.0821866308    2.1271909485
        C        0.6586062417   -0.0005495445    1.1771168299
        C       -0.7218673637   -0.0512289398    1.1263803686
        H       -1.2922480539   -0.1550274045    2.0377238497
        C       -1.3682069297    0.0262234602   -0.0957706434
        H       -2.4464775443   -0.0143946079   -0.1443980069
        C       -0.6287644469    0.1611958503   -1.2558786645
        H       -1.1305981043    0.2434344334   -2.2095647452
        C        0.7609221227    0.2202912208   -1.2185463116
        C        1.4153505362    0.1283553145    0.0165906712
        C        2.8860272242    0.1286042338    0.1134558978
        C        3.5300742505    0.9749334630    1.0113047892
        H        2.9434445471    1.6501346714    1.6184808256
        C        4.9092856982    0.9743426865    1.1118846260
        H        5.3976456758    1.6427000944    1.8064445777
        C        5.6631328744    0.1210917057    0.3242614036
        H        6.7404109854    0.1211236916    0.4025216543
        C        5.0300710653   -0.7344290167   -0.5621549172
        H        5.6132956512   -1.4075405705   -1.1735783420
        C        3.6522044767   -0.7315182035   -0.6685288459
        H        3.1615867319   -1.4065656074   -1.3543855146
        C        1.5064144323    0.4140453044   -2.4750979986
        C        2.4571879071    1.4243664184   -2.5915601389
        H        2.6500109290    2.0733627293   -1.7498760743
        C        3.1424769108    1.6078187675   -3.7775713848
        H        3.8769180470    2.3962173809   -3.8555811912
        C        2.8902741230    0.7856231452   -4.8633058376
        H        3.4303414002    0.9270458491   -5.7879590051
        C        1.9408138716   -0.2166304037   -4.7596363743
        H        1.7385781106   -0.8599502363   -5.6039491168
        C        1.2497622241   -0.3983505582   -3.5756958481
        H        0.5171222467   -1.1889186966   -3.4924774810
        """
        expected_atomtypes = {
            'H1': {
                'uff': 'H',
                'opls': 'opls_146',
                'gaff': 'ha'
            },
            'C2': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'C3': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H4': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C5': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H6': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C7': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H8': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C9': {
                'opls': 'opls_521',
                'gaff': 'cp',
                'uff': 'C'
            },
            'C10': {
                'opls': 'opls_521',
                'gaff': 'cq',
                'uff': 'C'
            },
            'C11': {
                'opls': 'opls_521',
                'gaff': 'cq',
                'uff': 'C'
            },
            'C12': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H13': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C14': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H15': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C16': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H17': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C18': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H19': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C20': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H21': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C22': {
                'opls': 'opls_521',
                'gaff': 'cp',
                'uff': 'C'
            },
            'C23': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H24': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C25': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H26': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C27': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H28': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C29': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H30': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C31': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H32': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            }
        }
        expected_equal_charges_list = [[1, 8], [2, 7], [3, 5], [4, 6], [9, 10],
                                       [11, 22], [12, 20, 23, 31],
                                       [13, 21, 24, 32], [14, 18, 25, 29],
                                       [15, 19, 26, 30], [16, 27], [17, 28]]
        expected_equiv_atoms = [
            'ha_00', 'ca_00', 'ca_01', 'ha_01', 'ca_01', 'ha_01', 'ca_00',
            'ha_00', 'cp_00', 'cp_00', 'cp_01', 'ca_02', 'ha_02', 'ca_03',
            'ha_03', 'ca_04', 'ha_04', 'ca_03', 'ha_03', 'ca_02', 'ha_02',
            'cp_01', 'ca_02', 'ha_02', 'ca_03', 'ha_03', 'ca_04', 'ha_04',
            'ca_03', 'ha_03', 'ca_02', 'ha_02'
        ]

        self.run_atomtypeidentifier(xyz_string, expected_atomtypes,
                                    expected_equal_charges_list,
                                    expected_equiv_atoms)

    def test_atomtypeidentifier_ca(self):

        xyz_string = """21
        xyz
        C        0.8897639463   -0.1087654484    1.0322914605
        C       -0.4774462285   -0.0844278058    1.2524483984
        H       -0.8695823487   -0.1506983070    2.2554761434
        C       -1.3175925426    0.0260732885    0.1614716518
        H       -2.3864598096    0.0467778064    0.3121273511
        C       -0.8011433781    0.1108588749   -1.1281276916
        H       -1.4761982453    0.1966148230   -1.9665051811
        C        0.5610040855    0.0864636795   -1.3455290887
        H        0.9586235234    0.1524139661   -2.3470643742
        C        1.4287346959   -0.0241758459   -0.2627176332
        C        2.8681996475   -0.0687984828   -0.2514214218
        C        3.3807381806   -0.1859946729    1.0517183078
        C        4.7432730758   -0.2461446845    1.2930569073
        H        5.1149144347   -0.3361434012    2.3020262208
        C        5.6055464484   -0.1881828550    0.2153118007
        H        6.6712242092   -0.2335276601    0.3822918333
        C        5.1153861239   -0.0720503492   -1.0820807374
        H        5.8072005735   -0.0283577995   -1.9099898215
        C        3.7578613784   -0.0122850155   -1.3204820193
        H        3.3806175601    0.0777846220   -2.3279786741
        S        2.1227346696   -0.2425447327    2.2672865678
        """
        expected_atomtypes = {
            'C1': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'C2': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H3': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C4': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H5': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C6': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H7': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C8': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H9': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C10': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'C11': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'C12': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'C13': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H14': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C15': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H16': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C17': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H18': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C19': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H20': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'S21': {
                'opls': 'opls_SS',
                'gaff': 'ss',
                'uff': 'S'
            }
        }
        expected_equal_charges_list = [[1, 12], [2, 13], [3, 14], [4, 15],
                                       [5, 16], [6, 17], [7, 18], [8, 19],
                                       [9, 20], [10, 11]]
        expected_equiv_atoms = [
            'ca_00', 'ca_01', 'ha_00', 'ca_02', 'ha_01', 'ca_03', 'ha_02',
            'ca_04', 'ha_03', 'ca_05', 'ca_05', 'ca_00', 'ca_01', 'ha_00',
            'ca_02', 'ha_01', 'ca_03', 'ha_02', 'ca_04', 'ha_03', 'ss_00'
        ]

        self.run_atomtypeidentifier(xyz_string, expected_atomtypes,
                                    expected_equal_charges_list,
                                    expected_equiv_atoms)

    def test_atomtypeidentifier_cp(self):

        xyz_string = """26
        xyz
        H        1.2137168024    0.0902122683    2.1092215444
        C        0.7223202047    0.0602901529    1.1483294823
        C       -0.6556876326    0.1559960961    1.0892283978
        H       -1.2282384130    0.2414466534    2.0013243712
        C       -1.2991715930    0.1463279158   -0.1359619808
        H       -2.3757173300    0.2177414304   -0.1863201165
        C       -0.5549456446    0.0597353253   -1.3003777877
        H       -1.0535894891    0.0664767739   -2.2597537177
        C        0.8258192837   -0.0296262095   -1.2530354598
        C        1.4799708311   -0.0487288848   -0.0145285957
        C        2.9422509352   -0.1673090514    0.0112283721
        C        3.6281289974   -0.6685800669    1.1141541991
        H        3.0788677991   -1.0142107532    1.9770269584
        C        5.0080976742   -0.7536042189    1.1112411030
        H        5.5244643908   -1.1464895305    1.9749754296
        C        5.7254749775   -0.3398829685    0.0024339116
        H        6.8037732184   -0.4017679549   -0.0025553806
        C        5.0533641524    0.1406531780   -1.1086962387
        H        5.6100474692    0.4511733059   -1.9820691787
        C        3.6714253373    0.2236465154   -1.1190913196
        C        1.6531074181   -0.1056740071   -2.5034558188
        C        2.9217684286    0.7239737633   -2.3196913633
        H        1.0796392306    0.2500540338   -3.3611283153
        H        1.9300727528   -1.1488414730   -2.6923782047
        H        3.5505781838    0.6708003671   -3.2100741635
        H        2.6458720151    1.7728973393   -2.1640061282
        """
        expected_atomtypes = {
            'H1': {
                'uff': 'H',
                'opls': 'opls_146',
                'gaff': 'ha'
            },
            'C2': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'C3': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H4': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C5': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H6': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C7': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H8': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C9': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'C10': {
                'opls': 'opls_521',
                'gaff': 'cp',
                'uff': 'C'
            },
            'C11': {
                'opls': 'opls_521',
                'gaff': 'cp',
                'uff': 'C'
            },
            'C12': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H13': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C14': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H15': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C16': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H17': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C18': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'H19': {
                'opls': 'opls_146',
                'gaff': 'ha',
                'uff': 'H'
            },
            'C20': {
                'opls': 'opls_145',
                'gaff': 'ca',
                'uff': 'C'
            },
            'C21': {
                'opls': 'opls_135',
                'gaff': 'c3',
                'uff': 'C'
            },
            'C22': {
                'opls': 'opls_135',
                'gaff': 'c3',
                'uff': 'C'
            },
            'H23': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'H24': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'H25': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            },
            'H26': {
                'opls': 'opls_140',
                'gaff': 'hc',
                'uff': 'H'
            }
        }
        expected_equal_charges_list = [[1, 13], [2, 12], [3, 14], [4, 15],
                                       [5, 16], [6, 17], [7, 18], [8, 19],
                                       [9, 20], [10, 11], [21, 22],
                                       [23, 24, 25, 26]]
        expected_equiv_atoms = [
            'ha_00', 'ca_00', 'ca_01', 'ha_01', 'ca_02', 'ha_02', 'ca_03',
            'ha_03', 'ca_04', 'cp_00', 'cp_00', 'ca_00', 'ha_00', 'ca_01',
            'ha_01', 'ca_02', 'ha_02', 'ca_03', 'ha_03', 'ca_04', 'c3_00',
            'c3_00', 'hc_00', 'hc_00', 'hc_00', 'hc_00'
        ]

        self.run_atomtypeidentifier(xyz_string, expected_atomtypes,
                                    expected_equal_charges_list,
                                    expected_equiv_atoms)
