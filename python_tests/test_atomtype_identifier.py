from veloxchem.atomtypeidentifier import AtomTypeIdentifier
from veloxchem.molecule import Molecule


class TestAtomTypeIdentifier:

    def run_atomtypeidentifier(self, xyz_string, expected_atomtypes):

        atomtypeidentifier = AtomTypeIdentifier()
        atomtypeidentifier.ostream.mute()

        molecule = Molecule.read_xyz_string(xyz_string)
        atomtypes = atomtypeidentifier.generate_gaff_atomtypes(molecule)

        assert atomtypes == expected_atomtypes

    def test_atomtypeidentifier_1(self):

        mol = """11
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
        expected_atomtypes = [
            'br', 'o', 'n2', 'c3', 'c2', 'c2', 'c1', 'h1', 'h1', 'ha', 'h4'
        ]

        self.run_atomtypeidentifier(mol, expected_atomtypes)

    def test_atomtypeidentifier_2(self):

        mol = """16
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
        expected_atomtypes = [
            'cx', 'cx', 'cx', 'cx', 'cx', 'cx', 'hc', 'hc', 'hc', 'hc', 'hc',
            'hc', 'hc', 'hc', 'hc', 'hc'
        ]

        self.run_atomtypeidentifier(mol, expected_atomtypes)

    def test_atomtypeidentifier_3(self):

        mol = """8
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
        expected_atomtypes = ['cl', 'ss', 's', 'cs', 'c3', 'h1', 'h1', 'h1']

        self.run_atomtypeidentifier(mol, expected_atomtypes)

    def test_atomtypeidentifier_4(self):

        mol = """8
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
        expected_atomtypes = ['c3', 'c3', 'cl', 'cl', 'cl', 'cl', 'cl', 'cl']

        self.run_atomtypeidentifier(mol, expected_atomtypes)

    def test_atomtypeidentifier_5(self):

        mol = """9
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
        expected_atomtypes = ['c', 'c', 'c3', 'o', 'o', 'h4', 'h4', 'hc', 'hc']

        self.run_atomtypeidentifier(mol, expected_atomtypes)

    def test_atomtypeidentifier_6(self):

        mol = """16
        1-methylcyclopentene
        C -0.756198 0.073934 -0.002520
        H -0.343604 2.171486 -0.066463
        C 0.023743 1.153783 -0.022376
        H -2.647844 -0.452806 0.859762
        H -2.606606 -0.553100 -0.889410
        H -2.679983 1.029419 -0.103430
        C -2.248786 0.030965 -0.036722
        H -0.043702 -1.617124 1.120770
        H -0.250301 -1.964412 -0.582053
        C 0.069948 -1.188503 0.117983
        H 1.780513 -0.916913 -1.185058
        H 2.249118 -1.218041 0.478619
        C 1.515599 -0.711239 -0.147209
        H 1.879408 1.089875 1.069321
        H 2.101884 1.357766 -0.645531
        C 1.489213 0.820035 0.081422
        """
        expected_atomtypes = [
            'c2', 'ha', 'c2', 'hc', 'hc', 'hc', 'c3', 'hc', 'hc', 'c5', 'hc',
            'hc', 'c5', 'hc', 'hc', 'c5'
        ]

        self.run_atomtypeidentifier(mol, expected_atomtypes)

    def test_atomtypeidentifier_7(self):

        mol = """10
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
        expected_atomtypes = [
            'br', 'o', 'n2', 'c3', 'c3', 'c1', 'h1', 'h1', 'h1', 'h1'
        ]

        self.run_atomtypeidentifier(mol, expected_atomtypes)

    def test_atomtypeidentifier_8(self):

        mol = """12
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
        expected_atomtypes = [
            'nb', 'ca', 'ca', 'ca', 'ca', 'ca', 'oh', 'ha', 'ha', 'ha', 'h4',
            'ho'
        ]

        self.run_atomtypeidentifier(mol, expected_atomtypes)

    def test_atomtypeidentifier_9(self):

        mol = """19
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
        expected_atomtypes = [
            'ca', 'ca', 'ca', 'ca', 'cc', 'ca', 'cd', 'ca', 'c3', 'na', 'ha',
            'ha', 'ha', 'ha', 'h4', 'hc', 'hc', 'hc', 'hn'
        ]

        self.run_atomtypeidentifier(mol, expected_atomtypes)

    def test_atomtypeidentifier_10(self):

        mol = """22
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
        expected_atomtypes = [
            'ca', 'ca', 'ca', 'ca', 'ca', 'ca', 'c3', 'c3', 'c3', 'oh', 'ha',
            'ha', 'ha', 'ha', 'hc', 'hc', 'hc', 'hc', 'hc', 'hc', 'hc', 'ho'
        ]

        self.run_atomtypeidentifier(mol, expected_atomtypes)

    def test_atomtypeidentifier_11(self):

        mol = """35
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
        expected_atomtypes = [
            'ca', 'ca', 'h4', 'nb', 'ca', 'h4', 'ca', 'ha', 'ca', 'cc', 'ha',
            'cd', 'h4', 'na', 'cc', 'ss', 'cd', 'ha', 'cd', 'ha', 'cc', 'cc',
            'cd', 'ha', 'cd', 'ha', 'ss', 'cc', 'c', 'o', 'os', 'c3', 'h1',
            'h1', 'h1'
        ]

        self.run_atomtypeidentifier(mol, expected_atomtypes)
