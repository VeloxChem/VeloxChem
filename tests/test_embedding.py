import numpy as np
import os
import pytest
from veloxchem.veloxchemlib import (Molecule, MolecularBasis, PolarizableEmbedding, ScfRestrictedDriver,
                                    LinearResponseSolver)


class TestPolarizableEmbedding:
    @staticmethod
    def get_molecule_and_basis(name):
        if name == 'acrolein':
            xyz_string = """8
                C3OH4
                C             -0.145335   -0.546770    0.000607
                C              1.274009   -0.912471   -0.000167
                C              1.630116   -2.207690   -0.000132
                O             -0.560104    0.608977    0.000534
                H             -0.871904   -1.386459    0.001253
                H              2.004448   -0.101417   -0.000710
                H              0.879028   -3.000685    0.000484
                H              2.675323   -2.516779   -0.000673
                """
            basis_label = 'sto-3g'
            mol = Molecule.read_xyz_string(xyz_string)
            bas = MolecularBasis.read(mol, basis_label, ostream=None)
            return mol, bas

    def run_scf_with_pe(self, name):
        mol, bas = self.get_molecule_and_basis(name)
        options = f'{os.path.dirname(__file__)}/data/{name}.json'
        scf_drv = ScfRestrictedDriver()
        scf_drv.conv_thresh = 1e-11
        scf_drv.embedding_options = {'settings': {'embedding_method': 'PE'},
                                     'inputs': {'json_file': options}
                                     }
        return scf_drv.compute(mol, bas)

    def run_lrs_with_pe(self, scf_results, name, freqs):
        mol, bas = self.get_molecule_and_basis(name)
        lrsolver = LinearResponseSolver()
        lrsolver.conv_thresh = 1e-10
        options = f'{os.path.dirname(__file__)}/data/{name}.json'
        lrsolver.embedding_options = {'settings': {'embedding_method': 'PE'},
                                      'inputs': {'json_file': options}
                                      }
        freqs_str = [str(x) for x in freqs]
        lrsolver.update_settings(
            {
                'frequencies': ','.join(freqs_str)
            })
        return lrsolver.compute(molecule=mol, basis=bas, scf_tensors=scf_results)

    def test_scf_with_pe(self):
        scf_results = self.run_scf_with_pe(name='acrolein')
        assert scf_results['scf_energy'] == pytest.approx(-188.314434428374, rel=1e-10)

    def test_lrs_with_pe(self):
        scf_results = self.run_scf_with_pe(name='acrolein')
        ref_freqs = [0.0, 0.04]
        lrs_results = self.run_lrs_with_pe(name='acrolein', scf_results=scf_results, freqs=ref_freqs)
        # scf 1e-11 lrs 1e-10
        reference_data = """
            2.097957814289E+01 -8.652608839280E+00 -7.686399537652E-03
            -8.652608839280E+00 3.353258603117E+01 -1.147335054308E-03
            -7.686399537652E-03 -1.147335054308E-03 4.586549080720E+00
            2.111014439306E+01 -8.811667762577E+00 -7.720196012987E-03
            -8.811667762577E+00 3.385544639017E+01 -1.118865578820E-03
            -7.720196012987E-03 -1.118865578820E-03 4.601670440902E+00
            """
        ref_prop = np.array([float(x) for x in reference_data.split()])
        prop = np.array([
            -lrs_results['response_functions'][(a, b, w)]
            for w in ref_freqs
            for a in 'xyz'
            for b in 'xyz'
        ])
        assert np.max(np.abs((prop - ref_prop) / prop)) < 6.0e-8
