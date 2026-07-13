"""QM force-field builder test (RESP charges computed on the fly).

This is the one EVB test that runs QM: RESP charges are computed for the small
reactant / product molecules. Reparameterisation is off by default and the
ethanol -> ethene + water reaction has only known GAFF parameters, so no Hessian
is computed.
"""

import pytest

from veloxchem.molecule import Molecule
from veloxchem.mmforcefieldgenerator import MMForceFieldGenerator
from veloxchem.reaffbuilder import ReactionForceFieldBuilder

from veloxchem.outputstream import OutputStream

from test_evb_helper import evb_compare_dict, evb_data_dir

pytestmark = [pytest.mark.timeconsuming]


class TestReactionForceFieldBuilder:

    def test_forcefield_builder(self):
        data_dir = evb_data_dir()
        # Build reactant and product force fields from unordered xyz inputs and
        # compare bonded parameters against the committed references.
        ffbuilder = ReactionForceFieldBuilder(ostream=OutputStream(None))
        ffbuilder.water_model = 'cspce'

        reactant = Molecule.read_xyz_file(str(data_dir / 'evb_ethanol.xyz'))
        product = [
            Molecule.read_xyz_file(str(data_dir / 'evb_ethene_H2O.xyz')),
        ]

        (reactant, product, formed_bonds, broken_bonds, reactants, products,
         mapping) = ffbuilder.build_force_fields(
             reactant=reactant,
             product=product,
         )

        reactant_ref = MMForceFieldGenerator.load_forcefield_from_json_file(
            str(data_dir / 'evb_ethanol_ff_data.json'),
            ostream=OutputStream(None))
        product_ref = MMForceFieldGenerator.load_forcefield_from_json_file(
            str(data_dir / 'evb_ethene_H2O_ff_data.json'),
            ostream=OutputStream(None))

        evb_compare_dict(reactant.bonds, reactant_ref.bonds)
        evb_compare_dict(reactant.angles, reactant_ref.angles)
        evb_compare_dict(reactant.dihedrals, reactant_ref.dihedrals)
        evb_compare_dict(reactant.impropers, reactant_ref.impropers)

        evb_compare_dict(product.bonds, product_ref.bonds)
        evb_compare_dict(product.angles, product_ref.angles)
        evb_compare_dict(product.dihedrals, product_ref.dihedrals)
        evb_compare_dict(product.impropers, product_ref.impropers)
