from pathlib import Path
import tempfile
from veloxchem import InputParser, Molecule, MolecularBasis
from veloxchem.molecularbasis import _basis_file_to_name
from veloxchem.environ import get_basis_path
from get_vlx_basis import get_vlx_basis_string


def test_get_vlx_basis(capsys):

    basis_path = Path(get_basis_path())
    basis_files = sorted((x for x in basis_path.iterdir() if x.is_file()))

    for x in basis_files:
        if x.stem.upper()[:6] in ['ANO-L-', 'ANO-S-']:
            continue
        if x.stem.upper() in ['STO-3G-OLD', 'AO-START-GUESS', 'CC-PVDZ-OPTGC']:
            continue

        with capsys.disabled():
            print(f'{x.stem.ljust(20)}', end='')

        basis_dict = InputParser(str(x)).input_dict

        elements = []
        for key in basis_dict:
            if key.startswith('atombasis_'):
                elements.append(key.replace('atombasis_', '').capitalize())

        xyzstr = ''
        for i, elem in enumerate(elements):
            xyzstr += f'{elem}  {0.0}  {0.0}  {i:.1f}\n'
        mol = Molecule.read_str(xyzstr)

        vlx_basis = MolecularBasis.read(mol, x.name, x.parent)

        bs_name = _basis_file_to_name(x.name)

        with tempfile.TemporaryDirectory() as temp_dir:
            bs_file = Path(temp_dir, bs_name.upper())
            vlx_basis_str = get_vlx_basis_string(bs_name,
                                                 str(bs_file),
                                                 use_gc_and_sp=False,
                                                 optimize_general=False)
            with Path(str(bs_file)).open('w') as f_out:
                f_out.write(vlx_basis_str)
            bse_basis = MolecularBasis.read(mol, bs_file.name, bs_file.parent)

        vlx_data = []
        for atombasis in vlx_basis:
            for basisfunc in atombasis:
                vlx_data.append((
                    atombasis.get_elemental_id(),
                    basisfunc.angular_momentum,
                    basisfunc.exponents,
                    basisfunc.normalization_factors,
                ))

        bse_data = []
        for atombasis in bse_basis:
            for basisfunc in atombasis:
                bse_data.append((
                    atombasis.get_elemental_id(),
                    basisfunc.angular_momentum,
                    basisfunc.exponents,
                    basisfunc.normalization_factors,
                ))

        assert sorted(vlx_data) == sorted(bse_data)

        with capsys.disabled():
            print('  [PASSED]')
