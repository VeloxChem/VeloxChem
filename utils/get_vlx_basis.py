import argparse
import json
import sys
import hashlib

import basis_set_exchange as bse

from pathlib import Path

def known_aliases_for_basis_sets():
    """
    Maps from common basis name to basis filename.

    :return:
        A dictionary that maps from common basis name to basis filename.
    """

    return {
        '6-31G*': '6-31G_D_',
        '6-31G**': '6-31G_D,P_',
        '6-31+G*': '6-31+G_D_',
        '6-31+G**': '6-31+G_D,P_',
        '6-31++G*': '6-31++G_D_',
        '6-31++G**': '6-31++G_D,P_',
        '6-311G*': '6-311G_D_',
        '6-311G**': '6-311G_D,P_',
        '6-311+G*': '6-311+G_D_',
        '6-311+G**': '6-311+G_D,P_',
        '6-311++G*': '6-311++G_D_',
        '6-311++G**': '6-311++G_D,P_',
        '6-31G(2DF,P)': '6-31G_2DF,P_',
        '6-31G(3DF,3PD)': '6-31G_3DF,3PD_',
        '6-311G(2DF,2PD)': '6-311G_2DF,2PD_',
        '6-311+G(2D,P)': '6-311+G_2D,P_',
        '6-311++G(2D,2P)': '6-311++G_2D,2P_',
        '6-311++G(3DF,3PD)': '6-311++G_3DF,3PD_',
        'DEF2-SV(P)': 'DEF2-SV_P_',
    }
    
def get_element_labels():

    """
    Gets chemical element labels for elements 1-86.

    :return:
        A list of chemical element labels.
    """
    
    return [
        '0', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg',
        'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr',
        'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
        'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
        'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
        'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf',
        'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po',
        'At', 'Rn'
    ]


def basis_name_to_file(name):
    """
    Determines basis set file name from common name.

    :param name:
        Common name of the basis set.

    :return:
        The file name of the basis set.
    """

    known_aliases = known_aliases_for_basis_sets()

    if name in known_aliases:
        return known_aliases[name]
    else:
        return name

def get_atom_basis(identifier,
                   gto_data,
                   ecp_data,
                   ecp_name):
    """
    Gets string for atom basis set.

    :param identifier:
        The identifier of chemical element.
    :param gto_data:
        The gto data of atomic basis.
    :param ecp_data:
        The ecp data of atomic basis.
    :param ecp_name:
        The name of ecp basis associated with gto basis.

    :return:
        The string with atomic basis group.
    """
    
    elem_labels = get_element_labels()
    
    ang_mom_labels = ['S', 'P', 'D', 'F', 'G', 'H', 'I']
    
    atom_basis_str = f'\n@ATOMBASIS {elem_labels[int(identifier)].upper()}\n'
    
    for gto in gto_data:
        # get basis function data
        ang_mom = gto['angular_momentum'][0]
        pexps = gto['exponents']
        pnorms = gto['coefficients'][0]
        # basis function header
        atom_basis_str += f'{ang_mom_labels[ang_mom]} {len(pexps)}\n'
        for pexp, pnorm in zip(pexps, pnorms):
            atom_basis_str += f'{float(pexp):18.12e} {float(pnorm):20.12e}\n'
        
    if ecp_data is not None:
        atom_basis_str += f'ECP {ecp_name.upper()}\n'
    atom_basis_str += '@END\n'
    
    return atom_basis_str
    
def get_atom_ecp(identifier,
                 ecp_elec,
                 ecp_data):
    """
    Gets string for atom ecp.

    :param identifier:
        The identifier of chemical element.
    :param ecp_elec:
        The number of core electrons.
    :param ecp_data:
        The ecp data of atomic basis.

    :return:
        The string with atomic ecp group.
    """
    
    elem_labels = get_element_labels()
    
    ang_mom_labels = ['S', 'P', 'D', 'F', 'G', 'H', 'I']
    
    atom_ecp_str = f'\n@ATOMECP {elem_labels[int(identifier)].upper()}\n'
    atom_ecp_str = f'NELEC {ecp_elec}\n'
    
    # local potential part
    
    ang_mom = ecp_data[0]['angular_momentum'][0]
    pexps = ecp_data[0]['gaussian_exponents']
    pnorms = ecp_data[0]['coefficients'][0]
    rfacts = ecp_data[0]['r_exponents']
    atom_ecp_str += f'UL {ang_mom_labels[ang_mom]} {len(pexps)}\n'
    for rfact, pexp, pnorm in zip(rfacts, pexps, pnorms):
        atom_ecp_str += f' {rfact} {float(pexp):18.12e} {float(pnorm):20.12e}\n'
    
    for ecp in ecp_data[1:]:
        ang_mom = ecp['angular_momentum'][0]
        pexps = ecp['gaussian_exponents']
        pnorms = ecp['coefficients'][0]
        rfacts = ecp['r_exponents']
        atom_ecp_str += f'UP {ang_mom_labels[ang_mom]} {len(pexps)}\n'
        for rfact, pexp, pnorm in zip(rfacts, pexps, pnorms):
            atom_ecp_str += f' {rfact} {float(pexp):18.12e} {float(pnorm):20.12e}\n'
        
    atom_ecp_str += '@END\n'
    
    return atom_ecp_str
    
def get_vlx_basis_string(basis_name,
                         basis_type,
                         ecp_name):
    """
    Gets string for VeloxChem basis set.

    :param basis_name:
        The name of basis set.
    :param basis_type:
        The type of basis set i.e. gto or ecp.
    :param ecp_name:
        The name of ecp basis associated with gto basis.

    :return:
        The string for VeloxChem basis set.
    """
    
    basis_label = basis_name

    # geometric augmentation

    geom_aug = 0
    if basis_name.upper().startswith('DAUG-'):
        basis_label = basis_name[1:]
        geom_aug = 1
    elif basis_name.upper().startswith('TAUG-'):
        basis_label = basis_name[1:]
        geom_aug = 2
    elif basis_name.upper().startswith('QAUG-'):
        basis_label = basis_name[1:]
        geom_aug = 3
        
    # convert basis set name to basis label
    if basis_name.upper() == 'SADLEJ-PVTZ':
        basis_label = 'SADLEJ PVTZ'
    
    # fetch basis set dictionary
    basis = bse.get_basis(basis_name,
                          header=False,
                          uncontract_general=True,
                          uncontract_spdf=True)
    
    # add geometrical augmentation
    if geom_aug > 0:
        basis = bse.manip.geometric_augmentation(basis, geom_aug)
        
    # write basis set string
    
    vlx_basis_str = None
    if basis_type == 'gto':
        vlx_basis_str = f'@BASIS_SET {basis_name}\n'
    if basis_type == 'ecp':
        vlx_basis_str = f'@ECP_SET {ecp_name}\n'
    
    basis_data = basis['elements']
    for elem, ao_basis in basis_data.items():
        if int(elem) > 86:
            continue
        gto_data = None
        if 'electron_shells' in ao_basis:
            gto_data = ao_basis['electron_shells']
        ecp_data = None
        if 'ecp_potentials' in ao_basis:
            ecp_data = ao_basis['ecp_potentials']
        ecp_elec = None
        if 'ecp_electrons' in ao_basis:
            ecp_elec = ao_basis['ecp_electrons']
        if basis_type == 'gto':
            vlx_basis_str += get_atom_basis(elem, gto_data, ecp_data, ecp_name)
        if basis_type == 'ecp' and ecp_data is not None:
            vlx_basis_str += get_atom_ecp(elem, ecp_elec, ecp_data)
            
    md5sum = hashlib.md5(vlx_basis_str.encode('utf-8')).hexdigest()
    
    return vlx_basis_str + md5sum
     
if __name__ == '__main__':

    # Process command line arguments

    usage = f'python {Path(__file__).name} '
    usage += '--basis <name_of_basis_set> '
    usage += '[--type <gto or ecp>] '
    usage += '[--ecp <name_of_ecp_set>] '
    
    parser = argparse.ArgumentParser(usage=usage)

    parser.add_argument('--basis',
                        dest='basis_name',
                        type=str,
                        default=None,
                        help='name of basis set')
                        
    parser.add_argument('--type',
                        dest='basis_type',
                        type=str,
                        default='gto',
                        help='type of basis set i.e. gto or ecp')
                        
    parser.add_argument('--ecp',
                        dest='ecp_name',
                        type=str,
                        default=None,
                        help='name of ecp set')

    args = parser.parse_args()

    if args.basis_name is None:
        parser.print_help()
        sys.exit(0)

    basis_name = args.basis_name.upper()
    basis_type = args.basis_type.lower()
    ecp_name = args.ecp_name
        
    vlx_basis_str = get_vlx_basis_string(basis_name,
                                         basis_type,
                                         ecp_name)
    
    output_name = basis_name_to_file(basis_name)
    with Path(output_name).open('w') as f_out:
        f_out.write(vlx_basis_str)
