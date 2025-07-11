#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from time import sleep
from urllib.request import urlopen
from urllib.error import HTTPError, URLError
import json
import re


def get_data_from_name(mol_name):
    """
    Accesses the PubChem database to retrieve data for a given molecule name.

    Note: This is seperate to the previous function because other data could
    also be retrieved. Some examples are charge, volume and 3D properties.
    At the time of implementation, relevant properties are unknown,
    but they are easily accessed by changing the url.
    A good idea would then be to make the output a dictionary instead.

    :param mol_name:
        Molecule name string.

    :return:
        Accessed data
        Current data includes:
        SMILES-string.
        Title-string.
        CID-string.
    """

    # sanity checks
    if not re.match(r'^[a-zA-Z0-9 \-(),\.]+$', mol_name):
        raise ValueError("Compound name contains invalid characters.")
    if not mol_name:
        raise ValueError("Compound name cannot be empty.")
    if len(mol_name) > 100:
        raise ValueError("Compound name is too long.")
    mol_name = mol_name.strip()
    mol_name = mol_name.replace(' ', '%20')
    mol_name = mol_name.replace(',', '_')

    url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compoundname}/property/title,SMILES/JSON'.format(
        compoundname=mol_name.lower())

    try:
        with urlopen(url) as page:
            data = json.loads(page.read())
            smiles_str = data['PropertyTable']['Properties'][0]['SMILES']
            title = data['PropertyTable']['Properties'][0]['Title']
            cid = data['PropertyTable']['Properties'][0]['CID']

            # avoid frequent request
            sleep(0.2)

            print(f"Reading {mol_name} from PubChem...\n")
            print("Reference: " + get_pubchem_reference() + "\n")
            print("Please double-check the compound since names may" +
                  " refer to more than one record.\n")

            return smiles_str, title, cid

    except HTTPError as e:
        print(f"HTTP Error: {e.code}. The compound may not exist.")

    except URLError as e:
        print(f"URL Error: {e.reason}.")


def get_all_conformer_IDs(mol_name):
    """
    Gets all conformer IDs for the PubChem database for a compound

    :param mol_name:
        The molecule name string.

    :return conformerID_list:
        List of conformer IDs
    """

    # sanity checks
    if not re.match(r'^[a-zA-Z0-9 \-(),\.]+$', mol_name):
        raise ValueError("Compound name contains invalid characters.")
    if not mol_name:
        raise ValueError("Compound name cannot be empty.")
    if len(mol_name) > 100:
        raise ValueError("Compound name is too long.")
    mol_name = mol_name.strip()
    mol_name = mol_name.replace(' ', '%20')
    mol_name = mol_name.replace(',', '_')

    url_conID = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compoundname}/conformers/JSON'.format(
        compoundname=mol_name.lower())

    try:
        with urlopen(url_conID) as page:
            data = json.loads(page.read())
            conformerID_list = data['InformationList']['Information'][0][
                'ConformerID']

            # avoid frequent request
            sleep(0.2)

            print(f"Reading {mol_name} from PubChem...\n")
            print("Reference: " + get_pubchem_reference() + "\n")
            print("Please double-check the compound since names may" +
                  " refer to more than one record.\n")

            return conformerID_list

    except HTTPError as e:
        print(f"HTTP Error: {e.code}. The compound may not have conformers, or it may not exist.")

    except URLError as e:
        print(f"URL Error: {e.reason}")


def get_conformer_data(conformer_ID):
    """
    Gets conformer coordinates from PubChem based on conformer ID.

    :param conformer_ID:
        Conformer ID for molecules compound

    :return xyz:
        The xyz-string of the molecule
    """

    # sanity checks

    if not re.match(r'^[a-zA-Z0-9 \-(),\.]+$', conformer_ID):
        raise ValueError("Conformer ID contains invalid characters.")
    if not conformer_ID:
        raise ValueError("Conformer ID cannot be empty.")
    if len(conformer_ID) > 100:
        raise ValueError("Conformer ID is too long.")
    conformer_ID = conformer_ID.strip()
    conformer_ID = conformer_ID.replace(' ', '%20')
    conformer_ID = conformer_ID.replace(',', '_')

    url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/conformers/{conformerID}/JSON'.format(
        conformerID=conformer_ID)

    try:
        with urlopen(url) as page:
            data = json.loads(page.read())
            elements = data['PC_Compounds'][0]['atoms']['element']
            xcoords = data['PC_Compounds'][0]['coords'][0]['conformers'][0]['x']
            ycoords = data['PC_Compounds'][0]['coords'][0]['conformers'][0]['y']
            zcoords = data['PC_Compounds'][0]['coords'][0]['conformers'][0]['z']

            elements = index_to_element(elements)

            xyz = '{numberofatoms}\n'.format(numberofatoms=len(elements))
            for j in range(len(elements)):
                newline = '\n{el}   {x}    {y}     {z}'.format(el=elements[j],
                                                               x=xcoords[j],
                                                               y=ycoords[j],
                                                               z=zcoords[j])
                xyz = xyz + newline

            # avoid frequent request
            sleep(0.2)

            print(f"Reading conformer ID {conformer_ID} from PubChem...")

            return xyz

    except HTTPError as e:
        print(f"HTTP Error: {e.code}. The compound may not have conformers, or it may not exist.")

    except URLError as e:
        print(f"URL Error: {e.reason}")


def index_to_element(indices):
    """Sends the numerical position of an element on the periodic table to their abbreviations

    :param indicies
        The numerical positions of a list of elements on the periodic table

    :return elements:
        The corresponding list of element abbreviations

    """

    periodic_table = [
        "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al",
        "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
        "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb",
        "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
        "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
        "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta",
        "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
        "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
        "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
        "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
    ]

    elements = []
    for i in indices:
        elements.append(periodic_table[i - 1])
    return elements


def get_pubchem_reference():
    """
    Gets PubChem reference.

    :return:
        The reference.
    """

    pubchem_ref = 'S. Kim, J. Chen, T. Cheng, A. Gindulyte, J. He, S. He,'
    pubchem_ref += ' Q. Li, B. A. Shoemaker, P. A. Thiessen, B. Yu,'
    pubchem_ref += ' L. Zaslavsky, J. Zhang, E. E. Bolton,'
    pubchem_ref += ' Nucleic Acids Res., 2025, 53, D1516-D1525.'

    return pubchem_ref
