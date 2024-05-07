#
#                              VELOXCHEM
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
#
#  SPDX-License-Identifier: LGPL-3.0-or-later
#
#  This file is part of VeloxChem.
#
#  VeloxChem is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Lesser General Public License as published by the Free
#  Software Foundation, either version 3 of the License, or (at your option)
#  any later version.
#
#  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

# C++ classes

from .veloxchemlib import ChemicalElement
from .veloxchemlib import BasisFunction
from .veloxchemlib import AtomBasis
from .veloxchemlib import GtoBlock
from .veloxchemlib import OverlapDriver

# C++ functions
from .veloxchemlib import bohr_in_angstroms
from .veloxchemlib import hartree_in_ev
from .veloxchemlib import cantor_index
from .veloxchemlib import cantor_pair
from .veloxchemlib import to_spherical_components
from .veloxchemlib import to_cartesian_components
from .veloxchemlib import angular_component_to_str

# Python functions
from .errorhandler import assert_msg_critical

# Python classes
from .inputparser import InputParser
from .outputstream import OutputStream
from .molecule import Molecule
from .molecularbasis import MolecularBasis
from .submatrix import SubMatrix
from .matrix import Matrix
from .sadguessdriver import SadGuessDriver
from .scfrestdriver import ScfRestrictedDriver

# Environment variable: basis set path
from .environment import set_vlxbasispath

set_vlxbasispath()

__version__ = "0.0"
