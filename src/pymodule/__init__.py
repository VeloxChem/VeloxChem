#
#                           VELOXCHEM 1.0-RC2
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
#  Contact: https://veloxchem.org/contact
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
from .veloxchemlib import AtomBasis
from .veloxchemlib import BasisFunction
from .veloxchemlib import OverlapIntegralsDriver
from .veloxchemlib import KineticEnergyIntegralsDriver
from .veloxchemlib import NuclearPotentialIntegralsDriver
from .veloxchemlib import ElectricDipoleIntegralsDriver
from .veloxchemlib import LinearMomentumIntegralsDriver
from .veloxchemlib import AngularMomentumIntegralsDriver
from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import GridDriver
from .veloxchemlib import DensityGridDriver
from .veloxchemlib import SADGuessDriver
from .veloxchemlib import DispersionModel
from .veloxchemlib import DenseMatrix
from .veloxchemlib import TwoIndexes
from .veloxchemlib import MOIntsBatch
from .veloxchemlib import ExcitationVector

# C++ functions
from .veloxchemlib import available_functionals
from .veloxchemlib import mpi_master
from .veloxchemlib import mpi_initialized
from .veloxchemlib import ao_matrix_to_veloxchem
from .veloxchemlib import ao_matrix_to_dalton
from .veloxchemlib import get_basis_function_indices_for_atom
from .veloxchemlib import bohr_in_angstroms
from .veloxchemlib import hartree_in_ev
from .veloxchemlib import hartree_in_kcalpermol
from .veloxchemlib import hartree_in_wavenumbers
from .veloxchemlib import dipole_in_debye
from .veloxchemlib import rotatory_strength_in_cgs
from .veloxchemlib import extinction_coefficient_from_beta
from .veloxchemlib import fine_structure_constant
from .veloxchemlib import mathconst_pi

# C++ enums
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .veloxchemlib import szblock
from .veloxchemlib import ericut
from .veloxchemlib import molorb
from .veloxchemlib import moints

# Python classes
from .inputparser import InputParser
from .outputstream import OutputStream
from .molecule import Molecule
from .molecularbasis import MolecularBasis
from .aodensitymatrix import AODensityMatrix
from .molecularorbitals import MolecularOrbitals
from .aofockmatrix import AOFockMatrix
from .scfrestdriver import ScfRestrictedDriver
from .scfunrestdriver import ScfUnrestrictedDriver
from .scfrestopendriver import ScfRestrictedOpenDriver
from .scfgradientdriver import ScfGradientDriver
from .xtbdriver import XTBDriver
from .xtbgradientdriver import XTBGradientDriver
from .optimizationdriver import OptimizationDriver
from .mointsdriver import MOIntegralsDriver
from .mp2driver import Mp2Driver
from .cubicgrid import CubicGrid
from .visualizationdriver import VisualizationDriver
from .excitondriver import ExcitonModelDriver
from .rspdriver import ResponseDriver
from .tdaexcidriver import TDAExciDriver
from .blockdavidson import BlockDavidsonSolver
from .lreigensolver import LinearResponseEigenSolver
from .lrsolver import LinearResponseSolver
from .cppsolver import ComplexResponse
from .c6solver import C6Solver
from .quadraticresponsedriver import QuadraticResponseDriver
from .cubicresponsedriver import CubicResponseDriver
from .shgdriver import SHGDriver
from .tpafulldriver import TPAFullDriver
from .tpareddriver import TPAReducedDriver
from .respchargesdriver import RespChargesDriver
from .rspproperty import ResponseProperty
from .rsplinabscross import LinearAbsorptionCrossSection
from .rspcdspec import CircularDichroismSpectrum
from .rsppolarizability import Polarizability
from .rspabsorption import Absorption
from .rspc6 import C6
from .rspshg import SHG
from .rsptpa import TPA
from .rspcustomproperty import CustomProperty
from .mpitask import MpiTask
from .subcommunicators import SubCommunicators
from .loprop import LoPropDriver
from .firstorderprop import FirstOrderProperties
from .orbitalresponse import OrbitalResponse
from .tdhfgradientdriver import TdhfGradientDriver
from .hessiandriver import HessianDriver
from .scfhessiandriver import ScfHessianDriver
from .cphfsolver import CphfSolver
from .tdhfhessiandriver import TdhfHessianDriver
from .impes import ImpesCoordinates
from .impes import ImpesDriver

# Python functions
from .errorhandler import assert_msg_critical
from .qqscheme import get_qq_type
from .qqscheme import get_qq_scheme
from .import_from_pyscf import overlap_deriv
from .import_from_pyscf import fock_deriv
from .import_from_pyscf import eri_deriv

# Environment variable: basis set path and number of OpenMP threads
from .environ import set_vlxbasispath, set_omp_num_threads, get_basis_path

set_vlxbasispath()
set_omp_num_threads()

__version__ = "1.0rc2"
