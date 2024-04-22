#
#                           VELOXCHEM 1.0-RC3
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2022 by VeloxChem developers. All rights reserved.
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
from .veloxchemlib import SADGuessDriver
from .veloxchemlib import DispersionModel
from .veloxchemlib import DenseMatrix
from .veloxchemlib import TwoIndexes
from .veloxchemlib import MOIntsBatch
from .veloxchemlib import ExcitationVector
from .veloxchemlib import XCIntegrator
from .veloxchemlib import XCFunctional

# C++ functions
from .veloxchemlib import available_functionals
from .veloxchemlib import mpi_master
from .veloxchemlib import mpi_initialized
from .veloxchemlib import ao_matrix_to_veloxchem
from .veloxchemlib import ao_matrix_to_dalton
from .veloxchemlib import get_basis_function_indices_for_atom
from .veloxchemlib import bohr_in_angstrom, bohr_in_angstroms
from .veloxchemlib import hartree_in_ev
from .veloxchemlib import hartree_in_kcalpermol
from .veloxchemlib import hartree_in_wavenumber, hartree_in_wavenumbers
from .veloxchemlib import dipole_in_debye
from .veloxchemlib import rotatory_strength_in_cgs
from .veloxchemlib import extinction_coefficient_from_beta
from .veloxchemlib import fine_structure_constant
from .veloxchemlib import mathconst_pi
from .veloxchemlib import parse_xc_func

# C++ enums
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .veloxchemlib import szblock
from .veloxchemlib import ericut
from .veloxchemlib import molorb
from .veloxchemlib import moints

# Python classes
from .atomtypeidentifier import AtomTypeIdentifier
from .seminario import Seminario
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
from .gradientdriver import GradientDriver
from .scfgradientdriver import ScfGradientDriver
from .tddftgradientdriver import TddftGradientDriver
from .xtbdriver import XtbDriver
from .xtbgradientdriver import XtbGradientDriver
from .xtbhessiandriver import XtbHessianDriver
from .optimizationdriver import OptimizationDriver
from .mointsdriver import MOIntegralsDriver
from .mp2driver import Mp2Driver
from .cubicgrid import CubicGrid
from .visualizationdriver import VisualizationDriver
from .excitondriver import ExcitonModelDriver
from .tdaeigensolver import TdaEigenSolver
from .blockdavidson import BlockDavidsonSolver
from .lreigensolver import LinearResponseEigenSolver
from .lrsolver import LinearResponseSolver
from .cppsolver import ComplexResponse
from .c6driver import C6Driver
from .quadraticresponsedriver import QuadraticResponseDriver
from .cubicresponsedriver import CubicResponseDriver
from .shgdriver import ShgDriver
from .tpatransitiondriver import TpaTransitionDriver
from .tpafulldriver import TpaFullDriver
from .tpareddriver import TpaReducedDriver
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
from .tddftorbitalresponse import TddftOrbitalResponse
from .tddftgradientdriver import TddftGradientDriver
from .hessiandriver import HessianDriver
from .scfhessiandriver import ScfHessianDriver
from .cphfsolver import CphfSolver
from .tdhfhessiandriver import TdhfHessianDriver
from .polorbitalresponse import PolOrbitalResponse
from .polarizabilitygradient import PolarizabilityGradient
from .vibrationalanalysis import VibrationalAnalysis
from .forcefieldgenerator import ForceFieldGenerator
from .openmmdriver import OpenMMDriver
from .openmmgradientdriver import OpenMMGradientDriver
from .orbitalviewer import OrbitalViewer
from .numerovdriver import NumerovDriver
from .mmdriver import MMDriver
from .mmgradientdriver import MMGradientDriver

# for backward compatibility
from .veloxchemlib import XCIntegrator as XCNewIntegrator
from .veloxchemlib import XCFunctional as XCNewFunctional
from .veloxchemlib import parse_xc_func as new_parse_xc_func
from .xtbdriver import XtbDriver as XTBDriver
from .xtbgradientdriver import XtbGradientDriver as XTBGradientDriver
from .xtbhessiandriver import XtbHessianDriver as XTBHessianDriver
from .tdaeigensolver import TdaEigenSolver as TDAExciDriver
from .shgdriver import ShgDriver as SHGDriver
from .tpafulldriver import TpaFullDriver as TPAFullDriver
from .tpareddriver import TpaReducedDriver as TPAReducedDriver

# Python functions
from .errorhandler import assert_msg_critical
from .qqscheme import get_qq_type
from .qqscheme import get_qq_scheme
from .features import print_features
from .import_from_pyscf import overlap_deriv
from .import_from_pyscf import fock_deriv
from .import_from_pyscf import eri_deriv

# Environment variable: basis set path, number of OpenMP threads, MKL linking
from .environment import (set_vlxbasispath, set_omp_num_threads, get_basis_path,
                          configure_mkl_rt)

set_vlxbasispath()
set_omp_num_threads()
configure_mkl_rt()

__version__ = "1.0rc3"
