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
from .veloxchemlib import AtomBasis
from .veloxchemlib import BasisFunction
from .veloxchemlib import OverlapDriver
from .veloxchemlib import KineticEnergyDriver
from .veloxchemlib import NuclearPotentialDriver
from .veloxchemlib import NuclearPotentialErfDriver
from .veloxchemlib import ElectricDipoleMomentDriver
from .veloxchemlib import OverlapGeom100Driver
from .veloxchemlib import OverlapGeom200Driver
from .veloxchemlib import OverlapGeom101Driver
from .veloxchemlib import KineticEnergyGeom100Driver
from .veloxchemlib import KineticEnergyGeom200Driver
from .veloxchemlib import KineticEnergyGeom101Driver
from .veloxchemlib import NuclearPotentialGeom010Driver
from .veloxchemlib import NuclearPotentialGeom100Driver
from .veloxchemlib import NuclearPotentialGeom200Driver
from .veloxchemlib import NuclearPotentialGeom110Driver
from .veloxchemlib import NuclearPotentialGeom101Driver
from .veloxchemlib import NuclearPotentialErfGeom100Driver
from .veloxchemlib import NuclearPotentialErfGeom010Driver
from .veloxchemlib import ElectricDipoleMomentGeom100Driver
from .veloxchemlib import T4CScreener
from .veloxchemlib import FockGeom1000Driver
from .veloxchemlib import FockGeom2000Driver
from .veloxchemlib import FockGeom1100Driver
from .veloxchemlib import FockGeom1010Driver
from .veloxchemlib import XCIntegrator
from .veloxchemlib import XCFunctional
from .veloxchemlib import DispersionModel
from .veloxchemlib import XCMolecularGradient
from .veloxchemlib import SubMatrix

# for backward compatibility only
from .veloxchemlib import ElectricDipoleIntegralsDriver

# C++ functions
from .veloxchemlib import is_chemical_element
from .veloxchemlib import chemical_element_name
from .veloxchemlib import chemical_element_label
from .veloxchemlib import chemical_element_identifier
from .veloxchemlib import chemical_element_mass
from .veloxchemlib import chemical_element_max_angular_momentum
from .veloxchemlib import available_functionals
from .veloxchemlib import available_pdft_functionals
from .veloxchemlib import mpi_master
from .veloxchemlib import bohr_in_angstrom
from .veloxchemlib import hartree_in_ev
from .veloxchemlib import hartree_in_kcalpermol
from .veloxchemlib import hartree_in_wavenumber, hartree_in_wavenumbers
from .veloxchemlib import dipole_in_debye
from .veloxchemlib import rotatory_strength_in_cgs
from .veloxchemlib import extinction_coefficient_from_beta
from .veloxchemlib import fine_structure_constant
from .veloxchemlib import parse_xc_func
from .veloxchemlib import make_matrix
from .veloxchemlib import make_matrices
from .veloxchemlib import partition_atoms

# C++ enums
from .veloxchemlib import mat_t
from .veloxchemlib import denmat

# Python enums
from .molecularorbitals import molorb

# Python classes
from .atomtypeidentifier import AtomTypeIdentifier
from .seminario import Seminario
from .inputparser import InputParser
from .outputstream import OutputStream
from .matrix import Matrix
from .matrices import Matrices
from .molecule import Molecule
from .molecularbasis import MolecularBasis
from .aodensitymatrix import AODensityMatrix
from .molecularorbitals import MolecularOrbitals
from .fockdriver import FockDriver
from .griddriver import GridDriver
from .scfrestdriver import ScfRestrictedDriver
from .scfunrestdriver import ScfUnrestrictedDriver
from .scfrestopendriver import ScfRestrictedOpenDriver
from .gradientdriver import GradientDriver
from .scfgradientdriver import ScfGradientDriver
from .xtbdriver import XtbDriver
from .xtbgradientdriver import XtbGradientDriver
from .xtbhessiandriver import XtbHessianDriver
from .optimizationdriver import OptimizationDriver
from .mointsdriver import MOIntegralsDriver
from .mp2driver import Mp2Driver
from .cpcmdriver import CpcmDriver
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
#from .rspcustomproperty import CustomProperty
from .mpitask import MpiTask
from .subcommunicators import SubCommunicators
from .loprop import LoPropDriver
from .firstorderprop import FirstOrderProperties
from .tddftorbitalresponse import TddftOrbitalResponse
from .tddftgradientdriver import TddftGradientDriver
from .hessiandriver import HessianDriver
from .scfhessiandriver import ScfHessianDriver
from .cphfsolver import CphfSolver
from .hessianorbitalresponse import HessianOrbitalResponse
from .tdhfhessiandriver import TdhfHessianDriver
from .polorbitalresponse import PolOrbitalResponse
from .polarizabilitygradient import PolarizabilityGradient
from .vibrationalanalysis import VibrationalAnalysis
from .forcefieldgenerator import ForceFieldGenerator
from .openmmdriver import OpenMMDriver
from .openmmgradientdriver import OpenMMGradientDriver
from .orbitalviewer import OrbitalViewer
from .densityviewer import DensityViewer
from .numerovdriver import NumerovDriver
from .mmdriver import MMDriver
from .mmgradientdriver import MMGradientDriver
from .symmetryanalyzer import SymmetryAnalyzer
from .solvationbuilder import SolvationBuilder
from .solvationfepdriver import SolvationFepDriver
from .openmmdynamics import OpenMMDynamics
from .evbdriver import EvbDriver
from .evbffbuilder import EvbForceFieldBuilder
from .evbsystembuilder import EvbSystemBuilder
from .evbfepdriver import FepDriver
from .evbdataprocessing import EvbDataProcessing


# Python functions
from .errorhandler import assert_msg_critical
from .features import print_features
from .oneeints import compute_overlap_integrals
from .oneeints import compute_kinetic_energy_integrals
from .oneeints import compute_nuclear_potential_integrals
from .oneeints import compute_electric_dipole_integrals
from .oneeints import compute_linear_momentum_integrals
from .oneeints import compute_angular_momentum_integrals

# Environment variable: basis set path, number of OpenMP threads, MKL linking
from .environment import (set_vlxbasispath, set_omp_num_threads, get_basis_path,
                          configure_mkl_rt)

set_vlxbasispath()
set_omp_num_threads()
configure_mkl_rt()

__version__ = "1.0rc3"
