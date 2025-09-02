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

# C++ classes
from .veloxchemlib import AtomBasis
from .veloxchemlib import BasisFunction
from .veloxchemlib import GtoBlock
from .veloxchemlib import GtoPairBlock
from .veloxchemlib import BlockedGtoPairBlock
from .veloxchemlib import T3FlatBuffer
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
from .veloxchemlib import ThreeCenterOverlapDriver
from .veloxchemlib import ThreeCenterElectronRepulsionDriver
from .veloxchemlib import ThreeCenterElectronRepulsionGeom100Driver
from .veloxchemlib import ThreeCenterElectronRepulsionGeom010Driver
from .veloxchemlib import TwoCenterElectronRepulsionDriver
from .veloxchemlib import TwoCenterElectronRepulsionGeom100Driver
from .veloxchemlib import T4CScreener
from .veloxchemlib import FockGeom1000Driver
from .veloxchemlib import FockGeom2000Driver
from .veloxchemlib import FockGeom1100Driver
from .veloxchemlib import FockGeom1010Driver
from .veloxchemlib import RIFockGradDriver
from .veloxchemlib import XCIntegrator
from .veloxchemlib import XCFunctional
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
from .veloxchemlib import hartree_in_kcalpermol, hartree_in_kjpermol
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
from .rifockdriver import RIFockDriver
from .fockdriver import FockDriver
from .griddriver import GridDriver
from .scfrestdriver import ScfRestrictedDriver
from .scfunrestdriver import ScfUnrestrictedDriver
from .scfrestopendriver import ScfRestrictedOpenDriver
from .gradientdriver import GradientDriver
from .scfgradientdriver import ScfGradientDriver
from .dispersionmodel import DispersionModel
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
from .threepatransitiondriver import ThreePATransitionDriver
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
from .peforcefieldgenerator import PEForceFieldGenerator
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
from .mmforcefieldgenerator import MMForceFieldGenerator
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
from .excitedstateanalysisdriver import ExcitedStateAnalysisDriver
from .evbdriver import EvbDriver
from .evbffbuilder import EvbForceFieldBuilder
from .evbsystembuilder import EvbSystemBuilder
from .evbsystembuilder import EvbForceGroup
from .evbfepdriver import EvbFepDriver
from .evbdataprocessing import EvbDataProcessing
from .evbreporter import EvbReporter

## interpolation section

from .atommapper import AtomMapper
from .imforcefieldgenerator import IMForceFieldGenerator
from .imdatabasepointcollecter import IMDatabasePointCollecter
from .interpolationdriver import InterpolationDriver
from .interpolationdatapoint import InterpolationDatapoint
from .alphaoptimizer import AlphaOptimizer
from .openmmimdynamics import OpenMMIMDynamics

## external interfaces for Interpolation module
from .clustermanager import ClusterManager
from .externalexcitedstatedriver import ExternalExcitedStatesScfDriver
from .externalexcitedstategradientdriver import ExternalExcitedStatesGradientDriver
from .externalexcitedstatehessiandriver import ExternalExcitedStatesHessianDriver
from .externalscfdriver import ExternalScfDriver
from .externalgradientdriver import ExternalGradientDriver
from .externalhessiandriver import ExternalHessianDriver
from .externaloptimdriver import ExternalOptimDriver
from .localbayesresidual import LocalBayesResidual



from .mofbuilder import MofBuilder
from .conformergenerator import ConformerGenerator
from .tsguesser import TransitionStateGuesser
from .smddriver import SmdDriver
# for backward compatibility only
from .peforcefieldgenerator import PEForceFieldGenerator as LoPropDriver

# Python functions
from .errorhandler import assert_msg_critical
from .features import print_features
from .oneeints import compute_overlap_integrals
from .oneeints import compute_kinetic_energy_integrals
from .oneeints import compute_nuclear_potential_integrals
from .oneeints import compute_electric_dipole_integrals
from .oneeints import compute_linear_momentum_integrals
from .oneeints import compute_angular_momentum_integrals
from .checkpoint import read_results

# Environment variable: basis set path, number of OpenMP threads, MKL linking
from .environment import (set_vlxbasispath, get_basis_path, set_vlxdatapath,
                          get_data_path, set_omp_num_threads, configure_mkl_rt)

set_vlxbasispath()
set_vlxdatapath()
set_omp_num_threads()
configure_mkl_rt()

__version__ = "1.0rc4"
