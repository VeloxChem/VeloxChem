# C++ classes
from .veloxchemlib import AtomBasis
from .veloxchemlib import BasisFunction
from .veloxchemlib import OverlapIntegralsDriver
from .veloxchemlib import KineticEnergyIntegralsDriver
from .veloxchemlib import NuclearPotentialIntegralsDriver
from .veloxchemlib import ElectricDipoleIntegralsDriver
from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import SADGuessDriver
from .veloxchemlib import TwoIndexes
from .veloxchemlib import MOIntsBatch
from .veloxchemlib import ExcitationVector

# C++ functions
from .veloxchemlib import mpi_master
from .veloxchemlib import mpi_initialized
from .veloxchemlib import bohr_in_angstroms
from .veloxchemlib import hartree_in_ev
from .veloxchemlib import mathconst_pi

# C++ enums
from .veloxchemlib import denmat
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
from .mointsdriver import MOIntegralsDriver
from .visualizationdriver import VisualizationDriver
from .rspdriver import ResponseDriver
from .tdaexcidriver import TDAExciDriver
from .blockdavidson import BlockDavidsonSolver
from .lrsolver import LinearResponseSolver
from .mpitask import MpiTask

# Python functions
from .errorhandler import assert_msg_critical
from .qqscheme import get_qq_type
from .qqscheme import get_qq_scheme
from .main import main

# Environment variable: basis set path

import os

if 'VLXBASISPATH' not in os.environ:
    module_path = os.path.dirname(os.path.abspath(__file__))
    os.environ['VLXBASISPATH'] = module_path + '/../../../basis'
