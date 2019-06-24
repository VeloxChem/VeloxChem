# Define names for sphinx-api

from unittest.mock import Mock

AtomBasis = Mock()
BasisFunction = Mock()
OverlapIntegralsDriver = Mock()
KineticEnergyIntegralsDriver = Mock()
NuclearPotentialIntegralsDriver = Mock()
ElectricDipoleIntegralsDriver = Mock()
LinearMomentumIntegralsDriver = Mock()
AngularMomentumIntegralsDriver = Mock()
ElectronRepulsionIntegralsDriver = Mock()
SADGuessDriver = Mock()
TwoIndexes = Mock()
MOIntsBatch = Mock()
ExcitationVector = Mock()


mpi_master = Mock()
mpi_initialized = Mock()
bohr_in_angstroms = Mock()
hartree_in_ev = Mock()
mathconst_pi = Mock()


denmat = Mock()
ericut = Mock()
molorb = Mock()
moints = Mock()

assert_msg_critical = Mock()
Molecule = Mock()
MolecularBasis = Mock()
ChemicalElement = Mock()
to_angular_momentum = Mock()
AODensityMatrix = Mock()
MolecularOrbitals = Mock()
AOFockMatrix = Mock()
fockmat = Mock()
VisualizationDriver = Mock()
CubicGrid = Mock()
get_dimer_ao_indices = Mock()
szblock = Mock()
TDASigmaVectorDriver = Mock()
rotatory_strength_in_cgs = Mock()
