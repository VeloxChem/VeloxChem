//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "XTBDriver.hpp"

#ifdef ENABLE_MKL
#include "xtb.h"
#endif

#include "MpiFunc.hpp"

CXTBDriver::CXTBDriver(MPI_Comm comm)
{
    _locRank = mpi::rank(comm);

    _locNodes = mpi::nodes(comm);

    mpi::duplicate(comm, &_locComm);
}

CXTBDriver::~CXTBDriver()
{
    mpi::destroy(&_locComm);
}

void 
CXTBDriver::compute()
{
  int    const natoms = 7;
  int    const attyp[7] = {6,6,6,1,1,1,1};
  double const charge = 0.0;
  int    const uhf = 0;
  double const coord[3*7] =
      {0.00000000000000, 0.00000000000000,-1.79755622305860,
       0.00000000000000, 0.00000000000000, 0.95338756106749,
       0.00000000000000, 0.00000000000000, 3.22281255790261,
      -0.96412815539807,-1.66991895015711,-2.53624948351102,
      -0.96412815539807, 1.66991895015711,-2.53624948351102,
       1.92825631079613, 0.00000000000000,-2.53624948351102,
       0.00000000000000, 0.00000000000000, 5.23010455462158};

  /*
   * All objects except for the molecular structure can be
   * constructued without other objects present.
   *
   * The construction of the molecular structure locks the
   * number of atoms, atomic number, total charge, multiplicity
   * and boundary conditions.
  **/
  xtb_TEnvironment env = xtb_newEnvironment();
  xtb_TCalculator calc = xtb_newCalculator();
  xtb_TResults res = xtb_newResults();
  xtb_TMolecule mol = xtb_newMolecule(
      env, &natoms, attyp, coord, &charge, &uhf, NULL, NULL);
  if (xtb_checkEnvironment(env)) {
    xtb_showEnvironment(env, NULL);
    return;
  }

  /*
   * Apply changes to the environment which will be respected
   * in all further API calls.
  **/
  xtb_setVerbosity(env, XTB_VERBOSITY_FULL);
  if (xtb_checkEnvironment(env)) {
    xtb_showEnvironment(env, NULL);
    return;
  }

  /*
   * Load a parametrisation, the last entry is a char* which can
   * be used to provide a particular parameter file.
   *
   * Otherwise the XTBPATH environment variable is used to find
   * parameter files.
   *
   * The calculator has to be reconstructed if the molecular
   * structure is reconstructed.
  **/
  xtb_loadGFN2xTB(env, mol, calc, NULL);
  if (xtb_checkEnvironment(env)) {
    xtb_showEnvironment(env, NULL);
    return;
  }

  /*
   * Actual calculation, will populate the results object,
   * the API can raise errors on failed SCF convergence or other
   * numerical problems.
   *
   * Not supported boundary conditions are usually raised here.
  **/
  xtb_singlepoint(env, mol, calc, res);
  if (xtb_checkEnvironment(env)) {
    xtb_showEnvironment(env, NULL);
    return;
  }

  /*
   * Query the environment for properties, an error in the environment
   * is not considered blocking for this calls and allows to query
   * for multiple entries before handling possible errors
  **/
  double energy = 0.0; 
  xtb_getEnergy(env, res, &energy);
  if (xtb_checkEnvironment(env)) {
    xtb_showEnvironment(env, NULL);
    return;
  }

  /*
   * deconstructor will deallocate the objects and overwrite the
   * pointer with NULL
  **/
  xtb_delResults(&res);
  xtb_delCalculator(&calc);
  xtb_delMolecule(&mol);
  xtb_delEnvironment(&env);
}

