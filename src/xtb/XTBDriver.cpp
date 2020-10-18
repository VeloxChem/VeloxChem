//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "XTBDriver.hpp"

#include <vector>

#include "MpiFunc.hpp"

CXTBDriver::CXTBDriver(MPI_Comm comm)
{
    _locRank = mpi::rank(comm);

    _locNodes = mpi::nodes(comm);

    _locComm = comm;

#ifdef ENABLE_XTB
    _environment = xtb_newEnvironment();
        
    _calculator = xtb_newCalculator();
          
    _results = xtb_newResults();
    
    xtb_setVerbosity(_environment, XTB_VERBOSITY_FULL);
#endif
}

CXTBDriver::~CXTBDriver()
{
#ifdef ENABLE_XTB
    xtb_delResults(&_results);
    
    xtb_delCalculator(&_calculator);
    
    xtb_delEnvironment(&_environment);
#endif
}

void 
CXTBDriver::compute(const CMolecule&   molecule,
                    const std::string& method)
{
#ifdef ENABLE_XTB
 
    // set up molecular data structure
    
    auto tmol = _set_molecule(molecule);

    // load DFT-B method parameters
    
    if (method == "gfn2-xtb")
    {
        xtb_loadGFN2xTB(_environment, tmol, _calculator, NULL);
    }
    else if (method == "gfn1-xtb")
    {
        xtb_loadGFN1xTB(_environment, tmol, _calculator, NULL);
    }
    else if (method == "gfn0-xtb")
    {
        xtb_loadGFN0xTB(_environment, tmol, _calculator, NULL);
    }
    else
    {
        return;
    }
  
    // perform single point calculation 

    xtb_singlepoint(_environment, tmol, _calculator, _results);

    // delete molecular data structure     

    xtb_delMolecule(&tmol);
#endif
}

bool 
CXTBDriver::isAvailable() const 
{
#ifdef ENABLE_XTB
    return true;
#endif
    return false; 
}

bool
CXTBDriver::getState()
{
#ifdef ENABLE_XTB
    return xtb_checkEnvironment(_environment) > 0;
#endif
    return false;
}

#ifdef ENABLE_XTB
xtb_TMolecule
CXTBDriver::_set_molecule(const CMolecule& molecule)
{
    int natoms = static_cast<int>(molecule.getNumberOfAtoms());
    
    double charge = molecule.getCharge();
    
    int uhf = 0;
    
    std::vector<int> atoms(natoms, 0);
    
    std::vector<double> coords(3 * natoms, 0.0);
    
    // reformat molecular data
    
    auto eleids = molecule.getIdsElemental();
    
    auto rx =  molecule.getCoordinatesX();
    
    auto ry =  molecule.getCoordinatesY();
    
    auto rz =  molecule.getCoordinatesZ();
    
    for (int i = 0; i < natoms; i++)
    {
        atoms.at(i) = static_cast<int>(eleids[i]);
        
        coords.at(3 * i) = rx[i];
        
        coords.at(3 * i + 1) = ry[i];
        
        coords.at(3 * i + 2) = rz[i];
    }
    
    return xtb_newMolecule(_environment, &natoms, atoms.data(),
                           coords.data(),
                           &charge, &uhf,
                           NULL, NULL);
}
#endif

