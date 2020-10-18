//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef XTBDriver_hpp
#define XTBDriver_hpp

#include <cstdint>
#include <string>

#include "mpi.h"
#include "Molecule.hpp"

#ifdef ENABLE_XTB
#include "xtb.h"
#endif

/**
 Class CXTBDriver enables DFT-B computations using XTB package from Grimme group.

 @author Z. Rinkevicius
 */
class CXTBDriver
{
    /**
     The rank of associated local MPI process.
     */
    int32_t _locRank;

    /**
     The total number of local MPI processes.
     */
    int32_t _locNodes;

    /**
     The MPI communicator.
     */
    MPI_Comm _locComm;
    
#ifdef ENABLE_XTB
    /**
     The XTB package environment.
    */
    xtb_TEnvironment _environment;
    
    /**
     The DFT-B calculator.
    */
    xtb_TCalculator _calculator;
    
    /**
     The DFT-B results data.
    */
    xtb_TResults _results;

    /**
     Converts VeloxChem molecule to XTB molecular system data structure.
      
     @param molecule the molecule.
     @returm the XTB  molecular system data structure.
    */
    xtb_TMolecule _set_molecule(const CMolecule& molecule);
#endif

   public:
    /**
     Creates a XTB driver object using MPI info.

     @param comm the MPI communicator.
     */
    CXTBDriver(MPI_Comm comm);

    /**
     Destroys a XTB driver object.
     */
    ~CXTBDriver();
    
    /**
     Computes DTB-B single point energy using XTB package.

     @param molecule the molecule.
     @param method the string with DFT-B method name supported by XTB package.
     */
    void compute(const CMolecule&   molecule,
                 const std::string& method);

    /**
     Checks if XTB package is available. 
     
     @return true if XTB package available, false otherwise.  
     */ 
    bool isAvailable() const;
    
    /**
     Gets state of XTB environment.
    
     @return the state of XTB environment.
    */
    bool getState();
};

#endif /* XTBDriver_hpp */
