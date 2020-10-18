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
     @param method the GTN-XT string method. 
     */
    void compute(const CMolecule&   molecule, 
		 const std::string& method); 

    /**
     Checks if XTB package is available. 
     
     @return true if XTB package available, false otherwise.  
     */ 
    bool isAvailable() const;
};

#endif /* XTBDriver_hpp */
