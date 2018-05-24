//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef SinglePointEnergy_hpp
#define SinglePointEnergy_hpp

#include <string>

#include "BaseJob.hpp"
#include "InputData.hpp"
#include "OutputStream.hpp"
#include "Molecule.hpp"
#include "MolecularBasis.hpp"

/**
 Class CSinglePointEnergy manages single point energy computation job for case
 of single molecule.
 
 @author Z. Rinkevicius
 */
class CSinglePointEnergy : public CBaseJob
{
    /**
     The molecular data.
     */
    CMolecule _molecule;      

    /**
     The AO basis set.
     */
    CMolecularBasis _aoBasis;

    /**
     The AO-RI/J(K) basis set.
     */
    CMolecularBasis _riBasis;

    //CAODensityMatrix _density;

    /**
     Prints start message for single point energy job to output stream.

     @param oStream the output stream.
     */
    void _startHeader(COutputStream& oStream) const;

public:
    
    /**
     Creates a single point energy computation job object.
     
     @param globRank the the rank of MPI process.
     @param globNodes the total number of MPI processes.
     @param runMode the execution mode of job.
     */
    CSinglePointEnergy(const int32_t  globRank,
                       const int32_t  globNodes,
                       const execmode runMode);

    /**
     Sets parameters of single point energy computation job.

     @param pathToBasisSets the path to basis sets library.
     @param pathToForceFields the path to force fields library.
     @param inputData the input data object.
     @param oStream the output stream.
     */
    void set(const std::string&   pathToBasisSets,
             const std::string&   pathToForceFields,
             const CInputData&    inputData,
                   COutputStream& oStream) override;

    /**
     Executes a single point energy computation job.

     @param comm the MPI communicator.
     @param oStream the output stream.
     */
    void run(COutputStream& oStream,
             MPI_Comm       comm) override;
};

#endif /* SinglePointEnergy_hpp */
