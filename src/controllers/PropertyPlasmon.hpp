//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef PropertyPlasmon_hpp
#define PropertyPlasmon_hpp

#include "BaseJob.hpp"

/**
 Class CPropertyPlasmon manages classical CMM plasmon job for case of
 single molecule.
 
 @author Z. Rinkevicius
 */
class CPropertyPlasmon : public CBaseJob
{
    /**
     Prints start message for classical CMM plasmon job to output stream.

     @param oStream the output stream.
     */
    void _startHeader(COutputStream& oStream) const;

public:
    
    /**
     Creates a classical CMM plasmon job object.
     
     @param globRank the the rank of MPI process.
     @param globNodes the total number of MPI processes.
     @param runMode the execution mode of job.
     */
    CPropertyPlasmon(const int32_t  globRank,
                     const int32_t  globNodes,
                     const execmode runMode);

    /**
     Sets parameters of classical CMM plasmon job.

     @param pathToBasisSets the path to basis set library.
     @param inputData the input data object.
     @param oStream the output stream.
     */
    void set(const std::string&   pathToBasisSets,
             const CInputData&    inputData,
                   COutputStream& oStream) override;

    /**
     Executes a classical CMM plasmon job.

     @param comm the MPI communicator.
     @param oStream the output stream.
     */
    void run(COutputStream& oStream,
             MPI_Comm       comm) override;
};

#endif /* PropertyPlasmony_hpp */
