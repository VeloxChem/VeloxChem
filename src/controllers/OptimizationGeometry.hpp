//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef OptimizationGeometry_hpp
#define OptimizationGeometry_hpp

#include "BaseJob.hpp"

/**
 Class COptimizationGeometry manages geometry optimization job for case of
 single molecule.
 
 @author Z. Rinkevicius
 */
class COptimizationGeometry : public CBaseJob
{
    /**
     Prints start message for geometry optimization job to output stream.

     @param oStream the output stream.
     */
    void _startHeader(COutputStream& oStream) const;

public:
    
    /**
     Creates a geometry optimization job object.
     
     @param globRank the the rank of MPI process.
     @param globNodes the total number of MPI processes.
     @param runMode the execution mode of job.
     */
    COptimizationGeometry(const int32_t globRank, const int32_t globNodes,
                          const execmode runMode);

    /**
     Sets parameters of geometry optimization job.

     @param pathToBasisSets the path to basis set library.
     @param inputData the input data object.
     @param oStream the output stream.
     */
    void set(const std::string& pathToBasisSets, const CInputData& inputData,
             COutputStream& oStream) override;

    /**
     Executes a geometry optimization job.

     @param comm the MPI communicator.
     @param oStream the output stream.
     */
    void run(COutputStream& oStream, MPI_Comm comm) override;
};

#endif /* OptimizationGeometry_hpp */
