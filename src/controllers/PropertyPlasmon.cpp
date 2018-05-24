//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "PropertyPlasmon.hpp"

#include "MpiFunc.hpp"

CPropertyPlasmon::CPropertyPlasmon(const int32_t  globRank,
                                   const int32_t  globNodes,
                                   const execmode runMode)

    : CBaseJob(globRank, globNodes, runMode)
{

}

void
CPropertyPlasmon::set(const std::string&   pathToBasisSets, 
                      const std::string&   pathToForceFields,
                      const CInputData&    inputData,
                            COutputStream& oStream)
{
    if (_globRank == mpi::master()) _startHeader(oStream);

    // FIX ME: read data
}

void
CPropertyPlasmon::run(COutputStream& oStream,
                      MPI_Comm       comm)
{
    
}

void
CPropertyPlasmon::_startHeader(COutputStream& oStream) const
{
    oStream << fmt::header;

    oStream << "=====================================";

    oStream << fmt::end;

    oStream << "= Classical CMM Plasmon Calculation =";

    oStream << fmt::end;

    oStream << "=====================================";

    oStream << fmt::end << fmt::blank;
}
