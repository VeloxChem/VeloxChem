//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "OptimizationGeometry.hpp"

#include "MpiFunc.hpp"

COptimizationGeometry::COptimizationGeometry(const int32_t globRank,
                                             const int32_t globNodes)

    : CBaseJob(globRank, globNodes)
{

}

void COptimizationGeometry::set(const std::string& pathToBasisSets,
                                const CInputData& inputData,
                                COutputStream& oStream)
{
    if (_globRank == mpi::master()) _startHeader(oStream);

    // FIX ME: read data
}

void COptimizationGeometry::run(COutputStream& oStream)
{
    // FIX ME: perform geometry optimization
}

void COptimizationGeometry::_startHeader(COutputStream& oStream) const
{
    oStream << fmt::header;

    oStream << "================================================";

    oStream << fmt::end;

    oStream << "=   Geometry Optimization For Single Molecule  =";

    oStream << fmt::end;

    oStream << "=================================================";

    oStream << fmt::end << fmt::blank;
}
