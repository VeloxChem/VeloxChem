//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "BaseJob.hpp"

CBaseJob::CBaseJob(const int32_t globRank, const int32_t globNodes,
                   const execmode runMode)

    : _state(true)

    , _globRank(globRank)

    , _globNodes(globNodes)

    , _runMode(runMode)
{

}

bool CBaseJob::getState() const
{
    return _state;
}
