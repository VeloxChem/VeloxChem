//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "GridDriver.hpp"

#include <array>

#include "MpiFunc.hpp"

CGridDriver::CGridDriver(const int32_t globRank, const int32_t globNodes,
                         execmode runMode, MPI_Comm comm)

    : _gridLevel(5)

    , _globRank(globRank)

    , _globNodes(globNodes)

    , _isLocalMode(false)

    , _thresholdOfWeight(1.0e-15)

    , _runMode(runMode)
{
    _locRank  = mpi::rank(comm);

    _locNodes = mpi::nodes(comm);

    _isLocalMode = !mpi::compare(comm, MPI_COMM_WORLD);
}

CGridDriver::~CGridDriver()
{

}

void CGridDriver::setLevel(const int32_t gridLevel, MPI_Comm comm)
{
    if (_globRank == mpi::master())
    {
        if ((gridLevel > 0) && (gridLevel < 7)) _gridLevel = gridLevel;
    }

    if (_isLocalMode)
    {
        // FIX ME: global master to local master transfer.
    }

    mpi::bcast(_gridLevel, comm);
}

// CMolecularGrid
//void CGridDriver::generate(const CMolecule& molecule, COutputStream& oStream,
//                           MPI_Comm comm)
//{
//}

int32_t CGridDriver::_getNumberOfRadialPoints(const int32_t idElemental) const
{
    // H, He atoms

    if ((idElemental == 1) || (idElemental == 2))
    {
        const std::array<int32_t, 6> nPoints{{20, 25, 30, 35, 45, 200}};

        return nPoints[_gridLevel - 1];
    }

    // Li-Ne atoms

    if ((idElemental > 2)  && (idElemental < 11))
    {
        const std::array<int32_t, 6> nPoints{{25, 30, 35, 40, 50, 200}};

        return nPoints[_gridLevel - 1];
    }

    // Na-Ar atoms

    if ((idElemental > 10)  && (idElemental < 19))
    {
        const std::array<int32_t, 6> nPoints{{30, 35, 40, 45, 55, 200}};

        return nPoints[_gridLevel - 1];
    }

    // K-Kr atoms

    if ((idElemental > 18)  && (idElemental < 37))
    {
        const std::array<int32_t, 6> nPoints{{35, 40, 45, 50, 60, 200}};

        return nPoints[_gridLevel - 1];
    }

    // Rb-Rn atoms

    if ((idElemental > 36)  && (idElemental < 87))
    {
        const std::array<int32_t, 6> nPoints{{40, 45, 50, 55, 65, 200}};

        return nPoints[_gridLevel - 1];
    }

    // unsupported chemical element

    return 0;
}

int32_t CGridDriver::_getNumberOfAngularPoints(const int32_t idElemental) const
{
    // H, He atoms

    if ((idElemental == 1) || (idElemental == 2))
    {
        const std::array<int32_t, 6> nPoints{{50, 110, 194, 302, 434, 2030}};

        return nPoints[_gridLevel - 1];
    }

    // Li-Kr atoms

    if ((idElemental > 2)  && (idElemental < 37))
    {
        const std::array<int32_t, 6> nPoints{{110, 194, 302, 434, 590, 2030}};

        return nPoints[_gridLevel - 1];
    }

    // Rb-Rn atoms

    if ((idElemental > 36)  && (idElemental < 87))
    {
        const std::array<int32_t, 6> nPoints{{194, 302, 434, 590, 770, 2030}};

        return nPoints[_gridLevel - 1];
    }

    // unsupported chemical element

    return 0;
}
