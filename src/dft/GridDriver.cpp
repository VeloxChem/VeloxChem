//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#include "GridDriver.hpp"

#include <array>
#include <cmath>
#include <sstream>

#include <mpi.h>

#include "ChemicalElement.hpp"
#include "LebedevLaikovQuadrature.hpp"
#include "Log3Quadrature.hpp"
#include "MathConst.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"
#include "MpiFunc.hpp"
#include "PartitionFunc.hpp"
#include "StringFormat.hpp"

CGridDriver::CGridDriver(MPI_Comm comm)
{
    _gridLevel = 5;

    _thresholdOfWeight = 1.0e-15;

    _runMode = execmode::cpu;

    _locRank = mpi::rank(comm);

    _locNodes = mpi::nodes(comm);

    _locComm = comm;
}

CGridDriver::~CGridDriver()
{
}

void
CGridDriver::setLevel(const int32_t gridLevel)
{
    if (_locRank == mpi::master())
    {
        if ((gridLevel > 0) && (gridLevel < 9)) _gridLevel = gridLevel;
    }

    mpi::bcast(_gridLevel, _locComm);
}

CMolecularGrid
CGridDriver::generate(const CMolecule& molecule) const
{
    // initialize molecular grid

    CMolecularGrid molgrid;

    // execution mode: CPU

    if (_runMode == execmode::cpu)
    {
        molgrid = _genGridPointsOnCPU(molecule);

        molgrid.partitionGridPoints();

        molgrid.distributeCountsAndDisplacements(_locRank, _locNodes, _locComm);
    }

    // execution mode: CPU/GPU

    if (_runMode == execmode::cpu_gpu)
    {
        // TODO: implement CPU/GPU code
    }

    return molgrid;
}

int32_t
CGridDriver::_getNumberOfRadialPoints(const int32_t idElemental) const
{
    // H, He atoms

    if ((idElemental == 1) || (idElemental == 2))
    {
        const std::array<int32_t, 8> nPoints{{20, 25, 30, 40, 45, 80, 150, 200}};

        return nPoints[_gridLevel - 1];
    }

    // Li-Ne atoms

    if ((idElemental > 2) && (idElemental < 11))
    {
        const std::array<int32_t, 8> nPoints{{25, 30, 35, 45, 50, 85, 155, 200}};

        return nPoints[_gridLevel - 1];
    }

    // Na-Ar atoms

    if ((idElemental > 10) && (idElemental < 19))
    {
        const std::array<int32_t, 8> nPoints{{30, 35, 40, 50, 55, 90, 160, 200}};

        return nPoints[_gridLevel - 1];
    }

    // K-Kr atoms

    if ((idElemental > 18) && (idElemental < 37))
    {
        const std::array<int32_t, 8> nPoints{{35, 40, 45, 55, 60, 95, 165, 200}};

        return nPoints[_gridLevel - 1];
    }

    // Rb-Rn atoms

    if ((idElemental > 36) && (idElemental < 87))
    {
        const std::array<int32_t, 8> nPoints{{40, 45, 50, 60, 65, 100, 170, 200}};

        return nPoints[_gridLevel - 1];
    }

    // unsupported chemical element

    return 0;
}

int32_t
CGridDriver::_getNumberOfAngularPoints(const int32_t idElemental) const
{
    // H, He atoms

    if ((idElemental == 1) || (idElemental == 2))
    {
        const std::array<int32_t, 8> nPoints{{50, 110, 194, 302, 434, 590, 770, 2030}};

        return nPoints[_gridLevel - 1];
    }

    // Li-Kr atoms

    if ((idElemental > 2) && (idElemental < 37))
    {
        const std::array<int32_t, 8> nPoints{{110, 194, 302, 434, 590, 770, 974, 2030}};

        return nPoints[_gridLevel - 1];
    }

    // Rb-Rn atoms

    if ((idElemental > 36) && (idElemental < 87))
    {
        const std::array<int32_t, 8> nPoints{{194, 302, 434, 590, 770, 974, 2030, 2030}};

        return nPoints[_gridLevel - 1];
    }

    // unsupported chemical element

    return 0;
}

std::string
CGridDriver::_startHeader(const CMolecule& molecule) const
{
    std::stringstream ss;

    ss << "Numerical Grid For DFT Integration"
       << "\n";

    ss << std::string(36, '=') << "\n\n";

    std::string str("Grid Level          : ");

    str.append(std::to_string(_gridLevel));

    ss << fstr::format(str, 54, fmt::left) << "\n";

    str.assign("Radial Quadrature   : Log3");

    ss << fstr::format(str, 54, fmt::left) << "\n";

    str.assign("Angular Quadrature  : Lebedev-Laikov");

    ss << fstr::format(str, 54, fmt::left) << "\n";

    str.assign("Partitioning Scheme : SSF");

    ss << fstr::format(str, 54, fmt::left) << "\n\n";

    ss << "  Atom ";

    ss << fstr::format(std::string("Angular Points"), 15, fmt::center);

    ss << fstr::format(std::string("Radial Points"), 15, fmt::center);

    ss << "\n\n";

    auto molcomp = molecule.getElementalComposition();

    int32_t npoints = 0;

    for (auto i = molcomp.cbegin(); i != molcomp.cend(); ++i)
    {
        std::string label("  ");

        CChemicalElement elem;

        elem.setAtomType(*i);

        label.append(elem.getName());

        ss << fstr::format(label, 6, fmt::left);

        auto apoints = _getNumberOfAngularPoints(*i);

        auto rpoints = _getNumberOfRadialPoints(*i);

        ss << fstr::format(std::to_string(apoints), 15, fmt::center);

        ss << fstr::format(std::to_string(rpoints), 15, fmt::center);

        ss << "\n";

        npoints += apoints * rpoints * molecule.getNumberOfAtoms(*i);
    }

    ss << "\n";

    str.assign("Full Grid       : ");

    str.append(std::to_string(npoints));

    ss << fstr::format(str, 54, fmt::left) << "\n";

    return ss.str();
}

std::string
CGridDriver::_finishHeader(const CMolecularGrid& molecularGrid) const
{
    std::stringstream ss;

    std::string str("Pruned Grid     : ");

    str.append(std::to_string(molecularGrid.getNumberOfGridPoints()));

    ss << fstr::format(str, 54, fmt::left) << "\n";

    return ss.str();
}

CMolecularGrid
CGridDriver::_genGridPointsOnCPU(const CMolecule& molecule) const
{
    // molecular data

    auto natoms = molecule.getNumberOfAtoms();

    auto idselem = molecule.getIdsElemental();

    auto molrx = molecule.getCoordinatesX();

    auto molry = molecule.getCoordinatesY();

    auto molrz = molecule.getCoordinatesZ();

    auto mdist = molecule.getMinDistances();

    auto molrm = mdist.data();

    // determine dimensions of atoms batch for each MPI node

    auto nodatm = mpi::batch_size(natoms, _locRank, _locNodes);

    auto nodoff = mpi::batch_offset(natoms, _locRank, _locNodes);

    // allocate raw grid

    int32_t bpoints = _getBatchSize(idselem, nodoff, nodatm);

    CMemBlock2D<double>* rawgrid = new CMemBlock2D<double>(bpoints, 4);

    // generate atomic grid for each atom in molecule

    #pragma omp parallel shared(rawgrid, idselem, molrx, molry, molrz, molrm,\
                                natoms, nodatm, nodoff)
    {
        #pragma omp single nowait
        {
            int32_t gridoff = 0;

            for (int32_t i = 0; i < nodatm; i++)
            {
                auto iatom = nodoff + i;

                auto ielem = idselem[iatom];

                auto minrad = molrm[iatom];

                #pragma omp task firstprivate(ielem, iatom, minrad, gridoff)
                {
                    _genAtomGridPoints(rawgrid, minrad, gridoff, molrx, molry, molrz, natoms, ielem, iatom);
                }

                gridoff += _getNumberOfRadialPoints(ielem)

                           * _getNumberOfAngularPoints(ielem);
            }
        }
    }

    // screen raw grid points & create prunned grid

    bpoints = _screenRawGridPoints(rawgrid);

    auto prngrid = rawgrid->slice(0, bpoints);

    delete rawgrid;

    // create molecular grid on master node

    return CMolecularGrid(prngrid.gather(_locRank, _locNodes, _locComm));
}

int32_t
CGridDriver::_getBatchSize(const int32_t* idsElemental, const int32_t offset, const int32_t nAtoms) const
{
    int32_t npoints = 0;

    for (int32_t i = 0; i < nAtoms; i++)
    {
        auto idx = idsElemental[offset + i];

        npoints += _getNumberOfRadialPoints(idx) * _getNumberOfAngularPoints(idx);
    }

    return npoints;
}

void
CGridDriver::_genAtomGridPoints(CMemBlock2D<double>* rawGridPoints,
                                const double         minDistance,
                                const int32_t        gridOffset,
                                const double*        atomCoordinatesX,
                                const double*        atomCoordinatesY,
                                const double*        atomCoordinatesZ,
                                const int32_t        nAtoms,
                                const int32_t        idElemental,
                                const int32_t        idAtomic) const
{
    // determine number of grid points

    auto nrpoints = _getNumberOfRadialPoints(idElemental);

    auto napoints = _getNumberOfAngularPoints(idElemental);

    // generate radial grid points

    CLog3Quadrature rquad(nrpoints, idElemental);

    auto rpoints = rquad.generate();

    auto rrx = rpoints.data(0);

    auto rrw = rpoints.data(1);

    auto rrf = 4.0 * mathconst::getPiValue();

    // generate angular grid points

    CLebedevLaikovQuadrature aquad(napoints);

    auto apoints = aquad.generate();

    auto rax = apoints.data(0);

    auto ray = apoints.data(1);

    auto raz = apoints.data(2);

    auto raw = apoints.data(3);

    // atom coordinates

    auto atmx = atomCoordinatesX[idAtomic];

    auto atmy = atomCoordinatesY[idAtomic];

    auto atmz = atomCoordinatesZ[idAtomic];

    // assemble atom grid points from quadratures

    for (int32_t i = 0; i < nrpoints; i++)
    {
        // radial quadrature point

        auto crx = rrx[i];

        auto crw = rrf * rrw[i];

        // set up pointers to raw grid

        auto gx = rawGridPoints->data(0, gridOffset + i * napoints);

        auto gy = rawGridPoints->data(1, gridOffset + i * napoints);

        auto gz = rawGridPoints->data(2, gridOffset + i * napoints);

        auto gw = rawGridPoints->data(3, gridOffset + i * napoints);

        // contract with angular quadrature

        #pragma omp simd
        for (int32_t j = 0; j < napoints; j++)
        {
            gx[j] = crx * rax[j] + atmx;

            gy[j] = crx * ray[j] + atmy;

            gz[j] = crx * raz[j] + atmz;

            gw[j] = crw * crx * crx * raw[j];
        }
    }

    // apply partitioning function

    partfunc::ssf(
        rawGridPoints, minDistance, gridOffset, nrpoints * napoints, atomCoordinatesX, atomCoordinatesY, atomCoordinatesZ, nAtoms, idAtomic);
}

int32_t
CGridDriver::_screenRawGridPoints(CMemBlock2D<double>* rawGridPoints) const
{
    // set up pointers to grid points data

    auto gx = rawGridPoints->data(0);

    auto gy = rawGridPoints->data(1);

    auto gz = rawGridPoints->data(2);

    auto gw = rawGridPoints->data(3);

    // loop over grid points

    int32_t cpoints = 0;

    for (int32_t i = 0; i < rawGridPoints->size(0); i++)
    {
        if (std::fabs(gw[i]) > _thresholdOfWeight)
        {
            gx[cpoints] = gx[i];

            gy[cpoints] = gy[i];

            gz[cpoints] = gz[i];

            gw[cpoints] = gw[i];

            cpoints++;
        }
    }

    return cpoints;
}
