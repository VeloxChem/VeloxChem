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
#include <iostream>
#include <sstream>

#include "ChemicalElement.hpp"
#include "GridPartitionFuncGPU.hpp"
#include "LebedevLaikovQuadrature.hpp"
#include "Log3Quadrature.hpp"
#include "MathConst.hpp"
#include "MathFunc.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"
#include "MpiFunc.hpp"
#include "StringFormat.hpp"
#include "Timer.hpp"

CGridDriver::CGridDriver(MPI_Comm comm)

    : _gridLevel(5)

    , _thresholdOfWeight(1.0e-15)
{
    _locComm = comm;
}

auto
CGridDriver::setLevel(const int64_t gridLevel) -> void
{
    if (mpi::rank(_locComm) == mpi::master())
    {
        if ((gridLevel > 0) && (gridLevel < 9)) _gridLevel = gridLevel;
    }

    _gridLevel = mpi::bcastScalar(_gridLevel, _locComm);
}

auto
CGridDriver::generate(const CMolecule& molecule, const int64_t numGpusPerNode) const -> CMolecularGrid
{
    CMolecularGrid molgrid;

    molgrid = _genGridPoints(molecule, numGpusPerNode);

    molgrid.partitionGridPoints();

    molgrid.distributeCountsAndDisplacements(_locComm);

    return molgrid;
}

auto
CGridDriver::_getNumberOfRadialPoints(const int64_t idElemental) const -> int64_t
{
    // H, He atoms

    if ((idElemental == 1) || (idElemental == 2))
    {
        const std::array<int64_t, 8> nPoints{{20, 25, 30, 40, 45, 80, 150, 200}};

        return nPoints[_gridLevel - 1];
    }

    // Li-Ne atoms

    if ((idElemental > 2) && (idElemental < 11))
    {
        const std::array<int64_t, 8> nPoints{{25, 30, 35, 45, 50, 85, 155, 200}};

        return nPoints[_gridLevel - 1];
    }

    // Na-Ar atoms

    if ((idElemental > 10) && (idElemental < 19))
    {
        const std::array<int64_t, 8> nPoints{{30, 35, 40, 50, 55, 90, 160, 200}};

        return nPoints[_gridLevel - 1];
    }

    // K-Kr atoms

    if ((idElemental > 18) && (idElemental < 37))
    {
        const std::array<int64_t, 8> nPoints{{35, 40, 45, 55, 60, 95, 165, 200}};

        return nPoints[_gridLevel - 1];
    }

    // Rb-Rn atoms

    if ((idElemental > 36) && (idElemental < 87))
    {
        const std::array<int64_t, 8> nPoints{{40, 45, 50, 60, 65, 100, 170, 200}};

        return nPoints[_gridLevel - 1];
    }

    // unsupported chemical element

    return 0;
}

auto
CGridDriver::_getNumberOfAngularPoints(const int64_t idElemental) const -> int64_t
{
    // H, He atoms

    if ((idElemental == 1) || (idElemental == 2))
    {
        const std::array<int64_t, 8> nPoints{{50, 110, 194, 302, 434, 590, 770, 2030}};

        return nPoints[_gridLevel - 1];
    }

    // Li-Kr atoms

    if ((idElemental > 2) && (idElemental < 37))
    {
        const std::array<int64_t, 8> nPoints{{110, 194, 302, 434, 590, 770, 974, 2030}};

        return nPoints[_gridLevel - 1];
    }

    // Rb-Rn atoms

    if ((idElemental > 36) && (idElemental < 87))
    {
        const std::array<int64_t, 8> nPoints{{194, 302, 434, 590, 770, 974, 2030, 2030}};

        return nPoints[_gridLevel - 1];
    }

    // unsupported chemical element

    return 0;
}

auto
CGridDriver::_startHeader(const CMolecule& molecule) const -> std::string
{
    std::stringstream ss;

    ss << "Numerical Grid For DFT Integration"
       << "\n";

    ss << std::string(36, '=') << "\n\n";

    std::string str("Grid Level          : ");

    str.append(std::to_string(_gridLevel));

    ss << fstr::format(str, 54, fmt_t::left) << "\n";

    str.assign("Radial Quadrature   : Log3");

    ss << fstr::format(str, 54, fmt_t::left) << "\n";

    str.assign("Angular Quadrature  : Lebedev-Laikov");

    ss << fstr::format(str, 54, fmt_t::left) << "\n";

    str.assign("Partitioning Scheme : SSF");

    ss << fstr::format(str, 54, fmt_t::left) << "\n\n";

    ss << "  Atom ";

    ss << fstr::format(std::string("Angular Points"), 15, fmt_t::center);

    ss << fstr::format(std::string("Radial Points"), 15, fmt_t::center);

    ss << "\n\n";

    auto molcomp = molecule.getElementalComposition();

    int64_t npoints = 0;

    for (auto i = molcomp.cbegin(); i != molcomp.cend(); ++i)
    {
        std::string label("  ");

        CChemicalElement elem;

        elem.setAtomType(*i);

        label.append(elem.getName());

        ss << fstr::format(label, 6, fmt_t::left);

        auto apoints = _getNumberOfAngularPoints(*i);

        auto rpoints = _getNumberOfRadialPoints(*i);

        ss << fstr::format(std::to_string(apoints), 15, fmt_t::center);

        ss << fstr::format(std::to_string(rpoints), 15, fmt_t::center);

        ss << "\n";

        npoints += apoints * rpoints * molecule.getNumberOfAtoms(*i);
    }

    ss << "\n";

    str.assign("Full Grid       : ");

    str.append(std::to_string(npoints));

    ss << fstr::format(str, 54, fmt_t::left) << "\n";

    return ss.str();
}

auto
CGridDriver::_finishHeader(const CMolecularGrid& molecularGrid) const -> std::string
{
    std::stringstream ss;

    std::string str("Pruned Grid     : ");

    str.append(std::to_string(molecularGrid.getNumberOfGridPoints()));

    ss << fstr::format(str, 54, fmt_t::left) << "\n";

    return ss.str();
}

auto
CGridDriver::_genGridPoints(const CMolecule& molecule, const int64_t numGpusPerNode) const -> CMolecularGrid
{
    CTimer timer;

    // molecular data

    auto natoms = molecule.getNumberOfAtoms();

    auto idselem_data = molecule.getIdsElemental();

    auto idselem = idselem_data.data();

    auto molcoords_data = molecule.getCoordinates(std::string("BOHR"));

    auto molcoords = molcoords_data.data();

    auto mdist = molecule.getMinDistances();

    auto molrm = mdist.data();

    // determine dimensions of atoms batch for each MPI process

    auto rank = mpi::rank(_locComm);

    auto nodes = mpi::nodes(_locComm);

    auto nodatm = mathfunc::batch_size(natoms, rank, nodes);

    auto nodoff = mathfunc::batch_offset(natoms, rank, nodes);

    // allocate raw grid

    int64_t bpoints = _getBatchSize(idselem, nodoff, nodatm);

    CDenseMatrix* rawgrid = new CDenseMatrix(4, bpoints);

    // generate atomic grid for each atom in molecule

    std::vector<int64_t> grid_offsets;

    for (int64_t i = 0, gridoff = 0; i < nodatm; i++)
    {
        grid_offsets.push_back(gridoff);

        auto ielem = idselem[nodoff + i];

        gridoff += _getNumberOfRadialPoints(ielem) * _getNumberOfAngularPoints(ielem);
    }

    std::vector<uint32_t> atom_ids_of_points(bpoints);
    std::vector<double> atom_min_distances(bpoints);

    std::cout << "*** Generating grid points for atoms ..." << std::endl;

    timer.reset();
    timer.start();

    for (int64_t i = 0; i < nodatm; i++)
    {
        auto iatom = nodoff + i;

        auto ielem = idselem[iatom];

        auto minrad = molrm[iatom];

        auto gridoff = grid_offsets[i];

        _genAtomGridPoints(rawgrid, minrad, gridoff, molcoords, natoms, ielem, iatom);

        auto nrpoints = _getNumberOfRadialPoints(ielem);

        auto napoints = _getNumberOfAngularPoints(ielem);

        // apply partitioning function

        std::fill(atom_ids_of_points.data() + gridoff, atom_ids_of_points.data() + gridoff + nrpoints * napoints, iatom);

        std::fill(atom_min_distances.data() + gridoff, atom_min_distances.data() + gridoff + nrpoints * napoints, minrad);
    }

    timer.stop();

    std::cout << "    Grid points generated in        " << timer.getElapsedTime() << std::endl;

    timer.reset();
    timer.start();

    gpu::applyGridPartitionFunc(rawgrid, atom_ids_of_points, atom_min_distances, bpoints, molcoords, natoms, numGpusPerNode);

    timer.stop();

    std::cout << "    Grid weights determined in      " << timer.getElapsedTime() << std::endl;

    // screen raw grid points & create prunned grid

    std::cout << "*** Screening grid points ";

    timer.reset();
    timer.start();

    bpoints = _screenRawGridPoints(rawgrid);

    auto prngrid = rawgrid->slice(0, bpoints);

    delete rawgrid;

    timer.stop();

    std::cout << "       " << timer.getElapsedTime() << std::endl;

    std::cout << "*** Communicating grid points ";

    timer.reset();
    timer.start();

    auto gathered_prngrid = mpi::gatherDenseMatricesByColumns(prngrid, _locComm);

    timer.stop();

    std::cout << "   " << timer.getElapsedTime() << std::endl;

    return CMolecularGrid(gathered_prngrid);
}

auto
CGridDriver::_getBatchSize(const int64_t* idsElemental, const int64_t offset, const int64_t nAtoms) const -> int64_t
{
    int64_t npoints = 0;

    for (int64_t i = 0; i < nAtoms; i++)
    {
        auto idx = idsElemental[offset + i];

        npoints += _getNumberOfRadialPoints(idx) * _getNumberOfAngularPoints(idx);
    }

    return npoints;
}

auto
CGridDriver::_genAtomGridPoints(CDenseMatrix*   rawGridPoints,
                                const double    minDistance,
                                const int64_t   gridOffset,
                                const TPoint3D* atomCoordinates,
                                const int64_t   nAtoms,
                                const int64_t   idElemental,
                                const int64_t   idAtomic) const -> void
{
    // determine number of grid points

    auto nrpoints = _getNumberOfRadialPoints(idElemental);

    auto napoints = _getNumberOfAngularPoints(idElemental);

    // generate radial grid points

    CLog3Quadrature rquad(nrpoints, idElemental);

    auto rpoints = rquad.generate();

    auto rrx = rpoints.row(0);

    auto rrw = rpoints.row(1);

    auto rrf = 4.0 * mathconst::getPiValue();

    // generate angular grid points

    CLebedevLaikovQuadrature aquad(napoints);

    auto apoints = aquad.generate();

    auto rax = apoints.row(0);

    auto ray = apoints.row(1);

    auto raz = apoints.row(2);

    auto raw = apoints.row(3);

    // atom coordinates

    auto atmx = atomCoordinates[idAtomic][0];

    auto atmy = atomCoordinates[idAtomic][1];

    auto atmz = atomCoordinates[idAtomic][2];

    // assemble atom grid points from quadratures

    for (int64_t i = 0; i < nrpoints; i++)
    {
        // radial quadrature point

        auto crx = rrx[i];

        auto crw = rrf * rrw[i];

        // set up pointers to raw grid

        auto gx = rawGridPoints->row(0) + gridOffset + i * napoints;

        auto gy = rawGridPoints->row(1) + gridOffset + i * napoints;

        auto gz = rawGridPoints->row(2) + gridOffset + i * napoints;

        auto gw = rawGridPoints->row(3) + gridOffset + i * napoints;

        // contract with angular quadrature

        for (int64_t j = 0; j < napoints; j++)
        {
            gx[j] = crx * rax[j] + atmx;

            gy[j] = crx * ray[j] + atmy;

            gz[j] = crx * raz[j] + atmz;

            gw[j] = crw * crx * crx * raw[j];
        }
    }
}

auto
CGridDriver::_screenRawGridPoints(CDenseMatrix* rawGridPoints) const -> int64_t
{
    // set up pointers to grid points data

    auto gx = rawGridPoints->row(0);

    auto gy = rawGridPoints->row(1);

    auto gz = rawGridPoints->row(2);

    auto gw = rawGridPoints->row(3);

    // loop over grid points

    int64_t cpoints = 0;

    for (int64_t i = 0; i < rawGridPoints->getNumberOfColumns(); i++)
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
