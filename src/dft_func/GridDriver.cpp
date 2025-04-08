//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "GridDriver.hpp"

#include <array>
#include <cmath>
#include <sstream>

#include "ChemicalElement.hpp"
#include "LebedevLaikovQuadrature.hpp"
#include "M4T2Quadrature.hpp"
#include "MathConst.hpp"
#include "MathFunc.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"
#include "PartitionFunc.hpp"
#include "StringFormat.hpp"

CGridDriver::CGridDriver()

    : _gridLevel(5)

    , _thresholdOfWeight(1.0e-15)
{
}

auto
CGridDriver::setLevel(const int gridLevel) -> void
{
    if ((gridLevel > 0) && (gridLevel < 9)) _gridLevel = gridLevel;
}

auto
CGridDriver::_getNumberOfRadialPoints(const int idElemental) const -> int
{
    // H, He atoms

    if ((idElemental == 1) || (idElemental == 2))
    {
        const std::array<int, 8> nPoints{{20, 25, 30, 40, 45, 80, 150, 200}};

        return nPoints[_gridLevel - 1];
    }

    // Li-Ne atoms

    if ((idElemental > 2) && (idElemental < 11))
    {
        const std::array<int, 8> nPoints{{25, 30, 35, 45, 50, 85, 155, 200}};

        return nPoints[_gridLevel - 1];
    }

    // Na-Ar atoms

    if ((idElemental > 10) && (idElemental < 19))
    {
        const std::array<int, 8> nPoints{{30, 35, 40, 50, 55, 90, 160, 200}};

        return nPoints[_gridLevel - 1];
    }

    // K-Kr atoms

    if ((idElemental > 18) && (idElemental < 37))
    {
        const std::array<int, 8> nPoints{{35, 40, 45, 55, 60, 95, 165, 200}};

        return nPoints[_gridLevel - 1];
    }

    // Rb-Rn atoms

    if ((idElemental > 36) && (idElemental < 87))
    {
        const std::array<int, 8> nPoints{{40, 45, 50, 60, 65, 100, 170, 200}};

        return nPoints[_gridLevel - 1];
    }

    // unsupported chemical element

    return 0;
}

auto
CGridDriver::_getNumberOfAngularPoints(const int idElemental) const -> int
{
    // H, He atoms

    if ((idElemental == 1) || (idElemental == 2))
    {
        const std::array<int, 8> nPoints{{50, 110, 194, 302, 434, 590, 770, 2030}};

        return nPoints[_gridLevel - 1];
    }

    // Li-Kr atoms

    if ((idElemental > 2) && (idElemental < 37))
    {
        const std::array<int, 8> nPoints{{110, 194, 302, 434, 590, 770, 974, 2030}};

        return nPoints[_gridLevel - 1];
    }

    // Rb-Rn atoms

    if ((idElemental > 36) && (idElemental < 87))
    {
        const std::array<int, 8> nPoints{{194, 302, 434, 590, 770, 974, 2030, 2030}};

        return nPoints[_gridLevel - 1];
    }

    // unsupported chemical element

    return 0;
}

auto
CGridDriver::generate_local_grid(const CMolecule& molecule, const int rank, const int nnodes) const -> CMolecularGrid
{
    // molecular data

    auto natoms = molecule.number_of_atoms();

    auto idselem_data = molecule.identifiers();

    auto idselem = idselem_data.data();

    auto molcoords_data = molecule.coordinates(std::string("BOHR"));

    auto molcoords = molcoords_data.data();

    auto mdist = molecule.min_distances();

    auto molrm = mdist.data();

    // determine dimensions of atoms batch for this MPI process

    auto nodatm = mathfunc::batch_size(natoms, rank, nnodes);

    auto nodoff = mathfunc::batch_offset(natoms, rank, nnodes);

    // allocate raw grid

    int bpoints = _getBatchSize(idselem, nodoff, nodatm);

    CDenseMatrix* rawgrid = new CDenseMatrix(4, bpoints);

    // generate atomic grid for each atom in molecule

#pragma omp parallel shared(rawgrid, idselem, molcoords, molrm, natoms, nodatm, nodoff)
    {
#pragma omp single nowait
        {
            int gridoff = 0;

            for (int i = 0; i < nodatm; i++)
            {
                auto iatom = nodoff + i;

                auto ielem = idselem[iatom];

                auto minrad = molrm[iatom];

#pragma omp task firstprivate(ielem, iatom, minrad, gridoff)
                {
                    _genAtomGridPoints(rawgrid, minrad, gridoff, molcoords, natoms, ielem, iatom);
                }

                gridoff += _getNumberOfRadialPoints(ielem) * _getNumberOfAngularPoints(ielem);
            }
        }
    }

    // screen raw grid points & create prunned grid

    bpoints = _screenRawGridPoints(rawgrid);

    auto prngrid = rawgrid->slice(0, bpoints);

    delete rawgrid;

    return CMolecularGrid(prngrid);
}

auto
CGridDriver::_getBatchSize(const int* idsElemental, const int offset, const int nAtoms) const -> int
{
    int npoints = 0;

    for (int i = 0; i < nAtoms; i++)
    {
        auto idx = idsElemental[offset + i];

        npoints += _getNumberOfRadialPoints(idx) * _getNumberOfAngularPoints(idx);
    }

    return npoints;
}

auto
CGridDriver::_genAtomGridPoints(CDenseMatrix*   rawGridPoints,
                                const double    minDistance,
                                const int   gridOffset,
                                const TPoint<double>* atomCoordinates,
                                const int   nAtoms,
                                const int   idElemental,
                                const int   idAtomic) const -> void
{
    // determine number of grid points

    auto nrpoints = _getNumberOfRadialPoints(idElemental);

    auto napoints = _getNumberOfAngularPoints(idElemental);

    // generate radial grid points

    CM4T2Quadrature rquad(nrpoints, idElemental);

    auto rpoints = rquad.generate();

    auto rrx = rpoints.row(0);

    auto rrw = rpoints.row(1);

    auto rrf = 4.0 * mathconst::pi_value();

    // generate angular grid points

    CLebedevLaikovQuadrature aquad(napoints);

    auto apoints = aquad.generate();

    auto rax = apoints.row(0);

    auto ray = apoints.row(1);

    auto raz = apoints.row(2);

    auto raw = apoints.row(3);

    // atom coordinates

    auto atmxyz = atomCoordinates[idAtomic].coordinates();

    auto atmx = atmxyz[0];

    auto atmy = atmxyz[1];

    auto atmz = atmxyz[2];

    // assemble atom grid points from quadratures

    for (int i = 0; i < nrpoints; i++)
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

#pragma omp simd
        for (int j = 0; j < napoints; j++)
        {
            gx[j] = crx * rax[j] + atmx;

            gy[j] = crx * ray[j] + atmy;

            gz[j] = crx * raz[j] + atmz;

            gw[j] = crw * crx * crx * raw[j];
        }
    }

    // apply partitioning function

    partfunc::ssf(rawGridPoints, minDistance, gridOffset, nrpoints * napoints, atomCoordinates, nAtoms, idAtomic);
}

auto
CGridDriver::_screenRawGridPoints(CDenseMatrix* rawGridPoints) const -> int
{
    // set up pointers to grid points data

    auto gx = rawGridPoints->row(0);

    auto gy = rawGridPoints->row(1);

    auto gz = rawGridPoints->row(2);

    auto gw = rawGridPoints->row(3);

    // loop over grid points

    int cpoints = 0;

    for (int i = 0; i < rawGridPoints->getNumberOfColumns(); i++)
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
