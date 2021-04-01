//
//                           VELOXCHEM 1.0-RC
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

#include "VisualizationDriver.hpp"

#include <cmath>
#include <cstring>

#include <mpi.h>

#include "BasisFunction.hpp"
#include "CubicGrid.hpp"
#include "ErrorHandler.hpp"
#include "MemBlock.hpp"
#include "MolecularBasis.hpp"
#include "MpiFunc.hpp"
#include "SphericalMomentum.hpp"
#include "StringFormat.hpp"

CVisualizationDriver::CVisualizationDriver(MPI_Comm comm)
{
    _locRank = mpi::rank(comm);

    _locNodes = mpi::nodes(comm);

    _locComm = comm;
}

CVisualizationDriver::~CVisualizationDriver()
{
}

std::vector<std::vector<int32_t>>
CVisualizationDriver::_buildCartesianAngularMomentum(int32_t angl) const
{
    std::vector<std::vector<int32_t>> lmn;

    // lexical order of Cartesian angular momentum
    // 1S: 0
    // 3P: x,y,z
    // 6D: xx,xy,xz,yy,yz,zz
    // 10F: ...

    for (int32_t i = 0; i <= angl; i++)
    {
        int32_t lx = angl - i;

        for (int32_t j = 0; j <= i; j++)
        {
            int32_t ly = i - j;

            int32_t lz = j;

            lmn.push_back(std::vector<int32_t>({lx, ly, lz}));
        }
    }

    return lmn;
}

std::vector<double>
CVisualizationDriver::_compPhiAtomicOrbitals(const CMolecule&       molecule,
                                             const CMolecularBasis& basis,
                                             const double           xp,
                                             const double           yp,
                                             const double           zp) const
{
    auto natoms = molecule.getNumberOfAtoms();

    auto max_angl = basis.getMolecularMaxAngularMomentum(molecule);

    std::vector<double> phi;

    // azimuthal quantum number: s,p,d,f,...

    for (int32_t aoidx = 0, angl = 0; angl <= max_angl; angl++)
    {
        CSphericalMomentum sphmom(angl);

        auto nsph = sphmom.getNumberOfComponents();

        auto lmn = _buildCartesianAngularMomentum(angl);

        // magnetic quantum number: s,p-1,p0,p+1,d-2,d-1,d0,d+1,d+2,...

        for (int32_t isph = 0; isph < nsph; isph++)
        {
            // prepare Cartesian components

            std::vector<double> fcarts, lx, ly, lz;

            auto ncomp = sphmom.getNumberOfFactors(isph);

            for (int32_t icomp = 0; icomp < ncomp; icomp++)
            {
                fcarts.push_back(sphmom.getFactors(isph)[icomp]);

                auto cartind = sphmom.getIndexes(isph)[icomp];

                lx.push_back(lmn[cartind][0]);

                ly.push_back(lmn[cartind][1]);

                lz.push_back(lmn[cartind][2]);
            }

            // go through atoms

            for (int32_t atomidx = 0; atomidx < natoms; atomidx++)
            {
                // process coordinates

                double rx = xp - molecule.getCoordinatesX()[atomidx];

                double ry = yp - molecule.getCoordinatesY()[atomidx];

                double rz = zp - molecule.getCoordinatesZ()[atomidx];

                double r2 = rx * rx + ry * ry + rz * rz;

                // process atomic orbitals

                auto idelem = molecule.getIdsElemental()[atomidx];

                auto nao = basis.getNumberOfBasisFunctions(idelem, angl);

                auto basisfunc = basis.getBasisFunctions(idelem, angl);

                for (int32_t i = 0; i < nao; i++, aoidx++)
                {
                    double phiao = 0.0;

                    // process primitives

                    auto nprims = basisfunc[i].getNumberOfPrimitiveFunctions();

                    auto exponents = basisfunc[i].getExponents();

                    auto normcoefs = basisfunc[i].getNormalizationFactors();

                    for (int32_t iprim = 0; iprim < nprims; iprim++)
                    {
                        double expon = std::exp(-exponents[iprim] * r2);

                        double coef1 = normcoefs[iprim];

                        // transform from Cartesian to spherical harmonics

                        for (int32_t icomp = 0; icomp < ncomp; icomp++)
                        {
                            double coef2 = coef1 * fcarts[icomp];

                            double powxyz = std::pow(rx, lx[icomp]) * std::pow(ry, ly[icomp]) * std::pow(rz, lz[icomp]);

                            phiao += coef2 * powxyz * expon;
                        }
                    }

                    phi.push_back(phiao);
                }
            }
        }
    }

    return phi;
}

int32_t
CVisualizationDriver::getRank() const
{
    return _locRank;
}

std::vector<std::vector<int32_t>>
CVisualizationDriver::getCountsAndDisplacements(int32_t nx) const
{
    int32_t ave = nx / _locNodes;

    int32_t res = nx % _locNodes;

    std::vector<int32_t> counts;

    for (int32_t i = 0; i < _locNodes; i++)
    {
        if (i < res)
        {
            counts.push_back(ave + 1);
        }
        else
        {
            counts.push_back(ave);
        }
    }

    std::vector<int32_t> displs;

    for (int32_t i = 0, disp = 0; i < _locNodes; i++)
    {
        displs.push_back(disp);

        disp += counts[i];
    }

    return std::vector<std::vector<int32_t>>({counts, displs});
}

void
CVisualizationDriver::compute(CCubicGrid&               grid,
                              const CMolecule&          molecule,
                              const CMolecularBasis&    basis,
                              const CMolecularOrbitals& mo,
                              const int32_t             moidx,
                              const std::string&        mospin) const
{
    // grid information

    auto x0 = grid.originX(), y0 = grid.originY(), z0 = grid.originZ();

    auto dx = grid.stepSizeX(), dy = grid.stepSizeY(), dz = grid.stepSizeZ();

    auto nx = grid.numPointsX(), ny = grid.numPointsY(), nz = grid.numPointsZ();

    // compute local grid on this MPI process

    auto xcntdsp = getCountsAndDisplacements(nx);

    auto xcounts = xcntdsp[0];

    auto xdispls = xcntdsp[1];

    std::vector<double> origin({x0 + dx * xdispls[_locRank], y0, z0});

    std::vector<double> stepsize({dx, dy, dz});

    std::vector<int32_t> npoints({xcounts[_locRank], ny, nz});

    CCubicGrid localgrid(origin, stepsize, npoints);

    compute_omp(localgrid, molecule, basis, mo, moidx, mospin);

    // gather local grids

    std::vector<int32_t> yzcounts, yzdispls;

    for (int32_t i = 0; i < static_cast<int32_t>(xcounts.size()); i++)
    {
        yzcounts.push_back(xcounts[i] * ny * nz);

        yzdispls.push_back(xdispls[i] * ny * nz);
    }

    MPI_Gatherv(localgrid.values(), yzcounts[_locRank], MPI_DOUBLE,
                grid.values(), yzcounts.data(), yzdispls.data(), MPI_DOUBLE,
                mpi::master(), _locComm);
}

void
CVisualizationDriver::compute(CCubicGrid&             grid,
                              const CMolecule&        molecule,
                              const CMolecularBasis&  basis,
                              const CAODensityMatrix& density,
                              const int32_t           denidx,
                              const std::string&      denspin) const
{
    // grid information

    auto x0 = grid.originX(), y0 = grid.originY(), z0 = grid.originZ();

    auto dx = grid.stepSizeX(), dy = grid.stepSizeY(), dz = grid.stepSizeZ();

    auto nx = grid.numPointsX(), ny = grid.numPointsY(), nz = grid.numPointsZ();

    // compute local grid on this MPI process

    auto xcntdsp = getCountsAndDisplacements(nx);

    auto xcounts = xcntdsp[0];

    auto xdispls = xcntdsp[1];

    std::vector<double> origin({x0 + dx * xdispls[_locRank], y0, z0});

    std::vector<double> stepsize({dx, dy, dz});

    std::vector<int32_t> npoints({xcounts[_locRank], ny, nz});

    CCubicGrid localgrid(origin, stepsize, npoints);

    compute_omp(localgrid, molecule, basis, density, denidx, denspin);

    // gather local grids

    std::vector<int32_t> yzcounts, yzdispls;

    for (int32_t i = 0; i < static_cast<int32_t>(xcounts.size()); i++)
    {
        yzcounts.push_back(xcounts[i] * ny * nz);

        yzdispls.push_back(xdispls[i] * ny * nz);
    }

    MPI_Gatherv(localgrid.values(), yzcounts[_locRank], MPI_DOUBLE,
                grid.values(), yzcounts.data(), yzdispls.data(), MPI_DOUBLE,
                mpi::master(), _locComm);
}

void
CVisualizationDriver::compute_omp(CCubicGrid&               grid,
                                  const CMolecule&          molecule,
                                  const CMolecularBasis&    basis,
                                  const CMolecularOrbitals& mo,
                                  const int32_t             moidx,
                                  const std::string&        mospin) const
{
    // grid information

    auto x0 = grid.originX(), y0 = grid.originY(), z0 = grid.originZ();

    auto dx = grid.stepSizeX(), dy = grid.stepSizeY(), dz = grid.stepSizeZ();

    auto nx = grid.numPointsX(), ny = grid.numPointsY(), nz = grid.numPointsZ();

    // sanity check

    std::string erridx("VisualizationDriver.compute: Invalid index of MO");

    std::string errspin("VisualizationDriver.compute: Invalid spin of MO");

    std::string errnao("VisualizationDriver.compute: Inconsistent number of AOs");

    auto morows = mo.getNumberOfRows();

    auto mocols = mo.getNumberOfColumns();

    errors::assertMsgCritical(0 <= moidx && moidx < mocols, erridx);

    bool alphaspin = (fstr::upcase(mospin) == std::string("ALPHA"));

    bool betaspin = (fstr::upcase(mospin) == std::string("BETA"));

    errors::assertMsgCritical(alphaspin || betaspin, errspin);

    auto phi0 = _compPhiAtomicOrbitals(molecule, basis, x0, y0, z0);

    auto nao = static_cast<int32_t>(phi0.size());

    errors::assertMsgCritical(morows == nao, errnao);

    // target MO

    auto mocoefs = alphaspin ? mo.alphaOrbitals() : mo.betaOrbitals();

    // calculate psi on grid points

    #pragma omp parallel for schedule(dynamic)
    for (int32_t ix = 0; ix < nx; ix++)
    {
        double xp = x0 + dx * ix;

        int32_t xstride = ix * ny * nz;

        for (int32_t iy = 0; iy < ny; iy++)
        {
            double yp = y0 + dy * iy;

            int32_t ystride = iy * nz;

            for (int32_t iz = 0; iz < nz; iz++)
            {
                double zp = z0 + dz * iz;

                int32_t zstride = iz;

                auto phi = _compPhiAtomicOrbitals(molecule, basis, xp, yp, zp);

                double psi = 0.0;

                for (int32_t aoidx = 0; aoidx < nao; aoidx++)
                {
                    double mocoef = mocoefs[aoidx * mocols + moidx];

                    psi += mocoef * phi[aoidx];
                }

                int32_t index = xstride + ystride + zstride;

                grid.values()[index] = psi;
            }
        }
    }
}

void
CVisualizationDriver::compute_omp(CCubicGrid&             grid,
                                  const CMolecule&        molecule,
                                  const CMolecularBasis&  basis,
                                  const CAODensityMatrix& density,
                                  const int32_t           denidx,
                                  const std::string&      denspin) const
{
    // grid information

    auto x0 = grid.originX(), y0 = grid.originY(), z0 = grid.originZ();

    auto dx = grid.stepSizeX(), dy = grid.stepSizeY(), dz = grid.stepSizeZ();

    auto nx = grid.numPointsX(), ny = grid.numPointsY(), nz = grid.numPointsZ();

    // sanity check

    std::string erridx("VisualizationDriver.compute: Invalid index of density matrix");

    std::string errspin("VisualizationDriver.compute: Invalid spin of density matrix");

    std::string errnao("VisualizationDriver.compute: Inconsistent number of AOs");

    auto numdens = density.getNumberOfDensityMatrices();

    errors::assertMsgCritical(0 <= denidx && denidx < numdens, erridx);

    bool alphaspin = (fstr::upcase(denspin) == std::string("ALPHA"));

    bool betaspin = (fstr::upcase(denspin) == std::string("BETA"));

    bool diffspin = (fstr::upcase(denspin) == std::string("SPIN"));

    errors::assertMsgCritical(alphaspin || betaspin || diffspin, errspin);

    auto phi0 = _compPhiAtomicOrbitals(molecule, basis, x0, y0, z0);

    const int32_t nao = static_cast<int32_t>(phi0.size());

    auto denrows = density.getNumberOfRows(denidx);

    auto dencols = density.getNumberOfColumns(denidx);

    errors::assertMsgCritical(denrows == nao && dencols == nao, errnao);

    // target density

    CMemBlock<double> rho(nao * nao);

    auto dens_alpha = density.alphaDensity(denidx);

    auto dens_beta = density.betaDensity(denidx);

    if (alphaspin)
    {
        std::memcpy(rho.data(), dens_alpha, nao * nao * sizeof(double));
    }
    else if (betaspin)
    {
        std::memcpy(rho.data(), dens_beta, nao * nao * sizeof(double));
    }
    else if (diffspin)
    {
        for (int32_t p = 0; p < nao * nao; p++)
        {
            rho.data()[p] = dens_alpha[p] - dens_beta[p];
        }
    }

    auto rho_data = rho.data();

    // calculate densities on grid points

    #pragma omp parallel for schedule(dynamic)
    for (int32_t ix = 0; ix < nx; ix++)
    {
        double xp = x0 + dx * ix;

        int32_t xstride = ix * ny * nz;

        for (int32_t iy = 0; iy < ny; iy++)
        {
            double yp = y0 + dy * iy;

            int32_t ystride = iy * nz;

            for (int32_t iz = 0; iz < nz; iz++)
            {
                double zp = z0 + dz * iz;

                int32_t zstride = iz;

                auto phi = _compPhiAtomicOrbitals(molecule, basis, xp, yp, zp);

                double psi = 0.0;

                for (int32_t iao = 0; iao < nao; iao++)
                {
                    for (int32_t jao = 0; jao < nao; jao++)
                    {
                        psi += phi[iao] * rho_data[iao * nao + jao] * phi[jao];
                    }
                }

                int32_t index = xstride + ystride + zstride;

                grid.values()[index] = psi;
            }
        }
    }
}

std::vector<double>
CVisualizationDriver::getMO(const std::vector<std::vector<double>>& coords,
                            const CMolecule&                        molecule,
                            const CMolecularBasis&                  basis,
                            const CMolecularOrbitals&               mo,
                            const int32_t                           moidx,
                            const std::string&                      mospin) const
{
    if (_locRank == mpi::master())
    {
        // sanity check

        std::string errspin("VisualizationDriver.get_mo: invalid spin");

        bool alphaspin = (fstr::upcase(mospin) == std::string("ALPHA"));

        bool betaspin = (fstr::upcase(mospin) == std::string("BETA"));

        errors::assertMsgCritical(alphaspin || betaspin, errspin);

        std::string erridx("VisualizationDriver.get_mo: invalid MO index");

        auto nao = mo.getNumberOfRows();

        auto nmo = mo.getNumberOfColumns();

        errors::assertMsgCritical(0 <= moidx && moidx < nmo, erridx);

        // compute MO

        auto mocoefs = alphaspin ? mo.alphaOrbitals() : mo.betaOrbitals();

        auto npoints = static_cast<int32_t>(coords.size());

        std::vector<double> psi(npoints, 0.0);

        for (int32_t p = 0; p < npoints; p++)
        {
            auto phi = _compPhiAtomicOrbitals(molecule, basis, coords[p][0], coords[p][1], coords[p][2]);

            for (int32_t aoidx = 0; aoidx < nao; aoidx++)
            {
                psi[p] += mocoefs[aoidx * nmo + moidx] * phi[aoidx];
            }
        }

        return psi;
    }

    return std::vector<double>();
}

std::vector<double>
CVisualizationDriver::getDensity(const std::vector<std::vector<double>>& coords,
                                 const CMolecule&                        molecule,
                                 const CMolecularBasis&                  basis,
                                 const CAODensityMatrix&                 density,
                                 const int32_t                           denidx,
                                 const std::string&                      denspin) const
{
    if (_locRank == mpi::master())
    {
        // sanity check

        std::string errspin("VisualizationDriver.get_density: invalid spin");

        bool alphaspin = (fstr::upcase(denspin) == std::string("ALPHA"));

        bool betaspin = (fstr::upcase(denspin) == std::string("BETA"));

        errors::assertMsgCritical(alphaspin || betaspin, errspin);

        std::string erridx("VisualizationDriver.get_density: invalid index for density matrix");

        auto numdens = density.getNumberOfDensityMatrices();

        errors::assertMsgCritical(0 <= denidx && denidx < numdens, erridx);

        // compute density

        auto nao = density.getNumberOfRows(denidx);

        auto dens = alphaspin ? density.alphaDensity(denidx) : density.betaDensity(denidx);

        auto npoints = static_cast<int32_t>(coords.size());

        std::vector<double> psi(npoints, 0.0);

        for (int32_t p = 0; p < npoints; p++)
        {
            auto phi = _compPhiAtomicOrbitals(molecule, basis, coords[p][0], coords[p][1], coords[p][2]);

            for (int32_t iao = 0; iao < nao; iao++)
            {
                for (int32_t jao = 0; jao < nao; jao++)
                {
                    psi[p] += phi[iao] * dens[iao * nao + jao] * phi[jao];
                }
            }
        }

        return psi;
    }

    return std::vector<double>();
}

std::vector<double>
CVisualizationDriver::getOneParticleDensity(const std::vector<std::vector<double>>& coords_1,
                                            const std::vector<std::vector<double>>& coords_2,
                                            const CMolecule&                        molecule,
                                            const CMolecularBasis&                  basis,
                                            const CAODensityMatrix&                 density,
                                            const int32_t                           denidx,
                                            const std::string&                      spin_1,
                                            const std::string&                      spin_2) const
{
    if (_locRank == mpi::master())
    {
        if (fstr::upcase(spin_1) != fstr::upcase(spin_2))
        {
            return std::vector<double>(coords_1.size(), 0.0);
        }

        bool alphaspin = (fstr::upcase(spin_1) == std::string("ALPHA"));

        // compute density

        auto nao = density.getNumberOfRows(denidx);

        auto dens = alphaspin ? density.alphaDensity(denidx) : density.betaDensity(denidx);

        auto npoints = static_cast<int32_t>(coords_1.size());

        std::vector<double> psi(npoints, 0.0);

        for (int32_t p = 0; p < npoints; p++)
        {
            auto phi_1 = _compPhiAtomicOrbitals(molecule, basis, coords_1[p][0], coords_1[p][1], coords_1[p][2]);

            auto phi_2 = _compPhiAtomicOrbitals(molecule, basis, coords_2[p][0], coords_2[p][1], coords_2[p][2]);

            for (int32_t iao = 0; iao < nao; iao++)
            {
                for (int32_t jao = 0; jao < nao; jao++)
                {
                    psi[p] += phi_1[iao] * dens[iao * nao + jao] * phi_2[jao];
                }
            }
        }

        return psi;
    }

    return std::vector<double>();
}

std::vector<double>
CVisualizationDriver::getTwoParticleDensity(const std::vector<std::vector<double>>& coords_1,
                                            const std::vector<std::vector<double>>& coords_2,
                                            const CMolecule&                        molecule,
                                            const CMolecularBasis&                  basis,
                                            const CAODensityMatrix&                 density,
                                            const int32_t                           denidx,
                                            const std::string&                      spin_1,
                                            const std::string&                      spin_2) const
{
    if (_locRank == mpi::master())
    {
        // sanity check

        std::string errspin("VisualizationDriver.get_two_particle_density: invalid spin");

        bool alphaspin_1 = (fstr::upcase(spin_1) == std::string("ALPHA"));

        bool betaspin_1 = (fstr::upcase(spin_1) == std::string("BETA"));

        bool alphaspin_2 = (fstr::upcase(spin_2) == std::string("ALPHA"));

        bool betaspin_2 = (fstr::upcase(spin_2) == std::string("BETA"));

        errors::assertMsgCritical(alphaspin_1 || betaspin_1, errspin);

        errors::assertMsgCritical(alphaspin_2 || betaspin_2, errspin);

        std::string erridx("VisualizationDriver.get_two_particle_density: invalid index for density matrix");

        auto numdens = density.getNumberOfDensityMatrices();

        errors::assertMsgCritical(0 <= denidx && denidx < numdens, erridx);

        // compute density

        auto g_11 = getOneParticleDensity(coords_1, coords_1, molecule, basis, density, denidx, spin_1, spin_1);
                                                                                                              
        auto g_22 = getOneParticleDensity(coords_2, coords_2, molecule, basis, density, denidx, spin_2, spin_2);
                                                                                                              
        auto g_12 = getOneParticleDensity(coords_1, coords_2, molecule, basis, density, denidx, spin_1, spin_2);
                                                                                                              
        auto g_21 = getOneParticleDensity(coords_2, coords_1, molecule, basis, density, denidx, spin_2, spin_1);

        auto npoints = static_cast<int32_t>(coords_1.size());

        std::vector<double> psi(npoints, 0.0);

        for (int32_t p = 0; p < npoints; p++)
        {
            psi[p] = g_11[p] * g_22[p] - g_12[p] * g_21[p];
        }

        return psi;
    }

    return std::vector<double>();
}
