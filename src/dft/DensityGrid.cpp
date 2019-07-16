//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "DensityGrid.hpp"

CDensityGrid::CDensityGrid()

    : _gridType(dengrid::undefined)

    , _nDensityMatrices(0)

    , _densityValues(CMemBlock2D<double>())
{
    
}

CDensityGrid::CDensityGrid(const CMemBlock2D<double>& densityValues, const dengrid gridType)

    : _gridType(gridType)

    , _nDensityMatrices(1)

    , _densityValues(densityValues)
{
    
}

CDensityGrid::CDensityGrid(const int32_t nGridPoints, const int32_t nDensityMatrices, const xcfun xcFuncType, const dengrid gridType)
{
    _gridType = gridType;
    
    _nDensityMatrices = nDensityMatrices;
    
    int32_t ncomp = 0;
    
    if (xcFuncType == xcfun::lda) ncomp = (_gridType == dengrid::ab) ? 2 : 1;
    
    if (xcFuncType == xcfun::gga) ncomp = (_gridType == dengrid::ab) ? 11 : 5;
    
    // NOTE: this needs to be checked with mgga functionals implementation
    
    if (xcFuncType == xcfun::mgga) ncomp = (_gridType == dengrid::ab) ? 13 : 6;
    
    _densityValues = CMemBlock2D<double>(nGridPoints, _nDensityMatrices * ncomp);
}

CDensityGrid::CDensityGrid(const CDensityGrid& source)

    : _gridType(source._gridType)

    , _nDensityMatrices(source._nDensityMatrices)

    , _densityValues(source._densityValues)
{
}

CDensityGrid::CDensityGrid(CDensityGrid&& source) noexcept

    : _gridType(std::move(source._gridType))

    , _nDensityMatrices(std::move(source._nDensityMatrices))

    , _densityValues(std::move(source._densityValues))
{
}

CDensityGrid::~CDensityGrid()
{
}

CDensityGrid&
CDensityGrid::operator=(const CDensityGrid& source)
{
    if (this == &source) return *this;
    
    _gridType = source._gridType;
    
    _nDensityMatrices = source._nDensityMatrices;
    
    _densityValues = source._densityValues;
    
    return *this;
}

CDensityGrid&
CDensityGrid::operator=(CDensityGrid&& source) noexcept
{
    if (this == &source) return *this;
    
    _gridType = std::move(source._gridType);
    
    _nDensityMatrices = std::move(source._nDensityMatrices);
    
    _densityValues = std::move(source._densityValues);
    
    return *this;
}

bool
CDensityGrid::operator==(const CDensityGrid& other) const
{
    if (_gridType != other._gridType) return false;
    
    if (_nDensityMatrices != other._nDensityMatrices) return false;
    
    if (_densityValues != other._densityValues) return false;
    
    return true;
}

bool
CDensityGrid::operator!=(const CDensityGrid& other) const
{
    return !(*this == other);
}

void
CDensityGrid::zero()
{
    _densityValues.zero(); 
}

void
CDensityGrid::slice(const int32_t nGridPoints)
{
    if (nGridPoints < getNumberOfGridPoints())
    {
        _densityValues = _densityValues.slice(0, nGridPoints);
    }
}

 void
CDensityGrid::updateBetaDensities()
{
    if (_gridType != dengrid::ab) return;
    
    auto ngpoints = getNumberOfGridPoints();
    
    if ((2 * _nDensityMatrices) == _densityValues.blocks())
    {
        for (int32_t i = 0; i < _nDensityMatrices; i++)
        {
            auto rhoa = alphaDensity(i);
            
            auto rhob = betaDensity(i);
            
            #pragma omp simd aligned(rhoa, rhob: VLX_ALIGN)
            for (int32_t j = 0; j < ngpoints; j++)
            {
                rhob[j] = rhoa[j];
            }
        }
    }
    
    if ((11 * _nDensityMatrices) == _densityValues.blocks())
    {
        for (int32_t i = 0; i < _nDensityMatrices; i++)
        {
            auto rhoa = alphaDensity(i);
            
            auto rhob = betaDensity(i);
            
            auto grada_x = alphaDensityGradientX(i);
            
            auto gradb_x = betaDensityGradientX(i);
            
            auto grada_y = alphaDensityGradientY(i);
            
            auto gradb_y = betaDensityGradientY(i);
            
            auto grada_z = alphaDensityGradientZ(i);
            
            auto gradb_z = betaDensityGradientZ(i);
            
            #pragma omp simd aligned(rhoa, rhob, grada_x, gradb_x, grada_y, gradb_y, grada_z, gradb_z: VLX_ALIGN)
            for (int32_t j = 0; j < ngpoints; j++)
            {
                rhob[j] = rhoa[j];
                
                gradb_x[j] = grada_x[j];
                
                gradb_y[j] = grada_y[j];
                
                gradb_z[j] = grada_z[j];
            }
        }
    }
}


void
CDensityGrid::computeDensityNorms()
{
    auto ngpoints = getNumberOfGridPoints();
    
    if (_gridType == dengrid::ab)
    {
        if ((11 * _nDensityMatrices) == _densityValues.blocks())
        {
            for (int32_t i = 0; i < _nDensityMatrices; i++)
            {
                auto grada_x = alphaDensityGradientX(i);
                
                auto gradb_x = betaDensityGradientX(i);
                
                auto grada_y = alphaDensityGradientY(i);
                
                auto gradb_y = betaDensityGradientY(i);
                
                auto grada_z = alphaDensityGradientZ(i);
                
                auto gradb_z = betaDensityGradientZ(i);
                
                auto grada = alphaDensityGradient(i);
                
                auto gradb = betaDensityGradient(i);
                
                auto gradab = mixedDensityGradient(i);
                
                #pragma omp simd aligned(grada_x, gradb_x, grada_y, gradb_y, grada_z, gradb_z: VLX_ALIGN)
                for (int32_t j = 0; j < ngpoints; j++)
                {
                    grada[j] = std::sqrt(grada_x[j] * grada_x[j] + grada_y[j] * grada_y[j] + grada_z[j] * grada_z[j]);
                    
                    gradb[j] = std::sqrt(gradb_x[j] * gradb_x[j] + gradb_y[j] * gradb_y[j] + gradb_z[j] * gradb_z[j]);
            
                    gradab[j] = grada_x[j] * gradb_x[j] + grada_y[j] * gradb_y[j] + grada_z[j] * gradb_z[j];
                }
            }
        }
    }
}

int32_t
CDensityGrid::getNumberOfGridPoints() const
{
    return _densityValues.size(0);
}

int32_t
CDensityGrid::getNumberOfDensityMatrices() const
{
    return _nDensityMatrices;
}

dengrid
CDensityGrid::getDensityGridType() const
{
    return _gridType; 
}

const double*
CDensityGrid::alphaDensity(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::lima) return nullptr;
    
    return _densityValues.data(iDensityMatrix);
}

double*
CDensityGrid::alphaDensity(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::lima) return nullptr;
    
    return _densityValues.data(iDensityMatrix);
}

const double*
CDensityGrid::betaDensity(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(_nDensityMatrices + iDensityMatrix);
    
    if (_gridType == dengrid::lima) return _densityValues.data(iDensityMatrix);
    
    return nullptr;
}

double*
CDensityGrid::betaDensity(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(_nDensityMatrices + iDensityMatrix);
    
    if (_gridType == dengrid::lima) return _densityValues.data(iDensityMatrix);
    
    return nullptr;
}

const double*
CDensityGrid::alphaDensityGradient(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(2 * _nDensityMatrices + iDensityMatrix);
    
    if (_gridType == dengrid::limb) return _densityValues.data(_nDensityMatrices + iDensityMatrix);
    
    return nullptr;
}

double*
CDensityGrid::alphaDensityGradient(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(2 * _nDensityMatrices + iDensityMatrix);
    
    if (_gridType == dengrid::limb) return _densityValues.data(_nDensityMatrices + iDensityMatrix);
    
    return nullptr;
}

const double*
CDensityGrid::betaDensityGradient(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(3 * _nDensityMatrices + iDensityMatrix);
    
    if (_gridType == dengrid::lima) return _densityValues.data(_nDensityMatrices + iDensityMatrix);
    
    return nullptr;
}

double*
CDensityGrid::betaDensityGradient(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(3 * _nDensityMatrices + iDensityMatrix);
    
    if (_gridType == dengrid::lima) return _densityValues.data(_nDensityMatrices + iDensityMatrix);
    
    return nullptr;
}

const double*
CDensityGrid::mixedDensityGradient(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(4 * _nDensityMatrices + iDensityMatrix);
    
    return nullptr;
}

double*
CDensityGrid::mixedDensityGradient(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(4 * _nDensityMatrices + iDensityMatrix);
    
    return nullptr;
}

const double*
CDensityGrid::alphaDensityGradientX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(5 * _nDensityMatrices + iDensityMatrix);
    
    if (_gridType == dengrid::limb) return _densityValues.data(2 * _nDensityMatrices + iDensityMatrix);
    
    return nullptr;
}

double*
CDensityGrid::alphaDensityGradientX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(5 * _nDensityMatrices + iDensityMatrix);
    
    if (_gridType == dengrid::limb) return _densityValues.data(2 * _nDensityMatrices + iDensityMatrix);
    
    return nullptr;
}

const double*
CDensityGrid::alphaDensityGradientY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(6 * _nDensityMatrices + iDensityMatrix);
    
    if (_gridType == dengrid::limb) return _densityValues.data(3 * _nDensityMatrices + iDensityMatrix);
    
    return nullptr;
}

double*
CDensityGrid::alphaDensityGradientY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(6 * _nDensityMatrices + iDensityMatrix);
    
    if (_gridType == dengrid::limb) return _densityValues.data(3 * _nDensityMatrices + iDensityMatrix);
    
    return nullptr;
}

const double*
CDensityGrid::alphaDensityGradientZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(7 * _nDensityMatrices + iDensityMatrix);
    
    if (_gridType == dengrid::limb) return _densityValues.data(4 * _nDensityMatrices + iDensityMatrix);
    
    return nullptr;
}

double*
CDensityGrid::alphaDensityGradientZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(7 * _nDensityMatrices + iDensityMatrix);
    
    if (_gridType == dengrid::limb) return _densityValues.data(4 * _nDensityMatrices + iDensityMatrix);
    
    return nullptr;
}

const double*
CDensityGrid::betaDensityGradientX(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(8 * _nDensityMatrices + iDensityMatrix);
    
    if (_gridType == dengrid::lima) return _densityValues.data(2 * _nDensityMatrices + iDensityMatrix);
    
    return nullptr;
}

double*
CDensityGrid::betaDensityGradientX(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(8 * _nDensityMatrices + iDensityMatrix);
    
    if (_gridType == dengrid::lima) return _densityValues.data(2 * _nDensityMatrices + iDensityMatrix);
    
    return nullptr;
}

const double*
CDensityGrid::betaDensityGradientY(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(9 * _nDensityMatrices + iDensityMatrix);
    
    if (_gridType == dengrid::lima) return _densityValues.data(3 * _nDensityMatrices + iDensityMatrix);
    
    return nullptr;
}

double*
CDensityGrid::betaDensityGradientY(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(9 * _nDensityMatrices + iDensityMatrix);
    
    if (_gridType == dengrid::lima) return _densityValues.data(3 * _nDensityMatrices + iDensityMatrix);
    
    return nullptr;
}

const double*
CDensityGrid::betaDensityGradientZ(const int32_t iDensityMatrix) const
{
    if (_gridType == dengrid::ab) return _densityValues.data(10 * _nDensityMatrices + iDensityMatrix);
    
    if (_gridType == dengrid::lima) return _densityValues.data(4 * _nDensityMatrices + iDensityMatrix);
    
    return nullptr;
}

double*
CDensityGrid::betaDensityGradientZ(const int32_t iDensityMatrix)
{
    if (_gridType == dengrid::ab) return _densityValues.data(10 * _nDensityMatrices + iDensityMatrix);
    
    if (_gridType == dengrid::lima) return _densityValues.data(4 * _nDensityMatrices + iDensityMatrix);
    
    return nullptr;
}

void
CDensityGrid::setScreenedGrids(      std::vector<CDensityGrid>&   densityGrids,
                                     std::vector<CMolecularGrid>& molecularGrids,
                               const double                       densityThreshold,
                               const xcfun                        xcFuncType) const
{
    // valid only for dengrid::ab
    
    if (_gridType != dengrid::ab) return;
    
    // set up number of grid points
    
    auto npoints = getNumberOfGridPoints();
    
    // local density approximation
    
    if (xcFuncType == xcfun::lda)
    {
        for (int32_t i = 0; i < _nDensityMatrices; i++)
        {
            // set up pointers to source density
            
            auto srhoa = alphaDensity(i);
            
            auto srhob = betaDensity(i);
            
            // set up pointers to destination density
            
            auto drhoa = densityGrids[i].alphaDensity(i);
            
            auto drhob = densityGrids[i].betaDensity(i);
            
            // set up pointers to molecular grid data
            
            auto gx = molecularGrids[i].getCoordinatesX();
            
            auto gy = molecularGrids[i].getCoordinatesY();
            
            auto gz = molecularGrids[i].getCoordinatesZ();
            
            auto gw = molecularGrids[i].getWeights();
            
            auto ipoint = 0;
            
            for (int32_t j = 0; j < npoints; j++)
            {
                if (std::fabs(srhoa[j] + srhob[i]) > densityThreshold)
                {
                    // density data
                    
                    drhoa[ipoint] = srhoa[j];
                    
                    drhob[ipoint] = srhob[j];
                    
                    // molecula grid data
                    
                    gx[ipoint] = gx[j];
                    
                    gy[ipoint] = gy[j];
                    
                    gz[ipoint] = gz[j];
                    
                    gw[ipoint] = gw[j];
                    
                    // update counter
                    
                    ipoint++;
                }
            }
            
            // slice grids
            
            densityGrids[i].slice(ipoint);
            
            molecularGrids[i].slice(ipoint);
        }
    }
    
    // general gradient density approximation
    
    if (xcFuncType == xcfun::gga)
    {
        for (int32_t i = 0; i < _nDensityMatrices; i++)
        {
            // set up pointers to source density
            
            auto srhoa = alphaDensity(i);
            
            auto srhob = betaDensity(i);
            
            auto sgrada = alphaDensityGradient(i);
            
            auto sgradb = betaDensityGradient(i);
            
            auto sgradab = mixedDensityGradient(i);
            
            // set up pointers to destination density
            
            auto drhoa = densityGrids[i].alphaDensity(i);
            
            auto drhob = densityGrids[i].betaDensity(i);
            
            auto dgrada = densityGrids[i].alphaDensityGradient(i);
            
            auto dgradb = densityGrids[i].betaDensityGradient(i);
            
            auto dgradab = densityGrids[i].mixedDensityGradient(i);
            
            // set up pointers to molecular grid data
            
            auto gx = molecularGrids[i].getCoordinatesX();
            
            auto gy = molecularGrids[i].getCoordinatesY();
            
            auto gz = molecularGrids[i].getCoordinatesZ();
            
            auto gw = molecularGrids[i].getWeights();
            
            auto ipoint = 0;
            
            for (int32_t j = 0; j < npoints; j++)
            {
                if (std::fabs(srhoa[j] + srhob[i]) > densityThreshold)
                {
                    // density data
                    
                    drhoa[ipoint] = srhoa[j];
                    
                    drhob[ipoint] = srhob[j];
                    
                    dgrada[ipoint] = sgrada[j];
                    
                    dgradb[ipoint] = sgradb[j];
                    
                    dgradab[ipoint] = sgradab[j];
                    
                    // molecula grid data
                    
                    gx[ipoint] = gx[j];
                    
                    gy[ipoint] = gy[j];
                    
                    gz[ipoint] = gz[j];
                    
                    gw[ipoint] = gw[j];
                    
                    // update counter
                    
                    ipoint++;
                }
            }
            
            // slice grids
            
            densityGrids[i].slice(ipoint);
            
            molecularGrids[i].slice(ipoint);
        }
    }
}

std::ostream&
operator<<(std::ostream& output, const CDensityGrid& source)
{
    output << std::endl;
    
    output << "[CDensityGrid (Object):" << &source << "]" << std::endl;
    
    output << "_gridType: " << to_string(source._gridType) << std::endl;
    
    output << "_densityValues: " << std::endl;
    
    output << source._densityValues << std::endl;
    
    return output;
}
