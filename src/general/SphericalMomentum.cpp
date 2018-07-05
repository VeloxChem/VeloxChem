//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "SphericalMomentum.hpp"

#include <cmath>

CSphericalMomentum::CSphericalMomentum()

    : _angularMomentum(-1)
{

}

CSphericalMomentum::CSphericalMomentum(const int32_t angularMomentum)

    : _angularMomentum(angularMomentum)
{
    // s-type angular momentum
    
    if (_angularMomentum == 0)
    {
        _factors = CMemBlock2D<double>(std::vector<double>({1.0}),
                                       std::vector<int32_t>({1}));
        
        _indexes = CMemBlock2D<int32_t>(std::vector<int32_t>({0}),
                                        std::vector<int32_t>({1}));
    }
    
    // p-type angular momentum

    if (_angularMomentum == 1)
    {
        // order: p_-1, p_0, p_1 i.e. p_y, p_z, p_x
        
        _factors = CMemBlock2D<double>({1.0, 1.0, 1.0}, {1, 1, 1});
        
        _indexes = CMemBlock2D<int32_t>({1, 2, 0}, {1, 1, 1});
    }

    // d-type angular momentum
    
    if (_angularMomentum == 2)
    {
        // order: d_-2, d_-1, d_0, d_1, d_2
        
         double f3 = 2.0 * std::sqrt(3.0);
        
        _factors = CMemBlock2D<double>({ f3, f3, -1.0, -1.0, 2.0, f3, 0.5 * f3,
                                        -0.5 * f3}, {1, 1, 3, 1, 2});
        
        _indexes = CMemBlock2D<int32_t>({1, 4, 0, 3, 5, 2, 0, 3},
                                        {1, 1, 3, 1, 2});
    }
    
    // f-type angular momentum

    if (_angularMomentum == 3)
    {
        // order: f_-3, f_-2, f_-1, f_0, f_1, f_2, f_3

        double f5 = std::sqrt(2.5);
        
        double f15 = 2.0 * std::sqrt(15.0);
        
        double f3 = std::sqrt(1.5);
        
        _factors = CMemBlock2D<double>({3.0 * f5, -f5, f15, 4.0 * f3, -f3, -f3,
                                        2.0, -3.0, -3.0, 4.0 * f3, -f3, -f3,
                                        0.5 * f15, -0.5 * f15, f5, -3.0 * f5},
                                       {2, 1, 3, 3, 3, 2, 2});
        
        _indexes = CMemBlock2D<int32_t>({1, 6, 4, 8, 1, 6, 9, 2, 7, 5, 0, 3,
                                         2, 7, 0, 3}, {2, 1, 3, 3, 3, 2, 2});
    }

    // g-type angular momentum
    
    if (_angularMomentum == 4)
    {
        // order: g_-4, g_-3, g_-2, g_-1, g_0, g_1, g_2, g_3, g_4

        double f35 = 4.0 * std::sqrt(35);
        
        double f17 = 4.0 * std::sqrt(17.5);
        
        double f5 = 4.0 * std::sqrt(5.0);
        
        double f2 = 4.0 * std::sqrt(2.5);
        
        _factors = CMemBlock2D<double>({ f35, -f35, 3.0 * f17, -f17, 6.0 * f5,
                                        -f5, -f5, 4.0 * f2, -3.0 * f2, -3.0 * f2,
                                         8.0, 3.0, 3.0, 6.0, -24.0, -24.0,
                                         4.0 * f2, -3.0 * f2, -3.0 * f2,
                                         3.0 * f5, -3.0 * f5, -0.5 * f5, 0.5 * f5,
                                         f17, -3.0 * f17, 0.25 * f35, 0.25 * f35,
                                         -1.50 * f35},
                                       {2, 2, 3, 3, 6, 3, 4, 2, 3});
        
        _indexes = CMemBlock2D<int32_t>({1, 6, 4, 11, 8, 1, 6, 13, 4, 11, 14, 0,
                                         10, 3, 5, 12, 9, 2, 7, 5, 12, 0, 10, 2,
                                         7, 0, 10, 3},
                                        {2, 2, 3, 3, 6, 3, 4, 2, 3});
    }

    // TODO: implement I, H type functions
}

CSphericalMomentum::CSphericalMomentum(const CSphericalMomentum& source)

    : _angularMomentum(source._angularMomentum)

    , _factors(source._factors)

    , _indexes(source._indexes)
{

}

CSphericalMomentum::CSphericalMomentum(CSphericalMomentum&& source) noexcept

    : _angularMomentum(std::move(source._angularMomentum))

    , _factors(std::move(source._factors))

    , _indexes(std::move(source._indexes))
{

}

CSphericalMomentum::~CSphericalMomentum()
{

}

CSphericalMomentum&
CSphericalMomentum::operator=(const CSphericalMomentum& source)
{
    if (this == &source) return *this;

    _angularMomentum = source._angularMomentum;

    _factors = source._factors;
    
    _indexes = source._indexes;

    return *this;
}

CSphericalMomentum&
CSphericalMomentum::operator=(CSphericalMomentum&& source) noexcept
{
    if (this == &source) return *this;

    _angularMomentum = std::move(source._angularMomentum);

    _factors = std::move(source._factors);
    
    _indexes = std::move(source._indexes);

    return *this;
}

bool
CSphericalMomentum::operator==(const CSphericalMomentum& other) const
{
    if (this == &other) return true;

    if (_angularMomentum != other._angularMomentum) return false;

    if (_factors != other._factors) return false;

    if (_indexes != other._indexes) return false;

    return true;
}

bool
CSphericalMomentum::operator!=(const CSphericalMomentum& other) const
{
    return !( (*this) == other);
}

int32_t
CSphericalMomentum::getAngularMomentum() const
{
    return _angularMomentum;
}

int32_t
CSphericalMomentum::getNumberOfComponents() const
{
    return _factors.blocks();
}

const double*
CSphericalMomentum::getFactors(const int32_t iComponent) const
{
    return _factors.data(iComponent);
}

int32_t
CSphericalMomentum::getNumberOfFactors(const int32_t iComponent) const
{
    return _factors.size(iComponent);
}

const int32_t*
CSphericalMomentum::getIndexes(const int32_t iComponent) const
{
    return _indexes.data(iComponent);
}

std::ostream&
operator<<(      std::ostream&       output,
           const CSphericalMomentum& source)
{
    output << std::endl;
    
    output << "[CSphericalMomentum (Object):" << &source << "]" << std::endl;

    output << "_angularMomentum: " << source._angularMomentum << std::endl;

    output << "_factors: " << source._factors << std::endl;
    
    output << "_indexes: " << source._indexes << std::endl;

    return output;
}
