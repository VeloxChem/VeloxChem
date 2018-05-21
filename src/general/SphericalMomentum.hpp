//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef SphericalMomentum_hpp
#define SphericalMomentum_hpp

#include <cstdint>

#include "MemBlock2D.hpp"

/**
 Class CSphericalMomentum stores spherical angular momentum data.
 
 @author Z. Rinkevicius
 */
class CSphericalMomentum
{
    /**
     The angular momentum.
     */
    int32_t _angularMomentum;              

    /**
     The tranformation factors for Carterian angular momentum.
     */
    CMemBlock2D<double> _factors; 

    /**
     The tranformation indexes map of Carterian angular momentum components.
     */
    CMemBlock2D<int32_t> _indexes;
    
public:

    /**
     Creates an empty spherical momentum object.
     */
    CSphericalMomentum();
    
    /**
     Creates a spherical momentum object.
     
     @param angularMomentum the angular momentum.
     */
    CSphericalMomentum(const int32_t angularMomentum);

    /**
     Creates a spherical momentum object by copying other spherical momentum
     object.
     
     @param source the spherical momentum object.
     */
    CSphericalMomentum(const CSphericalMomentum& source);

    /**
     Creates a spherical momentum object by moving other spherical momentum
     object.
     
     @param source the spherical momentum object.
     */
    CSphericalMomentum(CSphericalMomentum&& source) noexcept;

    /**
     Destroys a spherical momentum object.
     */
    ~CSphericalMomentum();

    /**
     Assigns a spherical momentum object by copying other spherical momentum
     object.
     
     @param source the spherical momentum object.
     */
    CSphericalMomentum& operator=(const CSphericalMomentum& source);

    /**
     Assigns a spherical momentum object by moving other spherical momentum
     object.
     
     @param source the spherical momentum object.
     */
    CSphericalMomentum& operator=(CSphericalMomentum&& source) noexcept;

    /**
     Compares spherical momentum object with other spherical momentum object.
     
     @param other the spherical momentum object.
     @return true if spherical momentum objects are equal, false otherwise.
     */
    bool operator==(const CSphericalMomentum& other) const;

    /**
     Compares spherical momentum object with other spherical momentum object.
     
     @param other the spherical momentum object.
     @return true if spherical momentum objects are not equal, false otherwise.
     */
    bool operator!=(const CSphericalMomentum& other) const;

    /**
     Gets angular momentum of spherical momentum.

     @return the angular momentum.
     */
    int32_t getAngularMomentum() const;

    /**
     Gets number of components of spherical momentum.

     @return the number of angulat momentum components.
     */
    int32_t getNumberOfComponents() const;

    /**
     Gets constant pointer to tranformation factors for Cartesian angular
     momentum component.

     @param iComponent the angular momentum Cartesian component.
     @return the vector of transformation factors.
     */
    const double* getFactors(const int32_t iComponent) const;
    
    /**
     Gets number of transformation factors for Cartesian angular momentum
     component.

     @param iComponent the angular momentum Cartesian component.
     @return the number of transformation factors.
     */
    int32_t getNumberOfFactors(const int32_t iComponent) const;
    
    
    /**
     Gets constant pointer to tranformation indexes map for Cartesian angular
     momentum component.
     
     @param iComponent the angular momentum Cartesian component.
     @return the vector of transformation indexes.
     */
    const int32_t* getIndexes(const int32_t iComponent) const;

    /**
     Converts spherical momentum object to text output and insert it into output
     text stream.
     
     @param output the output text stream.
     @param source the spherical momentum object.
     */
    friend std::ostream& operator<<(      std::ostream&       output,
                                    const CSphericalMomentum& source);
};

#endif /* SphericalMomentum_hpp */
