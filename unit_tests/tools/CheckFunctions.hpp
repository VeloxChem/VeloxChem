//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef CheckFunctions_hpp
#define CheckFunctions_hpp

#include <cstdint>
#include <vector>

namespace vlxtest {

    void compare(const std::vector<double>& aVector,
                 const double*              bVector);
    
    void compare(const std::vector<double>& aVector,
                 const double*              bVector,
                 const double               threshod);
    
    
    void compare(const std::vector<int32_t>& aVector,
                 const int32_t*              bVector);

    void compare(const int32_t* aVector,
                 const int32_t* bVector,
                 const int32_t  nElements);

    void compare(const double* aVector,
                 const double* bVector,
                 const int32_t nElements);

    void compare(const std::vector<double>& aVector,
                 const std::vector<double>& bVector);

    void compare(const std::vector<int32_t>& aVector,
                 const std::vector<int32_t>& bVector);

    void checkNorm(const double* aVector,
                   const int32_t nElements);
}

#endif /* CheckFunctions_hpp */
