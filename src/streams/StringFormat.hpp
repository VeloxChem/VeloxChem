//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef StringFormat_hpp
#define StringFormat_hpp

#include <string>
#include <cstdint>

#include "FmtType.hpp"

namespace fstr { // start of fstr namespace

    /**
     Creates uppercased string from string.

     @param source the string.
     @return the uppercased string.
     */
    std::string upcase(const std::string& source);

    /**
     Creates formatted string with requested width from string.

     @param source the string.
     @param width the width of formatted string.
     @param aligment the alignment of formatted string.
     @return the formatted string.
     */
    std::string format(const std::string& source, const size_t width,
                       const fmt aligment);

    /**
     Creates formatted string with requested width from C string.

     @param source the C string.
     @param width the width of formatted string.
     @param aligment the alignment of formatted string.
     @return the formatted string.
     */
    std::string to_string(const char* source, const size_t width,
                          const fmt aligment);

    /**
     Creates centered string with requested width from C string.

     @param source the C string.
     @param width the width of centered string.
     @return the centered string.
     */
    std::string to_string(const char* source, const size_t width);

    /**
     Creates formatted string with requested width from real number.

     @param source the real number.
     @param presicion the conversion precision in decimal places.
     @param width the width of formatted string.
     @param aligment the alignment of formatted string.
     @return the formatted string.
     */
    std::string to_string(const double source, const size_t presicion,
                          const size_t width, const fmt aligment);

    /**
     Creates formatted string with requested width from integer number.

     @param source the integer number.
     @param width the width of formatted string.
     @param aligment the alignment of formatted string.
     @return the formatted string.
     */
    std::string to_string(const int32_t source, const size_t width,
                          const fmt aligment);

    /**
     Creates formatted string with requested width from unsigned integer number.

     @param source the insigned integer number.
     @param width the width of formatted string.
     @param aligment the alignment of formatted string.
     @return the formatted string.
     */
    std::string to_string(const size_t source, const size_t width,
                          const fmt aligment);

    /**
     Converts angular momentum label to angular momentum quantum number.
     Supported angular momentum: from S to I.

     @param label the angular momentum label.
     @return the angular momentum quantum number.
     */
    int32_t to_AngularMomentum(const std::string& label);

    /**
     Converts angular momentum quantum number to angular momentum label.
     Supported angular momentum: from S to I.

     @param angularmomentum the angular momentum quantum number.
     @return the angular momentum label.
     */
    std::string to_AngularMomentum(const int32_t angularmomentum);

} // end of fstr namespace

#endif /* StringFormat_hpp */
