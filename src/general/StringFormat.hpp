//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef StringFormat_hpp
#define StringFormat_hpp

#include <array>
#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

#include "FmtType.hpp"

namespace fstr {  // fstr namespace

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
std::string format(const std::string& source, const size_t width, const fmt aligment);

/**
 Creates formatted string with requested width from C string.

 @param source the C string.
 @param width the width of formatted string.
 @param aligment the alignment of formatted string.
 @return the formatted string.
 */
std::string to_string(const char* source, const size_t width, const fmt aligment);

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
std::string to_string(const double source, const size_t presicion, const size_t width, const fmt aligment);

/**
 Creates string with requested precision from real number.

 @param source the real number.
 @param presicion the conversion precision in decimal places.
 @return the formatted string.
 */
std::string to_string(const double source, const size_t presicion);

/**
 Creates formatted string with requested width from integer number.

 @param source the integer number.
 @param width the width of formatted string.
 @param aligment the alignment of formatted string.
 @return the formatted string.
 */
std::string to_string(const int32_t source, const size_t width, const fmt aligment);

/**
 Creates formatted string with requested width from unsigned integer number.

 @param source the insigned integer number.
 @param width the width of formatted string.
 @param aligment the alignment of formatted string.
 @return the formatted string.
 */
std::string to_string(const size_t source, const size_t width, const fmt aligment);

/**
 Creates formatted string from boolean.

 @param source the boolean.
 @return the formatted string.
 */
std::string to_string(const bool source);

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

template <typename T>
auto
stream_collection(const T& coll) -> std::string
{
    std::ostringstream os;
    bool               first = true;
    os << "[";
    for (auto elem : coll)
    {
        if (!first) os << ", ";
        os << elem;
        first = false;
    }
    os << "]";
    return os.str();
}
}  // namespace fstr

template <typename T, size_t D>
auto
operator<<(std::ostream& os, const std::array<T, D>& coll) -> std::ostream&
{
    return (os << fstr::stream_collection(coll));
}

template <typename T>
auto
operator<<(std::ostream& os, const std::vector<T>& coll) -> std::ostream&
{
    return (os << fstr::stream_collection(coll));
}

#endif /* StringFormat_hpp */
