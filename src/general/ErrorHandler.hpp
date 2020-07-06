//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef ErrorHandler_hpp
#define ErrorHandler_hpp

#include <string>

namespace errors {  // errors namespace

/**
 Prints message and aborts in case of a critical error.
 */
void assertMsgCritical(const bool condition, const std::string& message);

}  // namespace errors

#endif /* ErrorHandler_hpp */
