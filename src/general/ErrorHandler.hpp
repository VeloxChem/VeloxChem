//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

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
