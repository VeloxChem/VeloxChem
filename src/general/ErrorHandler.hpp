//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef ErrorHandler_hpp
#define ErrorHandler_hpp

#include <string>

namespace errors { // errors namespace

/**
 Prints message and aborts in case of a critical error.
 */
void
assertMsgCritical(const bool         condition,
                  const std::string& message);

} // errors namespace

#endif /* ErrorHandler_hpp */
