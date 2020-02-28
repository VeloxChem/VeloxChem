//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef DummyFunctions_hpp
#define DummyFunctions_hpp

#include <vector>

#include "RecursionTerm.hpp"

namespace vlxtest {

std::vector<CRecursionTerm> dummy_func_10(const CRecursionTerm& recurcionTerm);

std::vector<CRecursionTerm> dummy_func_01(const CRecursionTerm& recurcionTerm);

std::vector<CRecursionTerm> dummy_func_11(const CRecursionTerm& recurcionTerm);
}  // namespace vlxtest

#endif /* DummyFunctions_hpp */
