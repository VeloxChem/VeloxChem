//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "DummyFunctions.hpp"

namespace vlxtest {

std::vector<CRecursionTerm>
dummy_func_10(const CRecursionTerm& recurcionTerm)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 7, 2, 3}, 1, 2, 5);

    return std::vector<CRecursionTerm>({rta});
}

std::vector<CRecursionTerm>
dummy_func_01(const CRecursionTerm& recurcionTerm)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 7, 2, 3}, 1, 2, 5);

    return std::vector<CRecursionTerm>({rta});
}

std::vector<CRecursionTerm>
dummy_func_11(const CRecursionTerm& recurcionTerm)
{
    auto rta = recurcionTerm;

    return std::vector<CRecursionTerm>({rta, rta.braShift(-1, 0)});
}

}  // namespace vlxtest
