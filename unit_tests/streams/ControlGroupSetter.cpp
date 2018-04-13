//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "ControlGroupSetter.hpp"

CControlGroup getMolXYZGroup()
{
    CControlGroup cg;

    cg.setHeader(CInputLine(std::string(" @molxyz ! molecular geometry")));

    cg.addCommand(CInputLine(std::string("0 1 ! charge multiplicity")));

    cg.addCommand(CInputLine(std::string("Li 0.0 0.0 0.0")));

    cg.addCommand(CInputLine(std::string("H  0.0 0.0 1.2")));

    return cg;
}

CControlGroup getSCFGroup()
{
    CControlGroup cg;

    cg.setHeader(CInputLine(std::string("@Scf")));

    cg.addCommand(CInputLine(std::string("MaxIterations 100")));

    cg.addCommand(CInputLine(std::string("ConvThresholds 1.0e-5 1.0e-6")));

    return cg;
}
