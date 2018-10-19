//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <boost/python.hpp>
#include <memory>
#include <string>

#include "AppManager.hpp"
#include "ExportControllers.hpp"

namespace bp = boost::python;

namespace bp_controllers { // bp_controllers namespace

// Exports classes/functions in src/controllers to python

void export_controllers()
{
    // CAppManager class

    bp::class_< CAppManager, std::shared_ptr<CAppManager> >
        (
            "AppManager",
            bp::init<
                const std::string&,
                const std::string&
                >()
        )
        .def("execute", &CAppManager::execute)
        .def("get_state", &CAppManager::getState)
    ;
}

} // bp_controllers namespace
