//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <boost/python.hpp>

#include <memory>
#include <vector>
#include <string>

#include "MpiFunc.hpp"
#include "AppManager.hpp"

namespace bp = boost::python;

// ==> boost python helper function <==
// for creating CAppManager object

static std::shared_ptr<CAppManager>
CAppManager_create(std::string input_string, std::string output_string)
{
    int argc = 3;

    char* argv[argc];

    std::vector<std::string> inputs ({std::string("exe"), input_string, output_string});

    for (int i = 0; i < argc; i++)
    {
        const char* text = inputs[i].c_str();

        argv[i] = (char*)malloc(sizeof(char) * (strlen(text) + 1));

        memset(argv[i], '\0', sizeof(char) * (strlen(text) + 1));

        memcpy(argv[i], text, sizeof(char) * strlen(text));
    }

    mpi::init(argc, argv);

    return std::shared_ptr<CAppManager>(new CAppManager (argc, argv));
}

// ==> boost python <==
// functions and classes

void export_controllers()
{
    // CAppManager class
    // Note: "create" is a static method that returns a CAppManager object

    bp::class_< CAppManager, std::shared_ptr<CAppManager> >
        (
            "CAppManager",
            bp::init<int, char**>()
        )
        .def("create",    &CAppManager_create)
        .def("execute",   &CAppManager::execute)
        .def("get_state", &CAppManager::getState)
        .staticmethod("create")
    ;
}
