//
//                           VELOXCHEM 1.0-RC
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#ifndef ExportOneInts_hpp
#define ExportOneInts_hpp

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>

#include <memory>
#include <string>

#include "CartesianComponents.hpp"
#include "ExportGeneral.hpp"
#include "ExportMath.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"

namespace py = pybind11;

namespace vlx_oneints {  // vlx_oneints namespace
/** Convert a VeloxChem matrix object with Cartesian components, i.e. CElectricDipoleMomentMatrix, to NumPy array.
 *
 * @tparam T type of the VeloxChem matrix object.
 * @tparam cart requested Cartesian component.
 * @param obj VeloxChem matrix object.
 * @return The NumPy array.
 */
template <typename T, cartesians cart>
inline py::array_t<double>
matrix_to_numpy(const T& obj)
{
    return vlx_general::pointer_to_numpy(obj.values(cart), obj.getNumberOfRows(), obj.getNumberOfColumns());
}

/** Convert a VeloxChem matrix object with Cartesian components, i.e. CElectricDipoleMomentMatrix, to NumPy array.
 *
 * @tparam T type of the VeloxChem matrix object.
 * @param obj VeloxChem matrix object.
 * @param cart requested Cartesian component.
 * @return The NumPy array.
 */
template <typename T>
inline py::array_t<double>
matrix_to_numpy(const T& obj, cartesians cart)
{
    return vlx_general::pointer_to_numpy(obj.values(cart), obj.getNumberOfRows(), obj.getNumberOfColumns());
}

/** Common binding code for matrix objects for the various one-electron operators
 *
 * @tparam T type of the VeloxChem matrix object.
 * @param m the pybind11 module object.
 * @param cls_name the name of the class to bind.
 * @return the bound class.
 * @note we return the `py::class_` object to allow for binding of additional
 * methods specific to the class.
 */
template <typename T>
PyClass<T>
bind_operator_matrix(py::module& m, const std::string& cls_name)
{
    return PyClass<T>(m, cls_name.c_str())
        .def(py::init<>())
        .def(py::init<const CDenseMatrix&>())
        .def(py::init<const T&>())
        .def(py::init([](const py::array_t<double>& np) {
            auto mp = vlx_math::CDenseMatrix_from_numpy(np);
            return std::make_shared<T>(*mp);
        }))
        .def("__str__", &T::getString)
        .def(
            "to_numpy",
            [](const T& obj) { return vlx_general::pointer_to_numpy(obj.values(), obj.getNumberOfRows(), obj.getNumberOfColumns()); },
            "Convert operator matrix object to NumPy array")
        .def(py::self == py::self);
}

/** Common binding code for integral driver objects for the various one-electron operators
 *
 * @tparam T type of the VeloxChem integral driver object.
 * @param m the pybind11 module object.
 * @param cls_name the name of the class to bind, e.g. "COverlapIntegralsDriver".
 * @param docstring the name of the class for documentation, e.g. "overlap integrals".
 * @return the bound class.
 * @note we return the `py::class_` object to allow for binding of additional
 * methods specific to the class.
 */
template <typename T>
PyClass<T>
bind_integrals_driver(py::module& m, const char* cls_name, const std::string& docstring)
{
    using namespace py::literals;
    auto docpreamble = "Compute AO-basis " + docstring;
    return PyClass<T>(m, cls_name)
        .def(py::init(&vlx_general::create<T>), "comm"_a = py::none())
        .def("compute",
             vlx_general::overload_cast_<const CMolecule&, const CMolecularBasis&>()(&T::compute, py::const_),
             (docpreamble + " integrals for given molecule and basis").c_str(),
             "molecule"_a,
             "basis"_a)
        .def("compute",
             vlx_general::overload_cast_<const CMolecule&, const CMolecularBasis&, const CMolecularBasis&>()(&T::compute, py::const_),
             (docpreamble + " integrals for given molecule in mixed basis").c_str(),
             "molecule"_a,
             "bra_basis"_a,
             "ket_basis"_a)
        .def("compute",
             vlx_general::overload_cast_<const CMolecule&, const CMolecule&, const CMolecularBasis&>()(&T::compute, py::const_),
             (docpreamble + " integrals for two molecules in given basis").c_str(),
             "bra_molecule"_a,
             "ket_molecule"_a,
             "basis"_a)
        .def("compute",
             vlx_general::overload_cast_<const CMolecule&, const CMolecule&, const CMolecularBasis&, const CMolecularBasis&>()(&T::compute,
                                                                                                                               py::const_),
             (docpreamble + " integrals for two molecules, each in its own basis").c_str(),
             "bra_molecule"_a,
             "ket_molecule"_a,
             "bra_basis"_a,
             "ket_basis"_a);
}

/**
 Exports classes/functions in src/oneints to python.
 */
void export_oneints(py::module& m);
}  // namespace vlx_oneints

#endif /* ExportOneInts_hpp */
