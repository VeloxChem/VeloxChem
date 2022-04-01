//
//                           VELOXCHEM 1.0-RC2
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

#ifndef ExportGeneral_hpp
#define ExportGeneral_hpp

#include <mpi.h>
// see here: https://github.com/mpi4py/mpi4py/issues/19#issuecomment-768143143
#ifdef MSMPI_VER
#define PyMPI_HAVE_MPI_Message 1
#endif
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "ErrorHandler.hpp"
#include "MemBlock.hpp"
#include "MemBlock2D.hpp"
#include "NumaPolicy.hpp"

namespace py = pybind11;

template <typename T>
using PyClass = py::class_<T, std::shared_ptr<T>>;

namespace vlx_general {  // vlx_general namespace

/**
 Gets MPI_Comm pointer from a mpi4py communicator object.

 @param py_comm mpi4py communicator object.
 @return the MPI_Comm pointer.
 */
MPI_Comm* get_mpi_comm(py::object py_comm);

/** Gets shape and strides parameters for array_t CTOR.
 *
 * @tparam T scalar type of array.
 * @tparam Container the object collecting dimensions, _e.g._ `std::vector` or `std::array`.
 * @param dimension the shape of numpy array.
 * @return shape and strides of numpy array.
 */
template <typename T, typename Container>
auto
array_t_params(const Container& dimension) -> std::tuple<std::vector<py::ssize_t>, std::vector<py::ssize_t>>
{
    std::vector<py::ssize_t> shape, strides;

    for (size_t i = 0; i < dimension.size(); i++)
    {
        shape.push_back(static_cast<py::ssize_t>(dimension[i]));

        size_t strd = 1;

        for (size_t j = i + 1; j < dimension.size(); j++)
        {
            strd *= dimension[j];
        }

        strides.push_back(static_cast<py::ssize_t>(strd * sizeof(T)));
    }

    return {shape, strides};
}

/** Gets numpy array from pointer and shape.
 *
 * @tparam T scalar type of array.
 * @tparam Container the object collecting dimensions, _e.g._ `std::vector` or `std::array`.
 * @param ptr pointer to data.
 * @param dimension the shape of numpy array.
 * @return numpy array.
 *
 * @note The return type is rather convoluted, but it boils down to "return an
 * object of `py::array<T>` only if the `dimension` argument implmentes a `size`
 * method and the type `T` is arithmetic".
 */
template <typename T, typename Container>
auto
pointer_to_numpy(const T* ptr, const Container& dimension) -> decltype((void)(dimension.size()), (void)(std::is_arithmetic_v<T>), py::array_t<T>())
{
    static_assert(std::is_arithmetic_v<T>, "NumPy array can only be instantiated with arithmetic types.");
    if (ptr == nullptr || dimension.size() == 0)
    {
        return py::array_t<T>();
    }
    else
    {
        const auto [shape, strides] = array_t_params<T>(dimension);
        return py::array_t<T>(shape, strides, ptr);
    }
}

/** Gets 1d numpy array from double pointer and int32_t dimension.
 *
 * @tparam T scalar type of array.
 * @param ptr pointer to data.
 * @param nElements number of elements.
 * @return numpy array.
 */
template <typename T>
auto
pointer_to_numpy(const T* ptr, int32_t nElements) -> py::array_t<T>
{
    return pointer_to_numpy(ptr, std::vector<int32_t>{nElements});
}

/** Gets 2d numpy array from double pointer and int32_t dimension.
 *
 * @tparam T scalar type of array.
 * @param ptr pointer to data.
 * @param nRows number of rows.
 * @param nColumns number of columns.
 * @return numpy array.
 */
template <typename T>
auto
pointer_to_numpy(const T* ptr, int32_t nRows, int32_t nColumns) -> py::array_t<T>
{
    return pointer_to_numpy(ptr, std::vector<int32_t>{nRows, nColumns});
}

/**
 Convert NumPy array to 1-dimensional memory block, i.e. contiguous array.

 @tparam T underlying scalar type.
 @param arr NumPy array.
 @return memory block object.
 */
template <typename T>
CMemBlock<T>
numpy_to_memblock(const py::array_t<T>& arr)
{
    std::string errmsg("numpy_to_memblock: Expecting a C-style contiguous numpy array");
    auto        c_style = py::detail::check_flags(arr.ptr(), py::array::c_style);
    errors::assertMsgCritical(c_style, errmsg);
    return CMemBlock<T>(arr.data(), arr.size(), numa::serial);
}

/**
 Convert Fortran-style NumPy array to 2-dimensional memory block, i.e. contiguous storage
 for 2-index quantity.

 @tparam T underlying scalar type.
 @param arr NumPy array.
 @return memory block object.
 */
template <typename T>
CMemBlock2D<T>
numpy_fstyle_to_memblock2d(const py::array_t<T, py::array::f_style>& arr)
{
    return CMemBlock2D<T>(arr.data(), arr.shape(0), arr.shape(1));
}

/** Wrapper for object constructors accepting an MPI communicator.
 *
 * @tparam T type of the object to wrap in a shared pointer.
 * @param py_comm Python object wrapping an MPI communicator.
 */
template <typename T, typename... Args>
inline std::shared_ptr<T>
create(py::object py_comm, Args&&... args)
{
    if (py_comm.is_none())
    {
        return std::make_shared<T>(MPI_COMM_WORLD, std::forward<Args>(args)...);
    }
    else
    {
        return std::make_shared<T>(*get_mpi_comm(py_comm), std::forward<Args>(args)...);
    }
}

/**
 Exports classes/functions in src/general to python.
 */
void export_general(py::module& m);

}  // namespace vlx_general

#endif /* ExportGeneral_hpp */
