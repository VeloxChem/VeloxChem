//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "ExportT4CIntegrals.hpp"

#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "ExportGeneral.hpp"
#include "FockDriver.hpp"
#include "FockGeomX000Driver.hpp"
#include "FockGeomXY00Driver.hpp"
#include "FockGeomX0Y0Driver.hpp"
#include "T4CScreener.hpp"
#include "T4CUtils.hpp"

namespace vlx_t4cintegrals {  // vlx_t4cintegrals namespace

// Exports classes/functions in src/t4c_* to python

void
export_t4cintegrals(py::module& m)
{
    // FockDriver class
    // Note: FockDriver is prefixed by an underscore and will be used in fockdriver.py

    PyClass<CFockDriver>(m, "_FockDriver")
        .def(py::init<>())
        .def("_set_block_size_factor", &CFockDriver::set_block_size_factor, "Sets block size factor.")
        .def(
            "compute_eri",
            [](const CFockDriver&  fock_drv,
               const CT4CScreener& screener,
               const int           nao,
               const int           ithreshold) -> py::array_t<double> {
                const auto eri_tensor = fock_drv.compute_eri(screener, nao, ithreshold);
                return vlx_general::pointer_to_numpy(eri_tensor.values(), {eri_tensor.getiIndex(), eri_tensor.getjIndex(), eri_tensor.getkIndex(), eri_tensor.getlIndex()});
            },
            "Computes single Fock matrix of requested type for two-electron integrals screener.")
        .def(
            "_compute_fock_omp",
            [](const CFockDriver&     fock_drv,
               const CMolecularBasis& basis,
               const CMolecule&       molecule,
               const CMatrix&         density,
               const std::string&     label,
               const double           exchange_factor,
               const double           omega) -> CMatrix {
                return fock_drv.compute(basis, molecule, density, label, exchange_factor, omega);
            },
            "Computes single Fock matrix of requested type for given molecule and basis.")
        .def(
            "_compute_fock_omp",
            [](const CFockDriver&  fock_drv,
               const CT4CScreener& screener,
               const CMatrix&      density,
               const std::string&  label,
               const double        exchange_factor,
               const double        omega,
               const int           ithreshold) -> CMatrix {
                return fock_drv.compute(screener, density, label, exchange_factor, omega, ithreshold);
            },
            "Computes single Fock matrix of requested type for two-electron integrals screener.")
        .def(
            "_compute_fock_omp",
            [](const CFockDriver&              fock_drv,
               const CT4CScreener&             screener,
               const CMatrices&                densities,
               const std::vector<std::string>& labels,
               const double                    exchange_factor,
               const double                    omega,
               const int                       ithreshold) -> CMatrices {
                return fock_drv.compute(screener, densities, labels, exchange_factor, omega, ithreshold);
            },
            "Computes Fock matrices of requested type for two-electron integrals screener.")
        .def(
            "_compute_local_fock",
            [](const CFockDriver&  fock_drv,
               const CT4CScreener& screener,
               const int           rank,
               const int           nodes,
               const CMatrix&      density,
               const std::string&  label,
               const double        exchange_factor,
               const double        omega,
               const int           ithreshold) -> CMatrix {
                return fock_drv.compute(screener, rank, nodes, density, label, exchange_factor, omega, ithreshold);
            },
            "Computes single Fock matrix of requested type for two-electron integrals screener.")
        .def(
            "_compute_local_fock",
            [](const CFockDriver&              fock_drv,
               const CT4CScreener&             screener,
               const int                       rank,
               const int                       nodes,
               const CMatrices&                densities,
               const std::vector<std::string>& labels,
               const double                    exchange_factor,
               const double                    omega,
               const int                       ithreshold) -> CMatrices {
                return fock_drv.compute(screener, rank, nodes, densities, labels, exchange_factor, omega, ithreshold);
            },
            "Computes Fock matrices of requested type for two-electron integrals screener.");

    // CT4CScreener class
    PyClass<CT4CScreener>(m, "T4CScreener")
        .def(py::init<>())
        .def(py::init<const CT4CScreener&>())
        .def(py::init<const std::vector<CBlockedGtoPairBlock>&>())
        .def(py::pickle([](const CT4CScreener& screener) { return py::make_tuple(screener.gto_pair_blocks()); },
                        [](py::tuple t) { return CT4CScreener(t[0].cast<std::vector<CBlockedGtoPairBlock>>()); }))
        .def("partition", &CT4CScreener::partition, "Partition basis funtion pairs blocks for given molecule and basis.")
        .def("partition_atom", &CT4CScreener::partition_atom, "Partition basis funtion pairs blocks for given molecule and basis.")
        .def("partition_atom_pair", &CT4CScreener::partition_atom_pair, "Partition basis funtion pairs blocks for given molecule and basis.")
        .def("gto_pair_blocks", &CT4CScreener::gto_pair_blocks, "Gets vector of blocked basis function pairs blocks.")
        .def("__eq__", [](const CT4CScreener& self, const CT4CScreener& other) { return self == other; })
        .def("__ne__", [](const CT4CScreener& self, const CT4CScreener& other) { return self != other; })
        .def("__copy__", [](const CT4CScreener& self) { return CT4CScreener(self); })
        .def("__deepcopy__", [](const CT4CScreener& self, py::dict) { return CT4CScreener(self); });

    // CFockGeom1000Driver class
    PyClass<CFockGeomX000Driver<1>>(m, "FockGeom1000Driver")
        .def(py::init<>())
        .def("_set_block_size_factor", &CFockGeomX000Driver<1>::set_block_size_factor, "Sets block size factor.")
        .def(
            "compute",
            [](const CFockGeomX000Driver<1>& self,
               const CMolecularBasis&        basis,
               const CMolecule&              molecule,
               const CMatrix&                density,
               const int                     iatom,
               const std::string&            label,
               const double                  exchange_factor,
               const double                  omega) -> CMatrices {
                return self.compute(basis, molecule, density, iatom, label, exchange_factor, omega);
            },
            "Computes gradient of Fock matrix of requested type for given molecule and basis.")
        .def(
            "compute",
            [](const CFockGeomX000Driver<1>& self,
               const CMolecularBasis&        basis,
               const CT4CScreener&           screener_atom,
               const CT4CScreener&           screener,
               const CMatrix&                density,
               const int                     iatom,
               const std::string&            label,
               const double                  exchange_factor,
               const double                  omega,
               const int                     ithreshold) -> CMatrices {
                return self.compute(basis, screener_atom, screener, density, iatom, label, exchange_factor, omega, ithreshold);
            },
            "Computes gradient of Fock matrix of requested type for given molecule and basis.")
        .def(
            "compute",
            [](const CFockGeomX000Driver<1>& self,
               const CMolecularBasis&        basis,
               const CT4CScreener&           screener_atom,
               const CT4CScreener&           screener,
               const CMatrix&                density,
               const CMatrix&                density2,
               const int                     iatom,
               const std::string&            label,
               const double                  exchange_factor,
               const double                  omega,
               const int                     ithreshold) -> std::vector<double> {
                return self.compute(basis, screener_atom, screener, density, density2, iatom, label, exchange_factor, omega, ithreshold);
            },
            "Computes gradient of Fock matrix of requested type for given molecule and basis.");
    
    // CFockGeom2000Driver class
    PyClass<CFockGeomX000Driver<2>>(m, "FockGeom2000Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CFockGeomX000Driver<2>& fock_drv,
               const CMolecularBasis&        basis,
               const CMolecule&              molecule,
               const CMatrix&                density,
               const int                     iatom,
               const std::string&            label,
               const double                  exchange_factor,
               const double                  omega) -> CMatrices {
                return fock_drv.compute(basis, molecule, density, iatom, label, exchange_factor, omega);
            },
            "Computes Hessian of Fock matrix of requested type for given molecule and basis.")
        .def(
            "compute",
            [](const CFockGeomX000Driver<2>& self,
               const CMolecularBasis&        basis,
               const CT4CScreener&           screener_atom,
               const CT4CScreener&           screener,
               const CMatrix&                density,
               const CMatrix&                density2,
               const int                     iatom,
               const std::string&            label,
               const double                  exchange_factor,
               const double                  omega,
               const int                     ithreshold) -> std::vector<double> {
                return self.compute(basis, screener_atom, screener, density, density2, iatom, label, exchange_factor, omega, ithreshold);
            },
            "Computes Hessian of Fock matrix of requested type for given molecule and basis.");
    
    // CFockGeom1100Driver class
    PyClass<CFockGeomXY00Driver<1,1>>(m, "FockGeom1100Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CFockGeomXY00Driver<1,1>& fock_drv,
               const CMolecularBasis&        basis,
               const CMolecule&              molecule,
               const CMatrix&                density,
               const int                     iatom,
               const int                     jatom,
               const std::string&            label,
               const double                  exchange_factor,
               const double                  omega) -> CMatrices {
                   return fock_drv.compute(basis, molecule, density, iatom, jatom, label, exchange_factor, omega);
            },
            "Computes Hessian of Fock matrix of requested type for given molecule and basis.")
        .def(
            "compute",
            [](const CFockGeomXY00Driver<1,1>& self,
               const CMolecularBasis&          basis,
               const CT4CScreener&             screener_atom_pair,
               const CT4CScreener&             screener,
               const CMatrix&                  density,
               const CMatrix&                  density2,
               const int                       iatom,
               const int                       jatom,
               const std::string&              label,
               const double                    exchange_factor,
               const double                    omega,
               const int                       ithreshold) -> std::vector<double> {
                return self.compute(basis, screener_atom_pair, screener, density, density2, iatom, jatom, label, exchange_factor, omega, ithreshold);
            },
            "Computes Hessian of Fock matrix of requested type for given molecule and basis.");
    
    // CFockGeom1010Driver class
    PyClass<CFockGeomX0Y0Driver<1,1>>(m, "FockGeom1010Driver")
        .def(py::init<>())
        .def(
            "compute",
            [](const CFockGeomX0Y0Driver<1,1>& fock_drv,
               const CMolecularBasis&        basis,
               const CMolecule&              molecule,
               const CMatrix&                density,
               const int                     iatom,
               const int                     jatom,
               const std::string&            label,
               const double                  exchange_factor,
               const double                  omega) -> CMatrices {
                   return fock_drv.compute(basis, molecule, density, iatom, jatom, label, exchange_factor, omega);
            },
            "Computes Hessian of Fock matrix of requested type for given molecule and basis.")
        .def(
            "compute",
            [](const CFockGeomX0Y0Driver<1,1>& self,
               const CMolecularBasis&          basis,
               const CT4CScreener&             screener_atom_i,
               const CT4CScreener&             screener_atom_j,
               const CMatrix&                  density,
               const CMatrix&                  density2,
               const int                       iatom,
               const int                       jatom,
               const std::string&              label,
               const double                    exchange_factor,
               const double                    omega,
               const int                       ithreshold) -> std::vector<double> {
                return self.compute(basis, screener_atom_i, screener_atom_j, density, density2, iatom, jatom, label, exchange_factor, omega, ithreshold);
            },
            "Computes Hessian of Fock matrix of requested type for given molecule and basis.");
}

}  // namespace vlx_t4cintegrals
