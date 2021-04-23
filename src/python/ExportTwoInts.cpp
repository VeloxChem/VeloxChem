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

#include "ExportTwoInts.hpp"

#include <mpi.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "AOFockMatrix.hpp"
#include "DenseMatrix.hpp"
#include "ElectronRepulsionIntegralsDriver.hpp"
#include "EriScreenerType.hpp"
#include "ErrorHandler.hpp"
#include "ExportGeneral.hpp"
#include "ExportMath.hpp"
#include "FockMatrixType.hpp"
#include "MOIntsBatch.hpp"
#include "MOIntsType.hpp"
#include "ScreeningContainer.hpp"
#include "ExportOrbData.hpp"

namespace py = pybind11;
using namespace py::literals;

namespace vlx_twoints {  // vlx_twoints namespace

// Helper function for converting CAOFockMatrix to numpy array

static py::array_t<double>
CAOFockMatrix_to_numpy(const CAOFockMatrix& self, const int32_t iFockMatrix)
{
    return vlx_general::pointer_to_numpy(self.getFock(iFockMatrix), self.getNumberOfRows(iFockMatrix), self.getNumberOfColumns(iFockMatrix));
}

static py::array_t<double>
CAOFockMatrix_alpha_to_numpy(const CAOFockMatrix& self, const int32_t iFockMatrix)
{
    return vlx_general::pointer_to_numpy(self.getFock(iFockMatrix, "alpha"), self.getNumberOfRows(iFockMatrix), self.getNumberOfColumns(iFockMatrix));
}

static py::array_t<double>
CAOFockMatrix_beta_to_numpy(const CAOFockMatrix& self, const int32_t iFockMatrix)
{
    return vlx_general::pointer_to_numpy(self.getFock(iFockMatrix, "beta"), self.getNumberOfRows(iFockMatrix), self.getNumberOfColumns(iFockMatrix));
}

// Helper function for CAOFockMatrix constructor

static std::shared_ptr<CAOFockMatrix>
CAOFockMatrix_from_numpy_list(const std::vector<py::array_t<double>>& arrays,
                              const std::vector<fockmat>&             types,
                              const std::vector<double>&              factors,
                              const std::vector<int32_t>&             ids)
{
    std::vector<CDenseMatrix> fmat;

    for (const auto& x : arrays)
    {
        auto mp = vlx_math::CDenseMatrix_from_numpy(x);
        fmat.push_back(*mp);
    }

    return std::make_shared<CAOFockMatrix>(fmat, types, factors, ids);
}

// Helper function for exporting CElectronRepulsionIntegralsDriver.computeInMemory

static void
CElectronRepulsionIntegralsDriver_compute_in_mem(const CElectronRepulsionIntegralsDriver& self,
                                                 const CMolecule&                         molecule,
                                                 const CMolecularBasis&                   basis,
                                                 py::array_t<double>&                     eri)
{
    std::string errsrc("ElectronRepulsionIntegralsDriver.compute_in_mem: Expecting a C-style contiguous numpy array");

    auto c_style = py::detail::check_flags(eri.ptr(), py::array::c_style);

    errors::assertMsgCritical(c_style, errsrc);

    std::string errshape("ElectronRepulsionIntegralsDriver.compute_in_mem: Invalid shape");

    errors::assertMsgCritical(eri.ndim() == 4, errshape);

    auto nao = vlx_orbdata::get_number_of_atomic_orbitals(molecule, basis);

    std::string errsize("ElectronRepulsionIntegralsDriver.compute_in_mem: Inconsistent size");

    errors::assertMsgCritical(eri.shape(0) == nao, errsize);

    errors::assertMsgCritical(eri.shape(1) == nao, errsize);

    errors::assertMsgCritical(eri.shape(2) == nao, errsize);

    errors::assertMsgCritical(eri.shape(3) == nao, errsize);

    self.computeInMemory(molecule, basis, eri.mutable_data());
}

// Helper function for reduce_sum CAOFockMatrix object

static void
CAOFockMatrix_reduce_sum(CAOFockMatrix& self, int32_t rank, int32_t nodes, py::object py_comm)
{
    auto comm = vlx_general::get_mpi_comm(py_comm);

    self.reduce_sum(rank, nodes, *comm);
}

// Helper function for converting CMOIntsBatch to numpy array

static py::array_t<double>
CMOIntsBatch_to_numpy(const CMOIntsBatch& self, const int32_t iBatch)
{
    return vlx_general::pointer_to_numpy(self.getBatch(iBatch), self.getNumberOfRows(), self.getNumberOfColumns());
}

static py::array_t<double>
CMOIntsBatch_to_numpy_2(const CMOIntsBatch& self, const CTwoIndexes& iGeneratorPair)
{
    return vlx_general::pointer_to_numpy(self.getBatch(iGeneratorPair), self.getNumberOfRows(), self.getNumberOfColumns());
}

static py::array_t<double>
CMOIntsBatchXY_to_numpy(const CMOIntsBatch& self, const int32_t iBatch)
{
    return vlx_general::pointer_to_numpy(self.getBatchXY(iBatch), self.getNumberOfRows(), self.getNumberOfColumns());
}

static py::array_t<double>
CMOIntsBatchYX_to_numpy(const CMOIntsBatch& self, const int32_t iBatch)
{
    return vlx_general::pointer_to_numpy(self.getBatchYX(iBatch), self.getNumberOfColumns(), self.getNumberOfRows());
}

// Helper function for collecting CMOIntsBatch on global master node

static void
CMOIntsBatch_collectBatches(CMOIntsBatch& self, int32_t cross_rank, int32_t cross_nodes, py::object py_cross_comm)
{
    auto cross_comm = vlx_general::get_mpi_comm(py_cross_comm);

    int32_t nrows = self.getNumberOfRows();

    int32_t ncols = self.getNumberOfColumns();

    // master: receive data

    if (cross_rank == mpi::master())
    {
        for (int32_t cross_id = 1; cross_id < cross_nodes; ++cross_id)
        {
            int32_t tag_id = cross_id;

            MPI_Status mstat;

            int32_t numbatches = 0;

            auto merror = MPI_Recv(&numbatches, 1, MPI_INT, cross_id, tag_id++, *cross_comm, &mstat);

            if (merror != MPI_SUCCESS) mpi::abort(merror, "collectBatches");

            std::vector<double> data(nrows * ncols);

            for (int32_t ibatch = 0; ibatch < numbatches; ibatch++)
            {
                int32_t first = -1, second = -1;

                merror = MPI_Recv(&first, 1, MPI_INT, cross_id, tag_id++, *cross_comm, &mstat);

                if (merror != MPI_SUCCESS) mpi::abort(merror, "collectBatches");

                merror = MPI_Recv(&second, 1, MPI_INT, cross_id, tag_id++, *cross_comm, &mstat);

                if (merror != MPI_SUCCESS) mpi::abort(merror, "collectBatches");

                merror = MPI_Recv(data.data(), nrows * ncols, MPI_DOUBLE, cross_id, tag_id++, *cross_comm, &mstat);

                if (merror != MPI_SUCCESS) mpi::abort(merror, "collectBatches");

                self.appendMOInts(CDenseMatrix(data, nrows, ncols), CTwoIndexes(first, second));
            }
        }
    }

    // worker process: send data

    else
    {
        int32_t tag_id = cross_rank;

        auto numbatches = self.getNumberOfBatches();

        auto merror = MPI_Send(&numbatches, 1, MPI_INT, mpi::master(), tag_id++, *cross_comm);

        if (merror != MPI_SUCCESS) mpi::abort(merror, "collectBatches");

        auto genpairs = self.getGeneratorPairs();

        for (int32_t ibatch = 0; ibatch < numbatches; ibatch++)
        {
            auto data = self.getBatch(ibatch);

            auto first = genpairs[ibatch].first();

            auto second = genpairs[ibatch].second();

            merror = MPI_Send(&first, 1, MPI_INT, mpi::master(), tag_id++, *cross_comm);

            if (merror != MPI_SUCCESS) mpi::abort(merror, "collectBatches");

            merror = MPI_Send(&second, 1, MPI_INT, mpi::master(), tag_id++, *cross_comm);

            if (merror != MPI_SUCCESS) mpi::abort(merror, "collectBatches");

            merror = MPI_Send(data, nrows * ncols, MPI_DOUBLE, mpi::master(), tag_id++, *cross_comm);

            if (merror != MPI_SUCCESS) mpi::abort(merror, "collectBatches");
        }
    }
}

// Exports classes/functions in src/twoints to python

void
export_twoints(py::module& m)
{
    // fockmat enum class

    py::enum_<fockmat>(m, "fockmat")
        .value("restjk", fockmat::restjk)
        .value("restjkx", fockmat::restjkx)
        .value("restj", fockmat::restj)
        .value("restk", fockmat::restk)
        .value("restkx", fockmat::restkx)
        .value("rgenjk", fockmat::rgenjk)
        .value("rgenjkx", fockmat::rgenjkx)
        .value("rgenj", fockmat::rgenj)
        .value("rgenk", fockmat::rgenk)
        .value("rgenkx", fockmat::rgenkx)
        .value("unrestjk", fockmat::unrestjk)
        .value("unrestj", fockmat::unrestj)
        .value("unrestjkx", fockmat::unrestjkx);

    // ericut enum class
    // clang-format off
    py::enum_<ericut>(m, "ericut")
        .value("qq", ericut::qq)
        .value("qqr", ericut::qqr)
        .value("qqden", ericut::qqden)
        .value("qqrden", ericut::qqrden);
    // clang-format on

    // moints enum class

    py::enum_<moints>(m, "moints")
        .value("oooo", moints::oooo)
        .value("ooov", moints::ooov)
        .value("oovv", moints::oovv)
        .value("ovov", moints::ovov)
        .value("ovvv", moints::ovvv)
        .value("vvvv", moints::vvvv)
        .value("asym_oooo", moints::asym_oooo)
        .value("asym_ooov", moints::asym_ooov)
        .value("asym_oovv", moints::asym_oovv)
        .value("asym_ovov", moints::asym_ovov)
        .value("asym_ovvv", moints::asym_ovvv)
        .value("asym_vvvv", moints::asym_vvvv);

    // CAOFockMatrix class

    PyClass<CAOFockMatrix>(m, "AOFockMatrix")
        .def(py::init<>())
        .def(py::init<const CAODensityMatrix&>())
        .def(py::init<const CAOFockMatrix&>())
        .def(py::init(&CAOFockMatrix_from_numpy_list))
        .def("__str__", &CAOFockMatrix::getString)
        .def("to_numpy", &CAOFockMatrix_to_numpy)
        .def("alpha_to_numpy", &CAOFockMatrix_alpha_to_numpy)
        .def("beta_to_numpy", &CAOFockMatrix_beta_to_numpy)
        .def("is_closed_shell", &CAOFockMatrix::isClosedShell)
        .def("number_of_fock_matrices", &CAOFockMatrix::getNumberOfFockMatrices)
        .def("set_fock_type",
             &CAOFockMatrix::setFockType,
             "Set type of Fock matrix at fock_idx with given spin",
             "fock_type"_a,
             "fock_idx"_a,
             "spin"_a = std::string("ALPHA"))
        .def("set_scale_factor",
             &CAOFockMatrix::setFockScaleFactor,
             "Set scaling factor of exact exchange for Fock matrix at fock_idx with given spin",
             "factor"_a,
             "fock_idx"_a,
             "spin"_a = std::string("ALPHA"))
        .def("get_fock_type",
             &CAOFockMatrix::getFockType,
             "Get type of Fock matrix at fock_idx with given spin",
             "fock_idx"_a,
             "spin"_a = std::string("ALPHA"))
        .def("get_scale_factor",
             &CAOFockMatrix::getScaleFactor,
             "Get scaling factor of exchange contribution for Fock matrix at fock_idx with given spin",
             "fock_idx"_a,
             "spin"_a = std::string("ALPHA"))
        .def("get_density_identifier",
             &CAOFockMatrix::getDensityIdentifier,
             "Get identifier of AO density matrix used to construct Fock matrix at fock_idx with given spin",
             "fock_idx"_a,
             "spin"_a = std::string("ALPHA"))
        .def("add_hcore", &CAOFockMatrix::addCoreHamiltonian)
        .def("add", &CAOFockMatrix::add)
        .def("add_matrix",
             &CAOFockMatrix::addOneElectronMatrix,
             "Add one-electron operator matrix to Fock matrix at fock_idx with given spin",
             "one_el"_a,
             "fock_idx"_a,
             "spin"_a = std::string("ALPHA"))
        .def("scale",
             &CAOFockMatrix::scale,
             "Scale Fock matrix at fock_idx with given spin by a factor",
             "factor"_a,
             "fock_idx"_a,
             "spin"_a = std::string("ALPHA"))
        .def("reduce_sum", &CAOFockMatrix_reduce_sum)
        .def("get_energy", &CAOFockMatrix::getElectronicEnergy)
        .def(py::self == py::self);

    // CCauchySchwarzScreener class

    PyClass<CCauchySchwarzScreener>(m, "CauchySchwarzScreener")
        .def(py::init<>())
        .def(py::init<const CCauchySchwarzScreener&>())
        .def("get_threshold", &CCauchySchwarzScreener::getThreshold)
        .def("get_screening_scheme", &CCauchySchwarzScreener::getScreeningScheme)
        .def(py::self == py::self);

    // CScreeningContainer class

    PyClass<CScreeningContainer>(m, "ScreeningContainer")
        .def(py::init<>())
        .def(py::init<const CScreeningContainer&>())
        .def("is_empty", &CScreeningContainer::isEmpty)
        .def("number_of_screeners", &CScreeningContainer::getNumberOfScreeners)
        .def("get_screener", &CScreeningContainer::getScreener)
        .def("set_threshold", &CScreeningContainer::setThreshold)
        .def(py::self == py::self);

    // CElectronRepulsionIntegralsDriver class

    PyClass<CElectronRepulsionIntegralsDriver>(m, "ElectronRepulsionIntegralsDriver")
        .def(py::init(&vlx_general::create<CElectronRepulsionIntegralsDriver>), "comm"_a = py::none())
        .def("compute",
             vlx_general::overload_cast_<const ericut, const double, const CMolecule&, const CMolecularBasis&>()(
                 &CElectronRepulsionIntegralsDriver::compute, py::const_),
             // the wonky format of the raw string literal is to get help(...) in Python to look nice
             R"pbdoc(
Computes Q values for electron repulsion integrals for molecule with specific AO
basis set and stores results in screening container object.

:param screening:
    The screening scheme for screening container object.
:param threshold:
    The screening threshold for screening container object.
:param molecule:
    The molecule.
:param ao_basis:
    The molecular AO basis.
:return: the screening container with Q values.
)pbdoc",
             "screening"_a,
             "threshold"_a,
             "molecule"_a,
             "ao_basis"_a)
        .def("compute",
             vlx_general::
                 overload_cast_<CAOFockMatrix&, const CAODensityMatrix&, const CMolecule&, const CMolecularBasis&, const CScreeningContainer&>()(
                     &CElectronRepulsionIntegralsDriver::compute, py::const_),
             // the wonky format of the raw string literal is to get help(...) in Python to look nice
             R"pbdoc(
Computes electron repulsion integrals and stores them in AO Fock matrix for
molecule with specific AO basis set. Performs screening according to
``screening`` parameters.

:param ao_fock:
    The AO Fock matrix.
:param ao_density:
    The AO density matrix.
:param molecule:
    The molecule.
:param ao_basis:
    The molecular AO basis.
:param screening:
    The screening container object.
)pbdoc",
             "ao_fock"_a,
             "ao_density"_a,
             "molecule"_a,
             "ao_basis"_a,
             "screening"_a)
        .def("compute_in_mem",
             &CElectronRepulsionIntegralsDriver_compute_in_mem,
             "Computes electron repulsion integrals as a full 4D array and stores them in memory.",
             "molecule"_a,
             "basis"_a,
             "eri_tensor"_a)
        .def(
            "compute_in_mem",
            [](const CElectronRepulsionIntegralsDriver& eridrv, const CMolecule& molecule, const CMolecularBasis& basis) {
                auto nao = vlx_orbdata::get_number_of_atomic_orbitals(molecule, basis);
                auto eri = py::array_t<double, py::array::c_style>({nao, nao, nao, nao});
                eridrv.computeInMemory(molecule, basis, eri.mutable_data());
                return eri;
            },
            "Computes electron repulsion integrals as a full 4D array and stores them in memory.",
            "molecule"_a,
            "basis"_a);

    // CMOIntsBatch class

    PyClass<CMOIntsBatch>(m, "MOIntsBatch")
        .def(py::init<>())
        .def("to_numpy", &CMOIntsBatch_to_numpy)
        .def("to_numpy", &CMOIntsBatch_to_numpy_2)
        .def("xy_to_numpy", &CMOIntsBatchXY_to_numpy)
        .def("yx_to_numpy", &CMOIntsBatchYX_to_numpy)
        .def("number_of_batches", &CMOIntsBatch::getNumberOfBatches)
        .def("number_of_rows", &CMOIntsBatch::getNumberOfRows)
        .def("number_of_columns", &CMOIntsBatch::getNumberOfColumns)
        .def("append", &CMOIntsBatch::append)
        .def("set_batch_type", &CMOIntsBatch::setBatchType)
        .def("set_ext_indexes", &CMOIntsBatch::setExternalIndexes)
        .def("get_batch_type", &CMOIntsBatch::getBatchType)
        .def("get_ext_indexes", &CMOIntsBatch::getExternalIndexes)
        .def("get_gen_pairs", &CMOIntsBatch::getGeneratorPairs)
        .def("collect_batches", &CMOIntsBatch_collectBatches);
}

}  // namespace vlx_twoints
