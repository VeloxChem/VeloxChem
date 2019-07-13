//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2019 by VeloxChem developers. All rights reserved.

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <mpi.h>

#include "AOFockMatrix.hpp"
#include "DenseMatrix.hpp"
#include "ElectronRepulsionIntegralsDriver.hpp"
#include "EriScreenerType.hpp"
#include "ExportGeneral.hpp"
#include "ExportMath.hpp"
#include "ExportTwoInts.hpp"
#include "FockMatrixType.hpp"
#include "MOIntsBatch.hpp"
#include "MOIntsType.hpp"
#include "ScreeningContainer.hpp"

namespace py = pybind11;

namespace vlx_twoints {  // vlx_twoints namespace

// Helper function for printing CAOFockMatrix

static std::string
CAOFockMatrix_str(const CAOFockMatrix& self)
{
    return self.getString();
}

// Helper function for converting CAOFockMatrix to numpy array

static py::array_t<double>
CAOFockMatrix_to_numpy(const CAOFockMatrix& self, const int32_t iFockMatrix)
{
    return vlx_general::pointer_to_numpy(
        self.getFock(iFockMatrix), self.getNumberOfRows(iFockMatrix), self.getNumberOfColumns(iFockMatrix));
}

static py::array_t<double>
CAOFockMatrix_alpha_to_numpy(const CAOFockMatrix& self, const int32_t iFockMatrix)
{
    return vlx_general::pointer_to_numpy(
        self.getFock(iFockMatrix, false), self.getNumberOfRows(iFockMatrix), self.getNumberOfColumns(iFockMatrix));
}

static py::array_t<double>
CAOFockMatrix_beta_to_numpy(const CAOFockMatrix& self, const int32_t iFockMatrix)
{
    return vlx_general::pointer_to_numpy(
        self.getFock(iFockMatrix, true), self.getNumberOfRows(iFockMatrix), self.getNumberOfColumns(iFockMatrix));
}

// Helper function for CAOFockMatrix constructor

static std::shared_ptr<CAOFockMatrix>
CAOFockMatrix_from_numpy_list(const std::vector<py::array_t<double>>& arrays,
                              const std::vector<fockmat>&             types,
                              const std::vector<double>&              factors,
                              const std::vector<int32_t>&             ids)
{
    std::vector<CDenseMatrix> fmat;

    for (size_t i = 0; i < arrays.size(); i++)
    {
        auto mp = vlx_math::CDenseMatrix_from_numpy(arrays[i]);

        fmat.push_back(*mp);
    }

    return std::shared_ptr<CAOFockMatrix>(new CAOFockMatrix(fmat, types, factors, ids));
}

// Helper function for CElectronRepulsionIntegralsDriver constructor

static std::shared_ptr<CElectronRepulsionIntegralsDriver>
CElectronRepulsionIntegralsDriver_create(py::object py_comm)
{
    MPI_Comm* comm_ptr = vlx_general::get_mpi_comm(py_comm);

    return std::shared_ptr<CElectronRepulsionIntegralsDriver>(new CElectronRepulsionIntegralsDriver(*comm_ptr));
}

// Helper function for reduce_sum CAOFockMatrix object

static void
CAOFockMatrix_reduce_sum(CAOFockMatrix& self, int32_t rank, int32_t nodes, py::object py_comm)
{
    MPI_Comm* comm_ptr = vlx_general::get_mpi_comm(py_comm);

    self.reduce_sum(rank, nodes, *comm_ptr);
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
    return vlx_general::pointer_to_numpy(
        self.getBatch(iGeneratorPair), self.getNumberOfRows(), self.getNumberOfColumns());
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
    MPI_Comm* cross_comm_ptr = vlx_general::get_mpi_comm(py_cross_comm);

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

            auto merror = MPI_Recv(&numbatches, 1, MPI_INT, cross_id, tag_id++, *cross_comm_ptr, &mstat);

            if (merror != MPI_SUCCESS) mpi::abort(merror, "collectBatches");

            std::vector<double> data(nrows * ncols);

            for (int32_t ibatch = 0; ibatch < numbatches; ibatch++)
            {
                int32_t first = -1, second = -1;

                merror = MPI_Recv(&first, 1, MPI_INT, cross_id, tag_id++, *cross_comm_ptr, &mstat);

                if (merror != MPI_SUCCESS) mpi::abort(merror, "collectBatches");

                merror = MPI_Recv(&second, 1, MPI_INT, cross_id, tag_id++, *cross_comm_ptr, &mstat);

                if (merror != MPI_SUCCESS) mpi::abort(merror, "collectBatches");

                merror = MPI_Recv(data.data(), nrows * ncols, MPI_DOUBLE, cross_id, tag_id++, *cross_comm_ptr, &mstat);

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

        auto merror = MPI_Send(&numbatches, 1, MPI_INT, mpi::master(), tag_id++, *cross_comm_ptr);

        if (merror != MPI_SUCCESS) mpi::abort(merror, "collectBatches");

        auto genpairs = self.getGeneratorPairs();

        for (int32_t ibatch = 0; ibatch < numbatches; ibatch++)
        {
            auto data = self.getBatch(ibatch);

            auto first = genpairs[ibatch].first();

            auto second = genpairs[ibatch].second();

            merror = MPI_Send(&first, 1, MPI_INT, mpi::master(), tag_id++, *cross_comm_ptr);

            if (merror != MPI_SUCCESS) mpi::abort(merror, "collectBatches");

            merror = MPI_Send(&second, 1, MPI_INT, mpi::master(), tag_id++, *cross_comm_ptr);

            if (merror != MPI_SUCCESS) mpi::abort(merror, "collectBatches");

            merror = MPI_Send(data, nrows * ncols, MPI_DOUBLE, mpi::master(), tag_id++, *cross_comm_ptr);

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
        .value("unrestjk", fockmat::unrestjk);

    // ericut enum class

    py::enum_<ericut>(m, "ericut")
        .value("qq", ericut::qq)
        .value("qqr", ericut::qqr)
        .value("qqden", ericut::qqden)
        .value("qqrden", ericut::qqrden);

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

    py::class_<CAOFockMatrix, std::shared_ptr<CAOFockMatrix>>(m, "AOFockMatrix")
        .def(py::init<>())
        .def(py::init<const CAODensityMatrix&>())
        .def(py::init<const CAOFockMatrix&>())
        .def(py::init(&CAOFockMatrix_from_numpy_list))
        .def("__str__", &CAOFockMatrix_str)
        .def("to_numpy", &CAOFockMatrix_to_numpy)
        .def("alpha_to_numpy", &CAOFockMatrix_alpha_to_numpy)
        .def("beta_to_numpy", &CAOFockMatrix_beta_to_numpy)
        .def("number_of_fock_matrices", &CAOFockMatrix::getNumberOfFockMatrices)
        .def("set_fock_type", &CAOFockMatrix::setFockType, py::arg(), py::arg(), py::arg("beta") = false)
        .def("set_scale_factor", &CAOFockMatrix::setFockScaleFactor, py::arg(), py::arg(), py::arg("beta") = false)
        .def("get_fock_type", &CAOFockMatrix::getFockType, py::arg(), py::arg("beta") = false)
        .def("get_scale_factor", &CAOFockMatrix::getScaleFactor, py::arg(), py::arg("beta") = false)
        .def("get_density_identifier", &CAOFockMatrix::getDensityIdentifier, py::arg(), py::arg("beta") = false)
        .def("add_hcore", &CAOFockMatrix::addCoreHamiltonian)
        .def("add", &CAOFockMatrix::add)
        .def("add_matrix", &CAOFockMatrix::addOneElectronMatrix)
        .def("scale", &CAOFockMatrix::scale)
        .def("reduce_sum", &CAOFockMatrix_reduce_sum)
        .def("get_energy", &CAOFockMatrix::getElectronicEnergy)
        .def(py::self == py::self);

    // CCauchySchwarzScreener class

    py::class_<CCauchySchwarzScreener, std::shared_ptr<CCauchySchwarzScreener>>(m, "CauchySchwarzScreener")
        .def(py::init<>())
        .def(py::init<const CCauchySchwarzScreener&>())
        .def("get_threshold", &CCauchySchwarzScreener::getThreshold)
        .def("get_screening_scheme", &CCauchySchwarzScreener::getScreeningScheme)
        .def(py::self == py::self);

    // CScreeningContainer class

    py::class_<CScreeningContainer, std::shared_ptr<CScreeningContainer>>(m, "ScreeningContainer")
        .def(py::init<>())
        .def(py::init<const CScreeningContainer&>())
        .def("is_empty", &CScreeningContainer::isEmpty)
        .def("number_of_screeners", &CScreeningContainer::getNumberOfScreeners)
        .def("get_screener", &CScreeningContainer::getScreener)
        .def("set_threshold", &CScreeningContainer::setThreshold)
        .def(py::self == py::self);

    // CElectronRepulsionIntegralsDriver class

    py::class_<CElectronRepulsionIntegralsDriver, std::shared_ptr<CElectronRepulsionIntegralsDriver>>(
        m, "ElectronRepulsionIntegralsDriver")
        .def(py::init(&CElectronRepulsionIntegralsDriver_create))
        .def("compute",
             (void (CElectronRepulsionIntegralsDriver::*)(CAOFockMatrix&,
                                                          const CAODensityMatrix&,
                                                          const CMolecule&,
                                                          const CMolecularBasis&,
                                                          const CScreeningContainer&) const) &
                 CElectronRepulsionIntegralsDriver::compute)
        .def("compute",
             (CScreeningContainer(CElectronRepulsionIntegralsDriver::*)(
                 const ericut, const double, const CMolecule&, const CMolecularBasis&) const) &
                 CElectronRepulsionIntegralsDriver::compute);

    // CMOIntsBatch class

    py::class_<CMOIntsBatch, std::shared_ptr<CMOIntsBatch>>(m, "MOIntsBatch")
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
