//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <mpi.h>

#include "DenseMatrix.hpp"
#include "AOFockMatrix.hpp"
#include "FockMatrixType.hpp"
#include "EriScreenerType.hpp"
#include "ScreeningContainer.hpp"
#include "ElectronRepulsionIntegralsDriver.hpp"
#include "ExportMath.hpp"
#include "ExportGeneral.hpp"
#include "ExportTwoInts.hpp"

namespace bp = boost::python;

namespace np = boost::python::numpy;

namespace bp_twoints { // bp_twoints namespace

// Helper function for printing CAOFockMatrix

static std::string
CAOFockMatrix_str (const CAOFockMatrix& self)
{
    return self.getString();
}

// Helper function for converting CAOFockMatrix to numpy array

static np::ndarray
CAOFockMatrix_to_numpy(const CAOFockMatrix& self,
                       const int32_t iFockMatrix)
{
    return bp_general::pointer_to_numpy(self.getFock(iFockMatrix),
                                        self.getNumberOfRows(iFockMatrix),
                                        self.getNumberOfColumns(iFockMatrix));
}

// Helper function for CAOFockMatrix constructor

static std::shared_ptr<CAOFockMatrix>
CAOFockMatrix_from_numpy_list(const bp::list& arr_list,
                              const bp::list& fock_type_list,
                              const bp::list& scale_fac_list,
                              const bp::list& id_dmat_list)
{
    std::vector<CDenseMatrix> fmat;
    std::vector<fockmat> types;
    std::vector<double> factors;
    std::vector<int32_t> ids;

    for (int i = 0; i < bp::len(arr_list); i++)
    {
        np::ndarray arr    = np::array(arr_list[i]);
        fockmat     type   = bp::extract<fockmat>(fock_type_list[i]);
        double      factor = bp::extract<double>(scale_fac_list[i]);
        int         id     = bp::extract<int>(id_dmat_list[i]);

        std::shared_ptr<CDenseMatrix> mp = bp_math::CDenseMatrix_from_numpy(arr);

        fmat.push_back(*mp);
        types.push_back(type);
        factors.push_back(factor);
        ids.push_back(id);
    }

    return std::shared_ptr<CAOFockMatrix>(
            new CAOFockMatrix(fmat, types, factors, ids)
            );
}

// Helper function for CElectronRepulsionIntegralsDriver constructor

static std::shared_ptr<CElectronRepulsionIntegralsDriver>
CElectronRepulsionIntegralsDriver_create(int32_t    globRank,
                                         int32_t    globNodes,
                                         bp::object py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    return std::shared_ptr<CElectronRepulsionIntegralsDriver>(
        new CElectronRepulsionIntegralsDriver(globRank, globNodes, *comm_ptr)
        );
}

// Helper functions for overloading CElectronRepulsionIntegralsDriver::compute

static void
CElectronRepulsionIntegralsDriver_compute_1(
          CElectronRepulsionIntegralsDriver& self,
          CAOFockMatrix&                     aoFockMatrix,
    const CAODensityMatrix&                  aoDensityMatrix,
    const CMolecule&                         molecule,
    const CMolecularBasis&                   aoBasis,
    const CScreeningContainer&               screeningContainer,
          bp::object                         py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    self.compute(aoFockMatrix, aoDensityMatrix, molecule, aoBasis,
                 screeningContainer, *comm_ptr);
}

static CScreeningContainer
CElectronRepulsionIntegralsDriver_compute_2(
          CElectronRepulsionIntegralsDriver& self,
    const ericut                             screeningScheme,
    const double                             threshold,
    const CMolecule&                         molecule,
    const CMolecularBasis&                   aoBasis)
{
    return self.compute(screeningScheme, threshold, molecule, aoBasis);
}

// Helper function for reduce_sum CAOFockMatrix object
    
static void
CAOFockMatrix_reduce_sum(CAOFockMatrix& self,
                         int32_t        rank,
                         int32_t        nodes,
                         bp::object     py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);
        
    self.reduce_sum(rank, nodes, *comm_ptr);
}
    
// Exports classes/functions in src/twoints to python

void export_twoints()
{
    // initialize numpy

    Py_Initialize();

    np::initialize();

    // fockmat enum class

    bp::enum_<fockmat> ("fockmat")
        .value("restjk",  fockmat::restjk )
        .value("restjkx", fockmat::restjkx)
        .value("restj",   fockmat::restj  )
        .value("restk",   fockmat::restk  )
        .value("restkx",  fockmat::restkx )
    ;

    // ericut enum class

    bp::enum_<ericut> ("ericut")
        .value("qq",     ericut::qq )
        .value("qqr",    ericut::qqr)
        .value("qqden",  ericut::qqden)
        .value("qqrden", ericut::qqrden)
    ;

    // CAOFockMatrix class

    bp::class_< CAOFockMatrix, std::shared_ptr<CAOFockMatrix> >
        (
            "AOFockMatrix",
            bp::init<>()
        )
        .def(bp::init<const CAODensityMatrix&>())
        .def(bp::init<const CAOFockMatrix&>())
        .def("__str__", &CAOFockMatrix_str)
        .def("to_numpy", &CAOFockMatrix_to_numpy)
        .def("from_numpy_list", &CAOFockMatrix_from_numpy_list)
        .staticmethod("from_numpy_list")
        .def("zero", &CAOFockMatrix::zero)
        .def("get_number_of_fock_matrices", &CAOFockMatrix::getNumberOfFockMatrices)
        .def("get_fock_type", &CAOFockMatrix::getFockType)
        .def("get_scale_factor", &CAOFockMatrix::getScaleFactor)
        .def("get_density_identifier", &CAOFockMatrix::getDensityIdentifier)
        .def("add_hcore", &CAOFockMatrix::addCoreHamiltonian)
        .def("reduce_sum", &CAOFockMatrix_reduce_sum)
        .def(bp::self == bp::other<CAOFockMatrix>())
    ;

    // CCauchySchwarzScreener class

    bp::class_< CCauchySchwarzScreener, std::shared_ptr<CCauchySchwarzScreener> >
        (
            "CauchySchwarzScreener",
            bp::init<>()
        )
        .def(bp::init<const CCauchySchwarzScreener&>())
        .def("get_threshold", &CCauchySchwarzScreener::getThreshold)
        .def("set_threshold", &CCauchySchwarzScreener::setThreshold)
        .def("get_screening_scheme", &CCauchySchwarzScreener::getScreeningScheme)
        .def(bp::self == bp::other<CCauchySchwarzScreener>())
    ;

    // CScreeningContainer class

    bp::class_< CScreeningContainer, std::shared_ptr<CScreeningContainer> >
        (
            "ScreeningContainer",
            bp::init<>()
        )
        .def(bp::init<const CScreeningContainer&>())
        .def("is_empty", &CScreeningContainer::isEmpty)
        .def("get_number_of_screeners", &CScreeningContainer::getNumberOfScreeners)
        .def("set_threshold", &CScreeningContainer::setThreshold)
        .def(bp::self == bp::other<CScreeningContainer>())
    ;

    // CElectronRepulsionIntegralsDriver class

    bp::class_< CElectronRepulsionIntegralsDriver,
                std::shared_ptr<CElectronRepulsionIntegralsDriver> >
        (
            "ElectronRepulsionIntegralsDriver",
            bp::init<const int32_t, const int32_t, MPI_Comm>()
        )
        .def("create", &CElectronRepulsionIntegralsDriver_create)
        .staticmethod("create")
        .def("compute", &CElectronRepulsionIntegralsDriver_compute_1)
        .def("compute", &CElectronRepulsionIntegralsDriver_compute_2)
    ;
}

} // bp_twoints namespace
