//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

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

std::string
CAOFockMatrix_str (const CAOFockMatrix& self)
{
    return self.getString();
}

// Helper function for converting CAOFockMatrix to numpy array

np::ndarray
CAOFockMatrix_to_numpy(const CAOFockMatrix& self,
                       const int32_t iFockMatrix)
{
    return bp_general::pointer_to_numpy(self.getFock(iFockMatrix),
                                        self.getNumberOfRows(iFockMatrix),
                                        self.getNumberOfColumns(iFockMatrix));
}

// Helper function for converting a list of numpy array to CAOFockMatrix

CAOFockMatrix
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

        fmat.push_back(bp_math::CDenseMatrix_from_numpy(arr));
        types.push_back(type);
        factors.push_back(factor);
        ids.push_back(id);
    }

    return CAOFockMatrix(fmat, types, factors, ids);
}

// Helper function for creating a CElectronRepulsionIntegralsDriver object

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

void
CElectronRepulsionIntegralsDriver_compute_1(
          CElectronRepulsionIntegralsDriver& self,
          CAOFockMatrix&                     aoFockMatrix,
    const CAODensityMatrix&                  aoDensityMatrix,
    const CMolecule&                         molecule,
    const CMolecularBasis&                   aoBasis,
    const CScreeningContainer&               screeningContainer,
          COutputStream&                     oStream,
          bp::object                         py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    self.compute(aoFockMatrix, aoDensityMatrix,
                 molecule, aoBasis, screeningContainer,
                 oStream, *comm_ptr);
}

CScreeningContainer
CElectronRepulsionIntegralsDriver_compute_2(
          CElectronRepulsionIntegralsDriver& self,
    const ericut                             screeningScheme,
    const double                             threshold,
    const CMolecule&                         molecule,
    const CMolecularBasis&                   aoBasis,
          COutputStream&                     oStream,
          bp::object                         py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    return self.compute(screeningScheme, threshold,
                        molecule, aoBasis,
                        oStream, *comm_ptr);
}

// Exports classes/functions in src/twoints to python

void export_twoints()
{
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
        .value("qq",  ericut::qq )
        .value("qqr", ericut::qqr)
    ;

    // CAOFockMatrix class

    bp::class_< CAOFockMatrix, std::shared_ptr<CAOFockMatrix> >
        (
            "AOFockMatrix",
            bp::init<
                const std::vector<CDenseMatrix>&,
                const std::vector<fockmat>&,
                const std::vector<double>&,
                const std::vector<int32_t>&
                >()
        )
        .def(bp::init<>())
        .def(bp::init<const CAODensityMatrix&>())
        .def(bp::init<const CAOFockMatrix&>())
        .def("__str__", &CAOFockMatrix_str)
        .def("to_numpy", &CAOFockMatrix_to_numpy)
        .def("from_numpy_list", &CAOFockMatrix_from_numpy_list)
        .staticmethod("from_numpy_list")
        .def("zero", &CAOFockMatrix::zero)
        .def("get_fock_type", &CAOFockMatrix::getFockType)
        .def("get_scale_factor", &CAOFockMatrix::getScaleFactor)
        .def("get_density_identifier", &CAOFockMatrix::getDensityIdentifier)
        .def(bp::self == bp::other<CAOFockMatrix>())
    ;

    // CScreeningContainer class

    bp::class_< CScreeningContainer, std::shared_ptr<CScreeningContainer> >
        (
            "ScreeningContainer",
            bp::init<
                const CVecMemBlock<double>&,
                const CVecMemBlock<double>&,
                const CGtoPairsContainer&,
                const CGtoPairsContainer&,
                const ericut,
                const double
                >()
        )
        .def(bp::init<>())
        .def(bp::init<const CScreeningContainer&>())
        .def("is_empty", &CScreeningContainer::isEmpty)
        .def("get_number_of_screeners", &CScreeningContainer::getNumberOfScreeners)
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
