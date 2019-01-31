//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

#include <mpi.h>
#include <string>
#include <fstream>

#include "MolecularBasis.hpp"
#include "StringFormat.hpp"
#include "ErrorHandler.hpp"
#include "DenseMatrix.hpp"
#include "AODensityMatrix.hpp"
#include "DensityMatrixType.hpp"
#include "MolecularOrbitals.hpp"
#include "MolecularOrbitalsType.hpp"
#include "ExportGeneral.hpp"
#include "ExportMath.hpp"
#include "ExportOrbData.hpp"

namespace py = pybind11;

namespace vlx_orbdata { // vlx_orbdata namespace

// Helper function for broadcasting CMolecularBasis object

static void
CMolecularBasis_broadcast(CMolecularBasis& self,
                          int32_t          rank,
                          py::object       py_comm)
{
    MPI_Comm* comm_ptr = vlx_general::get_mpi_comm(py_comm);

    self.broadcast(rank, *comm_ptr);
}
    
// Helper function for printing CAODensityMatrix

static std::string
CAODensityMatrix_str (const CAODensityMatrix& self)
{
    return self.getString();
}

// Helper function for converting CAODensityMatrix to numpy array

static py::array_t<double>
CAODensityMatrix_total_density_to_numpy(const CAODensityMatrix& self,
                                        const int32_t iDensityMatrix)
{
    return vlx_general::pointer_to_numpy(self.totalDensity(iDensityMatrix),
                                         self.getNumberOfRows(iDensityMatrix),
                                         self.getNumberOfColumns(iDensityMatrix));
}

static py::array_t<double>
CAODensityMatrix_alpha_density_to_numpy(const CAODensityMatrix& self,
                                        const int32_t iDensityMatrix)
{
    return vlx_general::pointer_to_numpy(self.alphaDensity(iDensityMatrix),
                                         self.getNumberOfRows(iDensityMatrix),
                                         self.getNumberOfColumns(iDensityMatrix));
}

static py::array_t<double>
CAODensityMatrix_beta_density_to_numpy(const CAODensityMatrix& self,
                                       const int32_t iDensityMatrix)
{
    return vlx_general::pointer_to_numpy(self.betaDensity(iDensityMatrix),
                                         self.getNumberOfRows(iDensityMatrix),
                                         self.getNumberOfColumns(iDensityMatrix));
}

// Helper function for CAODensityMatrix constructor

static std::shared_ptr<CAODensityMatrix>
CAODensityMatrix_from_numpy_list(const std::vector<py::array_t<double>>& arrays,
                                 const denmat                            den_type)
{
    std::vector<CDenseMatrix> dmat;

    for (size_t i = 0; i < arrays.size(); i++)
    {
        auto mp = vlx_math::CDenseMatrix_from_numpy(arrays[i]);

        dmat.push_back(*mp);
    }

    return std::shared_ptr<CAODensityMatrix>(new CAODensityMatrix(dmat, den_type));
}

// Helper function for broadcasting CAODensityMatrix object
    
static void
CAODensityMatrix_broadcast(CAODensityMatrix& self,
                           int32_t           rank,
                           py::object        py_comm)
{
    MPI_Comm* comm_ptr = vlx_general::get_mpi_comm(py_comm);
        
    self.broadcast(rank, *comm_ptr);
}
    
// Helper function for printing CMolecularOrbitals

static std::string
CMolecularOrbitals_str (const CMolecularOrbitals& self)
{
    return self.getString();
}

// Helper function for converting CMolecularOrbitals to numpy array

static py::array_t<double>
CMolecularOrbitals_alpha_orbitals_to_numpy(const CMolecularOrbitals& self)
{
    return vlx_general::pointer_to_numpy(self.alphaOrbitals(),
                                         self.getNumberOfRows(),
                                         self.getNumberOfColumns());
}

static py::array_t<double>
CMolecularOrbitals_beta_orbitals_to_numpy(const CMolecularOrbitals& self)
{
    return vlx_general::pointer_to_numpy(self.betaOrbitals(),
                                         self.getNumberOfRows(),
                                         self.getNumberOfColumns());
}

static py::array_t<double>
CMolecularOrbitals_alpha_energies_to_numpy(const CMolecularOrbitals& self)
{
    py::list ea;

    for (int32_t i = 0; i < self.getNumberOfColumns(); i++)
    {
        ea.append(self.alphaEnergies()[i]);
    }

    return py::array_t<double>(ea);
}

static py::array_t<double>
CMolecularOrbitals_beta_energies_to_numpy(const CMolecularOrbitals& self)
{
    py::list eb;

    for (int32_t i = 0; i < self.getNumberOfColumns(); i++)
    {
        eb.append(self.betaEnergies()[i]);
    }

    return py::array_t<double>(eb);
}

// Helper function for CMolecularOrbitals constructor
    
static std::shared_ptr<CMolecularOrbitals>
CMolecularOrbitals_from_numpy_list(const std::vector<py::array_t<double>>& mol_orbs,
                                   const std::vector<py::array_t<double>>& eig_vals,
                                   const molorb                            orbs_type)
{
    std::vector<CDenseMatrix> cmos;

    for (size_t i = 0; i < mol_orbs.size(); i++)
    {
        auto mp = vlx_math::CDenseMatrix_from_numpy(mol_orbs[i]);

        cmos.push_back(*mp);
    }

    std::vector<CMemBlock<double>> ceigs;

    for (size_t i = 0; i < eig_vals.size(); i++)
    {
        const py::array_t<double>& arr = eig_vals[i];

        std::string errdim("MolecularOrbitals eigenvalues: need 1D numpy arrays");

        errors::assertMsgCritical(arr.ndim() == 1, errdim);

        if (arr.data() == nullptr || arr.size() == 0)
        {
            return std::shared_ptr<CMolecularOrbitals>(new CMolecularOrbitals());
        }

        std::vector<double> vec(arr.data(), arr.data() + arr.size());

        ceigs.push_back(CMemBlock<double>(vec));
    }

    return std::shared_ptr<CMolecularOrbitals>(
            new CMolecularOrbitals(cmos, ceigs, orbs_type)
            );
}

// Exports classes/functions in src/orbdata to python

void export_orbdata(py::module& m)
{
    // CBasisFunction class

    py::class_< CBasisFunction, std::shared_ptr<CBasisFunction> >
        (
            m, "BasisFunction"
        )
        .def(py::init<>())
        .def(py::init<const std::vector<double>&,
                      const std::vector<double>&,
                      const int32_t>())
        .def("normalize", &CBasisFunction::normalize)
    ;

    // CAtomBasis class

    py::class_< CAtomBasis, std::shared_ptr<CAtomBasis> >
        (
            m, "AtomBasis"
        )
        .def(py::init<>())
        .def("add_basis_function", &CAtomBasis::addBasisFunction)
        .def("set_elemental_id", &CAtomBasis::setIdElemental)
    ;

    // CMolecularBasis class

    py::class_< CMolecularBasis, std::shared_ptr<CMolecularBasis> >
        (
            m, "MolecularBasis"
        )
        .def(py::init<>())
        .def("get_string", &CMolecularBasis::printBasis)
        .def("set_label", &CMolecularBasis::setLabel)
        .def("get_label", &CMolecularBasis::getLabel)
        .def("broadcast", &CMolecularBasis_broadcast)
        .def("get_valence_basis", &CMolecularBasis::reduceToValenceBasis)
        .def("add_atom_basis", &CMolecularBasis::addAtomBasis)
        .def(py::self == py::self)
    ;

    // denmat enum class

    py::enum_<denmat> (m, "denmat")
        .value("rest",   denmat::rest  )
        .value("unrest", denmat::unrest)
    ;

    // CAODensityMatrix class

    py::class_< CAODensityMatrix, std::shared_ptr<CAODensityMatrix> >
        (
            m, "AODensityMatrix"
        )
        .def(py::init<>())
        .def(py::init<const CAODensityMatrix&>())
        .def(py::init(&CAODensityMatrix_from_numpy_list))
        .def("__str__", &CAODensityMatrix_str)
        .def("total_to_numpy", &CAODensityMatrix_total_density_to_numpy)
        .def("alpha_to_numpy", &CAODensityMatrix_alpha_density_to_numpy)
        .def("beta_to_numpy", &CAODensityMatrix_beta_density_to_numpy)
        .def("number_of_density_matrices",
                &CAODensityMatrix::getNumberOfDensityMatrices)
        .def("get_density_type", &CAODensityMatrix::getDensityType)
        .def("sub", &CAODensityMatrix::sub)
        .def("broadcast", &CAODensityMatrix_broadcast)
        .def(py::self == py::self)
    ;

    // molorb enum class

    py::enum_<molorb> (m, "molorb")
        .value("rest",   molorb::rest  )
        .value("unrest", molorb::unrest)
    ;

    // CMolecularOrbitals class

    py::class_< CMolecularOrbitals, std::shared_ptr<CMolecularOrbitals> >
        (
            m, "MolecularOrbitals"
         )
        .def(py::init<>())
        .def(py::init<const CMolecularOrbitals&>())
        .def(py::init(&CMolecularOrbitals_from_numpy_list))
        .def("__str__", &CMolecularOrbitals_str)
        .def("alpha_to_numpy", &CMolecularOrbitals_alpha_orbitals_to_numpy)
        .def("beta_to_numpy", &CMolecularOrbitals_beta_orbitals_to_numpy)
        .def("ea_to_numpy", &CMolecularOrbitals_alpha_energies_to_numpy)
        .def("eb_to_numpy", &CMolecularOrbitals_beta_energies_to_numpy)
        .def("get_orbitals_type", &CMolecularOrbitals::getOrbitalsType)
        .def("get_ao_density",
             (CAODensityMatrix (CMolecularOrbitals::*)(const int32_t) const)
             &CMolecularOrbitals::getAODensity)
        .def("get_ao_density",
             (CAODensityMatrix (CMolecularOrbitals::*)(const int32_t,
                                                       const int32_t) const)
             &CMolecularOrbitals::getAODensity)
        .def("insert", &CMolecularOrbitals::insert)
        .def(py::self == py::self)
    ;
}

} // vlx_orbdata namespace
