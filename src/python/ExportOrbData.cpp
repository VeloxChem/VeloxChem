//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <boost/python.hpp>
#include <mpi.h>

#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "DenseMatrix.hpp"
#include "AODensityMatrix.hpp"
#include "DensityMatrixType.hpp"
#include "MolecularOrbitals.hpp"
#include "MolecularOrbitalsType.hpp"
#include "ExportGeneral.hpp"
#include "ExportMath.hpp"
#include "ExportOrbData.hpp"

namespace bp = boost::python;

namespace np = boost::python::numpy;

namespace bp_orbdata { // bp_orbdata namespace

// Helper function for broadcasting CMolecularBasis object

static void
CMolecularBasis_broadcast(CMolecularBasis& self,
                          int32_t          rank,
                          bp::object       py_comm)
{
    MPI_Comm* comm_ptr = bp_general::get_mpi_comm(py_comm);

    self.broadcast(rank, *comm_ptr);
}
    
// Helper function for printing CAODensityMatrix

static std::string
CAODensityMatrix_str (const CAODensityMatrix& self)
{
    return self.getString();
}

// Helper function for converting CAODensityMatrix to numpy array

static np::ndarray
CAODensityMatrix_total_density_to_numpy(const CAODensityMatrix& self,
                                        const int32_t iDensityMatrix)
{
    return bp_general::pointer_to_numpy(self.totalDensity(iDensityMatrix),
                                        self.getNumberOfRows(iDensityMatrix),
                                        self.getNumberOfColumns(iDensityMatrix));
}

static np::ndarray
CAODensityMatrix_alpha_density_to_numpy(const CAODensityMatrix& self,
                                        const int32_t iDensityMatrix)
{
    return bp_general::pointer_to_numpy(self.alphaDensity(iDensityMatrix),
                                        self.getNumberOfRows(iDensityMatrix),
                                        self.getNumberOfColumns(iDensityMatrix));
}

static np::ndarray
CAODensityMatrix_beta_density_to_numpy(const CAODensityMatrix& self,
                                       const int32_t iDensityMatrix)
{
    return bp_general::pointer_to_numpy(self.betaDensity(iDensityMatrix),
                                        self.getNumberOfRows(iDensityMatrix),
                                        self.getNumberOfColumns(iDensityMatrix));
}

// Helper function for CAODensityMatrix constructor

static std::shared_ptr<CAODensityMatrix>
CAODensityMatrix_from_numpy_list(const bp::list& arr_list,
                                 const denmat    den_type)
{
    std::vector<CDenseMatrix> dmat;

    for (int i = 0; i < bp::len(arr_list); i++)
    {
        np::ndarray arr = np::array(arr_list[i]);

        std::shared_ptr<CDenseMatrix> mp = bp_math::CDenseMatrix_from_numpy(arr);

        dmat.push_back(*mp);
    }

    return std::shared_ptr<CAODensityMatrix>(new CAODensityMatrix(dmat, den_type));
}
    
// Helper function for printing CMolecularOrbitals

static std::string
CMolecularOrbitals_str (const CMolecularOrbitals& self)
{
    return self.getString();
}

// Helper function for converting CMolecularOrbitals to numpy array

static np::ndarray
CMolecularOrbitals_alpha_orbitals_to_numpy(const CMolecularOrbitals& self)
{
    return bp_general::pointer_to_numpy(self.alphaOrbitals(),
                                        self.getNumberOfRows(),
                                        self.getNumberOfColumns());
}

static np::ndarray
CMolecularOrbitals_beta_orbitals_to_numpy(const CMolecularOrbitals& self)
{
    return bp_general::pointer_to_numpy(self.betaOrbitals(),
                                        self.getNumberOfRows(),
                                        self.getNumberOfColumns());
}

static np::ndarray
CMolecularOrbitals_alpha_energies_to_numpy(const CMolecularOrbitals& self)
{
    bp::list ea;

    for (int32_t i = 0; i < self.getNumberOfColumns(); i++)
    {
        ea.append(self.alphaEnergies()[i]);
    }

    return np::array(ea);
}

static np::ndarray
CMolecularOrbitals_beta_energies_to_numpy(const CMolecularOrbitals& self)
{
    bp::list eb;

    for (int32_t i = 0; i < self.getNumberOfColumns(); i++)
    {
        eb.append(self.betaEnergies()[i]);
    }

    return np::array(eb);
}

// Helper function for CMolecularOrbitals constructor
    
static std::shared_ptr<CMolecularOrbitals>
CMolecularOrbitals_from_numpy_list(const bp::list& orbs_list,
                                   const bp::list& eigs_list,
                                   const molorb    orbs_type)
{
    std::vector<CDenseMatrix> cmos;
    
    for (int i = 0; i < bp::len(orbs_list); i++)
    {
        np::ndarray arr = np::array(orbs_list[i]);

        std::shared_ptr<CDenseMatrix> mp = bp_math::CDenseMatrix_from_numpy(arr);
        
        cmos.push_back(*mp);
    }

    std::vector<CMemBlock<double>> ceigs;
    
    for (int i = 0; i < bp::len(eigs_list); i++)
    {
        np::ndarray arr = np::array(eigs_list[i]);
        
        const double* data = reinterpret_cast<double*>(arr.get_data());
        
        if (data == nullptr)
        {
            return std::shared_ptr<CMolecularOrbitals>(new CMolecularOrbitals());
        }
        
        auto size = static_cast<int32_t>(arr.shape(0));
        
        if (size == 0)
        {
            return std::shared_ptr<CMolecularOrbitals>(new CMolecularOrbitals());
        }
        
        std::vector<double> vec (data, data + size);
        
        ceigs.push_back(CMemBlock<double>(vec));
    }
    
    return std::shared_ptr<CMolecularOrbitals>(
            new CMolecularOrbitals(cmos, ceigs, orbs_type)
            );
}

// Helper function for restricted density matrix from CMolecularOrbitals
    
static CAODensityMatrix
CMolecularOrbitals_get_rest_density(const CMolecularOrbitals& self,
                                    const CMolecule&          molecule)
{
    auto nelec = molecule.getNumberOfElectrons();
    
    return self.getAODensity(nelec);
}

// Exports classes/functions in src/orbdata to python

void export_orbdata()
{
    // initialize numpy

    Py_Initialize();

    np::initialize();

    // CMolecularBasis class

    bp::class_< CMolecularBasis, std::shared_ptr<CMolecularBasis> >
        (
            "MolecularBasis",
            bp::init<>()
        )
        .def("get_label", &CMolecularBasis::getLabel)
        .def("broadcast", &CMolecularBasis_broadcast)
        .def("print_basis", &CMolecularBasis::printBasis)
        .def("get_valence_basis", &CMolecularBasis::reduceToValenceBasis)
    ;

    // denmat enum class

    bp::enum_<denmat> ("denmat")
        .value("rest",   denmat::rest  )
        .value("unrest", denmat::unrest)
    ;

    // CAODensityMatrix class

    bp::class_< CAODensityMatrix, std::shared_ptr<CAODensityMatrix> >
        (
            "AODensityMatrix",
            bp::init<>()
        )
        .def(bp::init<const CAODensityMatrix&>())
        .def("__str__", &CAODensityMatrix_str)
        .def("total_to_numpy", &CAODensityMatrix_total_density_to_numpy)
        .def("alpha_to_numpy", &CAODensityMatrix_alpha_density_to_numpy)
        .def("beta_to_numpy", &CAODensityMatrix_beta_density_to_numpy)
        .def("from_numpy_list", &CAODensityMatrix_from_numpy_list)
        .staticmethod("from_numpy_list")
        .def("get_number_of_density_matrices",
                &CAODensityMatrix::getNumberOfDensityMatrices)
        .def("get_density_type", &CAODensityMatrix::getDensityType)
        .def(bp::self == bp::other<CAODensityMatrix>())
    ;

    // molorb enum class

    bp::enum_<molorb> ("molorb")
        .value("rest",   molorb::rest  )
        .value("unrest", molorb::unrest)
    ;

    // CMolecularOrbitals class

    bp::class_< CMolecularOrbitals, std::shared_ptr<CMolecularOrbitals> >
        (
            "MolecularOrbitals",
            bp::init<>()
         )
        .def(bp::init<const CMolecularOrbitals&>())
        .def("__str__", &CMolecularOrbitals_str)
        .def("alpha_to_numpy", &CMolecularOrbitals_alpha_orbitals_to_numpy)
        .def("beta_to_numpy", &CMolecularOrbitals_beta_orbitals_to_numpy)
        .def("ea_to_numpy", &CMolecularOrbitals_alpha_energies_to_numpy)
        .def("eb_to_numpy", &CMolecularOrbitals_beta_energies_to_numpy)
        .def("from_numpy_list", &CMolecularOrbitals_from_numpy_list)
        .staticmethod("from_numpy_list")
        .def("get_orbitals_type", &CMolecularOrbitals::getOrbitalsType)
        .def(bp::self == bp::other<CMolecularOrbitals>())
        .def("get_rest_density", &CMolecularOrbitals_get_rest_density)
        .def("insert", &CMolecularOrbitals::insert)
    ;
}

} // bp_orbdata namespace
