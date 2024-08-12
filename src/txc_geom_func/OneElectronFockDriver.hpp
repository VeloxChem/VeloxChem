#ifndef OneElectronExpectationDriver_hpp
#define OneElectronExpectationDriver_hpp

#include <cstdint>
#include <vector>

#include "GeomJob.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "OneElectronExpectation.hpp"
#include "OverlapGeomX00Driver.hpp"
#include "OverlapGeomX0YDriver.hpp"

template <class T, class U, class W>
class COneElectronFockDriver
{
   public:
    COneElectronFockDriver() = default;

    auto
    compute(CGeomJob<T, U, W>& result, const CMolecularBasis& basis, const CMolecule& molecule) const -> void
    {
        task.allocate(basis);

        auto op_case = _set_op_case(task.get_operators());

        if (op_case == 1)
        {
            if (geom_order == T({1, 0, 0, 0}))
            {
                CKineticEnergyGeomX00Driver<1> kin_fck_drv;

                kin_fck_drv.compute(task.get_results(), basis, molecule);

                CNuclearPotentialGeomX00Driver<1> nuc_fck_drv;

                nuc_fck_drv.compute(task.get_results(), basis, molecule);
            }

            if (geom_order == T({2, 0, 0, 0}))
            {
                CKineticEnergyGeomX00Driver<2> kin_fck_drv;

                kin_fck_drv.compute(task.get_results(), basis, molecule);

                CNuclearPotentialGeomX00Driver<2> nuc_fck_drv;

                nuc_fck_drv.compute(task.get_results(), basis, molecule);
            }

            if (geom_order == T({1, 0, 1, 0}))
            {
                CKineticEnergyGeomX0YDriver<1, 1> kin_fck_drv;

                kin_fck_drv.compute(task.get_results(), basis, molecule);

                CNuclearPotentialGeomX0YDriver<1, 1> nuc_fck_drv;

                nuc_fck_drv.compute(task.get_results(), basis, molecule);
            }

            if (geom_order == T({0, 1, 0, 0}))
            {
                CNuclearPotentialGeom0X0Driver<1> nuc_fck_drv;

                nuc_fck_drv.compute(task.get_results(), basis, molecule);
            }

            if (geom_order == T({0, 2, 0, 0}))
            {
                CNuclearPotentialGeom0X0Driver<2> nuc_fck_drv;

                nuc_fck_drv.compute(task.get_results(), basis, molecule);
            }
        }

        if (op_case == 2)
        {
            for (const auto& tarch : task.get_tasks())
            {
                const auto geom_order = tarch.get_archetype();

                if (geom_order == T({1, 0, 0, 0}))
                {
                    COverlapGeomX00Driver<1> fck_drv;

                    exp_drv.compute(task.get_results(), basis, molecule);
                }

                if (geom_order == T({2, 0, 0, 0}))
                {
                    COverlapGeomX00Driver<2> fck_drv;

                    exp_drv.compute(task.get_results(), basis, molecule);
                }

                if (geom_order == T({1, 0, 1, 0}))
                {
                    COverlapGeomX0YDriver<1, 1> fck_drv;

                    exp_drv.compute(task.get_results(), basis, molecule);
                }
            }
        }

        if (op_case == 3)
        {
            for (const auto& tarch : task.get_tasks())
            {
                const auto geom_order = tarch.get_archetype();

                if (geom_order == T({1, 0, 0, 0}))
                {
                    CElectricDipoleMomentumGeomX00Driver<1> exp_drv;

                    exp_drv.compute(task.get_results(), densities, basis, molecule);
                }
            }
        }
    }

    auto
    _set_op_case(const std::vector<std::string>& op_labels) -> int
    {
        if (op_labels.size() == 1)
        {
            if (op_labels[0].compare("el_dipole") == 0)
            {
                return 3;
            }
            else if (op_labels[0].compare("overlap") == 0)
            {
                return 2;
            }
            else
            {
                return -1;
            }
        }
        else if (op_labels.size() == 2)
        {
            if (op_labels[0].compare("kinetic") == 0)
            {
                if (op_labels[1].compare("nucel") == 0)
                {
                    return 1;
                }
                else
                {
                    return -1;
                }
            }
            else if (op_labels[0].compare("nucel") == 0)
            {
                if (op_labels[1].compare("kinetic") == 0)
                {
                    return 1;
                }
                else
                {
                    return -1;
                }
            }
            else
            {
                return -1;
            }
        }
    }
};

#endif /* OneElectronExpectationDriver_hpp */
