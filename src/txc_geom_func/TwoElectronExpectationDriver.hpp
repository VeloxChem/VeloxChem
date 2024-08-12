#ifndef TwoElectronExpectationDriver_hpp
#define TwoElectronExpectationDriver_hpp

template <class T, class U, class W>
class CTwoElectronExpectationDriver
{
   public:
    CTwoElectronExpectationDriver() = default;

    auto
    compute(CGeomJob<T, U, W>&     result,
            CMatrices&             densities,
            const CMolecularBasis& basis,
            const CMolecule&       molecule) const -> void
    {
        task.allocate(basis);

        // set basis set partitioning into batches

        CT4CScreener screener;

        screener.partition(basis, molecule, "eri");

        for (const auto& tarch : task.getTasks())
        {
            const auto geom_order = tarch.getArchetype();

            if (geom_order == T({1, 0, 0, 0}))
            {
                CFockGeomX000Driver<1> fock_drv;

                fock_drv.compute_values(tarch, task.get_results(), densities, basis, molecule, screener);
            }

            if (geom_order == T({2, 0, 0, 0}))
            {
                CFockGeomX000Driver<2> fock_drv;

                fock_drv.compute_values(tarch, task.get_results(), densities, basis, molecule, screener);
            }

            if (geom_order == T({1, 1, 0, 0}))
            {
                CFockGeomXY00Driver<1, 1> fock_drv;

                fock_drv.compute_values(tarch, task.get_results(), densities, basis, molecule, screener);
            }

            if (geom_order == T({1, 0, 1, 0}))
            {
                CFockGeomX0Y0Driver<1, 1> fock_drv;

                fock_drv.compute_values(tarch, task.get_results(), densities, basis, molecule);
            }
        }
    }
};

#endif /* TwoElectronExpectationDriver_hpp */
