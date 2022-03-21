#include <mpi.h>

#include <cstdint>
#include <cstdio>
#include <iostream>
#include <vector>

#include "MemBlock.hpp"
#include "Molecule.hpp"
#include "MpiFunc.hpp"

int
main(int argc, char** argv)
{
    mpi::init(argc, argv);

    // CMemBlock<int32_t> ma({1, 2, 3, 9});

    // ma.broadcast(0, MPI_COMM_WORLD);

    std::cout << "My rank is " << mpi::rank(MPI_COMM_WORLD) << std::endl;
    auto mol = CMolecule();

    if (mpi::rank(MPI_COMM_WORLD) == mpi::master())
    {
        std::cout << "Rank " << mpi::rank(MPI_COMM_WORLD) << " is creating the molecule" << std::endl;
        std::vector coords{0.0, 0.0, 0.0, 0.0, 1.4, -1.4, 0.0, 1.1, 1.1};

        std::vector charges{8.0, 1.0, 1.0};

        std::vector masses{15.994915, 1.007825, 1.007825};

        std::vector<std::string> labels{{"O"}, {"H"}, {"H"}};

        std::vector<int32_t> idselem{8, 1, 1};

        CMolecule foo(coords, charges, masses, labels, idselem);

        mol = foo;
        std::cout << "Copy-assigned foo to mol" << std::endl;

        // create molecule
        mol.setCharge(0.0);
        mol.setMultiplicity(1);
    }

    mol.broadcast(mpi::rank(MPI_COMM_WORLD), MPI_COMM_WORLD);

    mpi::finalize();
}
