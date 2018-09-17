//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "ElectronRepulsionIntegralsDriverTest.hpp"

#include "mpi.h"

#include "ElectronRepulsionIntegralsDriver.hpp"
#include "MolecularBasisSetter.hpp"
#include "MoleculeSetter.hpp"
#include "CheckFunctions.hpp"

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeSSSSForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoPairsBlock bpairs(bgtos, 1.0e-13);
    
    CMemBlock<double> fints(1);
    
    eridrv.compute(fints.data(), bpairs.pick(0), bpairs.pick(0));
    
    std::vector<double> r0000vals{ 1.629232440467330};
    
    vlxtest::compare(r0000vals, fints.data());
    
    eridrv.compute(fints.data(), bpairs.pick(1), bpairs.pick(1));
    
    std::vector<double> r0101vals{ 1.296490921648330};
    
    vlxtest::compare(r0101vals, fints.data());
    
    eridrv.compute(fints.data(), bpairs.pick(1), bpairs.pick(8));
    
    std::vector<double> r0123vals{ 0.569573705357336};
    
    vlxtest::compare(r0123vals, fints.data());
    
    eridrv.compute(fints.data(), bpairs.pick(3), bpairs.pick(5));
    
    std::vector<double> r0312vals{ 0.091957062417462};
    
    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeSSSPForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 0);
    
    CGtoPairsBlock bpairs(agtos, 1.0e-13);
    
    CGtoBlock cgtos(mlih, mbas, 0);
    
    CGtoBlock dgtos(mlih, mbas, 1);
    
    CGtoPairsBlock kpairs(cgtos, dgtos, 1.0e-13);
    
    CMemBlock<double> fints(3);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));
    
    std::vector<double> r0000vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000};
    
    vlxtest::compare(r0000vals,fints.data());
    
    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));
    
    std::vector<double> r0101vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000};
    
    vlxtest::compare(r0101vals, fints.data());
    
    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(11));
    
    std::vector<double> r0123vals{-0.048499014664650, -0.088914860218526, -0.032332676443100};
    
    vlxtest::compare(r0123vals, fints.data());
    
    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(6));
    
    std::vector<double> r0312vals{-0.065342004894742, -0.119793675640361, -0.043561336596495};
    
    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeSPSSForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 0);
    
    CGtoBlock bgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(agtos, bgtos, 1.0e-13);
    
    CGtoBlock dgtos(mlih, mbas, 0);
    
    CGtoPairsBlock kpairs(dgtos, 1.0e-13);
    
    CMemBlock<double> fints(3);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));
    
    std::vector<double> r0000vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000};
    
    vlxtest::compare(r0000vals, fints.data());
    
    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));
    
    std::vector<double> r0101vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000};
    
    vlxtest::compare(r0101vals, fints.data());
    
    eridrv.compute(fints.data(), bpairs.pick(11), kpairs.pick(1));
    
    std::vector<double> r2301vals{-0.048499014664650, -0.088914860218526, -0.032332676443100};
    
    vlxtest::compare(r2301vals, fints.data());
    
    eridrv.compute(fints.data(), bpairs.pick(6), kpairs.pick(3));
    
    std::vector<double> r1203vals{-0.065342004894742, -0.119793675640361, -0.043561336596495};
    
    vlxtest::compare(r1203vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeSPSPForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 0);
    
    CGtoBlock bgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(agtos, bgtos, 1.0e-13);
    
    CMemBlock<double> fints(9);
    
    eridrv.compute(fints.data(), bpairs.pick(0), bpairs.pick(0));

    std::vector<double> r0000vals{ 0.271300324109232,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.271300324109232,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.271300324109232};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), bpairs.pick(1));

    std::vector<double> r0101vals{ 0.158211411411533,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.158211411411533,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.158211411411533};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), bpairs.pick(11));

    std::vector<double> r0123vals{ 0.014876488624304, -0.023320174063975, -0.008480063295991,
                                  -0.023320174063975, -0.015157068882331, -0.015546782709317,
                                  -0.008480063295991, -0.015546782709317,  0.021943208037630};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), bpairs.pick(6));

    std::vector<double> r0312vals{ 0.056736244044378,  0.081977329675453,  0.029809938063801,
                                   0.081977329675453,  0.162313108020340,  0.054651553116969,
                                   0.029809938063801,  0.054651553116969,  0.031894628991210};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeSSPPForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 0);
    
    CGtoPairsBlock bpairs(agtos, 1.0e-13);
    
    CGtoBlock dgtos(mlih, mbas, 1);
    
    CGtoPairsBlock kpairs(dgtos, 1.0e-13);
    
    CMemBlock<double> fints(9);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 1.357416654829900,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  1.357416654829900,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  1.357416654829900};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 1.043481505026390,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  1.043481505026390,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  1.043481505026390};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(8));

    std::vector<double> r0123vals{ 0.472688274973429,  0.033593694100454,  0.012215888763801,
                                   0.033593694100454,  0.515952881011892,  0.022395796066969,
                                   0.012215888763801,  0.022395796066969,  0.462508367670261};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(5));

    std::vector<double> r0312vals{ 0.023339095978464, -0.095784363951802, -0.034830677800655,
                                  -0.095784363951802, -0.100019554565523, -0.063856242634535,
                                  -0.034830677800655, -0.063856242634535,  0.052364660812343};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputePPSSForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(agtos, 1.0e-13);
    
    CGtoBlock dgtos(mlih, mbas, 0);
    
    CGtoPairsBlock kpairs(dgtos, 1.0e-13);
    
    CMemBlock<double> fints(9);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 1.357416654829900,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  1.357416654829900,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  1.357416654829900};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 1.043481505026390,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  1.043481505026390,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  1.043481505026390};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(8));

    std::vector<double> r0123vals{ 0.527902875740139,  0.037278021172033,  0.013555644062557,
                                   0.037278021172033,  0.575912448461696,  0.024852014114689,
                                   0.013555644062557,  0.024852014114689,  0.516606505688008};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(5));

    std::vector<double> r0312vals{ 0.032590929327399, -0.067864620590098, -0.024678043850945,
                                  -0.067864620590098, -0.054810475978031, -0.045243080393399,
                                  -0.024678043850945, -0.045243080393399,  0.053155965869853};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputePPSPForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(agtos, 1.0e-13);
    
    CGtoBlock cgtos(mlih, mbas, 0);
    
    CGtoBlock dgtos(mlih, mbas, 1);
    
    CGtoPairsBlock kpairs(cgtos, dgtos, 1.0e-13);
    
    CMemBlock<double> fints(27);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(11));

    std::vector<double> r0123vals{-0.024740456719223, -0.079866951735308, -0.029042527903748,
                                   0.006120732901895, -0.011000649215887, -0.007422660648759,
                                   0.002225721055234, -0.007422660648759,  0.006712518241378,
                                   0.006120732901895, -0.011000649215887, -0.007422660648759,
                                  -0.057903022654361, -0.071646093782930, -0.038602015102907,
                                  -0.007422660648759, -0.007333766143925,  0.012306283442527,
                                   0.002225721055234, -0.007422660648759,  0.006712518241378,
                                  -0.007422660648759, -0.007333766143925,  0.012306283442527,
                                  -0.040189855197096, -0.073681401194676, -0.014244346707131};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(6));

    std::vector<double> r0312vals{-0.025187805950051, -0.039734523702261, -0.014448917709913,
                                   0.027147577114811,  0.097741908029784,  0.032527520363146,
                                   0.009871846223568,  0.032527520363146,  0.020119416254094,
                                   0.063991863435128,  0.077645024582338,  0.032527520363146,
                                   0.041163878682117,  0.069023990377715,  0.027442585788078,
                                   0.032527520363146,  0.051763349721559,  0.036885596465839,
                                   0.023269768521865,  0.032527520363146,  0.000022532806649,
                                   0.032527520363146,  0.065161272019856,  0.000041310145522,
                                  -0.036458613093572, -0.066840790671549, -0.026648694985836};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeSPPPForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 0);
    
    CGtoBlock bgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(agtos, bgtos, 1.0e-13);
    
    CGtoBlock dgtos(mlih, mbas, 1);
    
    CGtoPairsBlock kpairs(dgtos, 1.0e-13);
    
    CMemBlock<double> fints(27);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));

    std::vector<double> r0000vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));

    std::vector<double> r0101vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(8));

    std::vector<double> r0123vals{ 0.027005659425976, -0.006710881773728, -0.002440320644992,
                                  -0.006710881773728,  0.063368692728713,  0.008139179236825,
                                  -0.002440320644992,  0.008139179236825,  0.043945651368109,
                                   0.087349676872220,  0.012062933467287,  0.008139179236825,
                                   0.012062933467287,  0.078336635411377,  0.008041955644858,
                                   0.008139179236825,  0.008041955644858,  0.080567027508199,
                                   0.031763518862625,  0.008139179236825, -0.007360107893318,
                                   0.008139179236825,  0.042245795152476, -0.013493531137749,
                                  -0.007360107893318, -0.013493531137749,  0.015537355000097};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(5));

    std::vector<double> r0312vals{-0.014522101231219,  0.053494871659473,  0.019452680603445,
                                   0.082466044299379,  0.067897378051834,  0.043882287079406,
                                   0.029987652472501,  0.043882287079406, -0.036821716114931,
                                  -0.030937906977868,  0.129754078748514,  0.043882287079406,
                                   0.113951620944929,  0.128792581148996,  0.075967747296619,
                                   0.043882287079406,  0.086502719165676, -0.067506479544040,
                                  -0.011250147991952,  0.043882287079406,  0.025034984581749,
                                   0.043882287079406,  0.045264918701223,  0.045897471733207,
                                   0.009232526778164,  0.016926299093301, -0.022979063572148};

    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputePPPPForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(agtos, 1.0e-13);
    
    CMemBlock<double> fints(81);
    
    eridrv.compute(fints.data(), bpairs.pick(0), bpairs.pick(0));

    std::vector<double> r0000vals{ 1.324120036848930,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  1.165219477049000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  1.165219477049000,
                                   0.000000000000000,  0.079450279899964,  0.000000000000000,
                                   0.079450279899964,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.079450279899964,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.079450279899964,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.079450279899964,  0.000000000000000,
                                   0.079450279899964,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   1.165219477049000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  1.324120036848930,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  1.165219477049000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.079450279899964,
                                   0.000000000000000,  0.079450279899964,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.079450279899964,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.079450279899964,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.079450279899964,
                                   0.000000000000000,  0.079450279899964,  0.000000000000000,
                                   1.165219477049000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  1.165219477049000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  1.324120036848930};

    vlxtest::compare(r0000vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), bpairs.pick(1));

    std::vector<double> r0101vals{ 0.985817209017402,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.866952226906758,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.866952226906758,
                                   0.000000000000000,  0.059432491055322,  0.000000000000000,
                                   0.059432491055322,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.059432491055322,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.059432491055322,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.059432491055322,  0.000000000000000,
                                   0.059432491055322,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.866952226906758,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.985817209017402,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.866952226906758,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.059432491055322,
                                   0.000000000000000,  0.059432491055322,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.059432491055322,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.059432491055322,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.059432491055322,
                                   0.000000000000000,  0.059432491055322,  0.000000000000000,
                                   0.866952226906758,  0.000000000000000,  0.000000000000000,
                                   0.000000000000000,  0.866952226906758,  0.000000000000000,
                                   0.000000000000000,  0.000000000000000,  0.985817209017402};

    vlxtest::compare(r0101vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(1), bpairs.pick(8));

    std::vector<double> r0123vals{ 0.434400987263440,  0.012306427508151,  0.004475064548418,
                                   0.012306427508151,  0.477839565430902,  0.020216204162471,
                                   0.004475064548418,  0.020216204162471,  0.429596350952278,
                                   0.012099176064835, -0.001117051089860, -0.001395907024805,
                                  -0.001117051089860,  0.028426445525904,  0.005175785063223,
                                  -0.001395907024805,  0.005175785063223,  0.026275344338962,
                                   0.004399700387213, -0.001395907024805,  0.002214090673878,
                                  -0.001395907024805,  0.016888845185986, -0.005935570998635,
                                   0.002214090673878, -0.005935570998635,  0.003002714764875,
                                   0.012099176064835, -0.001117051089860, -0.001395907024805,
                                  -0.001117051089860,  0.028426445525904,  0.005175785063223,
                                  -0.001395907024805,  0.005175785063223,  0.026275344338962,
                                   0.477572650693298,  0.028633696969220,  0.016964209347191,
                                   0.028633696969220,  0.486859990413360,  0.019089131312813,
                                   0.016964209347191,  0.019089131312813,  0.463435809570639,
                                   0.020078036533594,  0.005175785063223, -0.005935570998635,
                                   0.005175785063223,  0.018950963683936, -0.005430205309213,
                                  -0.005935570998635, -0.005430205309213,  0.005504977068938,
                                   0.004399700387213, -0.001395907024805,  0.002214090673878,
                                  -0.001395907024805,  0.016888845185986, -0.005935570998635,
                                   0.002214090673878, -0.005935570998635,  0.003002714764875,
                                   0.020078036533594,  0.005175785063223, -0.005935570998635,
                                   0.005175785063223,  0.018950963683936, -0.005430205309213,
                                  -0.005935570998635, -0.005430205309213,  0.005504977068938,
                                   0.429659154419950,  0.026482595782278,  0.003078078926081,
                                   0.026482595782278,  0.463765527775914,  0.005643144697815,
                                   0.003078078926081,  0.005643144697815,  0.428169504502362};

    vlxtest::compare(r0123vals, fints.data());

    eridrv.compute(fints.data(), bpairs.pick(3), bpairs.pick(5));

    std::vector<double> r0312vals{ 0.016566258811693, -0.029593027022679, -0.010761100735520,
                                  -0.038062728927181, -0.032163055222041, -0.022143397895054,
                                  -0.013840992337157, -0.022143397895054,  0.020679144300248,
                                  -0.020486981860483,  0.102859587179471,  0.033939425136731,
                                   0.039799609976381,  0.077734864017329,  0.044178971788063,
                                   0.014971833782310,  0.051445039107019, -0.037952465807687,
                                  -0.007449811585630,  0.033939425136731,  0.021867777194091,
                                   0.014971833782310,  0.027245589113936,  0.028596021673293,
                                   0.004071370268595,  0.008614336546164, -0.012779262492248,
                                  -0.013642431687775,  0.049661060177614,  0.020076584474713,
                                   0.101813087958240,  0.070890313844621,  0.051740521240073,
                                   0.033397707892799,  0.041394562319480, -0.037952465807687,
                                  -0.010574122614948,  0.053314566777922,  0.021149890495462,
                                   0.061784268682424,  0.068705659285387,  0.041189512454949,
                                   0.021149890495462,  0.035543044518615, -0.028199031361167,
                                  -0.010967950427853,  0.041394562319480,  0.028596021673293,
                                   0.051740521240073,  0.047260209229747,  0.058695986924846,
                                   0.008614336546164,  0.015165591578048, -0.023428647902454,
                                  -0.004960884250100,  0.020076584474713,  0.001751029044775,
                                   0.033397707892799,  0.027245589113936,  0.000957210507559,
                                   0.022114012304969,  0.029408597539190, -0.015268189827778,
                                  -0.010967950427853,  0.051445039107019,  0.000957210507559,
                                   0.044178971788063,  0.051823242678219,  0.002983800152995,
                                   0.029408597539190,  0.059988721256955, -0.027991681350926,
                                   0.015599395451521, -0.054715637008488, -0.021659370579305,
                                  -0.054715637008488, -0.054867712816987, -0.039708846062059,
                                  -0.018579478977668, -0.034062378125724,  0.038503569391985};

    vlxtest::compare(r0312vals, fints.data());
}
