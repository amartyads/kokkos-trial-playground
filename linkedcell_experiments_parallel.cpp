#include <iostream>
#include <random>

#include <Kokkos_Core.hpp>

#include <moleculecontainer_parallel_benching.hpp>

int main(int argc, char* argv[])
{
    Kokkos::ScopeGuard guard(argc, argv);
    //std::random_device rd;
    //std::mt19937 gen(rd());

    std::mt19937 gen(1984);
    std::uniform_int_distribution<> dis(0, RAND_MAX);

    MoleculeContainer container(10,12,gen,dis);
    container.populateRandomly();
    container.printData();
    container.makeRandomHoles();
    container.printData();
    //container.compactor_onesweep();
    //container.compactor_twosweep();
    //container.compactor_pullback();
    container.compactor_onesweep_noOrder();
    return 0;
}