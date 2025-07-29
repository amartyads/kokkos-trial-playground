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

    int num_benchmarks = 1000;
    long onesweepTimer, twosweepTimer, pullbackTimer, onesweep_noOrderTimer = 0;
    for(int i = 0; i < num_benchmarks; i++)
    {
        int numcells = dis(gen) % 50000 + 2;
        int cellSize = dis(gen) % 30 + 2;
        MoleculeContainer container(numcells,cellSize,gen,dis);
        //std::cout << "Container made" << std::endl;
        container.populateRandomly();
        //std::cout << "Container populated" << std::endl;
        container.makeRandomHoles();
        //std::cout << "Holes made" << std::endl;
        onesweepTimer += container.compactor_onesweep();
        //std::cout << "oneSweep after: " << onesweepTimer << std::endl;
        twosweepTimer += container.compactor_twosweep();
        //std::cout << "twoSweep after: " << twosweepTimer << std::endl;
        pullbackTimer += container.compactor_pullback();
        //std::cout << "pullback after: " << pullbackTimer << std::endl;
        onesweep_noOrderTimer += container.compactor_onesweep_noOrder();
        //std::cout << "oneSweep_noOrder after: " << onesweep_noOrderTimer << std::endl;
        //std::cout << "----------------------------------------------" << std::endl;
    }
    std::cout << "oneSweep: " << onesweepTimer << std::endl;
    std::cout << "twoSweep: " << twosweepTimer << std::endl;
    std::cout << "pullback: " << pullbackTimer << std::endl;
    std::cout << "oneSweep_noOrder: " << onesweep_noOrderTimer << std::endl;
    return 0;
}