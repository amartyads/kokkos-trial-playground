#include <iostream>
#include <random>

#include <Kokkos_Core.hpp>

#include <molecule_container.hpp>
#include <linked_cells.hpp>
#include <molecule.hpp>
#include <index_converter.hpp>

int main(int argc, char* argv[])
{
    Kokkos::ScopeGuard guard(argc, argv);
    //std::random_device rd;
    //std::mt19937 gen(rd());

    std::mt19937 gen(1984);
    std::uniform_int_distribution<> dis(0, RAND_MAX);

    int domainSizeVolume = 2;
    double cellSizeVolume = 1;
    int numCellsPerDim = static_cast<int>(domainSizeVolume/cellSizeVolume);
    int totNumCells = numCellsPerDim * numCellsPerDim * numCellsPerDim;
    int cellSizeMolecules = 2;
    int extraCellSpaceFactor = 3;

    IndexConverter indexConverter(domainSizeVolume, numCellsPerDim);
    MoleculeContainer container(totNumCells, cellSizeMolecules, gen, dis);

    container.populateRandomly(domainSizeVolume); 
    container.printData();
    container.grow(cellSizeMolecules * extraCellSpaceFactor);
    for(int i = 0 ; i < totNumCells; i++)
        container.sort(i, indexConverter);
    
    container.printData();

    std::cout << "Iteration: " << std::endl;
    LinkedCell cell = container[0];
    for(auto x: cell)
        std::cout << x.to_string() << " ";
    std::cout << std::endl;
    //container.compact();
    //container.printData();

    return 0;
}