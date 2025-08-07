#pragma once

#include <random>
#include <iostream>
#include <chrono>
#include <vector>
#include <sstream>
#include <algorithm>
#include <numeric>

#include <Kokkos_Core.hpp>

#include <linked_cells.hpp>
#include <molecule.hpp>
#include <index_converter.hpp>

class MoleculeContainer
{
public:
    MoleculeContainer(int numCellsPerDim, int cellSize, std::mt19937 gen, std::uniform_int_distribution<> dis) : _numCellsPerDim(numCellsPerDim), _numCells(numCellsPerDim*numCellsPerDim*numCellsPerDim), _cellSize(cellSize), _gen(gen), 
        _dis(dis), moleculeData("moleculeData", numCellsPerDim*numCellsPerDim*numCellsPerDim, cellSize), linkedCellNumMolecules("linkedCellNumMolecules", numCellsPerDim*numCellsPerDim*numCellsPerDim), linkedCells("linkedCells", numCellsPerDim*numCellsPerDim*numCellsPerDim)
        {}

    void grow(int cellSize)
    {
        assert(_cellSize <= cellSize);
        _cellSize = cellSize;
        Kokkos::resize(moleculeData, _numCells, _cellSize);
        // new space created is filled with garbage data, so size of _linkedCell does not change
    }

    KOKKOS_INLINE_FUNCTION void insert(int cellIdx, Molecule& molecule)
    {
        moleculeData(cellIdx, linkedCellNumMolecules(cellIdx)) = molecule;
        linkedCellNumMolecules(cellIdx) += 1;
    }

    KOKKOS_INLINE_FUNCTION void remove(int cellIdx, int moleculeIdx)
    {
        moleculeData(cellIdx, moleculeIdx) = moleculeData(cellIdx, linkedCellNumMolecules(cellIdx) - 1);
        linkedCellNumMolecules(cellIdx) -= 1;
    }

    KOKKOS_FUNCTION void clearLinkedCell(int cellIdx)
    {
        linkedCellNumMolecules(cellIdx) = 0;
    }

    void sort(const IndexConverter& indexConverter)
    {
        Kokkos::fence();
        //find red-black cells
        for (int z = 0; z < 2; z++)
        {
            for (int y = 0; y < 2; y++)
            {
                for (int x = 0; x < 2; x++)
                {
                    const int lengthVector[3] = {(_numCellsPerDim + (_numCellsPerDim % 2) * (x == 0)) / 2, (_numCellsPerDim + (_numCellsPerDim % 2) * (y == 0)) / 2, (_numCellsPerDim + (_numCellsPerDim % 2) * (z == 0)) / 2};
                    const int length = lengthVector[0] * lengthVector[1] * lengthVector[2];
                    const int numCellsPerDim = _numCellsPerDim;
                    auto linkedCellLocal(linkedCellNumMolecules);
                    auto moleculeDataLocal(moleculeData);
                    Kokkos::parallel_for(length, KOKKOS_LAMBDA(const unsigned int j) {
                        // compute index of the current cell
                        int index = 0;
                        int helpIndex1 = j;
                        int helpIndex2 = 0;
                        // determine plane within traversed block
                        helpIndex2 = helpIndex1 / (lengthVector[0] * lengthVector[1]);
                        // save rest of index in helpIndex1
                        helpIndex1 = helpIndex1 - helpIndex2 * (lengthVector[0] * lengthVector[1]);
                        // compute contribution to index
                        index += (0+ 2 * helpIndex2 + z) * numCellsPerDim * numCellsPerDim;
                        // determine plane within traversed block
                        helpIndex2 = helpIndex1 / lengthVector[0];
                        // save rest of index in helpIndex1
                        helpIndex1 = helpIndex1 - helpIndex2 * lengthVector[0];
                        // compute contribution to index
                        index += (0 + 2 * helpIndex2 + y) * numCellsPerDim;
                        // compute contribution for last dimension
                        index += (0 + 2 * helpIndex1 + x);

                        for (size_t i = 0; i < linkedCellLocal(index); i++)
                        {
                            int curMolIdx = indexConverter.getIndex(moleculeDataLocal(index, i).pos);
                            if(curMolIdx != index) // if molecule does not belong to current cell anymore
                            {
                                // write data to target end
                                moleculeDataLocal(curMolIdx, linkedCellLocal(curMolIdx)) = moleculeDataLocal(index, i);
                                // increment target end
                                linkedCellLocal(curMolIdx)++;
                                // delete molecule at own position
                                moleculeDataLocal(index, i) = moleculeDataLocal(index, linkedCellLocal(index) - 1);
                                linkedCellLocal(index) -= 1;
                                // decrement iterator as the molecule at position i is now new
                                i--;
                            }
                        }
                    });
                }
            }
            
        }
        Kokkos::fence();
    }

    void printData() const
    {
        std::cout << "Container contents: " << std::endl;
        for (size_t i = 0; i < _numCells; i++)
        {
            std::cout << std::endl;
            std::cout << "Cell #" << i << ": " ;
            for (size_t j = 0; j < linkedCellNumMolecules(i); j++)
            {
                std::cout << std::endl;
                std::cout << moleculeData(i,j).to_string() << " ";
            }
            std::cout << std::endl;
        }
    }

    void populateRandomly(int domainSize)
    {
        //not in parallel to make sure the same rows have same data for some seed
        for (size_t i = 0; i < _numCells; i++)
        {
            for (size_t j = 0; j < _cellSize; j++)
            {
                Molecule m((i*_cellSize + j), _dis(_gen) % domainSize, _dis(_gen) % domainSize, _dis(_gen) % domainSize, _dis(_gen) % domainSize, _dis(_gen) % domainSize,_dis(_gen) % domainSize, _dis(_gen) % domainSize, _dis(_gen) % domainSize, _dis(_gen) % domainSize);
                moleculeData(i,j) = m;
            }
            linkedCellNumMolecules(i)= _cellSize;
        }
    }

    void makeRandomHoles()
    {
        int numHoles = _dis(_gen) % (_numCells*_cellSize);
        std::vector<int> allCoords(_numCells*_cellSize);
        std::iota(allCoords.begin(), allCoords.end(), 0);
        std::shuffle(allCoords.begin(), allCoords.end(),_gen);
        for (size_t i = 0; i < numHoles; i++)
        {
            Molecule m;
            moleculeData(allCoords[i]/_cellSize, allCoords[i]%_cellSize) = m;
        }
    }

    KOKKOS_INLINE_FUNCTION Molecule& getMoleculeAt(int i, int j) const { return moleculeData(i,j); }

    KOKKOS_FUNCTION LinkedCell& operator[](unsigned int idx) const
    {
        // Kokkos::View<Molecule*, Kokkos::LayoutRight, Kokkos::SharedSpace> lcMoleculeSlice(moleculeData, idx, Kokkos::ALL);
        // Kokkos::View<int, Kokkos::LayoutRight, Kokkos::SharedSpace> lcSizeSlice(linkedCellNumMolecules, idx);
        // return LinkedCell(lcSizeSlice, lcMoleculeSlice);

        // LinkedCell* cell = new (&linkedCells(idx)) LinkedCell(&linkedCellNumMolecules, &moleculeData, idx);
        // return *cell;

        // return LinkedCell(&linkedCellNumMolecules, &moleculeData, idx);

        linkedCells(idx) = LinkedCell(&linkedCellNumMolecules, &moleculeData, idx);
        return linkedCells(idx);
    }

    KOKKOS_FUNCTION int getNumCells() const { return _numCells; }

    
    void testTestData() {
        MoleculeContainer container = (*this);
        auto linkedCells2 = linkedCells;
        Kokkos::parallel_for(linkedCellNumMolecules.size(), KOKKOS_LAMBDA(const unsigned int i)
        {
            LinkedCell curCell = container[i];
            // linkedCells2(i) = LinkedCell(&linkedCellNumMolecules, &moleculeData, i);
            // auto curCell = linkedCells2(i);

            for (auto it = curCell.begin(); it != curCell.end(); ++it)
                (*it).f[0] = 69420;
        });
    }

    Kokkos::View<Molecule**, Kokkos::LayoutRight, Kokkos::SharedSpace> moleculeData;
    Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::SharedSpace> linkedCellNumMolecules;
    Kokkos::View<LinkedCell*, Kokkos::LayoutRight, Kokkos::SharedSpace> linkedCells;


private:
    int _numCells;
    int _numCellsPerDim;
    int _cellSize;
    std::mt19937 _gen;
    std::uniform_int_distribution<> _dis;
};