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
    MoleculeContainer(int numCells, int cellSize, std::mt19937 gen, std::uniform_int_distribution<> dis) : _numCells(numCells), _cellSize(cellSize), _gen(gen), 
        _dis(dis), moleculeData("moleculeData", numCells, cellSize), linkedCellNumMolecules("linkedCellNumMolecules", numCells)
        {}

    void grow(int cellSize)
    {
        assert(_cellSize <= cellSize);
        _cellSize = cellSize;
        Kokkos::resize(moleculeData, _numCells, _cellSize);
        // new space created is filled with garbage data, so size of _linkedCell does not change
    }

    KOKKOS_FUNCTION void insert(int cellIdx, Molecule& molecule)
    {
        moleculeData(cellIdx, linkedCellNumMolecules(cellIdx)) = molecule;
        linkedCellNumMolecules(cellIdx) += 1;
    }

    KOKKOS_FUNCTION void remove(int cellIdx, int moleculeIdx)
    {
        moleculeData(cellIdx, moleculeIdx) = moleculeData(cellIdx, linkedCellNumMolecules(cellIdx) - 1);
        linkedCellNumMolecules(cellIdx) -= 1;
    }

    KOKKOS_FUNCTION void clearLinkedCell(int cellIdx)
    {
        linkedCellNumMolecules(cellIdx) = 0;
    }

    void sort(int cellIdx, IndexConverter& indexConverter) //set all outgoing molecules
    {
        for (size_t i = 0; i < linkedCellNumMolecules(cellIdx); i++)
        {
            int curMolIdx = indexConverter.getIndex(moleculeData(cellIdx, i).pos);
            if(curMolIdx != cellIdx) // if molecule does not belong to current cell anymore
            {
                // write data to target end
                moleculeData(curMolIdx, linkedCellNumMolecules(curMolIdx)) = moleculeData(cellIdx, i);
                // increment target end
                linkedCellNumMolecules(curMolIdx)++;
                // delete molecule at own position
                remove(cellIdx, i);
                // decrement iterator as the molecule at position i is now new
                i--;
            }
        }
    }

    void sort(IndexConverter& indexConverter)
    {
        //find red-black cells
        for (size_t z = 0; z < 2; z++)
        {
            for (size_t y = 0; y < 2; y++)
            {
                for (size_t x = 0; x < 2; x++)
                {
                    int lengthVector[] = {};
                }
                
            }
            
        }
        
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

    LinkedCell operator[](int idx) const
    {
        Kokkos::View<Molecule*, Kokkos::LayoutRight, Kokkos::SharedSpace> lcMoleculeSlice(moleculeData, idx, Kokkos::ALL);
        Kokkos::View<int, Kokkos::LayoutRight, Kokkos::SharedSpace> lcSizeSlice(linkedCellNumMolecules, idx);
        return LinkedCell(lcSizeSlice, lcMoleculeSlice);
    }

    int getNumCells() const { return _numCells; }

    Kokkos::View<Molecule**, Kokkos::LayoutRight, Kokkos::SharedSpace> moleculeData;
    Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::SharedSpace> linkedCellNumMolecules;

private:
    int _numCells;
    int _cellSize;
    std::mt19937 _gen;
    std::uniform_int_distribution<> _dis;
};