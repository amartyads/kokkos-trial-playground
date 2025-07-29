#include <random>
#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <numeric>

#include <Kokkos_Core.hpp>

class MoleculeContainer
{
public:
    MoleculeContainer(int numCells, int cellSize, std::mt19937 gen, std::uniform_int_distribution<> dis) : _numCells(numCells), _cellSize(cellSize), _gen(gen), 
        _dis(dis), moleculeData("moleculeData", numCells, cellSize) {}

    void grow(int cellSize)
    {
        assert(_cellSize <= cellSize);
        _cellSize = cellSize;
        Kokkos::resize(moleculeData,_numCells, _cellSize);
    }

    void printData() const
    {
        std::cout << "Container contents: " << std::endl;
        for (size_t i = 0; i < _numCells; i++)
        {
            std::cout << "Cell #" << i << ": " ;
            for (size_t j = 0; j < _cellSize; j++)
            {
                std::cout << moleculeData(i,j) << " ";
            }
            std::cout << std::endl;
        }
    }

    void populateRandomly()
    {
        //not in parallel to make sure the same rows have same data
        for (size_t i = 0; i < _numCells; i++)
        {
            for (size_t j = 0; j < _cellSize; j++)
            {
                moleculeData(i,j) = _dis(_gen) % 9 + 1;
            }
        }
    }

    void makeRandomHoles(bool display = false)
    {
        int numHoles = _dis(_gen) % (_numCells*_cellSize);
        std::vector<int> allCoords(_numCells*_cellSize);
        std::iota(allCoords.begin(), allCoords.end(), 0);
        std::shuffle(allCoords.begin(), allCoords.end(),_gen);
        for (size_t i = 0; i < numHoles; i++)
        {
            moleculeData(allCoords[i]/_cellSize, allCoords[i]%_cellSize) = 0;
        }
        if(display) std::cout << numHoles << " holes created!" << std::endl;
    }

    long compactor_onesweep(bool display = false)
    {
        Kokkos::View<int**> containerCopy("copy", _numCells, _cellSize);
        Kokkos::deep_copy(containerCopy, moleculeData);
        auto t1 = std::chrono::high_resolution_clock::now();
        Kokkos::parallel_for(_numCells, KOKKOS_LAMBDA(const unsigned int i)
        {
            int j = 0;
            int lastHole = -1;
            //find first hole
            for (; j < _cellSize; j++)
            {
                if(containerCopy(i,j) == 0)
                    break;
            }
            lastHole = j;
            //iterate first hole onwards
            for (int j = lastHole; j < _cellSize; j++)
            {
                if(containerCopy(i,j) != 0)
                {
                    containerCopy(i, lastHole) = containerCopy(i,j);
                    containerCopy(i,j) = 0;
                    while (containerCopy(i,lastHole) != 0 && lastHole < j) lastHole++;
                }
            }
        }); //kokkos parallel for
        auto t2 = std::chrono::high_resolution_clock::now();
        if(display)
        {
            std::cout << "Data compacted with onesweep! Data:" << std::endl;
            for (size_t i = 0; i < _numCells; i++)
            {
                std::cout << "Cell #" << i << ": " ;
                for (size_t j = 0; j < _cellSize; j++)
                {
                    std::cout << containerCopy(i,j) << " ";
                }
                std::cout << std::endl;
            }
        }
        return std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
    }

    long compactor_twosweep(bool display = false)
    {
        Kokkos::View<int**> containerCopy("copy", _numCells, _cellSize);
        Kokkos::deep_copy(containerCopy, moleculeData);
        auto t1 = std::chrono::high_resolution_clock::now();
        Kokkos::parallel_for(_numCells, KOKKOS_LAMBDA(const unsigned int i)
        {
            int shiftAmt[_cellSize];
            int curShift = 0;
            for (int j = 0; j < _cellSize; j++)
            {
                shiftAmt[j] = curShift;
                if(containerCopy(i,j) == 0)
                {
                    curShift++;
                }
            }
            for (int j = 0; j < _cellSize; j++)
            {
                if(containerCopy(i,j) != 0 && shiftAmt[j] != 0)
                {
                    containerCopy(i, j-shiftAmt[j]) = containerCopy(i,j);
                    containerCopy(i,j) = 0;
                }
            }
        }); //kokkos parallel for
        auto t2 = std::chrono::high_resolution_clock::now();
        if(display)
        {
            std::cout << "Data compacted with twosweep! Data:" << std::endl;
            for (size_t i = 0; i < _numCells; i++)
            {
                std::cout << "Cell #" << i << ": " ;
                for (size_t j = 0; j < _cellSize; j++)
                {
                    std::cout << containerCopy(i,j) << " ";
                }
                std::cout << std::endl;
            }
        }
        return std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
    }

    long compactor_pullback(bool display = false)
    {
        Kokkos::View<int**> containerCopy("copy", _numCells, _cellSize);
        Kokkos::deep_copy(containerCopy, moleculeData);
        auto t1 = std::chrono::high_resolution_clock::now();
        Kokkos::parallel_for(_numCells, KOKKOS_LAMBDA(const unsigned int i)
        {
            int sourceIdx[_cellSize];
            int sourceIdxIdx = 0;
            for (int j = 0; j < _cellSize; j++)
            {
                if(containerCopy(i,j) != 0)
                {
                    sourceIdx[sourceIdxIdx++] = j;
                }
            }
            for (int j = 0; j < _cellSize; j++)
            {
                if(j < sourceIdxIdx)
                    containerCopy(i, j) = containerCopy(i,sourceIdx[j]);
                else
                    containerCopy(i, j) = 0;
            }
        }); //kokkos parallel for
        auto t2 = std::chrono::high_resolution_clock::now();
        if(display)
        {
            std::cout << "Data compacted with pullback! Data:" << std::endl;
            for (size_t i = 0; i < _numCells; i++)
            {
                std::cout << "Cell #" << i << ": " ;
                for (size_t j = 0; j < _cellSize; j++)
                {
                    std::cout << containerCopy(i,j) << " ";
                }
                std::cout << std::endl;
            }
        }
        return std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
    }

    long compactor_onesweep_noOrder(bool display = false)
    {
        Kokkos::View<int**> containerCopy("copy", _numCells, _cellSize);
        Kokkos::deep_copy(containerCopy, moleculeData);
        auto t1 = std::chrono::high_resolution_clock::now();
        //Kokkos::RangePolicy<Kokkos::Schedule<Kokkos::Static>> rp(0, _numCells, Kokkos::ChunkSize(_cellSize));
        Kokkos::parallel_for(_numCells, KOKKOS_LAMBDA(const unsigned int i)
        {
            int j = 0, k = _cellSize - 1;
            while(j < k)
            {
                //find first hole on left side
                while(j < _cellSize && containerCopy(i,j) != 0) j++;
                //find first data on right side
                while(k > -1 && containerCopy(i,k) == 0) k--;
                if(k <= j) break;
                containerCopy(i,j) = containerCopy(i,k);
                containerCopy(i,k) = 0;
            }
        }); //kokkos parallel for
        auto t2 = std::chrono::high_resolution_clock::now();
        if(display)
        {
            std::cout << "Data compacted with onesweep_noOrder! Data:" << std::endl;
            for (size_t i = 0; i < _numCells; i++)
            {
                std::cout << "Cell #" << i << ": " ;
                for (size_t j = 0; j < _cellSize; j++)
                {
                    std::cout << containerCopy(i,j) << " ";
                }
                std::cout << std::endl;
            }
        }
        return std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
    }


private:
    int _numCells;
    int _cellSize;
    std::mt19937 _gen;
    std::uniform_int_distribution<> _dis;
    Kokkos::View<int**> moleculeData;
};