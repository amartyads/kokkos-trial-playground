#pragma once

#include <iostream>
#include <sstream>
#include <iterator>
#include <cstddef>

#include <Kokkos_Core.hpp>

#include <molecule.hpp>

class LinkedCell
{
public:
    LinkedCell(Kokkos::View<int, Kokkos::LayoutRight, Kokkos::SharedSpace> nMolecules, Kokkos::View<Molecule*, Kokkos::LayoutRight, Kokkos::SharedSpace> moleculeSlice) 
        : numMolecules(nMolecules), moleculeData(moleculeSlice) {}
    
    class Iterator
    {
        using iterator_category = std::bidirectional_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = Molecule;
        using pointer = value_type*;
        using reference = value_type&;
public:
        Iterator(pointer ptr) : _myPtr(ptr), _idx(0) {}

        reference operator*() const { return *_myPtr;} 
        pointer operator->() { return _myPtr; }
        Iterator& operator++() { _myPtr++; _idx++; return *this; }
        Iterator operator++(int) { Iterator temp = *this; ++(*this); _idx++; return temp; }
        Iterator& operator--() { _myPtr--; _idx--; return *this; }
        Iterator operator--(int) { Iterator temp = *this; --(*this); _idx--; return temp; }
        
        friend bool operator== (const Iterator& a, const Iterator& b) { return a._myPtr == b._myPtr; }
        friend bool operator!= (const Iterator& a, const Iterator& b) { return a._myPtr != b._myPtr; }

        unsigned int getIndex() const { return _idx; }

private:
        pointer _myPtr;
        unsigned int _idx;
    };

    Iterator begin() { return Iterator(&moleculeData(0)); }
    Iterator end() { return Iterator(&moleculeData(numMolecules())); }

    void insert(Molecule& molecule)
    {
        moleculeData(numMolecules()) = molecule;
        numMolecules() += 1;
    }
    void remove(int moleculeIdx)
    {
        moleculeData(moleculeIdx) = moleculeData(numMolecules() - 1);
        numMolecules() -= 1;
    }
    void clear()
    {
        numMolecules() = 0;
    }

    std::string to_string() const
    {
        std::stringstream to_ret;
        to_ret << " numMol: " << numMolecules() << std::endl;
        return to_ret.str();
    }

    Kokkos::View<int, Kokkos::LayoutRight, Kokkos::SharedSpace> numMolecules;
    Kokkos::View<Molecule*, Kokkos::LayoutRight, Kokkos::SharedSpace> moleculeData;
};
