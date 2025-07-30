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
    LinkedCell(Kokkos::View<int> nMolecules, Kokkos::View<Molecule*> moleculeSlice) 
        : numMolecules(nMolecules), molecules(moleculeSlice) {}
    
    struct Iterator
    {
        using iterator_category = std::forward_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = Molecule;
        using pointer = value_type*;
        using reference = value_type&;
public:
        Iterator(pointer ptr) : _myPtr(ptr) {}

        reference operator*() const { return *_myPtr;} 
        pointer operator->() { return _myPtr; }
        Iterator& operator++() { _myPtr++; return *this; }
        Iterator operator++(int) { Iterator temp = *this; ++(*this); return temp; }
        
        friend bool operator== (const Iterator& a, const Iterator& b) { return a._myPtr == b._myPtr; }
        friend bool operator!= (const Iterator& a, const Iterator& b) { return a._myPtr != b._myPtr; }

private:
        pointer _myPtr;
    };

    Iterator begin() { return Iterator(&molecules(0)); }
    Iterator end() { return Iterator(&molecules(numMolecules())); }
    std::string to_string() const
    {
        std::stringstream to_ret;
        to_ret << " numMol: " << numMolecules() << std::endl;
        return to_ret.str();
    }
    Kokkos::View<int> numMolecules;
    Kokkos::View<Molecule*> molecules;
};
