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
    KOKKOS_FUNCTION LinkedCell() : linkedCellNumMolecules(NULL), moleculeData(NULL), linkedCellIndex(0) {}

    KOKKOS_FUNCTION LinkedCell(const Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::SharedSpace>* nMolecules, const Kokkos::View<Molecule**, Kokkos::LayoutRight, Kokkos::SharedSpace>* moleculeSlice, unsigned int cellIndex) 
        : linkedCellNumMolecules(nMolecules), moleculeData(moleculeSlice), linkedCellIndex(cellIndex) {}
    
    class Iterator
    {
        using iterator_category = std::bidirectional_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = Molecule;
        using pointer = value_type*;
        using reference = value_type&;
public:
        KOKKOS_FUNCTION Iterator(pointer ptr) : _myPtr(ptr), _idx(0) {}

        KOKKOS_INLINE_FUNCTION reference operator*() const { return *_myPtr;} 
        KOKKOS_INLINE_FUNCTION pointer operator->() { return _myPtr; }
        KOKKOS_INLINE_FUNCTION Iterator& operator++() { _myPtr++; _idx++; return *this; }
        KOKKOS_INLINE_FUNCTION Iterator operator++(int) { Iterator temp = *this; ++(*this); _idx++; return temp; }
        KOKKOS_INLINE_FUNCTION Iterator& operator--() { _myPtr--; _idx--; return *this; }
        KOKKOS_INLINE_FUNCTION Iterator operator--(int) { Iterator temp = *this; --(*this); _idx--; return temp; }
        
        KOKKOS_INLINE_FUNCTION friend bool operator== (const Iterator& a, const Iterator& b) { return a._myPtr == b._myPtr; }
        KOKKOS_INLINE_FUNCTION friend bool operator!= (const Iterator& a, const Iterator& b) { return a._myPtr != b._myPtr; }

        KOKKOS_INLINE_FUNCTION unsigned int getIndex() const { return _idx; }

private:
        pointer _myPtr;
        unsigned int _idx;
    };

    KOKKOS_FUNCTION Iterator begin() { return Iterator(&(*moleculeData)(linkedCellIndex, 0)); }
    KOKKOS_FUNCTION Iterator end() { return Iterator(&(*moleculeData)(linkedCellIndex, numMolecules())); }

    KOKKOS_FUNCTION unsigned int numMolecules() const {
        return (*linkedCellNumMolecules)(linkedCellIndex);
    }

    KOKKOS_FUNCTION void changeMoleculeCount(int by) {
        (*linkedCellNumMolecules)(linkedCellIndex) += by;
    }

    KOKKOS_FUNCTION void insert(Molecule& molecule)
    {
        (*moleculeData)(linkedCellIndex, numMolecules()) = molecule;
        changeMoleculeCount(+1);
    }
    KOKKOS_FUNCTION void remove(int moleculeIdx)
    {
        (*moleculeData)(linkedCellIndex, moleculeIdx) = (*moleculeData)(linkedCellIndex, numMolecules() - 1);
        changeMoleculeCount(-1);
    }
    KOKKOS_FUNCTION void clear()
    {
        (*linkedCellNumMolecules)(linkedCellIndex) = 0;
    }

    std::string to_string() const
    {
        std::stringstream to_ret;
        to_ret << " numMol: " << numMolecules() << std::endl;
        return to_ret.str();
    }

    const Kokkos::View<Molecule**, Kokkos::LayoutRight, Kokkos::SharedSpace>* moleculeData;
    const Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::SharedSpace>* linkedCellNumMolecules;
    unsigned int linkedCellIndex;
};
