#pragma once

#include <iostream>
#include <sstream>

#include <Kokkos_Core.hpp>

class LinkedCell
{
public:
    LinkedCell() : numMolecules(0) {}
    std::string to_string() const
    {
        std::stringstream to_ret;
        to_ret << " numMol: " << numMolecules << std::endl;
        return to_ret.str();
    }
    int numMolecules;
};
