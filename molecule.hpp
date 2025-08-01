#pragma once

#include <iostream>
#include <vector>
#include <sstream>

class Molecule
{
public:
    Molecule(int id, double px, double py, double pz, double vx, double vy, double vz,double fx, double fy, double fz) : dirty(false) 
    {
        this->id = id;
        pos[0] = px; pos[1] = py; pos[2] = pz;
        vel[0] = vx; vel[1] = vy; vel[2] = vz;
        f[0] = fx; f[1] = fy; f[2] = fz;
    }
    Molecule(int id, double px, double py, double pz) : dirty(false)
    {
        this->id = id;
        pos[0] = px; pos[1] = py; pos[2] = pz;
        vel[0] = 0; vel[1] = 0; vel[2] = 0;
        f[0] = 0; f[1] = 0; f[2] = 0;
    }
    Molecule() : id(-1), dirty(true) {}
    std::string to_string() const
    {
        std::stringstream to_ret;
        to_ret << "ID: " << id << std::endl;
        to_ret << "pos: " << pos[0] << ", " << pos[1] << ", " << pos[2];
        to_ret << " vel: " << vel[0] << ", " << vel[1] << ", " << vel[2];
        to_ret << " f: " << f[0] << ", " << f[1] << ", " << f[2];
        return to_ret.str();
    }
    bool dirty;
    int id;
    double pos[3], vel[3], f[3];
};