#pragma once

#include <iostream>

class IndexConverter
{
public:
    IndexConverter(int domainSize, int numCellsPerDim) : domainSize(domainSize), numCellsPerDim(numCellsPerDim) {}
    IndexConverter() : domainSize(0), numCellsPerDim(0) {}
    void reset(int domainSize, int numCellsPerDim)
    {
        this->domainSize = domainSize;
        this->numCellsPerDim = numCellsPerDim;
    }
    int getIndex(double posx, double posy, double posz) const
    {
        int xid = -1, yid = -1, zid = -1;
        posx -= (domainSize/numCellsPerDim);
        posy -= (domainSize/numCellsPerDim);
        posz -= (domainSize/numCellsPerDim);
        int i = 0;
        while(xid == -1 || yid == -1 || zid == -1)
        {
            if(xid == -1 && posx < 0) xid = i; else posx -= (domainSize/numCellsPerDim);
            if(yid == -1 && posy < 0) yid = i; else posy -= (domainSize/numCellsPerDim);
            if(zid == -1 && posz < 0) zid = i; else posz -= (domainSize/numCellsPerDim);
            i++;
        }
        return xid + yid*numCellsPerDim + zid*numCellsPerDim*numCellsPerDim;
    }
    KOKKOS_FUNCTION int getIndex(double pos[3]) const
    {
        return getIndex(pos[0],pos[1], pos[2]);
    }
    KOKKOS_FUNCTION bool isInIndex(double posx, double posy, double posz, int index) const
    {
        return getIndex(posx, posy, posz) == index;
    }
    KOKKOS_FUNCTION bool isInIndex(double pos[3], int index) const
    {
        return isInIndex(pos[0], pos[1], pos[2], index);
    }
    int domainSize, numCellsPerDim;
};