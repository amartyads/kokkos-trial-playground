#pragma once

#include <vector>

class IndexConverter
{
public:
    IndexConverter(int domainSize, int numCellsPerDim) : domainSize(domainSize), numCellsPerDim(numCellsPerDim) 
    {
        for(double i = domainSize/numCellsPerDim; i <= domainSize; i+= domainSize/numCellsPerDim)
        {
            xspans.push_back(i);
            yspans.push_back(i);
            zspans.push_back(i);
        }
    }
    IndexConverter() : domainSize(0), numCellsPerDim(0) {}
    void reset(int domainSize, int numCellsPerDim)
    {
        this->domainSize = domainSize;
        this->numCellsPerDim = numCellsPerDim;
        xspans.clear();
        yspans.clear();
        zspans.clear();
        for(double i = domainSize/numCellsPerDim; i <= domainSize; i+= domainSize/numCellsPerDim)
        {
            xspans.push_back(i);
            yspans.push_back(i);
            zspans.push_back(i);
        }
    }
    int getIndex(double posx, double posy, double posz)
    {
        int toRet = 0;
        int xid = -1, yid = -1, zid = -1;
        for(int i = 0; i < xspans.size(); i++)
        {
            if(xid == -1 && posx < xspans[i]) xid = i;
            if(yid == -1 && posy < yspans[i]) yid = i;
            if(zid == -1 && posz < zspans[i]) zid = i;
            if(xid != -1 && yid != -1 && zid != -1) break;
        }
        return xid + yid*numCellsPerDim + zid*numCellsPerDim*numCellsPerDim;
    }
    int getIndex(double pos[3])
    {
        return getIndex(pos[0],pos[1], pos[2]);
    }
    bool isInIndex(double posx, double posy, double posz, int index)
    {
        return getIndex(posx, posy, posz) == index;
    }
    bool isInIndex(double pos[3], int index)
    {
        return isInIndex(pos[0], pos[1], pos[2], index);
    }
    int domainSize, numCellsPerDim;
    std::vector<double> xspans, yspans, zspans;
};