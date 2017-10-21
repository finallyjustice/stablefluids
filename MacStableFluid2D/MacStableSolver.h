/** File:    MacStableSolver.h
 ** Author:  Dongli Zhang
 ** Contact: dongli.zhang0129@gmail.com
 **
 ** Copyright (C) Dongli Zhang 2013
 **
 ** This program is free software;  you can redistribute it and/or modify
 ** it under the terms of the GNU General Public License as published by
 ** the Free Software Foundation; either version 2 of the License, or
 ** (at your option) any later version.
 **
 ** This program is distributed in the hope that it will be useful,
 ** but WITHOUT ANY WARRANTY;  without even the implied warranty of
 ** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
 ** the GNU General Public License for more details.
 **
 ** You should have received a copy of the GNU General Public License
 ** along with this program;  if not, write to the Free Software 
 ** Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#ifndef __MACSTABLESOLVER_H__
#define __MACSTABLESOLVER_H__

#include "Vector2f.h"
#include <stdio.h>

class StableSolver
{
public:
    StableSolver();
    ~StableSolver();
    void init();
    void reset();
    void cleanBuffer();
    void start(){ running=1; }
    void stop(){ running=0; }
    int isRunning(){ return running; }

    //animation
    void setVelBoundary(int flag);
    void setCellBoundary(float *value);
    void projection();
    void advectVel();
    void advectCell(float *value, float *value0);
    void diffuseVel();
    void diffuseCell(float *value, float *value0);
    void addSource();
    void animVel();
    void animDen();

    //getter
    int getRowCell(){ return rowCell; }
    int getColCell(){ return colCell; }
    int getTotCell(){ return totCell; }
    int getRowVelX(){ return rowVelX; }
    int getcolVelX(){ return colVelX; }
    int getTotVelX(){ return totVelX; }
    int getRowVelY(){ return rowVelY; }
    int getColVelY(){ return colVelY; }
    int getTotVelY(){ return totVelY; }
    float* getVX(){ return vx; }
    float* getVY(){ return vy; }
    float* getD(){ return d;}
    Vec2f* getPVX(){ return pvx; }
    Vec2f* getPVY(){ return pvy; }
    int vxIdx(int i, int j){ return j*rowVelX+i; }
    int vyIdx(int i, int j){ return j*rowVelY+i; }
    int cIdx(int i, int j){ return j*rowCell+i; }
    Vec2f getCellVel(int i, int j){ return Vec2f((vx[vxIdx(i, j)]+vx[vxIdx(i+1, j)])/2, (vy[vyIdx(i, j)]+vy[vyIdx(i, j+1)])/2); }
    float getDens(int i, int j){ return (d[cIdx(i-1, j-1)]+d[cIdx(i, j-1)]+d[cIdx(i-1, j)]+d[cIdx(i, j)])/4.0f; }

    //setter
    void setVel0(int i, int j, float _vx0, float _vy0)
    { 
        vx0[vxIdx(i, j)] += _vx0;
        vx0[vxIdx(i+1, j)] += _vx0;
        vy0[vyIdx(i, j)] += _vy0;
        vy0[vyIdx(i, j+1)] += _vy0;
    }
    void setD0(int i, int j, float _d0){ d0[cIdx(i, j)]=_d0; }

private:
    int rowCell;
    int colCell;
    int totCell;
    int rowVelX;
    int colVelX;
    int totVelX;
    int rowVelY;
    int colVelY;
    int totVelY;
    float minX;
    float maxX;
    float minY;
    float maxY;

    //params
    int running;
    float timeStep;
    float diff;
    float visc;

    float *vx;
    float *vy;
    float *vx0;
    float *vy0;
    float *d;
    float *d0;
    float *div;
    float *p;
    Vec2f *pvx;
    Vec2f *pvy;
};

#endif
