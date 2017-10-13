/** File:    StableSolver2D.h
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

#ifndef __STABLESOLVER2D_H__
#define __StABLESOLVER2D_H__

class StableSolver2D
{
public:
    StableSolver2D();
    ~StableSolver2D();

    void start(){ running = 1; }
    void stop(){ running = 0; }
    int isRunning(){ return running; }

    void reset(int _rowSize, int _colSize);
    void clear();
    void addSource();
    void anim_vel();
    void anim_den();
    void anim_tex();

    float getForce(){ return force; }
    float getSource(){ return source; }

    int getRowSize(){ return rowSize; }
    int getColSize(){ return colSize; }
    int getTotSize(){ return totSize; }
    int getIndex(int i, int j){ return j*(rowSize+2)+i; }

    float* getPX(){ return px; }
    float* getPY(){ return py; }
    float* getVX(){ return vx; }
    float* getVY(){ return vy; }
    float* getD(){ return d; }
    float* getTX(){ return tx; }
    float* getTY(){ return ty;}
    float getDens(int i, int j){ return (d[getIndex(i-1, j-1)] +  d[getIndex(i, j-1)] + d[getIndex(i-1, j)] + d[getIndex(i, j)])/4.0f; }

    void setVX0(int i, int j, float _vx0){ vx0[getIndex(i, j)] = _vx0; }
    void setVY0(int i, int j, float _vy0){ vy0[getIndex(i, j)] = _vy0; }
    void setD0(int i, int j, float _d0){ d0[getIndex(i, j)] = _d0; }

    void cleanBuffer()
    {
        for(int i=0; i<totSize; i++) 
        {
            vx0[i] = 0.0f;
            vy0[i] = 0.0f;
            d0[i] = 0.0f;
        }
    }

private:

    //animtate
    void setBoundary(float *value, int flag);
    void lin_solve(float *value, float * value0, float a, float c, int flag);
    void advection(float *value, float *value0, float *u, float *v, int flag);
    void diffusion(float *value, float *value0, float diff, int flag);
    void projection();

private:
    int running;
    float time_step;
    float diff;
    float visc;
    float force;
    float source;

    int rowSize;
    int colSize;
    int totSize;

    float minX;
    float minY;
    float maxX;
    float maxY;

    float *px;
    float *py;

    float *vx;
    float *vy;

    float *vx0;
    float *vy0;

    float *d;
    float *d0;

    float *tx;
    float *ty;

    float *tx0;
    float *ty0;

    float *p;
    float *div;
};

#endif
