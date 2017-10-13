/** File:    StableSolver2D.cpp
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

#include "StableSolver2D.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define SWAP(x0,x) {float *tmp=x0; x0=x; x=tmp;}

StableSolver2D::StableSolver2D()
{
    running = 1;
    time_step = 1.0f;
    diff = 0.0f;
    visc = 0.0f;
    force = 5.0f;
    source = 2.0f;
}

StableSolver2D::~StableSolver2D()
{
    free(px);
    free(py);
    free(vx);
    free(vy);
    free(vx0);
    free(vy0);
    free(d);
    free(d0);
    free(tx);
    free(ty);
    free(tx0);
    free(ty0);
    free(p);
    free(div);
}

void StableSolver2D::reset(int _rowSize, int _colSize)
{
    rowSize = _rowSize;
    colSize = _colSize;
    totSize = (rowSize+2)*(colSize+2);

    minX = 1.0f;
    minY = 1.0f;
    maxX = (float)(rowSize+1);
    maxY = (float)(colSize+1);

    px = (float *)malloc(sizeof(float)*totSize);
    py = (float *)malloc(sizeof(float)*totSize);

    vx = (float *)malloc(sizeof(float)*totSize);
    vy = (float *)malloc(sizeof(float)*totSize);

    vx0 = (float *)malloc(sizeof(float)*totSize);
    vy0 = (float *)malloc(sizeof(float)*totSize);

    d  = (float *)malloc(sizeof(float)*totSize);
    d0 = (float *)malloc(sizeof(float)*totSize);

    tx = (float *)malloc(sizeof(float)*totSize);
    ty = (float *)malloc(sizeof(float)*totSize);

    tx0 = (float *)malloc(sizeof(float)*totSize);
    ty0 = (float *)malloc(sizeof(float)*totSize);

    p   = (float *)malloc(sizeof(float)*totSize);
    div = (float *)malloc(sizeof(float)*totSize);

    clear();
}

void StableSolver2D::clear()
{
    int index;

    for(int i=0; i<rowSize+2; i++)
    {
        for(int j=0; j<colSize+2; j++)
        {
            index = getIndex(i, j);

            px[index] = 0.0f;
            py[index] = 0.0f;
            vx[index] = 0.0f;
            vy[index] = 0.0f;
            vx0[index] = 0.0f;
            vy0[index] = 0.0f;
            d[index] = 0.0f;
            d0[index] = 0.0f;
            tx[index] = 0.0f;
            ty[index] = 0.0f;
            tx0[index] = 0.0f;
            ty0[index] = 0.0f;
            p[index] = 0.0f;
            div[index] = 0.0f;
        }
    }

    for(int i=0; i<rowSize+2; i++)
    {
        for(int j=0; j<colSize+2; j++)
        {
            index = getIndex(i, j);

            px[index] = i + 0.5f;
            py[index] = j + 0.5f;

            tx[index] = i + 0.5f;
            ty[index] = j + 0.5f;
        }
    }
}

void StableSolver2D::addSource()
{
    if(running == 0) return;

    for(int i=0; i<totSize; i++)
    {
        vx[i] += vx0[i];
        vy[i] += vy0[i];
        d[i]  += d0[i];
    }

    setBoundary(vx, 1);
    setBoundary(vy, 2);
    setBoundary(d, 0);
}

//Animating Velocity
void StableSolver2D::anim_vel()
{
    if(running == 0) return;

    SWAP(vx0, vx); 
    SWAP(vy0, vy); 
    diffusion(vx, vx0, visc, 1);
    diffusion(vy, vy0, visc, 2);

    projection();

    SWAP(vx0, vx);
    SWAP(vy0, vy);
    advection(vx, vx0, vx0, vy0, 1);
    advection(vy, vy0, vx0, vy0, 2);

    projection();
}

//Animating Density
void StableSolver2D::anim_den()
{
    if(running == 0) return;

    SWAP(d0, d); 
    advection(d, d0, vx, vy, 0);

    SWAP(d0, d); 
    diffusion(d, d0, diff, 0);
}

void StableSolver2D::anim_tex()
{
    if(running == 0) return;

    SWAP(tx0, tx); 
    SWAP(ty0, ty); 
    advection(tx, tx0, vx, vy, 0);
    advection(ty, ty0, vx, vy, 0);

    SWAP(tx0, tx); 
    SWAP(ty0, ty);  
    diffusion(tx, tx0, diff, 0);
    diffusion(ty, ty0, diff, 0);
}

void StableSolver2D::setBoundary(float *value, int flag)
{
    int dim = rowSize;

    for(int i=1; i<=dim; i++) 
    {
        if(flag == 1)
        {
            value[getIndex(0, i)]       = -value[getIndex(1,i)];
            value[getIndex(dim+1,i)]    = -value[getIndex(dim,i)];
            value[getIndex(i,0  )]      = value[getIndex(i,1)];
            value[getIndex(i,dim+1)]    = value[getIndex(i,dim)];

        }

        if(flag == 2)
        {
            value[getIndex(0, i)]       = value[getIndex(1,i)];
            value[getIndex(dim+1,i)]    = value[getIndex(dim,i)];
            value[getIndex(i,0  )]      = -value[getIndex(i,1)];
            value[getIndex(i,dim+1)]    = -value[getIndex(i,dim)];
        }

        if(flag == 0)
        {
            value[getIndex(0, i)]       = value[getIndex(1,i)];
            value[getIndex(dim+1,i)]    = value[getIndex(dim,i)];
            value[getIndex(i,0  )]      = value[getIndex(i,1)];
            value[getIndex(i,dim+1)]    = value[getIndex(i,dim)];
        }
    }

    value[getIndex(0, 0)]            = 0.5f*(value[getIndex(1, 0)]+value[getIndex(0, 1)]);
    value[getIndex(0, dim+1)]        = 0.5f*(value[getIndex(1, dim+1)]+value[getIndex(0, dim)]);
    value[getIndex(dim+1, 0)]        = 0.5f*(value[getIndex(dim, 0)]+value[getIndex(dim+1, 1)]);
    value[getIndex(dim+1, dim+1)]    = 0.5f*(value[getIndex(dim, dim+1)]+value[getIndex(dim+1, dim)]);
}

void StableSolver2D::lin_solve(float *value, float * value0, float a, float c, int flag)
{
    for(int iteration=0; iteration<20; iteration++) 
    {
        for(int i=1; i<=rowSize; i++) 
        { 
            for(int j=1; j<=colSize; j++) 
            {
                value[getIndex(i, j)] = (value0[getIndex(i, j)] + a*(value[getIndex(i-1, j)]+value[getIndex(i+1, j)]+value[getIndex(i, j-1)]+value[getIndex(i, j+1)]))/c;
            } 
        }

        setBoundary(value, flag);
    }
}

void StableSolver2D::advection(float *value, float *value0,  float *u, float *v, int flag)
{
    int idxNow;
    float oldX;
    float oldY;

    int i0;
    int j0;
    int i1;
    int j1;

    float iL;
    float iR;
    float iT;
    float iB;

    for(int i=1; i<=rowSize; i++)
    {
        for(int j=1; j<=colSize; j++)
        {
            idxNow = getIndex(i, j);

            //implicit method, trace the position back to old position
            oldX = px[idxNow] - u[idxNow] * time_step;
            oldY = py[idxNow] - v[idxNow] * time_step;

            if(oldX < minX) oldX = minX;
            if(oldX > maxX) oldX = maxX;
            if(oldY < minY) oldY = minY;
            if(oldY > maxY) oldY = maxY;

            i0 = int(oldX - 0.5f);
            j0 = int(oldY - 0.5f);

            i1 = i0 + 1;
            j1 = j0 + 1;

            iR = oldX - px[getIndex(i0, j0)];
            iT = oldY - py[getIndex(i0, j0)];
            iL = 1.0f - iR;
            iB = 1.0f - iT;

            value[idxNow] = iB * (iL*value0[getIndex(i0, j0)] + iR*value0[getIndex(i1, j0)]) +
                            iT * (iL*value0[getIndex(i0, j1)] + iR*value0[getIndex(i1, j1)]);
        }
    }

    setBoundary(value, flag);
}

void StableSolver2D::diffusion(float *value, float *value0, float diff, int flag)
{
    float a=time_step*diff; 
    lin_solve(value, value0, a, 1+4*a, flag);
}

void StableSolver2D::projection()
{
    for(int i=1; i<=rowSize; i++) 
    { 
        for(int j=1; j<=colSize; j++) 
        {
            div[getIndex(i, j)] = -0.5f*(vx[getIndex(i+1,j)]-vx[getIndex(i-1,j)]+vy[getIndex(i,j+1)]-vy[getIndex(i,j-1)]);
            p[getIndex(i,j)] = 0;
        }
    }
    setBoundary(div, 0); 
    setBoundary(p, 0);

    lin_solve(p, div, 1.0, 4.0, 0);

    for(int i=1; i<=rowSize; i++) 
    { 
        for(int j=1; j<=colSize; j++) 
        {
            vx[getIndex(i,j)] -= 0.5f*(p[getIndex(i+1,j)]-p[getIndex(i-1,j)]);
            vy[getIndex(i,j)] -= 0.5f*(p[getIndex(i,j+1)]-p[getIndex(i,j-1)]);
        }
    }
    setBoundary(vx, 1); 
    setBoundary(vy, 2);
}
