/** File:		GridStableSolver.cpp
 ** Author:		Dongli Zhang
 ** Contact:	dongli.zhang0129@gmail.com
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

#include "GridStableSolver.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SWAP(value0,value) {float *tmp=value0;value0=value;value=tmp;}

StableSolver::StableSolver()
{
}

StableSolver::~StableSolver()
{
	free(vx);
	free(vy);
	free(vx0);
	free(vy0);
	free(d);
	free(d0);
	free(px);
	free(py);
	free(div);
	free(p);

	//vorticity confinement
	free(vort);
	free(absVort);
	free(gradVortX);
	free(gradVortY);
	free(lenGrad);
	free(vcfx);
	free(vcfy);
}

void StableSolver::init()
{
	rowSize = 128;
	colSize = 128;
	totSize = rowSize*colSize;
	h = 1.0f;
	simSizeX = (float)rowSize;
	simSizeY = (float)colSize;
	minX = 1.0f;
	maxX = rowSize-1.0f;
	minY = 1.0f;
	maxY = colSize-1.0f;

	running = 1;
	visc = 0.0f;
	diff = 0.0f;
	vorticity = 0.0f;
	timeStep = 1.0f;

	vx = (float *)malloc(sizeof(float)*totSize);
	vy = (float *)malloc(sizeof(float)*totSize);
	vx0 = (float *)malloc(sizeof(float)*totSize);
	vy0 = (float *)malloc(sizeof(float)*totSize);
	d = (float *)malloc(sizeof(float)*totSize);
	d0 = (float *)malloc(sizeof(float)*totSize);
	px = (float *)malloc(sizeof(float)*totSize);
	py = (float *)malloc(sizeof(float)*totSize);
	div = (float *)malloc(sizeof(float)*totSize);
	p = (float *)malloc(sizeof(float)*totSize);

	//vorticity confinement
	vort = (float *)malloc(sizeof(float)*totSize);
	absVort = (float *)malloc(sizeof(float)*totSize);
	gradVortX = (float *)malloc(sizeof(float)*totSize);
	gradVortY = (float *)malloc(sizeof(float)*totSize);
	lenGrad = (float *)malloc(sizeof(float)*totSize);
	vcfx = (float *)malloc(sizeof(float)*totSize);
	vcfy = (float *)malloc(sizeof(float)*totSize);

	for(int i=0; i<rowSize; i++)
	{
		for(int j=0; j<colSize; j++)
		{
			px[cIdx(i, j)] = (float)i+0.5f;
			py[cIdx(i, j)] = (float)j+0.5f;
		}
	}
}

void StableSolver::reset()
{
	for(int i=0; i<totSize; i++)
	{
		vx[i] = 0.0f;
		vy[i] = 0.0f;
		d[i] = 0.0f;
	}
}

void StableSolver::cleanBuffer()
{
	memset(vx0, 0, sizeof(float)*totSize);
	memset(vy0, 0, sizeof(float)*totSize);
	memset(d0, 0, sizeof(float)*totSize);
}

void StableSolver::setBoundary(float *value, int flag)
{
	//for velocity along x-axis
	if(flag == 1)
	{
		for(int i=1; i<=rowSize-2; i++)
		{
			value[cIdx(i, 0)] = value[cIdx(i, 1)];
			value[cIdx(i, colSize-1)] = value[cIdx(i, colSize-2)];
		}
		for(int j=1; j<=colSize-2; j++)
		{
			value[cIdx(0, j)] = -value[cIdx(1, j)];
			value[cIdx(rowSize-1, j)] = -value[cIdx(rowSize-2, j)];
		}
	}

	//for velocity along y-axis
	if(flag == 2)
	{
		for(int i=1; i<=rowSize-2; i++)
		{
			value[cIdx(i, 0)] = -value[cIdx(i, 1)];
			value[cIdx(i, colSize-1)] = -value[cIdx(i, colSize-2)];
		}
		for(int j=1; j<=colSize-2; j++)
		{
			value[cIdx(0, j)] = value[cIdx(1, j)];
			value[cIdx(rowSize-1, j)] = value[cIdx(rowSize-2, j)];
		}
	}

	//density
	if(flag == 0)
	{
		for(int i=1; i<=rowSize-2; i++)
		{
			value[cIdx(i, 0)] = value[cIdx(i, 1)];
			value[cIdx(i, colSize-1)] = value[cIdx(i, colSize-2)];
		}
		for(int j=1; j<=colSize-2; j++)
		{
			value[cIdx(0, j)] = value[cIdx(1, j)];
			value[cIdx(rowSize-1, j)] = value[cIdx(rowSize-2, j)];
		}
	}

	value[cIdx(0, 0)] = (value[cIdx(0, 1)]+value[cIdx(1, 0)])/2;
	value[cIdx(rowSize-1, 0)] = (value[cIdx(rowSize-2, 0)]+value[cIdx(rowSize-1, 1)])/2;
	value[cIdx(0, colSize-1)] = (value[cIdx(0, colSize-2)]+value[cIdx(1, colSize-1)])/2;
	value[cIdx(rowSize-1, colSize-1)] = (value[cIdx(rowSize-2, colSize-1)]+value[cIdx(rowSize-1, colSize-2)])/2;
}

void StableSolver::projection()
{
	for(int i=1; i<=rowSize-2; i++)
	{
		for(int j=1; j<=colSize-2; j++)
		{
			div[cIdx(i, j)] = 0.5f * (vx[cIdx(i+1, j)]-vx[cIdx(i-1, j)]+vy[cIdx(i, j+1)]-vy[cIdx(i, j-1)]);
			p[cIdx(i, j)] = 0.0f;;
		}
	}
	setBoundary(div, 0);
	setBoundary(p, 0);

	//projection iteration
	for(int k=0; k<20; k++)
	{
		for(int i=1; i<=rowSize-2; i++)
		{
			for(int j=1; j<=colSize-2; j++)
			{
				p[cIdx(i, j)] = (p[cIdx(i+1, j)]+p[cIdx(i-1, j)]+p[cIdx(i, j+1)]+p[cIdx(i, j-1)]-div[cIdx(i, j)])/4.0f;
			}
		}
		setBoundary(p, 0);
	}

	//velocity minus grad of Pressure
	for(int i=1; i<=rowSize-2; i++)
	{
		for(int j=1; j<=colSize-2; j++)
		{
			vx[cIdx(i, j)] -= 0.5f*(p[cIdx(i+1, j)]-p[cIdx(i-1, j)]);
			vy[cIdx(i, j)] -= 0.5f*(p[cIdx(i, j+1)]-p[cIdx(i, j-1)]);
		}
	}
	setBoundary(vx, 1);
	setBoundary(vy, 2);
}

void StableSolver::advection(float *value, float *value0, float *u, float *v, int flag)
{
	float oldX;
	float oldY;
	int i0;
	int i1;
	int j0;
	int j1;
	float wL;
	float wR;
	float wB;
	float wT;

	for(int i=1; i<=rowSize-2; i++)
	{
		for(int j=1; j<=colSize-2; j++)
		{
			oldX = px[cIdx(i, j)] - u[cIdx(i, j)]*timeStep;
			oldY = py[cIdx(i, j)] - v[cIdx(i, j)]*timeStep;

			if(oldX < minX) oldX = minX;
			if(oldX > maxX) oldX = maxX;
			if(oldY < minY) oldY = minY;
			if(oldY > maxY) oldY = maxY;

			i0 = (int)(oldX-0.5f);
			j0 = (int)(oldY-0.5f);
			i1 = i0+1;
			j1 = j0+1;
			
			wL = px[cIdx(i1, j0)]-oldX;
			wR = 1.0f-wL;
			wB = py[cIdx(i0, j1)]-oldY;
			wT = 1.0f-wB;

			value[cIdx(i, j)] = wB*(wL*value0[cIdx(i0, j0)]+wR*value0[cIdx(i1, j0)])+
								wT*(wL*value0[cIdx(i0, j1)]+wR*value0[cIdx(i1, j1)]);
		}
	}
	
	setBoundary(value, flag);
}

void StableSolver::diffusion(float *value, float *value0, float rate, int flag)
{
	for(int i=0; i<totSize; i++) value[i] = 0.0f;
	float a = rate*timeStep;

	for(int k=0; k<20; k++)
	{
		for(int i=1; i<=rowSize-2; i++)
		{
			for(int j=1; j<=colSize-2; j++)
			{
				value[cIdx(i, j)] = (value0[cIdx(i, j)]+a*(value[cIdx(i+1, j)]+value[cIdx(i-1, j)]+value[cIdx(i, j+1)]+value[cIdx(i, j-1)])) / (4.0f*a+1.0f);
			}
		}
		setBoundary(value, flag);
	}
}

void StableSolver::vortConfinement()
{
	for(int i=1; i<=rowSize-2; i++)
	{
		for(int j=1; j<=colSize-2; j++)
		{
			vort[cIdx(i, j)] = 0.5f*(vy[cIdx(i+1, j)]-vy[cIdx(i-1, j)]-vx[cIdx(i, j+1)]+vx[cIdx(i, j-1)]);
			if(vort[cIdx(i, j)] >= 0.0f) absVort[cIdx(i, j)] = vort[cIdx(i, j)];
			else absVort[cIdx(i, j)] = -vort[cIdx(i, j)];
		}
	}
	setBoundary(vort, 0);
	setBoundary(absVort, 0);

	for(int i=1; i<=rowSize-2; i++)
	{
		for(int j=1; j<=colSize-2; j++)
		{
			gradVortX[cIdx(i, j)] = 0.5f*(absVort[cIdx(i+1, j)]-absVort[cIdx(i-1, j)]);
			gradVortY[cIdx(i, j)] = 0.5f*(absVort[cIdx(i, j+1)]-absVort[cIdx(i, j-1)]);
			lenGrad[cIdx(i, j)] = sqrt(gradVortX[cIdx(i, j)]*gradVortX[cIdx(i, j)]+gradVortY[cIdx(i, j)]*gradVortY[cIdx(i, j)]);
			if(lenGrad[cIdx(i, j)] < 0.01f)
			{
				vcfx[cIdx(i, j)] = 0.0f;
				vcfy[cIdx(i, j)] = 0.0f;
			}
			else
			{
				vcfx[cIdx(i, j)] = gradVortX[cIdx(i, j)] / lenGrad[cIdx(i, j)];
				vcfy[cIdx(i, j)] = gradVortY[cIdx(i, j)] / lenGrad[cIdx(i, j)];
			}
		}
	}
	setBoundary(vcfx, 0);
	setBoundary(vcfy, 0);

	for(int i=1; i<=rowSize-2; i++)
	{
		for(int j=1; j<=colSize-2; j++)
		{
			vx[cIdx(i, j)] += vorticity * (vcfy[cIdx(i, j)] * vort[cIdx(i, j)]);
			vy[cIdx(i, j)] += vorticity * (-vcfx[cIdx(i, j)] * vort[cIdx(i, j)]);
		}
	}

	setBoundary(vx, 1);
	setBoundary(vy, 2);
}

void StableSolver::addSource()
{
	int index;
	for(int i=1; i<=rowSize-2; i++)
	{
		for(int j=1; j<=colSize-2; j++)
		{
			index = cIdx(i, j);
			vx[index] += vx0[index];
			vy[index] += vy0[index];
			d[index] += d0[index];
		}
	}

	setBoundary(vx, 1);
	setBoundary(vy, 2);
	setBoundary(d, 0);
}

void StableSolver::animVel()
{
	if(diff > 0.0f)
	{
		SWAP(vx0, vx);
		SWAP(vy0, vy);
		diffusion(vx, vx0, diff, 1);
		diffusion(vy, vy0, diff, 2);
	}

	projection();

	SWAP(vx0, vx);
	SWAP(vy0, vy);
	advection(vx, vx0, vx0, vy0, 1);
	advection(vy, vy0, vx0, vy0, 2);

	projection();
}

void StableSolver::animDen()
{
	if(visc > 0.0f)
	{
		SWAP(d0, d);
		diffusion(d, d0, visc, 0);
	}
	SWAP(d0, d);
	advection(d, d0, vx, vy, 0);
}
