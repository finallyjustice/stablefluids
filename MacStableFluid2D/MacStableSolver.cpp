/** File:		MacStableSolver.cpp
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

#include "MacStableSolver.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SWAP(value0,value) {float *tmp=value0;value0=value;value=tmp;}

StableSolver::StableSolver()
{
}

StableSolver::~StableSolver()
{

}

void StableSolver::init()
{
	rowCell = 128;
	colCell = 128;
	totCell = rowCell*colCell;
	rowVelX = rowCell+1;
	colVelX = colCell;
	totVelX = rowVelX*colVelX;
	rowVelY = rowCell;
	colVelY = colCell+1;
	totVelY = rowVelY*colVelY;
	minX = 0.0f;
	maxX = (float)rowCell;
	minY = 0.0f;
	maxY = (float)colCell;

	//params
	running = 1;
	timeStep = 1.0f;
	diff = 0.0f;
	visc = 0.0f;

	vx = (float *)malloc(sizeof(float)*totVelX);
	vy = (float *)malloc(sizeof(float)*totVelY);
	vx0 = (float *)malloc(sizeof(float)*totVelX);
	vy0 = (float *)malloc(sizeof(float)*totVelY);
	d = (float *)malloc(sizeof(float)*totCell);
	d0 = (float *)malloc(sizeof(float)*totCell);
	div = (float *)malloc(sizeof(float)*totCell);
	p = (float *)malloc(sizeof(float)*totCell);
	pvx = (Vec2f *)malloc(sizeof(Vec2f)*totVelX);
	pvy = (Vec2f *)malloc(sizeof(Vec2f)*totVelY);

	for(int i=0; i<rowVelX; i++)
	{
		for(int j=0; j<colVelX; j++)
		{
			pvx[vxIdx(i, j)].x = (float)i;
			pvx[vxIdx(i, j)].y = (float)j+0.5f;
		}
	}

	for(int i=0; i<rowVelY; i++)
	{
		for(int j=0; j<colVelY; j++)
		{
			pvy[vyIdx(i, j)].x = (float)i+0.5f;
			pvy[vyIdx(i, j)].y = (float)j;
		}
	}
}

void StableSolver::reset()
{
	for(int i=0; i<totCell; i++) d[i] = 0.0f;
	for(int i=0; i<totVelX; i++) vx[i] = 0.0f;
	for(int i=0; i<totVelY; i++) vy[i] = 0.0f;
}

void StableSolver::cleanBuffer()
{
	for(int i=0; i<totCell; i++) d0[i] = 0.0f;
	for(int i=0; i<totVelX; i++) vx0[i] =0.0f;
	for(int i=0; i<totVelY; i++) vy0[i] = 0.0f;
}

void StableSolver::setVelBoundary(int flag)
{
	//x-axis
	if(flag == 1)
	{
		for(int i=1; i<=rowVelX-2; i++)
		{
			vx[vxIdx(i, 0)] = vx[vxIdx(i, 1)];
			vx[vxIdx(i, colVelX-1)] = vx[vxIdx(i, colVelX-2)];
		}
		for(int j=1; j<=colVelX-2; j++)
		{
			vx[vxIdx(0, j)] = -vx[vxIdx(1, j)];
			vx[vxIdx(rowVelX-1, j)] = -vx[vxIdx(rowVelX-2, j)];
		}
		vx[vxIdx(0, 0)] = (vx[vxIdx(1, 0)]+vx[vxIdx(0, 1)])/2;
		vx[vxIdx(rowVelX-1, 0)] = (vx[vxIdx(rowVelX-2, 0)]+vx[vxIdx(rowVelX-1, 1)])/2;
		vx[vxIdx(0, colVelX-1)] = (vx[vxIdx(1, colVelX-1)]+vx[vxIdx(0, colVelX-2)])/2;
		vx[vxIdx(rowVelX-1, colVelX-1)] = (vx[vxIdx(rowVelX-2, colVelX-1)]+vx[vxIdx(rowVelX-1, colVelX-2)])/2;
	}

	//y-axis
	if(flag == 2)
	{
		for(int i=1; i<=rowVelY-2; i++)
		{
			vy[vyIdx(i, 0)] = -vy[vyIdx(i, 1)];
			vy[vyIdx(i, colVelY-1)] = -vy[vyIdx(i, colVelY-2)];
		}
		for(int j=1; j<=colVelY-2; j++)
		{
			vy[vyIdx(0, j)] = vy[vyIdx(1, j)];
			vy[vyIdx(rowVelY-1, j)] = vy[vyIdx(rowVelY-2, j)];
		}
		vy[vyIdx(0, 0)] = (vy[vyIdx(1, 0)]+vy[vyIdx(0, 1)])/2;
		vy[vyIdx(rowVelY-1, 0)] = (vy[vyIdx(rowVelY-2, 0)]+vy[vyIdx(rowVelY-1, 1)])/2;
		vy[vyIdx(0, colVelY-1)] = (vy[vyIdx(1, colVelY-1)]+vy[vyIdx(0, colVelY-2)])/2;
		vy[vyIdx(rowVelY-1, colVelY-1)] = (vy[vyIdx(rowVelY-2, colVelY-1)]+vy[vyIdx(rowVelY-1, colVelY-2)])/2;
	}
}

void StableSolver::setCellBoundary(float *value)
{
	for(int i=1; i<=rowCell-2; i++)
	{
		value[cIdx(i, 0)] = value[cIdx(i, 1)];
		value[cIdx(i, colCell-1)] = value[cIdx(i, colCell-2)];
	}
	for(int j=1; j<=colCell-2; j++)
	{
		value[cIdx(0, j)] = value[cIdx(1, j)];
		value[cIdx(rowCell-1, j)] = value[cIdx(rowCell-2, j)];
	}
	value[cIdx(0, 0)] = (value[cIdx(1, 0)]+value[cIdx(0, 1)])/2;
	value[cIdx(rowCell-1, 0)] = (value[cIdx(rowCell-2, 0)]+value[cIdx(rowCell-1, 1)])/2;
	value[cIdx(0, colCell-1)] = (value[cIdx(1, colCell-1)]+value[cIdx(0, colCell-2)])/2;
	value[cIdx(rowCell-1, colCell-1)] = (value[cIdx(rowCell-1, colCell-2)]+value[cIdx(rowCell-1, colCell-2)])/2;
}

void StableSolver::projection()
{
	int static count=0;
	for(int i=1; i<=rowCell-2; i++)
	{
		for(int j=1; j<=colCell-2; j++)
		{
			div[cIdx(i, j)] = (vx[vxIdx(i+1, j)]-vx[vxIdx(i, j)]+vy[vyIdx(i, j+1)]-vy[vyIdx(i, j)]);
			p[cIdx(i, j)] = 0.0f;
		}
	}
	count++;
	setCellBoundary(p);
	setCellBoundary(div);

	//projection iteration
	for(int k=0; k<20; k++)
	{
		for(int i=1; i<=rowCell-2; i++)
		{
			for(int j=1; j<=colCell-2; j++)
			{
				p[cIdx(i, j)] = (p[cIdx(i+1, j)]+p[cIdx(i-1, j)]+p[cIdx(i, j+1)]+p[cIdx(i, j-1)]-div[cIdx(i, j)])/4.0f;
			}
		}
		setCellBoundary(p);
	}

	//velocity minus grad of Pressure
	for(int i=1; i<=rowVelX-2; i++)
	{
		for(int j=1; j<=colVelX-2; j++)
		{
			vx[vxIdx(i, j)] -= (p[cIdx(i, j)]-p[cIdx(i-1, j)]);
		}
	}
	for(int i=1; i<=rowVelY-2; i++)
	{
		for(int j=1; j<=colVelY-2; j++)
		{
			vy[vyIdx(i, j)] -= (p[cIdx(i, j)]-p[cIdx(i, j-1)]);
		}
	}
	setVelBoundary(1);
	setVelBoundary(2);
}

void StableSolver::advectVel()
{
	for(int i=1; i<=rowVelX-2; i++)
	{
		for(int j=1; j<=colVelX-2; j++)
		{
			float nvx = vx0[vxIdx(i, j)];
			float nvy = (vy0[vyIdx(i-1, j)]+vy0[vyIdx(i-1, j+1)]+vy0[vyIdx(i, j)]+vy0[vyIdx(i, j+1)])/4;

			float oldX = pvx[vxIdx(i, j)].x - nvx*timeStep;
			float oldY = pvx[vxIdx(i, j)].y - nvy*timeStep;

			if(oldX < 0.5f) oldX = 0.5f;
			if(oldX > maxX-0.5f) oldX = maxX-0.5f;
			if(oldY < 1.0f) oldY = 1.0f;
			if(oldY > maxY-1.0f) oldY = maxY-1.0f;

			int i0 = (int)oldX;
			int j0 = (int)(oldY-0.5f);
			int i1 = i0+1;
			int j1 = j0+1;

			float wL = pvx[vxIdx(i1, j0)].x-oldX;
			float wR = 1.0f-wL;
			float wB = pvx[vxIdx(i0, j1)].y-oldY;
			float wT = 1.0f-wB;

			//printf("%f, %f, %f, %f\n", wL, wR, wB, wT);

			vx[vxIdx(i, j)] = wB*(wL*vx0[vxIdx(i0, j0)]+wR*vx0[vxIdx(i1, j0)])+
							  wT*(wL*vx0[vxIdx(i0, j1)]+wR*vx0[vxIdx(i1, j1)]);
		}
	}

	for(int i=1; i<=rowVelY-2; i++)
	{
		for(int j=1; j<=colVelY-2; j++)
		{
			float nvx = (vx0[vxIdx(i, j-1)]+vx0[vxIdx(i+1, j-1)]+vx0[vxIdx(i, j)]+vx0[vxIdx(i+1, j)])/4;
			float nvy = vy0[vyIdx(i, j)];

			float oldX = pvy[vyIdx(i, j)].x - nvx*timeStep;
			float oldY = pvy[vyIdx(i, j)].y - nvy*timeStep;

			if(oldX < 1.0f) oldX = 1.0f;
			if(oldX > maxX-1.0f) oldX = maxX-1.0f;
			if(oldY < 0.5f) oldY = 0.5f;
			if(oldY > maxY-0.5f) oldY = maxY-0.5f;

			int i0 = (int)(oldX-0.5f);
			int j0 = (int)oldY;
			int i1 = i0+1;
			int j1 = j0+1;

			float wL = pvy[vyIdx(i1, j0)].x-oldX;
			float wR = 1.0f-wL;
			float wB = pvy[vyIdx(i0, j1)].y-oldY;
			float wT = 1.0f-wB;

			vy[vyIdx(i, j)] = wB*(wL*vy0[vyIdx(i0, j0)]+wR*vy0[vyIdx(i1, j0)])+
							  wT*(wL*vy0[vyIdx(i0, j1)]+wR*vy0[vyIdx(i1, j1)]);
		}
	}

	setVelBoundary(1);
	setVelBoundary(2);
}

void StableSolver::advectCell(float *value, float *value0)
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

	for(int i=1; i<=rowCell-2; i++)
	{
		for(int j=1; j<=colCell-2; j++)
		{
			float cvx = getCellVel(i, j).x;
			float cvy = getCellVel(i, j).y;

			oldX = (float)i+0.5f - cvx*timeStep;
			oldY = (float)j+0.5f - cvy*timeStep;

			if(oldX < 1.0f) oldX = 1.0f;
			if(oldX > rowCell-1.0f) oldX = rowCell-1.0f;
			if(oldY < 1.0f) oldY = 1.0f;
			if(oldY > colCell-1.0f) oldY = colCell-1.0f;

			i0 = (int)(oldX-0.5f);
			j0 = (int)(oldY-0.5f);
			i1 = i0+1;
			j1 = j0+1;
			
			wL = (float)i1+0.5f-oldX;
			wR = 1.0f-wL;
			wB = (float)j1+0.5f-oldY;
			wT = 1.0f-wB;

			value[cIdx(i, j)] = wB*(wL*value0[cIdx(i0, j0)]+wR*value0[cIdx(i1, j0)])+
								wT*(wL*value0[cIdx(i0, j1)]+wR*value0[cIdx(i1, j1)]);
		}
	}
	
	setCellBoundary(d);
}

void StableSolver::diffuseVel()
{
	for(int i=0; i<totVelX; i++) vx[i] = 0.0f;
	for(int i=0; i<totVelY; i++) vy[i] = 0.0f;
	float a = diff*timeStep;

	for(int k=0; k<20; k++)
	{
		//diffuse velX
		for(int i=1; i<=rowVelX-2; i++)
		{
			for(int j=1; j<=colVelX-2; j++)
			{
				vx[vxIdx(i, j)] = (vx0[vxIdx(i, j)]+a*(vx[vxIdx(i+1, j)]+vx[vxIdx(i-1, j)]+vx[vxIdx(i, j+1)]+vx[vxIdx(i, j-1)])) / (4.0f*a+1.0f);
			}
		}
		//diffuse velY
		for(int i=1; i<=rowVelY-2; i++)
		{
			for(int j=1; j<=colVelY-2; j++)
			{
				vy[vyIdx(i, j)] = (vy0[vyIdx(i, j)]+a*(vy[vyIdx(i+1, j)]+vy[vyIdx(i-1, j)]+vy[vyIdx(i, j+1)]+vy[vyIdx(i, j-1)])) / (4.0f*a+1.0f);
			}
		}

		//boundary
		setVelBoundary(1);
		setVelBoundary(2);
	}
}

void StableSolver::diffuseCell(float *value, float *value0)
{
	for(int i=0; i<totCell; i++) value[i] = 0.0f;
	float a = visc*timeStep;

	for(int k=0; k<20; k++)
	{
		for(int i=1; i<=rowCell-2; i++)
		{
			for(int j=1; j<=colCell-2; j++)
			{
				value[cIdx(i, j)] = (value0[cIdx(i, j)]+a*(value[cIdx(i+1, j)]+value[cIdx(i-1, j)]+value[cIdx(i, j+1)]+value[cIdx(i, j-1)])) / (4.0f*a+1.0f);
			}
		}
		setCellBoundary(value);
	}
}

void StableSolver::addSource()
{
	for(int i=0; i<totCell; i++) d[i] += d0[i];
	for(int i=0; i<totVelX; i++) vx[i] += vx0[i];
	for(int i=0; i<totVelY; i++) vy[i] += vy0[i];

	setVelBoundary(1);
	setVelBoundary(2);
	setCellBoundary(d);
}

void StableSolver::animVel()
{
	projection();

	if(diff > 0.0f)
	{
		SWAP(vx0, vx);
		SWAP(vy0, vy);
		diffuseVel();
	}

	SWAP(vx0, vx);
	SWAP(vy0, vy);
	advectVel();

	projection();
}

void StableSolver::animDen()
{
	if(visc > 0.0f)
	{
		SWAP(d0, d);
		diffuseCell(d, d0);
	}

	SWAP(d0, d);
	advectCell(d, d0);
}
