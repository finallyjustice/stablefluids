/** File:		Main.cpp
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

#include <GL\glut.h>
#include "GridStableSolver.h"

StableSolver *solver;

int disp_type=0;

int win_x=800;
int win_y=800;

int mouse_down[3];
int omx;
int omy;
int mx;
int my;

void draw_velocity()
{
	float *px = solver->getPX();
	float *py = solver->getPY();
	float *vx = solver->getVX();
	float *vy = solver->getVY();

	glColor3f(0.0f, 1.0f, 0.0f);
	glLineWidth(1.0f);

	glBegin(GL_LINES);
		for(int i=0; i<solver->getTotSize(); i++)
		{
			glVertex2f(px[i], py[i]);
			glVertex2f(px[i]+vx[i]*10.0f, py[i]+vy[i]*10.0f);
		}
	glEnd ();
}

void draw_density()
{
	float x;
	float y;
	float d00;
	float d01;
	float d10;
	float d11;

	int rowSize = solver->getRowSize();
	int colSize = solver->getColSize();

	glBegin(GL_QUADS);
		for(int i=1; i<=rowSize-2; i++)
		{
			x = (float)i;
			for(int j=1; j<=colSize-2; j++)
			{
				y = (float)j;

				d00 = solver->getDens(i, j);
				d01 = solver->getDens(i, j+1);
				d10 = solver->getDens(i+1, j);
				d11 = solver->getDens(i+1, j+1);

				glColor3f(1.0f-d00, 1.0f, 1.0f-d00); glVertex2f(x, y);
				glColor3f(1.0f-d10, 1.0f, 1.0f-d10); glVertex2f(x+1.0f, y);
				glColor3f(1.0f-d11, 1.0f, 1.0f-d11); glVertex2f(x+1.0f, y+1.0f);
				glColor3f(1.0f-d01, 1.0f, 1.0f-d01); glVertex2f(x, y+1.0f);
			}
		}
	glEnd();
}

void get_input()
{
	solver->cleanBuffer();

	int totSize = solver->getTotSize();
	int rowSize = solver->getRowSize();
	int colSize = solver->getColSize();

	int xPos;
    int yPos;

    if(mouse_down[0] || mouse_down[2])
    {
		xPos = (int)((float)(omx)/win_x*(rowSize));
		yPos = (int)((float)(win_y - omy)/win_y*(colSize));

        if(xPos > 0 && xPos < rowSize-1 && yPos > 0 && yPos < colSize-1)
        {
            if(mouse_down[0])
            {
				solver->setVX0(xPos, yPos, 1.0f * (mx - omx));
				solver->setVY0(xPos, yPos, 1.0f * (omy - my));
            }

            if(mouse_down[2])
            {
				solver->setD0(xPos, yPos, 10.0f);
            }

            omx = mx;
            omy = my;
        }

		solver->addSource();
    }
}

void key_func(unsigned char key, int x, int y)
{
	switch(key)
	{
		case 'v':
		case 'V':
			disp_type = (disp_type+1) % 2;
			break;

		case ' ':
			if(solver->isRunning() == 1)
			{
				solver->stop();
			}
			else
			{
				solver->start();
			}
			break;
		case 'c':
		case 'C':
			solver->reset();
			break;
	}
}

void mouse_func(int button, int state, int x, int y)
{
	omx = x;
	omy = y;

	mx = x;
	my = y;

	mouse_down[button] = state == GLUT_DOWN;
}

void motion_func(int x, int y)
{
	mx = x;
	my = y;
}

void reshape_func (int width, int height)
{
	glutReshapeWindow (width, height);

	win_x = width;
	win_y = height;
}

void idle_func()
{
	glutPostRedisplay ();
}

 void display_func()
{
	get_input();
	solver->vortConfinement();
	solver->animVel();
	solver->animDen();

	glViewport(0, 0, win_x, win_y);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity ();
	gluOrtho2D(0.0f, (float)(solver->getRowSize()), 0.0f, (float)(solver->getColSize()));

	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	if(disp_type == 0) draw_density();
	if(disp_type == 1) draw_velocity();
	
	glColor3f(1.0f, 0.0f, 0.0f);
	glPointSize(1.0f);
	glBegin(GL_POINTS);
		for(int i=0; i<solver->getTotSize(); i++)
		{
			glVertex2f(solver->getPX()[i], solver->getPY()[i]);
		}
	glEnd();

	glutSwapBuffers ();
}

int main(int argc, char** argv)
{
	solver=new StableSolver();
	solver->init();
	solver->reset();

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(win_x, win_y);
	glutCreateWindow("StableFluid2D");

	glutKeyboardFunc(key_func);
	glutMouseFunc(mouse_func);
	glutMotionFunc(motion_func);
	glutReshapeFunc(reshape_func);
	glutIdleFunc(idle_func);
	glutDisplayFunc(display_func);

	glutMainLoop ();

	return 0;
}
