/** File:		main.cpp
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

#include <GL\glaux.h>
#include <GL\glut.h>
#include "StableSolver2D.h"

#pragma comment(lib, "glaux.lib")

StableSolver2D *solver;

int disp_type=0;

int win_x=800;
int win_y=800;

int mouse_down[3];
int omx;
int omy;
int mx;
int my;

GLuint tex;

int LoadGLTextures(GLuint& unTexture, const char* chFileName)                
{
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    AUX_RGBImageRec *TextureImage;                  
    TextureImage = auxDIBImageLoad(chFileName); 

    glGenTextures(1, &unTexture);                    
    glBindTexture(GL_TEXTURE_2D, unTexture);

    glTexImage2D(GL_TEXTURE_2D, 0, 3, TextureImage->sizeX, TextureImage->sizeY, 0, GL_RGB, GL_UNSIGNED_BYTE, TextureImage->data);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    if(TextureImage)           
    {
        if(TextureImage->data)
        {
            free(TextureImage->data);
        }

        free(TextureImage);
    }

    return 1;
}

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
		for(int i=1; i<=rowSize; i++)
		{
			x = (float)i;
			for(int j=1; j<=colSize; j++)
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

void draw_texture()
{
	glBindTexture(GL_TEXTURE_2D, tex);
	glColor3f(1.0f, 1.0f, 1.0f);

	float rowSize = (float)(solver->getRowSize());
	float colSize = (float)(solver->getColSize());

	float x;
	float y;

	float *tx = solver->getTX();
	float *ty = solver->getTY();

	glBegin(GL_QUADS); 
		for(int i=0; i<rowSize; i++)
		{
			x = (float)i;

			for(int j=0; j<colSize; j++)
			{
				y = (float)j;

				glTexCoord2f((tx[solver->getIndex(i, j)] - 0.5f)/((float)rowSize), (ty[solver->getIndex(i, j)] - 0.5f)/((float)rowSize)); glVertex2f(x+1.0f, y+1.0f);
				glTexCoord2f((tx[solver->getIndex(i+1, j)] - 0.5f)/((float)rowSize), (ty[solver->getIndex(i+1, j)] - 0.5f)/((float)rowSize)); glVertex2f(x+2.0f, y+1.0f);
				glTexCoord2f((tx[solver->getIndex(i+1, j+1)] - 0.5f)/((float)rowSize), (ty[solver->getIndex(i+1, j+1)] - 0.5f)/((float)rowSize)); glVertex2f(x+2.0f, y+2.0f);
				glTexCoord2f((tx[solver->getIndex(i, j+1)] - 0.5f)/((float)rowSize), (ty[solver->getIndex(i, j+1)] - 0.5f)/((float)rowSize)); glVertex2f(x+1.0f, y+2.0f);
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
		xPos = (int)((float)(omx)/win_x*(rowSize+2));
		yPos = (int)((float)(win_y - omy)/win_y*(colSize+2));

        if(xPos > 0 && xPos < rowSize+1 && yPos > 0 && yPos < colSize+1)
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
			disp_type = (disp_type+1) % 3;
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
			solver->clear();
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
	solver->anim_vel();
	solver->anim_tex();
	solver->anim_den();

	glViewport(0, 0, win_x, win_y);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity ();
	gluOrtho2D(0.0f, (float)(solver->getRowSize()+2), 0.0f, (float)(solver->getColSize()+2));

	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	if(disp_type == 2) 
	{
		draw_density ();
	}

	if(disp_type == 1) 
	{
		draw_velocity();
	}
		
	if(disp_type == 0)
	{
		draw_texture();
	}

	/*glColor3f(1.0f, 0.0f, 0.0f);
	glPointSize(1.0f);
	glBegin(GL_POINTS);
		for(int i=0; i<solver->getTotSize(); i++)
		{
			glVertex2f(solver->getPX()[i], solver->getPY()[i]);
		}
	glEnd();*/

	glutSwapBuffers ();
}

int main(int argc, char** argv)
{
	solver=new StableSolver2D();
	solver->reset(128, 128);

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(win_x, win_y);
	glutCreateWindow("StableFluid2D");

	glEnable(GL_TEXTURE_2D);
	LoadGLTextures(tex, "ali.bmp");

	glutKeyboardFunc(key_func);
	glutMouseFunc(mouse_func);
	glutMotionFunc(motion_func);
	glutReshapeFunc(reshape_func);
	glutIdleFunc(idle_func);
	glutDisplayFunc(display_func);

	glutMainLoop ();

	return 0;
}
