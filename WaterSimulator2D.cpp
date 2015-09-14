#include "stdafx.h"
#include "WaterSimulator2D.h"

#include <gl/glut.h>
#include <gl/GL.h>
#include "LevelSet2D.h"
#include "Particle2D.h"
#include "ParticleSet2D.h"
#include "StableFluid2D.h"
#include "MarchingSquares.h"

#include "Timer.h"
using namespace std;

WaterSimulator2D::WaterSimulator2D(int gridSizeX, int gridSizeY)
	:lset(gridSizeX, gridSizeY),
	pset(gridSizeX, gridSizeY),
	grid0(gridSizeX, gridSizeY, gridSizeX, gridSizeY),
	grid1(gridSizeX, gridSizeY, gridSizeX, gridSizeY),
	grid2(gridSizeX, gridSizeY, gridSizeX, gridSizeY),
	gridSizeX(gridSizeX),
	gridSizeY(gridSizeY),
	frame(0),
	particle_mode(0),
	bfill(true),
	pause(true),
	step(false),
	debugMarchingCubes(false),
	render_pressure(false),
	render_surface(true),
	OUTPUT_FILE(false),
	writePOVRAYFile(false)
{
}


WaterSimulator2D::~WaterSimulator2D()
{
}

void WaterSimulator2D::DisplayFunc()
{
	frame++;

	static double time = 0.0;
	static int frames = 0;

	time += timer.GetElapsedTime();
	timer.Reset();
	frames++;

	if (time >= 1.0)
	{
		cout << "Framerate:" << (frames / time) << endl;
		frames = 0;
		time = 0;
	}

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	if (bfill)
	{
		glPolygonMode(GL_FRONT, GL_FILL);
	}
	else
	{
		glPolygonMode(GL_FRONT, GL_LINE);
	}

	int verts[8][3];

	for (int i = 0; i < 8; i++)
	{
		verts[i][0] = i & 1 ? 1 : -1;
		verts[i][1] = i & 2 ? 1 : -1;
		verts[i][2] = i & 4 ? 1 : -1;
	}

	glBegin(GL_LINES);
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j<3; j++)
		{
			glVertex3iv(verts[i]);
			glVertex3iv(verts[i ^ (1 << j)]);
		}
	}
	glEnd();

	if (!pause || step)
	{

		lset.Update(grid0, grid1, grid2, TIMESTEP);
		//lset.UpdateThreaded(grid0, grid1, grid2, TIMESTEP);
		//pset.Update(grid0, grid1, grid2, TIMESTEP);
		pset.UpdateThreaded(grid0, grid1, grid2, TIMESTEP);

		//		lset.Fix(pset);
		lset.Reinitialize();
		//	lset.Reinitialize();
		//lset.CheckGrid();

		static int count = 0;
		count++;

		grid0 = grid2;
		grid1.step(TIMESTEP);
		grid2.step(TIMESTEP);

		if (count % 20 == 0)
		{
			//			pset.Reseed(lset);
		}
		if (count % 10 == 0)
		{
			//	lset.Reinitialize();
		}
		step = false;
	}

	//glEnable(GL_LIGHTING);
	if (render_surface)
	{
		marchingSquares.Update();
		marchingSquares.Draw();
	}

	if (writePOVRAYFile)
	{
		ostringstream outs;
		outs << "povray/water" << frame << ".pov";

		ofstream out(outs.str().c_str());;
		//	out << *surface << endl;
		out.close();



		ostringstream cmds;
		cmds << "\"C:\\Users\\alanga\\AppData\\Roaming\\POV-Ray\\v3.6\\bin\\pvengine64.exe\" /exit -w1600 -h1200 +a0.3 -d +Ipovray/water";
		cmds << frame;
		cmds << ".pov +f +opovray/water";
		cmds << frame;
		cmds << ".bmp -p";
		system(cmds.str().c_str());
	}

	if (particle_mode)
	{
		glDisable(GL_LIGHTING);
		glBegin(GL_POINTS);
		for (int i = 0; i<gridSizeX; i++)
		{
			for (int j = 0; j<gridSizeY; j++)
			{
				ParticleSet2D::Iterator it = pset.begin(i, j);
				while (it != pset.end(i, j))
				{
					double x, y;

					(*it)->GetPosition(x, y);
					x = x * 2 / gridSizeX - 1;
					y = y * 2 / gridSizeY - 1;

					if ((*it)->Sign() < 0)
					{
						glColor3f(1, 0, 0);
						if (particle_mode & 1)
							glVertex3d(x, y, 0);
					}
					else
					{
						glColor3f(0, 0, 1);
						if (particle_mode & 2)
							glVertex3d(x, y, 0);
					}


					it++;
				}
				
			}
		}
		glEnd();
	}
	if (render_pressure)
	{
		glDisable(GL_LIGHTING);
		glBegin(GL_POINTS);
		for (int i = 0; i<gridSizeX; i++)
		{
			for (int j = 0; j<gridSizeY; j++)
			{
				double x = i, y = j;
				x = x * 2 / gridSizeX - 1;
				y = y * 2 / gridSizeY - 1;

				double pressureVal = grid0.pressure[grid0.makeIndex(i, j)];
				//cout << i << " " << j << " " << k << " " << pressureVal << endl;
				double r = pressureVal / 20;
				double g = 0;
				double b = (20 - pressureVal) / 20;
				glColor3f(r, g, b);
				glVertex3d(x, y, 0);
				
			}

		}
		glEnd();
	}

	glColor3f(1, 1, 1);
	glutSwapBuffers();
	/*
	if (OUTPUT_FILE)
	{
	static int frame = 0;
	ostringstream out;
	out << "movie/water" << frame << ".ppm";
	frame++;
	FILE* fp;
	fopen_s(&fp, out.str().c_str(), "w");

	DumpPPM(fp, 0, v_width, v_height);
	fclose(fp);
	}
	*/
}

void WaterSimulator2D::IdleFunc()
{
	glutPostRedisplay();

}

void WaterSimulator2D::ReshapeFunc(int width, int height)
{
	if (height == 0)
		height = 1;

	v_height = height;
	v_width = width;

	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(-1, 1, -1, 1);
	glMatrixMode(GL_MODELVIEW);
}

void WaterSimulator2D::KeyboardFunc(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'a':
		break;
	case 'd':
		break;
	case 'w':
		break;
	case 's':
		break;
	case 'r':
		break;
	case 'f':
		break;
	case 'q':
		break;
	case 'e':
		break;
	case 't':
		break;
	case 'g':
		break;
	case 'n':
		step = true;
		break;
	case ' ':
		pause = !pause;
		break;
	case 'p':
		particle_mode++;
		particle_mode %= 4;
		break;
	case 'z':
		bfill = !bfill;
		break;

	default:
		break;
	}

}


void WaterSimulator2D::MouseFunc(int button, int state, int x, int y)
{
	y = v_height - y;

	if (state == 0)  //button down
	{
		if (button == GLUT_LEFT_BUTTON)
		{

		}
		else if (button == GLUT_MIDDLE_BUTTON)
		{

		}
		else if (button == GLUT_RIGHT_BUTTON)
		{

		}
	}
	else if (state == 1) //button up
	{
		if (button == GLUT_LEFT_BUTTON)
		{

		}
		else if (button == GLUT_MIDDLE_BUTTON)
		{

		}
		else if (button == GLUT_RIGHT_BUTTON)
		{

		}
	}
}
WaterSimulator2D* simulator;

void GlutDisplayFunc2D()
{
	simulator->DisplayFunc();
}

void GlutKeyboardFunc2D(unsigned char key, int x, int y)
{
	simulator->KeyboardFunc(key, x, y);
}

void GlutReshapeFunc2D(int width, int height)
{
	simulator->ReshapeFunc(width, height);
}

void GlutIdleFunc2D()
{
	simulator->IdleFunc();
}

void GlutMouseFunc2D(int button, int state, int x, int y)
{
	simulator->MouseFunc(button, state, x, y);
}

double implicit(const vec2d& o)
{
	double x = (o[0] + 1) * simulator->gridSizeX / 2.0;
	double y = (o[1] + 1) * simulator->gridSizeY / 2.0;
	double d = simulator->lset.LinearSample(x, y);

	double scale = double(simulator->gridSizeX) / 2.0;
	return d / scale;
}

void WaterSimulator2D::RunSimulation()
{
	simulator = this;
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA); // set display mode
	glutInitWindowSize(400, 300); // set window size
	glutInitWindowPosition(0, 0); // set window position on screen
	glutCreateWindow("Glut Window"); // set window title

	glutMouseFunc(GlutMouseFunc2D); // register the mouse action function
	glutKeyboardFunc(GlutKeyboardFunc2D); // register the keyboard action function
	glutDisplayFunc(GlutDisplayFunc2D); //register the redraw function
	glutReshapeFunc(GlutReshapeFunc2D);
	glutIdleFunc(GlutIdleFunc2D);

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

	lset.MakeCircle(.3*gridSizeX, .5*gridSizeY, .15*gridSizeX);
	lset.Reinitialize();

//	grid1.step(TIMESTEP * 0.5);
//	grid2.step(TIMESTEP);

	//surface = new IsoSurface(&lset);

	marchingSquares.setImplicitFunction(implicit);
	marchingSquares.createTable();

	pset.Reseed(lset);
	timer.Reset();
	glutMainLoop();
}
