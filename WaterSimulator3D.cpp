#include "stdafx.h"
#include "WaterSimulator3D.h"

#include <gl/glut.h>
#include <gl/GL.h>
#include "LevelSet3D.h"
#include "Particle3D.h"
#include "ParticleSet3D.h"
#include "StableFluid3D.h"
#include "MarchingCubes.h"
#include "PerspectiveCamera.h"
#include "Timer.h"
using namespace std;

WaterSimulator3D::WaterSimulator3D(int gridSizeX, int gridSizeY, int gridSizeZ)
	:lset(gridSizeX, gridSizeY, gridSizeZ), 
	pset(gridSizeX, gridSizeY, gridSizeZ), 
	grid0(gridSizeX, gridSizeY, gridSizeZ, gridSizeX, gridSizeY, gridSizeZ),
	grid1(gridSizeX, gridSizeY, gridSizeZ, gridSizeX, gridSizeY, gridSizeZ),
	grid2(gridSizeX, gridSizeY, gridSizeZ, gridSizeX, gridSizeY, gridSizeZ),
	gridSizeX(gridSizeX), 
	gridSizeY(gridSizeY),
	gridSizeZ(gridSizeZ),
	cam(0,0,5),
	frame(0),
	particle_mode(0),
	bfill(true),
	pause(true),
	step(false),
	debugMarchingCubes(false),
	render_pressure(true),
	render_surface(true),
	OUTPUT_FILE(false),
	writePOVRAYFile(false)
{
}


WaterSimulator3D::~WaterSimulator3D()
{
}

bool vec3flt(vec3d& lhs, vec3d& rhs)
{
	if (lhs[0] < rhs[0]) return true;
	if (lhs[0] > rhs[0]) return false;

	if (lhs[1] < rhs[1]) return true;
	if (lhs[1] > rhs[1]) return false;

	if (lhs[2] < rhs[2]) return true;
	if (lhs[2] > rhs[2]) return false;

	return false;
}

bool tri3flt(triangle3d& lhs, triangle3d& rhs)
{
	if (vec3flt(lhs.v[0], rhs.v[0])) return true;
	if (vec3flt(rhs.v[0], lhs.v[0])) return false;

	if (vec3flt(lhs.v[1], rhs.v[1])) return true;
	if (vec3flt(rhs.v[1], lhs.v[1])) return false;

	if (vec3flt(lhs.v[2], rhs.v[2])) return true;
	if (vec3flt(rhs.v[2], lhs.v[2])) return false;

	return false;
}

void WaterSimulator3D::DisplayFunc()
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

	cam.Draw();
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
		
		//lset.Update(grid0, grid1, grid2, TIMESTEP);
		lset.UpdateThreaded(grid0, grid1, grid2, TIMESTEP);
		//pset.Update(grid0, grid1, grid2, TIMESTEP);
		pset.UpdateThreaded(grid0, grid1, grid2, TIMESTEP);
		
		//		lset.Fix(pset);
		lset.Reinitialize();
		//	lset.Reinitialize();
		//		lset.CheckGrid();

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

	glEnable(GL_LIGHTING);
	if (render_surface)
	{
		
		marchingCubes.Update();
			
		marchingCubes.Draw();
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
				for (int k = 0; k<gridSizeZ; k++)
				{
					ParticleSet3D::Iterator it = pset.begin(i, j, k);
					while (it != pset.end(i, j, k))
					{
						double x, y, z;

						(*it)->GetPosition(x, y, z);
						x = x * 2 / gridSizeX - 1;
						y = y * 2 / gridSizeY - 1;
						z = z * 2 / gridSizeZ - 1;

						if ((*it)->Sign() < 0)
						{
							glColor3f(1, 0, 0);
							if (particle_mode & 1)
								glVertex3d(x, y, z);
						}
						else
						{
							glColor3f(0, 0, 1);
							if (particle_mode & 2)
								glVertex3d(x, y, z);
						}


						it++;
					}
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
				for (int k = 0; k<gridSizeZ; k++)
				{
					double x = i, y = j, z = k;
					x = x * 2 / gridSizeX - 1;
					y = y * 2 / gridSizeY - 1;
					z = z * 2 / gridSizeZ - 1;

					double pressureVal = grid0.pressure[grid0.makeIndex(i, j, k)];
					//cout << i << " " << j << " " << k << " " << pressureVal << endl;
					double r = pressureVal / 20;
					double g = 0;
					double b = (20 - pressureVal) / 20;
					glColor3f(r, g, b);
					glVertex3d(x, y, z);
				}
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

void WaterSimulator3D::IdleFunc()
{
	glutPostRedisplay();

}

void WaterSimulator3D::ReshapeFunc(int width, int height)
{
	if (height == 0)
		height = 1;

	v_height = height;
	v_width = width;

	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45, float(width) / float(height), 0.01, 1000);
	glMatrixMode(GL_MODELVIEW);
}

void WaterSimulator3D::KeyboardFunc(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'a':
		cam.MoveLeft(0.2f);
		break;
	case 'd':
		cam.MoveRight(0.2f);
		break;
	case 'w':
		cam.MoveForward(0.2f);
		break;
	case 's':
		cam.MoveBackward(0.2f);
		break;
	case 'r':
		cam.MoveUp(0.2f);
		break;
	case 'f':
		cam.MoveDown(0.2f);
		break;
	case 'q':
		cam.RotateY(-3);
		break;
	case 'e':
		cam.RotateY(3);
		break;
	case 't':
		cam.RotateX(3);
		break;
	case 'g':
		cam.RotateX(-3);
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
	case 'o':
		marchingCubes.level++;
		cout << "level is " << marchingCubes.level << endl;
		break;
	case 'i':
		marchingCubes.level--;
		cout << "level is " << marchingCubes.level << endl;
		break;
	case 'z':
		bfill = !bfill;
		break;

	default:
		break;
	}

}


void WaterSimulator3D::MouseFunc(int button, int state, int x, int y)
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
WaterSimulator3D* simulator;

void GlutDisplayFunc()
{
	simulator->DisplayFunc();
}

void GlutKeyboardFunc(unsigned char key, int x, int y)
{
	simulator->KeyboardFunc(key, x, y);
}

void GlutReshapeFunc(int width, int height)
{
	simulator->ReshapeFunc(width, height);
}

void GlutIdleFunc()
{
	simulator->IdleFunc();
}

void GlutMouseFunc(int button, int state, int x, int y)
{
	simulator->MouseFunc(button, state, x, y);
}

double implicit(const vec3d& o)
{
	double x = (o[0] + 1) * simulator->gridSizeX / 2.0;
	double y = (o[1] + 1) * simulator->gridSizeY / 2.0;
	double z = (o[2] + 1) * simulator->gridSizeZ / 2.0;
	double d = simulator->lset.LinearSample(x, y, z);

	double scale = double(simulator->gridSizeX) / 2.0;
	return d / scale;
}

void WaterSimulator3D::RunSimulation()
{
	simulator = this;
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA); // set display mode
	glutInitWindowSize(400, 300); // set window size
	glutInitWindowPosition(0, 0); // set window position on screen
	glutCreateWindow("Glut Window"); // set window title

	glutMouseFunc(GlutMouseFunc); // register the mouse action function
	glutKeyboardFunc(GlutKeyboardFunc); // register the keyboard action function
	glutDisplayFunc(GlutDisplayFunc); //register the redraw function
	glutReshapeFunc(GlutReshapeFunc);
	glutIdleFunc(GlutIdleFunc);

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST);
	cam.RotateY(180);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

	lset.MakeSphere(.5*gridSizeX, .9*gridSizeY, .5*gridSizeZ, .1*gridSizeX);

	//lset.MakeSphere(.65*grid_sizex,.65*grid_sizey,.65*grid_sizez,.15*grid_sizex );
	lset.Reinitialize();
	//lset.CheckGrid();
	for (int i = 0; i<gridSizeX; i++)
	{
		for (int j = 0; j<gridSizeY; j++)
		{
			for (int k = 0; k<gridSizeZ; k++)
			{
				vec3d center(.45*gridSizeX, .45*gridSizeY, .45*gridSizeZ);
				vec3d here(i, j, k);

				if ((center - here).Magnitude() < .45*gridSizeX)
				{
					int index = grid0.makeIndex(i, j, k);
					grid0.pressure[index] = 1.0;
					grid1.pressure[index] = 1.0;
					grid2.pressure[index] = 1.0;
				}
			}
		}
	}

	
	grid1.step(TIMESTEP * 0.5);
	grid2.step(TIMESTEP);
	
	//surface = new IsoSurface(&lset);

	marchingCubes.setImplicitFunction(implicit);
	marchingCubes.createTable();

	pset.Reseed(lset);
	timer.Reset();
	glutMainLoop();
}
