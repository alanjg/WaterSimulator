#pragma once
#include "LevelSet3D.h"
#include "ParticleSet3D.h"
#include "StableFluid3D.h"
#include "MarchingCubes.h"
#include "Timer.h"
#include "PerspectiveCamera.h"

class WaterSimulator3D
{
	bool bfill;
	bool pause;
	bool step;
	int gridSizeX;
	int gridSizeY;
	int gridSizeZ;
	int particle_mode;
	
	bool debugMarchingCubes;
	bool render_pressure;
	bool render_surface;
	bool OUTPUT_FILE;
	bool writePOVRAYFile;

	const double TIMESTEP = 0.1;

	Timer timer;

	PerspectiveCamera cam;
	int frame = 0;
	int v_height, v_width;


	LevelSet3D lset;
	ParticleSet3D pset;
	StableFluid3D grid0, grid1, grid2;
	MarchingCubes marchingCubes;

	friend double implicit(const vec3d& o);
public:
	WaterSimulator3D(int gridSizeX, int gridSizeY, int gridSizeZ);
	~WaterSimulator3D();

	void RunSimulation();

	// GLUT callbacks
	void DisplayFunc();
	void KeyboardFunc(unsigned char key, int x, int y);
	void ReshapeFunc(int width, int height);
	void IdleFunc();
	void MouseFunc(int button, int state, int x, int y);
};

