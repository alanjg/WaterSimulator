#pragma once
#include "LevelSet2D.h"
#include "ParticleSet2D.h"
#include "StableFluid2D.h"
#include "MarchingSquares.h"
#include "Timer.h"

class WaterSimulator2D
{
	bool bfill;
	bool pause;
	bool step;
	int gridSizeX;
	int gridSizeY;
	int particle_mode;

	bool debugMarchingCubes;
	bool render_pressure;
	bool render_surface;
	bool OUTPUT_FILE;
	bool writePOVRAYFile;

	const double TIMESTEP = 0.1;

	Timer timer;

	int frame = 0;
	int v_height, v_width;


	LevelSet2D lset;
	ParticleSet2D pset;
	StableFluid2D grid0, grid1, grid2;
	MarchingSquares marchingSquares;

	friend double implicit(const vec2d& o);
public:
	WaterSimulator2D(int gridSizeX, int gridSizeY);
	~WaterSimulator2D();

	void RunSimulation();

	// GLUT callbacks
	void DisplayFunc();
	void KeyboardFunc(unsigned char key, int x, int y);
	void ReshapeFunc(int width, int height);
	void IdleFunc();
	void MouseFunc(int button, int state, int x, int y);
};

