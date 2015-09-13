// WaterSimulator.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include"WaterSimulator2D.h"
#include"WaterSimulator3D.h"
#include "Tests.h"
#include <gl/glut.h>
#include <gl/GL.h>
#pragma comment(lib,"opengl32.lib")
#pragma comment(lib,"glu32.lib")
#pragma comment(lib,"glut32.lib")

int main(int argc, char** argv)
{
	bool testResults = RunTests();
	if (!testResults)
	{
		return 0;
	}
	glutInit(&argc, argv);

	//WaterSimulator3D simulator(10, 10, 10);
	WaterSimulator2D simulator(100, 100);
	simulator.RunSimulation();
    return 0;
}

