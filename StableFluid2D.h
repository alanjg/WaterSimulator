#pragma once
#include "vec.h"
#include "Vector.h"
#include "SparseVector.h"
#include "SparseMatrix.h"

static int xOffset2[2] = { 1, 0 };
static int yOffset2[2] = { 0, 1 };


class StableFluid2D
{
	const double GRAV_ACCEL = 9.8;
	const double PCG_EPS = 1e-10;
	const double PCG_MAXITER = 500;
public:
	StableFluid2D();
	~StableFluid2D();
	StableFluid2D(double width_, double height_, int resX_, int resY_);
	StableFluid2D(const StableFluid2D& toCopy);

	StableFluid2D& operator =(const StableFluid2D& toCopy);

	void step(double dt);

	void getvalue(double i, double j, double& ux, double& uy) const;
	vec2d linearSamp(vec2d& pos, double*** v) const;
	vec2d linearSamp(double x, double y, double*** v) const;

	double linearSamp(vec2d& pos, Vector& v) const;
	double linearSamp(double x, double y, Vector& v) const;

	void setVisc(double viscosity_);
	void setDensity(double density_);
	void setSize(int width_, int height_);
	void setRes(int resX_, int resY_);

protected:
	void initGrids();
	void clearGrids();
	void initPoisson();
	void clearPoisson();
public:
	int makeIndex(int x, int y) const;
	void splitIndex(int index, int& x, int& y);
	void setBoundary();
	void addForcesToVelocity(double dt);
	void advectVelocity(double dt);
	void projectVelocity(double dt);

	SparseMatrix laplacian;
	SparseMatrix preconditioner;

	Vector vel;
	Vector	pressure, tempPressure;

	std::vector<Vector>	r;
	std::vector<Vector>	z;
	Vector	p;

	Vector	ones;

	double*** velocities;
	double*** tempVelocities;

	double viscosity;
	double density;
	double width;
	double height;
	int	resX;
	int resY;
};
