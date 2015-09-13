#pragma once
#include "vec.h"
#include "Vector.h"
#include "SparseVector.h"
#include "SparseMatrix.h"

static int xOffset3[3] = { 1, 0, 0 };
static int yOffset3[3] = { 0, 1, 0 };
static int zOffset3[3] = { 0, 0, 1 };

class StableFluid3D
{
	const double GRAV_ACCEL = 9.8;
	const double PCG_EPS = 1e-10;
	const double PCG_MAXITER = 500;
public:
	StableFluid3D();
	~StableFluid3D();
	StableFluid3D(double width_, double height_, double depth_, int resX_, int resY_, int resZ_);
	StableFluid3D(const StableFluid3D& toCopy);

	StableFluid3D& operator =(const StableFluid3D& toCopy);

	void step(double dt);

	void getvalue(double i, double j, double k, double& ux, double& uy, double& uz) const;
	vec3d linearSamp(vec3d& pos, double**** v) const;
	vec3d linearSamp(double x, double y, double z, double**** v) const;

	double linearSamp(vec3d& pos, Vector& v) const;
	double linearSamp(double x, double y, double z, Vector& v) const;

	void setVisc(double viscosity_);
	void setDensity(double density_);
	void setSize(int width_, int height_, int depth_);
	void setRes(int resX_, int resY_, int resZ_);

protected:
	void initGrids();
	void clearGrids();
	void initPoisson();
	void clearPoisson();
public:
	int makeIndex(int x, int y, int z) const;
	void splitIndex(int index, int& x, int& y, int& z);
	void setBoundary();
	void addForcesToVelocity(double dt);
	void advectVelocity(double dt);
	void advectPressure(double dt);
	void projectVelocity(double dt);

	SparseMatrix laplacian;
	SparseMatrix preconditioner;

	Vector vel;
	Vector	pressure, tempPressure;

	std::vector<Vector>	r;
	std::vector<Vector>	z;
	Vector	p;

	Vector	ones;

	double**** velocities;
	double**** tempVelocities;

	double viscosity;
	double density;
	double width;
	double height;
	double depth;
	int	resX;
	int resY;
	int resZ;
};
