#pragma once
const double PARTICLE_RADIUS_MIN = 0.1;
const double PARTICLE_RADIUS_MAX = 0.5;

class Particle3D
{
	double x, y, z;
	double sign;
	double radius;

public:
	Particle3D(double x, double y, double z, double phi);
	double phi(double x, double y, double z);
	double Sign();
	double Radius(double phi);
	void Update(double ux0, double uy0, double uz0, double ux1, double uy1, double uz1, double ux2, double uy2, double uz2, double timestep);
	void GetPosition(double& x, double& y, double& z);
};
