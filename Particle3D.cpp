#include "stdafx.h"
#include "Particle3D.h"

Particle3D::Particle3D(double px, double py, double pz, double phi)
{
	sign = phi < 0 ? -1 : 1;
	radius = sign * phi;
	if (radius > PARTICLE_RADIUS_MAX)
		radius = PARTICLE_RADIUS_MAX;
	else if (radius < PARTICLE_RADIUS_MIN)
		radius = PARTICLE_RADIUS_MIN;

	x = px;
	y = py;
	z = pz;
}

double Particle3D::phi(double px, double py, double pz)
{
	return sign * (radius - std::sqrt((px - x)*(px - x) + (py - y)*(py - y) + (pz - z)*(pz - z)));
}

void Particle3D::Update(double ux0, double uy0, double uz0, double ux1, double uy1, double uz1, double ux2, double uy2, double uz2, double timestep)
{
	/*
	//do an euler integration
	x += ux * timestep;
	y += uy * timestep;
	z += uz * timestep;
	*/
	//do an RK3 integration
	double x1 = x + ux0 * timestep;
	double y1 = y + uy0 * timestep;
	double z1 = z + uz0 * timestep;

	double x2 = x1 + ux2 * timestep;
	double y2 = y1 + uy2 * timestep;
	double z2 = z1 + uz2 * timestep;

	double x05 = 0.75 * x + 0.25 * x2;
	double y05 = 0.75 * y + 0.25 * y2;
	double z05 = 0.75 * z + 0.25 * z2;

	double x15 = x05 + ux1 * timestep;
	double y15 = y05 + uy1 * timestep;
	double z15 = z05 + uz1 * timestep;

	x = x / 3.0 + x15 * 2 / 3.0;
	y = y / 3.0 + y15 * 2 / 3.0;
	z = z / 3.0 + z15 * 2 / 3.0;
}

void Particle3D::GetPosition(double& px, double& py, double& pz)
{
	px = x;
	py = y;
	pz = z;
}

double Particle3D::Sign()
{
	return sign;
}

double Particle3D::Radius(double phi)
{
	radius = fabs(phi);
	if (radius > PARTICLE_RADIUS_MAX)
		radius = PARTICLE_RADIUS_MAX;
	else if (radius < PARTICLE_RADIUS_MIN)
		radius = PARTICLE_RADIUS_MIN;
	return radius;
}

