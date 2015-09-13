#include "stdafx.h"
#include "Particle2D.h"

Particle2D::Particle2D(double px, double py, double phi)
{
	sign = phi < 0 ? -1 : 1;
	radius = sign * phi;
	if (radius > PARTICLE_RADIUS_MAX)
		radius = PARTICLE_RADIUS_MAX;
	else if (radius < PARTICLE_RADIUS_MIN)
		radius = PARTICLE_RADIUS_MIN;

	x = px;
	y = py;
}

double Particle2D::phi(double px, double py)
{
	return sign * (radius - std::sqrt((px - x)*(px - x) + (py - y)*(py - y)));
}

void Particle2D::Update(double ux0, double uy0, double ux1, double uy1, double ux2, double uy2, double timestep)
{
	/*
	//do an euler integration
	x += ux * timestep;
	y += uy * timestep;
	*/
	//do an RK3 integration
	double x1 = x + ux0 * timestep;
	double y1 = y + uy0 * timestep;

	double x2 = x1 + ux2 * timestep;
	double y2 = y1 + uy2 * timestep;

	double x05 = 0.75 * x + 0.25 * x2;
	double y05 = 0.75 * y + 0.25 * y2;

	double x15 = x05 + ux1 * timestep;
	double y15 = y05 + uy1 * timestep;

	x = x / 3.0 + x15 * 2 / 3.0;
	y = y / 3.0 + y15 * 2 / 3.0;
}

void Particle2D::GetPosition(double& px, double& py)
{
	px = x;
	py = y;
}

double Particle2D::Sign()
{
	return sign;
}

double Particle2D::Radius(double phi)
{
	radius = fabs(phi);
	if (radius > PARTICLE_RADIUS_MAX)
		radius = PARTICLE_RADIUS_MAX;
	else if (radius < PARTICLE_RADIUS_MIN)
		radius = PARTICLE_RADIUS_MIN;
	return radius;
}

