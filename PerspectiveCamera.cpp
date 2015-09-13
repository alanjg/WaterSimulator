#include "stdafx.h"
#include "PerspectiveCamera.h"
#include <gl/glut.h>
#include "math.h"

using std::max;
using std::min;

PerspectiveCamera::PerspectiveCamera(float x, float y, float z) :
	theta(90), phi(0), dirty(false), dir(0, 0, -1), pos(x, y, z)
{
}

void PerspectiveCamera::MoveForward(float dist)
{
	CalcDir();
	pos += dir*dist;
}
void PerspectiveCamera::MoveBackward(float dist)
{
	CalcDir();
	pos -= dir*dist;
}
void PerspectiveCamera::MoveLeft(float dist)
{
	CalcDir();
	pos -= right*dist;
}
void PerspectiveCamera::MoveRight(float dist)
{
	CalcDir();
	pos += right*dist;
}
void PerspectiveCamera::MoveUp(float dist)
{
	CalcDir();
	pos += up*dist;
}
void PerspectiveCamera::MoveDown(float dist)
{
	CalcDir();
	pos -= up*dist;
}

void PerspectiveCamera::RotateY(float degrees)
{
	dirty = true;
	theta += degrees;
}
void PerspectiveCamera::RotateX(float degrees)
{
	dirty = true;
	phi += degrees;
	phi = max(-85.0f, phi);
	phi = min(85.0f, phi);
}


vec3d PerspectiveCamera::GetDir()
{
	CalcDir();
	return dir;
}

void PerspectiveCamera::Draw()
{
	CalcDir();
	vec3d at = pos + dir * 100;
	gluLookAt(pos[0], pos[1], pos[2], at[0], at[1], at[2], up[0], up[1], up[2]);
}

vec3d PerspectiveCamera::GetPos()
{
	CalcDir();
	return pos;
}

void PerspectiveCamera::CalcDir()
{
	if (!dirty)
		return;
	float rphi = Deg2Rad(phi);
	float rtheta = Deg2Rad(theta);
	float r = cos(rphi);

	dir[0] = r*cos(rtheta);
	dir[1] = sin(rphi);
	dir[2] = r*sin(rtheta);

	vec3d vup(0, 1, 0);
	right = dir.Cross(vup);

	dir.Normalize();
	right.Normalize();

	up = right.Cross(dir);
	up.Normalize();

	dirty = false;

}