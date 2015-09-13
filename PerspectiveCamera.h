#pragma once

#include "vec.h"

class PerspectiveCamera
{
public:
	PerspectiveCamera(float x, float y, float z);
	void MoveForward(float dist);
	void MoveBackward(float dist);
	void MoveLeft(float dist);
	void MoveRight(float dist);
	void MoveUp(float dist);
	void MoveDown(float dist);

	void RotateY(float degrees);
	void RotateX(float degrees);

	vec3d GetDir();
	vec3d GetPos();

	void Draw();
private:
	void CalcDir();
	float theta, phi;
	vec3d pos;
	vec3d dir, up, right;
	bool dirty;

};