#include "stdafx.h"

#include "MarchingCubes.h"
#include "vec.h"
#include <GL/glut.h>

using namespace std;
// these are the original 0-14 cube patterns from the marching cubes
// algorithm
patInfo patInfo::patInitTable[] = {
	patInfo(0x0),
	patInfo(0x1, triangle(0, 3, 4)), // 1
	patInfo(0x3, triangle(5, 3, 4), triangle(5, 1, 3)), // 2
	patInfo(0x21, triangle(0, 3, 4), triangle(5, 9, 10)),
	patInfo(0x81, triangle(0, 3, 4), triangle(6, 10, 11)),
	patInfo(0xe, triangle(0, 5, 3), triangle(3, 5, 7), triangle(5, 6, 7)),
	patInfo(0x83, triangle(5, 3, 4), triangle(1, 3, 5), triangle(6, 10, 11)),
	patInfo(0x92, triangle(4, 8, 9), triangle(0, 5, 1), triangle(6, 10, 11)),
	patInfo(0xf, triangle(4, 5, 6), triangle(4, 6, 7)),
	patInfo(0x4d, triangle(4, 11, 8), triangle(0, 11, 4), triangle(0, 6, 11),
		triangle(0, 1, 6)),
	patInfo(0x99, triangle(0, 3, 9), triangle(9, 3, 8), triangle(2, 1, 10),
		triangle(2, 10, 11)),
	patInfo(0x8d, triangle(0, 1, 10), triangle(0, 10, 7), triangle(0, 7, 4),
		triangle(7, 10, 11)),
	patInfo(0x1e, triangle(4, 8, 9), triangle(3, 0, 5), triangle(3, 5, 7),
		triangle(5, 6, 7)),
	patInfo(0x69, triangle(0, 3, 4), triangle(1, 6, 2), triangle(5, 9, 10),
		triangle(7, 11, 8)),
	patInfo(0x4e, triangle(0, 8, 3), triangle(0, 6, 8), triangle(0, 5, 6),
		triangle(6, 11, 8))
};


int patInfo::edgeToCorner[12][2] = {
	{ 0, 1 },
	{ 1, 3 },
	{ 2, 3 },
	{ 0, 2 },
	{ 0, 4 },
	{ 1, 5 },
	{ 3, 7 },
	{ 2, 6 },
	{ 4, 6 },
	{ 4, 5 },
	{ 5, 7 },
	{ 6, 7 },
};

// given edge, axis, will translate newEdge = edgeTable[edge][axis] 
vec3d patInfo::edgeTable[] = {
	vec3d(2, 1, 5),
	vec3d(6, 2, 10),
	vec3d(11, 3, 6),
	vec3d(7, 0, 1),
	vec3d(3, 5, 0),
	vec3d(1, 6, 9),
	vec3d(10, 7, 11),
	vec3d(8, 4, 2),
	vec3d(4, 9, 3),
	vec3d(0, 10, 4),
	vec3d(5, 11, 8),
	vec3d(9, 8, 7)
};

patInfo MarchingCubes::cubeLut[0x100];

patInfo::patInfo(unsigned short p) {

	vertPattern = p;
}

patInfo::patInfo(unsigned short p, const triangle &a) {

	vertPattern = p;

	triList[0] = a;
}

patInfo::patInfo(unsigned short p, const triangle &a,
	const triangle &b) {

	vertPattern = p;

	triList[0] = a;
	triList[1] = b;
}

patInfo::patInfo(unsigned short p, const triangle &a,
	const triangle &b, const triangle &c) {

	vertPattern = p;

	triList[0] = a;
	triList[1] = b;
	triList[2] = c;
}

patInfo::patInfo(unsigned short p, const triangle &a,
	const triangle &b, const triangle &c, const triangle &d) {

	vertPattern = p;

	triList[0] = a;
	triList[1] = b;
	triList[2] = c;
	triList[3] = d;

}

void patInfo::invert() {

	int i;

	vertPattern ^= 0xff;

	for (i = 0; i < 4; i++) {
		if (!triList[i].valid) {
			break;
		}
		triList[i].invert();
	}
}

void patInfo::getVertex(const vec3d &center, const vec3d &size, int vertex, vec3d &pos)
{
	double dx = vertex % 2 ? 1 : -1;
	double dy = (vertex % 4) / 2 ? -1 : 1;
	double dz = vertex / 4 ? 1 : -1;

	vec3d d(dx*size[0], dz*size[2], dy*size[1]);
	pos = center + d;
	/*
	switch(vertex) {
	case 0:
	pos.x() = center.x()-size.x();
	pos.y() = center.z()-size.z();
	pos.z() = center.y()+size.y();
	break;
	case 1:
	pos.x() = center.x()+size.x();
	pos.y() = center.z()-size.z();
	pos.z() = center.y()+size.y();
	break;
	case 2:
	pos.x() = center.x()-size.x();
	pos.y() = center.z()-size.z();
	pos.z() = center.y()-size.y();
	break;
	case 3:
	pos.x() = center.x()+size.x();
	pos.y() = center.z()-size.z();
	pos.z() = center.y()-size.y();
	break;
	case 4:
	pos.x() = center.x()-size.x();
	pos.y() = center.z()+size.z();
	pos.z() = center.y()+size.y();
	break;
	case 5:
	pos.x() = center.x()+size.x();
	pos.y() = center.z()+size.z();
	pos.z() = center.y()+size.y();
	break;
	case 6:
	pos.x() = center.x()-size.x();
	pos.y() = center.z()+size.z();
	pos.z() = center.y()-size.y();
	break;
	case 7:
	pos.x() = center.x()+size.x();
	pos.y() = center.z()+size.z();
	pos.z() = center.y()-size.y();
	break;
	}
	*/
}


MarchingCubes::MarchingCubes()
{
	level = 7;

	newway = false;
	int i;

	testCube = patInfo::patInitTable[0];

	for (i = 0; i <= 14; i++) {
		patInfo::patInitTable[i].origIndex = i;
	}
	iFunc = 0;

}

MarchingCubes::~MarchingCubes() {

}

void MarchingCubes::setImplicitFunction(implicitFuncType func)
{
	iFunc = func;
}

void MarchingCubes::Draw()
{
	float white[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
	glEnable(GL_LIGHTING);
	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, white);

	if (vertices.size() > 0)
	{
		glVertexPointer(3, GL_DOUBLE, 0, &vertices[0]);
		glNormalPointer(GL_DOUBLE, 0, &normals[0]);
		glDrawArrays(GL_TRIANGLES, 0, int(vertices.size()));
	}
}

void MarchingCubes::GetVertices(std::vector<vec3d>& _vertices)
{
	for each(auto v in vertices)
	{
		_vertices.push_back(v);
	}
}


void patInfo::rotPattern(int axis) {

	int i;
	unsigned short mask;
	bool vert[8], tmp;

	mask = 0x1;
	for (i = 0; i < 8; i++) {
		//vert[i] = (bool)(mask & pat);
		vert[i] = (mask & vertPattern) ? true : false;

		mask = (mask << 1);
	}

	switch (axis) {
	case X:
		tmp = vert[1];
		vert[1] = vert[5];
		vert[5] = vert[7];
		vert[7] = vert[3];
		vert[3] = tmp;
		tmp = vert[0];
		vert[0] = vert[4];
		vert[4] = vert[6];
		vert[6] = vert[2];
		vert[2] = tmp;
		break;
	case Y:
		tmp = vert[1];
		vert[1] = vert[0];
		vert[0] = vert[2];
		vert[2] = vert[3];
		vert[3] = tmp;
		tmp = vert[5];
		vert[5] = vert[4];
		vert[4] = vert[6];
		vert[6] = vert[7];
		vert[7] = tmp;
		break;
	case Z:
		tmp = vert[1];
		vert[1] = vert[0];
		vert[0] = vert[4];
		vert[4] = vert[5];
		vert[5] = tmp;
		tmp = vert[3];
		vert[3] = vert[2];
		vert[2] = vert[6];
		vert[6] = vert[7];
		vert[7] = tmp;
		break;
	}

	mask = 0x1;
	//pat = 0;
	vertPattern = 0;
	for (i = 0; i < 8; i++) {
		if (vert[i]) {
			//pat |= mask;
			vertPattern |= mask;
		}
		mask = (mask << 1);
	}
}

void patInfo::rotate(int axis) {

	int i, j;
	int vert;
	orderInfo oInfo;

	oInfo.pattern = vertPattern;
	oInfo.rotation = axis;
	rotOrder.push_back(oInfo);

	rotPattern(axis);

	for (i = 0; i < 4; i++) {
		if (!triList[i].valid) {
			break;
		}

		for (j = 0; j < 3; j++) {
			vert = triList[i].vert[j];
			triList[i].vert[j] = patInfo::edgeTable[vert][axis];
		}
	}

}

void invertList(std::set<patInfo> &patList) {

	std::set<patInfo>::iterator it;
	std::deque<patInfo> list;

	for (it = patList.begin(); it != patList.end(); it++) {
		list.push_back(*it);
	}

	for (unsigned int i = 0; i < list.size(); i++) {
		list[i].invert();
		patList.insert(list[i]);
	}

}
void invertList(std::set<unsigned short> &patList) {

	std::set<unsigned short>::iterator it;
	std::deque<unsigned char> list;


	for (it = patList.begin(); it != patList.end(); it++) {
		list.push_back(*it);
	}

	for (unsigned int i = 0; i < list.size(); i++) {
		list[i] ^= 0xff;
		patList.insert(list[i]);
	}

}

void makeList(std::set<unsigned short> &patList, unsigned short pat) {

	unsigned short P;

	patList.clear();
	P = pat;

	// X positive
	patList.insert(pat);
	//rotPattern(X, pat);
	patList.insert(pat);
	//rotPattern(X, pat);
	patList.insert(pat);
	//rotPattern(X, pat);
	patList.insert(pat);
	//rotPattern(X, pat);

	//rotPattern(Y, pat);
	//rotPattern(Y, pat);

	// X negative
	patList.insert(pat);
	//rotPattern(X, pat);
	patList.insert(pat);
	//rotPattern(X, pat);
	patList.insert(pat);
	//rotPattern(X, pat);
	patList.insert(pat);
	//rotPattern(X, pat);

	pat = P;
	//rotPattern(Z, pat);

	// y positive
	patList.insert(pat);
	//rotPattern(X, pat);
	patList.insert(pat);
	//rotPattern(X, pat);
	patList.insert(pat);
	//rotPattern(X, pat);
	patList.insert(pat);
	//rotPattern(X, pat);

	//rotPattern(Y, pat);
	//rotPattern(Y, pat);

	// y negative
	patList.insert(pat);
	//rotPattern(X, pat);
	patList.insert(pat);
	//rotPattern(X, pat);
	patList.insert(pat);
	//rotPattern(X, pat);
	patList.insert(pat);
	//rotPattern(X, pat);

	pat = P;
	//rotPattern(Y, pat);
	//rotPattern(Y, pat);
	//rotPattern(Y, pat);

	// z positive
	patList.insert(pat);
	//rotPattern(X, pat);
	patList.insert(pat);
	//rotPattern(X, pat);
	patList.insert(pat);
	//rotPattern(X, pat);
	patList.insert(pat);
	//rotPattern(X, pat);

	//rotPattern(Y, pat);
	//rotPattern(Y, pat);

	// z negative
	patList.insert(pat);
	//rotPattern(X, pat);
	patList.insert(pat);
	//rotPattern(X, pat);
	patList.insert(pat);
	//rotPattern(X, pat);
	patList.insert(pat);
	//rotPattern(X, pat);

}

void patInfo::clearRotations() {
	rotOrder.clear();
}

void makeList(std::set<patInfo> &patList, patInfo &pat) {

	patInfo P;

	pat.clearRotations();
	patList.clear();
	P = pat;


	// X positive
	patList.insert(pat);
	pat.rotate(patInfo::X);
	patList.insert(pat);
	pat.rotate(patInfo::X);
	patList.insert(pat);
	pat.rotate(patInfo::X);
	patList.insert(pat);
	pat.rotate(patInfo::X);

	pat.rotate(patInfo::Y);
	pat.rotate(patInfo::Y);

	// patInfo::X negative
	patList.insert(pat);
	pat.rotate(patInfo::X);
	patList.insert(pat);
	pat.rotate(patInfo::X);
	patList.insert(pat);
	pat.rotate(patInfo::X);
	patList.insert(pat);
	pat.rotate(patInfo::X);

	pat = P;
	pat.rotate(patInfo::Z);

	// y positive
	patList.insert(pat);
	pat.rotate(patInfo::X);
	patList.insert(pat);
	pat.rotate(patInfo::X);
	patList.insert(pat);
	pat.rotate(patInfo::X);
	patList.insert(pat);
	pat.rotate(patInfo::X);

	pat.rotate(patInfo::Y);
	pat.rotate(patInfo::Y);

	// y negative
	patList.insert(pat);
	pat.rotate(patInfo::X);
	patList.insert(pat);
	pat.rotate(patInfo::X);
	patList.insert(pat);
	pat.rotate(patInfo::X);
	patList.insert(pat);
	pat.rotate(patInfo::X);

	pat = P;
	pat.rotate(patInfo::Y);
	pat.rotate(patInfo::Y);
	pat.rotate(patInfo::Y);

	// z positive
	patList.insert(pat);
	pat.rotate(patInfo::X);
	patList.insert(pat);
	pat.rotate(patInfo::X);
	patList.insert(pat);
	pat.rotate(patInfo::X);
	patList.insert(pat);
	pat.rotate(patInfo::X);

	pat.rotate(patInfo::Y);
	pat.rotate(patInfo::Y);

	// z negative
	patList.insert(pat);
	pat.rotate(patInfo::X);
	patList.insert(pat);
	pat.rotate(patInfo::X);
	patList.insert(pat);
	pat.rotate(patInfo::X);
	patList.insert(pat);
	pat.rotate(patInfo::X);

}

void MarchingCubes::createTable() {

	std::set<patInfo> cubeSet;
	std::set<patInfo>::iterator it;
	patInfo info;
	int i, total;

	total = 0;
	for (i = 0; i <= 14; i++) {
		info = patInfo::patInitTable[i];

		makeList(cubeSet, info);
		invertList(cubeSet);

		total += int(cubeSet.size());

		for (it = cubeSet.begin(); it != cubeSet.end(); it++) {
			//patInfo& entry = it.;

			(*it).origIndex = i;
			cubeLut[(*it).vertPattern] = *it;
		}
	}

}

static unsigned int bitIndex[8] = {
	0x1,
	0x2,
	0x4,
	0x8,
	0x10,
	0x20,
	0x40,
	0x80
};

void MarchingCubes::March(double lx, double ly, double lz, double rx, double ry, double rz)
{
	unsigned int index;
	vec3d positions[8];
	double values[8];

	double x[] = { lx,rx };
	double y[] = { ly,ry };
	double z[] = { lz,rz };

	for (int i = 0; i<8; i++)
	{
		positions[i].Set(x[i % 2], y[i / 4], z[1 - (i % 4) / 2]);
		values[i] = iFunc(positions[i]);
	}
	/*
	positions[0].Set(lx,ly,rz);
	values[0] = iFunc(positions[0]);
	positions[1].Set(rx,ly,rz);
	values[1] = iFunc(positions[1]);
	positions[2].Set(lx,ly,lz);
	values[2] = iFunc(positions[2]);
	positions[3].Set(rx,ly,lz);
	values[3] = iFunc(positions[3]);
	positions[4].Set(lx,ry,rz);
	values[4] = iFunc(positions[4]);
	positions[5].Set(rx,ry,rz);
	values[5] = iFunc(positions[5]);
	positions[6].Set(lx,ry,lz);
	values[6] = iFunc(positions[6]);
	positions[7].Set(rx,ry,lz);
	values[7] = iFunc(positions[7]);
	*/
	index = 0;


	// calc our index
	for (int l = 0; l < 8; l++)
		if (values[l] > 0)
			index |= bitIndex[l];

	// create our triangles
	for (int l = 0; l < 4; l++)
	{
		if (!cubeLut[index].triList[l].valid)
		{
			break;
		}
		vec3d verts[3], norms[3];

		// calc the three vertices
		for (int m = 0; m < 3; m++)
		{

			const double GRAD_EPSILON = 1e-5;
			const double GRAD_EPS_INV = 1.0 / GRAD_EPSILON;
			int edge = cubeLut[index].triList[l].vert[m];
			vec3d& aPos = positions[patInfo::edgeToCorner[edge][0]];
			vec3d& bPos = positions[patInfo::edgeToCorner[edge][1]];

			double aVal = values[patInfo::edgeToCorner[edge][0]];
			double bVal = values[patInfo::edgeToCorner[edge][1]];

			double t = fabs(aVal) / fabs(bVal - aVal);
			//double t = fabs(aVal)/(fabs(bVal) + fabs(aVal));
			verts[m] = aPos * (1 - t) + bPos * t;
			vec3d pos = verts[m];
			double dx = iFunc(vec3d(pos[0] + GRAD_EPSILON,
				pos[1],
				pos[2])) -
				iFunc(vec3d(pos[0] - GRAD_EPSILON,
					pos[1],
					pos[2]));
			double dy = iFunc(vec3d(pos[0],
				pos[1] + GRAD_EPSILON,
				pos[2])) -
				iFunc(vec3d(pos[0],
					pos[1] - GRAD_EPSILON,
					pos[2]));
			double dz = iFunc(vec3d(pos[0],
				pos[1],
				pos[2] + GRAD_EPSILON)) -
				iFunc(vec3d(pos[0],
					pos[1],
					pos[2] - GRAD_EPSILON));
			vec3d g(dx, dy, dz);
			g *= GRAD_EPS_INV * 0.5;
			g.Normalize();
			norms[m] = g;
		}

		vertices.push_back(verts[0]);
		normals.push_back(norms[0]);

		vertices.push_back(verts[2]);
		normals.push_back(norms[2]);

		vertices.push_back(verts[1]);
		normals.push_back(norms[1]);
	}
}

void MarchingCubes::recurse(double lx, double ly, double lz, double rx, double ry, double rz, int subdivisions)
{
	//	cout << "Recurse at "<<subdivisions<<endl;
	if (subdivisions == 0)
	{
		March(lx, ly, lz, rx, ry, rz);
		return;
	}

	double cx = (lx + rx) / 2;
	double cy = (ly + ry) / 2;
	double cz = (lz + rz) / 2;
	double ex = (rx - lx) / 2;
	double ey = (ry - ly) / 2;
	double ez = (rz - lz) / 2;

	double x[] = { lx,cx,rx };
	double y[] = { ly,cy,ry };
	double z[] = { lz,cz,rz };

	double ext = sqrt(ex*ex + ey*ey + ez*ez);
	vec3d val(cx, cy, cz);
	double phi = iFunc(val);
	bool isOuter = false;
	if (fabs(phi) > 2 * ext)
	{
		//empty cell
		isOuter = true;
		return;
	}

	bool sawPositive = false, sawNegative = false;
	for (int i = 0; i<3; i++)
		for (int j = 0; j<3; j++)
			for (int k = 0; k<3; k++)
			{
				double d = iFunc(vec3d(x[i], y[j], z[k]));
				if (d < 0) sawNegative = true;
				if (d > 0) sawPositive = true;
			}
	if (isOuter && sawNegative && sawPositive)
	{
		cout << "Error!" << endl;
	}
	for (int i = 0; i<2; i++)
		for (int j = 0; j<2; j++)
			for (int k = 0; k<2; k++)
				recurse(x[i], y[j], z[k], x[i + 1], y[j + 1], z[k + 1], subdivisions - 1);
}

void MarchingCubes::Update()
{
	vertices.resize(0);
	normals.resize(0);
	//	vertices.clear();
	//	normals.clear();

	recurse(-1, -1, -1, 1, 1, 1, level);
	return;
}
