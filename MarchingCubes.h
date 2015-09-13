#pragma once

#ifdef WIN32
#pragma warning(disable:4251 4786)
#endif

#include "vec.h"

class triangle {

public:
	bool valid;
	int vert[3];

	triangle(int a, int b, int c) {
		vert[0] = a;
		vert[1] = b;
		vert[2] = c;
	}

	triangle() { valid = false; }

	void invert() {
		int tmp;
		tmp = vert[0];
		vert[0] = vert[1];
		vert[1] = tmp;
	}
};

class patInfo {
protected:

	typedef struct {
		unsigned short pattern;
		int rotation;
	} orderInfo;

	std::deque<orderInfo> rotOrder;

	static vec3d edgeTable[];

	static void getVertex(const vec3d &center, const vec3d &size,
		int vertex, vec3d &pos);
	void rotPattern(int axis);

public:

	static int edgeToCorner[12][2];
	static patInfo patInitTable[];

	mutable int origIndex;
	unsigned short vertPattern;
	triangle triList[4];

	void clearRotations();

	patInfo() { vertPattern = 0x0; }
	patInfo(unsigned short p);
	patInfo(unsigned short p, const triangle &a);
	patInfo(unsigned short p, const triangle &a,
		const triangle &b);
	patInfo(unsigned short p, const triangle &a,
		const triangle &b, const triangle &c);
	patInfo(unsigned short p, const triangle &a,
		const triangle &b, const triangle &c, const triangle &d);

	void clear() {
		triList[0].valid = false;
		triList[1].valid = false;
		triList[2].valid = false;
		triList[3].valid = false;
	}

	bool operator<(const patInfo &tr) const {
		return vertPattern<tr.vertPattern;
	}

	enum { X, Y, Z };

	void rotate(int axis);
	void invert();
};

typedef double(*implicitFuncType)(const vec3d& vec);

class MarchingCubes 
{

protected:



	implicitFuncType iFunc;

	patInfo testCube;

	// void quatRoot(double w, double x, double y, double z);

	static patInfo cubeLut[0x100];



	void recurse(double lx, double ly, double lz, double rx, double ry, double rz, int subdivisions);
	void March(double lx, double ly, double lz, double rx, double ry, double rz);

public:

	bool newway;
	int level;
	std::vector<vec3d> vertices;
	std::vector<vec3d> normals;

	MarchingCubes();
	~MarchingCubes();

	void Draw();
	void GetVertices(std::vector<vec3d>& vertices);

	virtual void setImplicitFunction(implicitFuncType func);

	static void createTable();

	virtual void Update();
};
