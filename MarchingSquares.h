#pragma once

#include "vec.h"
typedef double(*implicitFuncType2D)(const vec2d& vec);

struct MarchingSquaresPatternEntry
{
	MarchingSquaresPatternEntry(int s1, int s2, int e1, int e2);
	int indexes[2][2];
};

struct MarchingSquaresPattern
{
	std::vector<MarchingSquaresPatternEntry> entries;
};

class MarchingSquares
{
	static int xOffsets[4], yOffsets[4];
	implicitFuncType2D iFunc;
	std::vector<line2d> lines;
	static std::vector<MarchingSquaresPattern> pattern;
public:
	MarchingSquares();
	~MarchingSquares();

	void Draw();
	virtual void setImplicitFunction(implicitFuncType2D func);
	static void createTable();

	void Update();
};

