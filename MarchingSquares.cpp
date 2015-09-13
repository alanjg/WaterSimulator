#include "stdafx.h"
#include "MarchingSquares.h"

#include <GL/glut.h>

std::vector<MarchingSquaresPattern> MarchingSquares::pattern;
int MarchingSquares::xOffsets[4];
int MarchingSquares::yOffsets[4];


MarchingSquaresPatternEntry::MarchingSquaresPatternEntry(int s1, int s2, int e1, int e2)
{
	indexes[0][0] = s1;
	indexes[0][1] = s2;
	indexes[1][0] = e1;
	indexes[1][1] = e2;
}

MarchingSquares::MarchingSquares()
{
}

MarchingSquares::~MarchingSquares()
{
}

void MarchingSquares::Update()
{
	lines.clear();
	int xDiv = 24, yDiv = 24;
	double xBegin = -1;
	double yBegin = -1;
	double xEnd = 1;
	double yEnd = 1;
	double xDelta = (xEnd - xBegin) / xDiv;
	double yDelta = (yEnd - yBegin) / yDiv;

	for (int i = 0; i < xDiv; i++)
	{		
		for (int j = 0; j < yDiv; j++)
		{
			vec2d points[4];
			points[0] = vec2d(xBegin + xDelta * i, yBegin + yDelta * (j + 1));
			points[1] = vec2d(xBegin + xDelta * (i + 1), yBegin + yDelta * (j + 1));
			points[2] = vec2d(xBegin + xDelta * (i + 1), yBegin + yDelta * j);
			points[3] = vec2d(xBegin + xDelta * i, yBegin + yDelta * j);
			
			double phi[4];
			phi[0] = iFunc(points[0]);
			phi[1] = iFunc(points[1]);
			phi[2] = iFunc(points[2]);
			phi[3] = iFunc(points[3]);
			
			int index = 0;
			for (int k = 0; k < 4; k++)
			{
				if (phi[k] >= 0) index += 1 << k;
			}
			
			for each(auto entry in pattern[index].entries)
			{
				line2d line;
				for (int k = 0; k < 2; k++)
				{
					double length = abs(phi[entry.indexes[k][0]]) + abs(phi[entry.indexes[k][1]]);
					double lerp = abs(phi[entry.indexes[k][0]]) / length;
						
					line.v[k] = points[entry.indexes[k][0]] * (1-lerp) + points[entry.indexes[k][1]] * lerp;
				}
				
				lines.push_back(line);
			}
		}
	}
}


void MarchingSquares::setImplicitFunction(implicitFuncType2D func)
{
	iFunc = func;
}


void MarchingSquares::Draw()
{
	float white[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
	glEnable(GL_LIGHTING);
	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, white);

	if (lines.size() > 0)
	{
		glPushAttrib(GL_NORMAL_ARRAY);
		glDisable(GL_NORMAL_ARRAY);
		glVertexPointer(2, GL_DOUBLE, 0, &lines[0].v[0]);
		glDrawArrays(GL_LINES, 0, int(lines.size()*2));
		glPopAttrib();
	}	
}

void MarchingSquares::createTable()
{
	pattern.resize(16);
	// 0 has no entries
	pattern[1].entries.push_back(MarchingSquaresPatternEntry(0, 3, 0, 1));
	pattern[2].entries.push_back(MarchingSquaresPatternEntry(0, 1, 1, 2));
	pattern[3].entries.push_back(MarchingSquaresPatternEntry(0, 3, 1, 2));
	pattern[4].entries.push_back(MarchingSquaresPatternEntry(3, 2, 1, 2));
	pattern[5].entries.push_back(MarchingSquaresPatternEntry(0, 3, 3, 2));
	pattern[5].entries.push_back(MarchingSquaresPatternEntry(0, 1, 1, 2));
	pattern[6].entries.push_back(MarchingSquaresPatternEntry(0, 1, 3, 2));
	pattern[7].entries.push_back(MarchingSquaresPatternEntry(0, 3, 3, 2));
	pattern[8].entries.push_back(MarchingSquaresPatternEntry(0, 3, 2, 3));
	pattern[9].entries.push_back(MarchingSquaresPatternEntry(0, 1, 2, 3));
	pattern[10].entries.push_back(MarchingSquaresPatternEntry(0, 1, 0, 3));
	pattern[10].entries.push_back(MarchingSquaresPatternEntry(2, 3, 1, 2));
	pattern[11].entries.push_back(MarchingSquaresPatternEntry(2, 3, 1, 2));
	pattern[12].entries.push_back(MarchingSquaresPatternEntry(0, 3, 1, 2));
	pattern[13].entries.push_back(MarchingSquaresPatternEntry(0, 1, 1, 2));
	pattern[14].entries.push_back(MarchingSquaresPatternEntry(0, 1, 0, 3));
	// 15 has no entries
}