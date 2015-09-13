#include "stdafx.h"
#include "LevelSetGrid3D.h"

LevelSetGrid3D::LevelSetGrid3D(int xResolution, int yResolution, int zResolution)
{
	size[0] = xResolution + 2 * bufferSize;
	size[1] = yResolution + 2 * bufferSize;
	size[2] = zResolution + 2 * bufferSize;

	gridLength = size[0] * size[1] * size[2];
	grid = new double[gridLength];
	std::fill(grid, grid + gridLength, 0);
}

LevelSetGrid3D::~LevelSetGrid3D()
{
	delete[] grid;
}

double& LevelSetGrid3D::get(int x, int y, int z)
{
	x += bufferSize;
	y += bufferSize;
	z += bufferSize;
	return grid[getindex(x, y, z)];
}

void LevelSetGrid3D::UpdateBorders()
{
	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{
			for (int k = 0; k < size[2]; k++)
			{
				int index[3] = { i,j,k };
				double dist[3] = { 0,0,0 };
				for (int m = 0; m < 3; m++)
				{
					if (index[m] < bufferSize)
					{
						dist[m] = abs(index[m] - bufferSize);
						index[m] = bufferSize;
					}
					else if (size[m] - index[m] <= bufferSize)
					{
						dist[m] = abs(size[m] - 1 - index[m] - bufferSize);
						index[m] = size[m] - 1 - bufferSize;
					}
				}
				double distance = sqrt(dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2]);
				double& phi = grid[getindex(i, j, k)];
				phi = grid[getindex(index[0], index[1], index[2])];

				//added this so that phi values beyond the border are more accurate
				if (phi < 0)
				{
					phi -= distance;
				}
				else
				{
					phi += distance;
				}
			}
		}
	}
}

int LevelSetGrid3D::getindex(int x, int y, int z)
{
	return x * size[1] * size[2] + y * size[2] + z;
}