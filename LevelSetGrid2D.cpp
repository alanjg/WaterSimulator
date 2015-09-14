#include "stdafx.h"
#include "LevelSetGrid2D.h"

LevelSetGrid2D::LevelSetGrid2D(int xResolution, int yResolution)
{
	size[0] = xResolution + 2 * bufferSize;
	size[1] = yResolution + 2 * bufferSize;

	gridLength = size[0] * size[1];
	grid = new double[gridLength];
	std::fill(grid, grid + gridLength, 0);
}

LevelSetGrid2D::~LevelSetGrid2D()
{
	delete[] grid;
}

double& LevelSetGrid2D::get(int x, int y)
{
	x += bufferSize;
	y += bufferSize;
	if (x < 0 || x >= size[0] || y < 0 || y >= size[1]) throw std::exception();
	return grid[getindex(x, y)];
}

void LevelSetGrid2D::UpdateBorders()
{
	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{
			int index[2] = { i, j };
			double dist[2] = { 0, 0 };
			for (int m = 0; m < 2; m++)
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
			double distance = sqrt(dist[0] * dist[0] + dist[1] * dist[1]);
			double& phi = grid[getindex(i, j)];
			phi = grid[getindex(index[0], index[1])];

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

int LevelSetGrid2D::getindex(int x, int y)
{
	return x * size[1] + y;
}