#pragma once

class LevelSetGrid2D
{
public:
	LevelSetGrid2D(int xResolution, int yResolution);
	~LevelSetGrid2D();
	double& get(int x, int y);
	void UpdateBorders();
	void Write(std::ostream& out)
	{
		out << (size[0] - 2 * bufferSize) << " " << (size[1] - 2 * bufferSize) << std::endl;
		for (int i = bufferSize; i<size[0] - bufferSize; i++)
		{
			for (int j = bufferSize; j<size[1] - bufferSize; j++)
			{
				out << grid[getindex(i, j)] << " ";
			}
		}
	}
private:
	static const int bufferSize = 3;
	int size[2];
	int gridLength;
	double* grid;

	int getindex(int x, int y);
};