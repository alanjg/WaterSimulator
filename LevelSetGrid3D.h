#pragma once

class LevelSetGrid3D
{
public:
	LevelSetGrid3D(int xResolution, int yResolution, int zResolution);
	~LevelSetGrid3D();
	double& get(int x, int y, int z);
	void UpdateBorders();
	void Write(std::ostream& out)
	{
		out << (size[0] - 2 * bufferSize) << " " << (size[1] - 2 * bufferSize) << " " << (size[2] - 2 * bufferSize) << std::endl;
		for (int i = bufferSize; i<size[0] - bufferSize; i++)
		{
			for (int j = bufferSize; j<size[1] - bufferSize; j++)
			{
				for (int k = bufferSize; k<size[2] - bufferSize; k++)
				{
					out << grid[getindex(i, j, k)] << " ";
				}
			}
		}
	}
private:
	static const int bufferSize = 3;
	int size[3];
	int gridLength;
	double* grid;

	int getindex(int x, int y, int z);
};