#include "stdafx.h"
#include "LevelSetFMGrid2D.h"
#include "math.h"

LevelSetFMGrid2D::LevelSetFMGrid2D(int xResolution, int yResolution)
{
	size[0] = xResolution;
	size[1] = yResolution;

	gridLength = size[0] * size[1];
	FMHeap.reserve(gridLength);
	ClosePoints.reserve(gridLength);
	grid = new FMContainer[gridLength];
}

LevelSetFMGrid2D::~LevelSetFMGrid2D()
{
	delete[] grid;
}

void LevelSetFMGrid2D::AddClose(int index)
{
	if (grid[index].HeapPosition == -1)
	{
		AddToHeap(index);
	}
	else
	{
		UpdateHeap(index);
	}
}

void LevelSetFMGrid2D::AddToHeap(int index) 
{
	pushed.push_back(index);

	FMHeap.push_back(index);

	grid[index].HeapPosition = heapSize;
	int j, i = heapSize;
	for (i; i > 0; i = j) 
	{
		j = (i - 1) / 2;
		if (grid[(FMHeap[i])].value < grid[(FMHeap[j])].value) 
		{
			FMHeap[i] = FMHeap[j];
			grid[(FMHeap[j])].HeapPosition = i;
			FMHeap[j] = index;
			grid[index].HeapPosition = j;
		}
		else
		{
			break;
		}
	}
	heapSize++;
}

void LevelSetFMGrid2D::UpdateHeap(int index) 
{
	int j, i = grid[index].HeapPosition;
	for (i; i > 0; i = j) 
	{
		j = (i - 1) / 2;
		if (grid[(FMHeap[i])].value < grid[(FMHeap[j])].value) 
		{
			FMHeap[i] = FMHeap[j];
			grid[(FMHeap[j])].HeapPosition = i;
			FMHeap[j] = index;
			grid[index].HeapPosition = j;
		}
		else
		{
			break;
		}
	}
}

int LevelSetFMGrid2D::PopHeap()
{
	if (heapSize == 0)
	{
		return -1;
	}

	int j, index = FMHeap[0];
	grid[index].DoneFlag = 1;
	heapSize--;
	FMHeap[0] = FMHeap[heapSize];
	grid[FMHeap[heapSize]].HeapPosition = 0;
	for (int i = 0; i < (heapSize - 1); i = j) 
	{
		int lc = 2 * i + 1;
		int rc = 2 * i + 2;
		double current = grid[(FMHeap[i])].value;
		double lv, rv;
		if (lc < heapSize) 
		{
			lv = grid[(FMHeap[lc])].value;
			if (rc < heapSize) 
			{
				rv = grid[(FMHeap[rc])].value;
				if (lv > rv) 
				{
					lc = rc;
					lv = rv;
				}
			}
			if (current > lv) 
			{
				FMHeap[i] = FMHeap[lc];
				grid[FMHeap[i]].HeapPosition = i;
				FMHeap[lc] = FMHeap[heapSize];
				grid[FMHeap[heapSize]].HeapPosition = lc;
				j = lc;
			}
			else
			{
				break;
			}
		}
		else
		{
			break;
		}
	}
	FMHeap.pop_back();
	popped.push_back(index);
	return index;
}

int LevelSetFMGrid2D::getindex(int x, int y) 
{
	return x * size[1] + y;
}

double LevelSetFMGrid2D::Get(int x, int y)
{
	return grid[getindex(x, y)].value;
}

void LevelSetFMGrid2D::Reinitialize()
{
	pushed.clear();
	popped.clear();
	heapSize = 0;
	closeSize = 0;
	FMHeap.clear();
	ClosePoints.clear();
	Initialize();
	InitHeap();
	FastMarch();

	return;
	std::sort(pushed.begin(), pushed.end());
	std::sort(popped.begin(), popped.end());
	if (pushed.size() != popped.size())
	{
		std::cout << "bad heap size" << std::endl;
	}
	std::vector<int> diff1, diff2;
	std::set_difference(pushed.begin(), pushed.end(), popped.begin(), popped.end(), std::back_inserter(diff1));
	std::set_difference(popped.begin(), popped.end(), pushed.begin(), pushed.end(), std::back_inserter(diff2));
	if (diff1.size() > 0 || diff2.size() > 0)
	{
		std::cout << "bad heap contents" << std::endl;
	}
	std::cout << "pushed nodes: " << pushed.size() << std::endl;
	int done = 0;
	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{
			if (grid[getindex(i, j)].DoneFlag == 1)
			{
				done++;
			}
		}
	}
	std::cout << "done: " << done << std::endl;
}

void LevelSetFMGrid2D::getxyz(int index, int& x, int& y) 
{
	x = index / size[1];
	y = index % size[1];
}

void LevelSetFMGrid2D::FindPhi(int x, int y)
{

	if (x == 15 && y == 59)
	{
		int p = 5;
	}

	double phiX = 0;
	double phiY = 0;
	int a = 0;
	bool flagX = 0;
	bool flagY = 0;
	//Find The phiS
	CheckFront(phiX, a, flagX, x + 1, size[0], x + 1, y);
	CheckBehind(phiX, a, flagX, x - 1, x - 1, y);
	CheckFront(phiY, a, flagY, y + 1, size[1], x, y + 1);
	CheckBehind(phiY, a, flagY, y - 1, x, y - 1);

	if (a == 2)
	{		
		if (phiX >= phiY)
		{
			CheckMax2(a, phiX, phiY);
		}
		else
		{
			CheckMax2(a, phiY, phiX);
		}
	}

	if (a == 0)
	{
		throw std::exception();
	}

	double b = phiX + phiY;
	double quotient = square(b) - double(a) * (square(phiX) + square(phiY) - 1.0);
	if (quotient < 0)
	{
		std::cout << "0x ";
	}
	else
	{
		//cout << "1 ";
		double phi = b + std::sqrt(quotient);
		phi /= double(a);
		int index = getindex(x, y);

		if (phiX != 0 && abs(phi - phiX) > 1.001 || phiY != 0 && abs(phi - phiY) > 1.001)
		{
			std::cout << "bad phi " << phi << " " << phiX << " " << phiY << std::endl;
		}
		//if(phi < 2)
		//cout << x << " " << y << " " << z << " " << a << " " << 
		//	phiX << " " << phiY << " " << phiZ << " " << 
		//	grid[index].value << " " << phi << endl;
		grid[index].value = phi;
		AddClose(index);
	}
}

void LevelSetFMGrid2D::CheckMax2(int& a, double& phi1, double phi2) 
{
	if (square(phi1 - phi2) > 1) 
	{
		phi1 = 0;
		a = 1;
	}
}

void LevelSetFMGrid2D::CheckFront(double& phi, int& a, bool& flag, int check, int size, int x, int y) 
{
	int index = getindex(x, y);
	if (check < size) 
	{
		if (grid[index].DoneFlag == 1)
		{
			phi = grid[index].value;
			flag = 1;
			a++;
		}
	}
}

void LevelSetFMGrid2D::CheckBehind(double& phi, int& a, bool& flag, int check, int x, int y) 
{
	int index = getindex(x, y);
	if (check >= 0)
	{
		if (grid[index].DoneFlag == 1)
		{
			if (!flag) 
			{
				phi = grid[index].value;
				a++;
				flag = 1;
			}
			else
			{
				phi = std::min(grid[index].value, phi);
			}
		}
	}
}

void LevelSetFMGrid2D::FastMarch() 
{
	int x, y;
	for (int index = PopHeap(); index != -1; index = PopHeap()) 
	{
		getxyz(index, x, y);
		CheckPosExtent(x + 1, size[0], x + 1, y);
		CheckNegExtent(x - 1, x - 1, y);
		CheckPosExtent(y + 1, size[1], x, y + 1);
		CheckNegExtent(y - 1, x, y - 1);
	}
	if (heapSize != 0)
	{
		std::cout << "bad heap!!!" << std::endl;
	}
}

void LevelSetFMGrid2D::CheckPosExtent(int index, int size, int x, int y) 
{
	if (index < size) 
	{
		if (grid[getindex(x, y)].DoneFlag == 0) 
		{
			FindPhi(x, y);
		}
	}
}

void LevelSetFMGrid2D::CheckNegExtent(int index, int x, int y)
{
	if (index >= 0) 
	{
		if (grid[getindex(x, y)].DoneFlag == 0) 
		{
			FindPhi(x, y);
		}
	}
}

void LevelSetFMGrid2D::Set(int x, int y, double value) 
{
	int index = getindex(x, y);
	grid[index].value = value;
	grid[index].HeapPosition = -1;
	if (value < 0)
		grid[index].DoneFlag = -1;
	else
		grid[index].DoneFlag = 0;
}

void LevelSetFMGrid2D::InitHeap() 
{
	for (int i = 0; i < closeSize; i++) 
	{
		if (grid[ClosePoints[i]].HeapPosition == -1 && grid[ClosePoints[i]].DoneFlag == 0) 
		{
			//if(grid[ClosePoints[i]].HeapPosition == -1) {
			int x, y;
			getxyz(ClosePoints[i], x, y);
			FindPhi(x, y);
		}
	}
}

void LevelSetFMGrid2D::CloseBehind(int index, int x, int y) 
{
	if (index >= 0) 
	{
		//if(grid[getindex(x,y,z)].value >= 0) {
		if (grid[getindex(x, y)].DoneFlag == 0) 
		{
			//Add index to list for initialization of close band
			ClosePoints.push_back(getindex(x, y));
			closeSize++;
		}
	}
}

void LevelSetFMGrid2D::CloseFront(int size, int index, int x, int y) 
{
	if (index < size) 
	{
		//if(grid[getindex(x,y,z)].value >= 0) {
		if (grid[getindex(x, y)].DoneFlag == 0) 
		{
			//Add index to list for initialization of close band
			ClosePoints.push_back(getindex(x, y));
			closeSize++;
		}
	}
}

void LevelSetFMGrid2D::Initialize()
{
	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1] - 1; j++)
		{
			int cindex = getindex(i, j);
			int nindex = getindex(i, j + 1);
			double current = grid[cindex].value;
			double next = grid[nindex].value;

			if (current * next <= 0) 
			{
				if (abs(current) - abs(next) > 1.001)
				{
			//		std::cout << "too far:" << (abs(current) - abs(next)) << " " << i << " " << j << std::endl;
				}
				if (current >= 0) 
				{
					grid[cindex].DoneFlag = 1;
					CloseBehind(j - 1, i, j - 1);
				}
				else 
				{
					grid[nindex].DoneFlag = 1;
					CloseFront(size[1], j + 2, i, j + 2);
				}
			}
		}		
	}

	for (int j = 0; j < size[1]; j++)
	{
		for (int i = 0; i < size[0] - 1; i++)
		{
			int cindex = getindex(i, j);
			int nindex = getindex(i + 1, j);
			double current = grid[cindex].value;
			double next = grid[nindex].value;

			if (current * next <= 0) 
			{
				if (current >= 0) 
				{
					grid[cindex].DoneFlag = 1;
					CloseBehind(i - 1, i - 1, j);
				}
				else 
				{
					grid[nindex].DoneFlag = 1;
					CloseFront(size[0], i + 2, i + 2, j);
				}
			}
		}
	}
}