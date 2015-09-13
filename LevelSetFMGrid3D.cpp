#include "stdafx.h"
#include "LevelSetFMGrid3D.h"

LevelSetFMGrid3D::LevelSetFMGrid3D(int xResolution, int yResolution, int zResolution)
{
	size[0] = xResolution;
	size[1] = yResolution;
	size[2] = zResolution;

	gridLength = size[0] * size[1] * size[2];
	FMHeap.reserve(gridLength);
	ClosePoints.reserve(gridLength);
	grid = new FMContainer[gridLength];
}

LevelSetFMGrid3D::~LevelSetFMGrid3D()
{
	delete[] grid;
}

void LevelSetFMGrid3D::AddClose(int index) {
	if (grid[index].HeapPosition == -1)
		AddToHeap(index);
	else
		UpdateHeap(index);
}

void LevelSetFMGrid3D::AddToHeap(int index) {
	FMHeap.push_back(index);

	grid[index].HeapPosition = heapSize;
	int j, i = heapSize;
	for (i; i > 0; i = j) {
		j = (i - 1) / 2;
		if (grid[(FMHeap[i])].value < grid[(FMHeap[j])].value) {
			FMHeap[i] = FMHeap[j];
			grid[(FMHeap[j])].HeapPosition = i;
			FMHeap[j] = index;
			grid[index].HeapPosition = j;
		}
		else
			break;
	}
	heapSize++;
}

void LevelSetFMGrid3D::UpdateHeap(int index) {
	int j, i = grid[index].HeapPosition;
	for (i; i > 0; i = j) {
		j = (i - 1) / 2;
		if (grid[(FMHeap[i])].value < grid[(FMHeap[j])].value) {
			FMHeap[i] = FMHeap[j];
			grid[(FMHeap[j])].HeapPosition = i;
			FMHeap[j] = index;
			grid[index].HeapPosition = j;
		}
		else
			break;
	}
}

int LevelSetFMGrid3D::PopHeap() {
	if (heapSize == 0)
		return -1;
	int j, index = FMHeap[0];
	grid[index].DoneFlag = 1;
	heapSize--;
	FMHeap[0] = FMHeap[heapSize];
	grid[FMHeap[heapSize]].HeapPosition = 0;
	for (int i = 0; i < (heapSize - 1); i = j) {
		int lc = 2 * i + 1;
		int rc = 2 * i + 2;
		double current = grid[(FMHeap[i])].value;
		double lv, rv;
		if (lc < heapSize) {
			lv = grid[(FMHeap[lc])].value;
			if (rc < heapSize) {
				rv = grid[(FMHeap[rc])].value;
				if (lv > rv) {
					lc = rc;
					lv = rv;
				}
			}
			if (current > lv) {
				FMHeap[i] = FMHeap[lc];
				grid[FMHeap[i]].HeapPosition = i;
				FMHeap[lc] = FMHeap[heapSize];
				grid[FMHeap[heapSize]].HeapPosition = lc;
				j = lc;
			}
			else
				break;
		}
		else
			break;
	}
	FMHeap.pop_back();
	return index;
}

int LevelSetFMGrid3D::getindex(int x, int y, int z) {
	return x * size[1] * size[2] + y * size[2] + z;
}

double LevelSetFMGrid3D::Get(int x, int y, int z) {
	return grid[getindex(x, y, z)].value;
}

void LevelSetFMGrid3D::Reinitialize()
{
	heapSize = 0;
	closeSize = 0;
	FMHeap.clear();
	ClosePoints.clear();
	Initialize();
	InitHeap();
	FastMarch();
}

void LevelSetFMGrid3D::getxyz(int index, int& x, int& y, int&z) {
	x = index / (size[1] * size[2]);
	y = index % (size[1] * size[2]);
	z = y % size[2];
	y = y / size[2];
}

void LevelSetFMGrid3D::FindPhi(int x, int y, int z) {
	double phiX = 0;
	double phiY = 0;
	double phiZ = 0;
	int a = 0;
	bool flagX = 0;
	bool flagY = 0;
	bool flagZ = 0;
	//Find The phiS
	CheckFront(phiX, a, flagX, x + 1, size[0], x + 1, y, z);
	CheckBehind(phiX, a, flagX, x - 1, x - 1, y, z);
	CheckFront(phiY, a, flagY, y + 1, size[1], x, y + 1, z);
	CheckBehind(phiY, a, flagY, y - 1, x, y - 1, z);
	CheckFront(phiZ, a, flagZ, z + 1, size[2], x, y, z + 1);
	CheckBehind(phiZ, a, flagZ, z - 1, x, y, z - 1);

	//Max Tests
	if (a == 3)
	{
		if ((phiX >= phiY) && (phiX >= phiZ))
			CheckMax3(a, flagX, phiX, phiY, phiZ);
		else if ((phiY >= phiX) && (phiY >= phiZ))
			CheckMax3(a, flagY, phiY, phiX, phiZ);
		else
			CheckMax3(a, flagZ, phiZ, phiX, phiY);
	}
	if (a == 2)
	{
		if (!flagX)
		{
			if (phiY >= phiZ)
				CheckMax2(a, phiY, phiZ);
			else
				CheckMax2(a, phiZ, phiY);
		}
		else if (!flagY)
		{
			if (phiX >= phiZ)
				CheckMax2(a, phiX, phiZ);
			else
				CheckMax2(a, phiZ, phiX);
		}
		else
		{
			if (phiX >= phiY)
				CheckMax2(a, phiX, phiY);
			else
				CheckMax2(a, phiY, phiX);
		}
	}

	double b = phiX + phiY + phiZ;
	double quotient = square(b) - double(a) * (square(phiX) + square(phiY) + square(phiZ) - 1.0);
	if (quotient < 0)
	{
		std::cout << "0x ";
	}
	else
	{
		//cout << "1 ";
		double phi = b + std::sqrt(quotient);
		phi /= double(a);
		int index = getindex(x, y, z);

		//if(phi < 2)
		//cout << x << " " << y << " " << z << " " << a << " " << 
		//	phiX << " " << phiY << " " << phiZ << " " << 
		//	grid[index].value << " " << phi << endl;
		grid[index].value = phi;
		AddClose(index); //------------------------------------------------------------
	}
}

void LevelSetFMGrid3D::CheckMax2(int& a, double& phi1, double phi2) {
	if (square(phi1 - phi2) > 1) {
		phi1 = 0;
		a = 1;
	}
}

void LevelSetFMGrid3D::CheckMax3(int& a, bool& flag, double& phi1, double phi2, double phi3) {
	if ((square(phi1 - phi2) + square(phi1 - phi3)) > 1) {
		phi1 = 0;
		a = 2;
		flag = 0;
	}
}

void LevelSetFMGrid3D::CheckFront(double& phi, int& a, bool& flag, int check,
	int size, int x, int y, int z) {
	int index = getindex(x, y, z);
	if (check < size) {
		if (grid[index].DoneFlag == 1) {
			phi = grid[index].value;
			flag = 1;
			a++;
		}
	}
}

void LevelSetFMGrid3D::CheckBehind(double& phi, int& a, bool& flag, int check,
	int x, int y, int z) {
	int index = getindex(x, y, z);
	if (check >= 0) {
		if (grid[index].DoneFlag == 1) {
			if (!flag) {
				phi = grid[index].value;
				a++;
				flag = 1;
			}
			else
				phi = std::min(grid[index].value, phi);
		}
	}
}

void LevelSetFMGrid3D::FastMarch() {
	int x, y, z;
	for (int index = PopHeap(); index != -1; index = PopHeap()) {
		getxyz(index, x, y, z);
		//cout << index << "   " << x << "   " << y << "   " << z << "   "
		//	<< getindex(x,y,z) << endl;
		//cout << endl << grid[index].DoneFlag << "    " << heapSize << "       " 
		//	<< index << "      ";
		CheckPosExtent(x + 1, size[0], x + 1, y, z);
		CheckNegExtent(x - 1, x - 1, y, z);
		CheckPosExtent(y + 1, size[1], x, y + 1, z);
		CheckNegExtent(y - 1, x, y - 1, z);
		CheckPosExtent(z + 1, size[2], x, y, z + 1);
		CheckNegExtent(z - 1, x, y, z - 1);
		//cout << "     " << heapSize;
	}
	if (heapSize != 0)
	{
		std::cout << "bad heap!!!" << std::endl;
	}
}

void LevelSetFMGrid3D::CheckPosExtent(int index, int size, int x, int y, int z) {
	if (index < size) {
		if (grid[getindex(x, y, z)].DoneFlag == 0) {
			FindPhi(x, y, z);
			//if(getindex(x,y,z) >= 1331) 
			//cout << "   "  << getindex(x,y,z);
		}
	}
}

void LevelSetFMGrid3D::CheckNegExtent(int index, int x, int y, int z) {
	if (index >= 0) {
		if (grid[getindex(x, y, z)].DoneFlag == 0) {
			FindPhi(x, y, z);
			//if(getindex(x,y,z) >= 1331)
			//cout <<  "     " << getindex(x,y,z);
		}
	}
}

void LevelSetFMGrid3D::Set(int x, int y, int z, double value) {
	int index = getindex(x, y, z);
	grid[index].value = value;
	grid[index].HeapPosition = -1;
	if (value < 0)
		grid[index].DoneFlag = -1;
	else
		grid[index].DoneFlag = 0;
}

void LevelSetFMGrid3D::InitHeap() {
	for (int i = 0; i < closeSize; i++) {
		if (grid[ClosePoints[i]].HeapPosition == -1 && grid[ClosePoints[i]].DoneFlag == 0) {
			//if(grid[ClosePoints[i]].HeapPosition == -1) {
			int x, y, z;
			getxyz(ClosePoints[i], x, y, z);
			FindPhi(x, y, z);
		}
	}
}

void LevelSetFMGrid3D::CloseBehind(int index, int x, int y, int z) {
	if (index >= 0) {
		//if(grid[getindex(x,y,z)].value >= 0) {
		if (grid[getindex(x, y, z)].DoneFlag == 0) {
			//Add index to list for initialization of close band
			ClosePoints.push_back(getindex(x, y, z));
			closeSize++;
		}
	}
}

void LevelSetFMGrid3D::CloseFront(int size, int index, int x, int y, int z) {
	if (index < size) {
		//if(grid[getindex(x,y,z)].value >= 0) {
		if (grid[getindex(x, y, z)].DoneFlag == 0) {
			//Add index to list for initialization of close band
			ClosePoints.push_back(getindex(x, y, z));
			closeSize++;
		}
	}
}

void LevelSetFMGrid3D::Initialize()
{
	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{
			for (int k = 0; k < size[2] - 1; k++)
			{
				int cindex = getindex(i, j, k);
				int nindex = getindex(i, j, k + 1);
				double current = grid[cindex].value;
				double next = grid[nindex].value;

				if (current * next <= 0) {
					if (current >= 0) {
						grid[cindex].DoneFlag = 1;
						CloseBehind(k - 1, i, j, k - 1);
					}
					else {
						grid[nindex].DoneFlag = 1;
						CloseFront(size[2], k + 2, i, j, k + 2);
					}
				}
			}
		}
	}
	for (int i = 0; i < size[0]; i++)
	{
		for (int k = 0; k < size[2]; k++)
		{
			for (int j = 0; j < size[1] - 1; j++)
			{
				int cindex = getindex(i, j, k);
				int nindex = getindex(i, j + 1, k);
				double current = grid[cindex].value;
				double next = grid[nindex].value;

				if (current * next <= 0) {
					if (current >= 0) {
						grid[cindex].DoneFlag = 1;
						CloseBehind(j - 1, i, j - 1, k);
					}
					else {
						grid[nindex].DoneFlag = 1;
						CloseFront(size[1], j + 2, i, j + 2, k);
					}
				}
			}
		}
	}
	for (int j = 0; j < size[1]; j++)
	{
		for (int k = 0; k < size[2]; k++)
		{
			for (int i = 0; i < size[0] - 1; i++)
			{
				int cindex = getindex(i, j, k);
				int nindex = getindex(i + 1, j, k);
				double current = grid[cindex].value;
				double next = grid[nindex].value;

				if (current * next <= 0) {
					if (current >= 0) {
						grid[cindex].DoneFlag = 1;
						CloseBehind(i - 1, i - 1, j, k);
					}
					else {
						grid[nindex].DoneFlag = 1;
						CloseFront(size[0], i + 2, i + 2, j, k);
					}
				}
			}
		}
	}
}

/*
//SECOND-ORDER ACCURATE
void LevelSetFMGrid::FindPhi(int x, int y, int z) {
double phiXf, phiXb, phiYf, phiYb, phiZf, phiZb;
double quadCoef[3] = {0};
bool flagX[4] = {0};
bool flagY[4] = {0};
bool flagZ[4] = {0};

CheckFront(phiXf, flagX, x+1, size[0], x+1, y, z, x+2, y, z);
CheckBehind(phiXb, flagX, x-1, x-1, y, z, x-2, y, z);
CheckFront(phiYf, flagY, y+1, size[1], x, y+1, z, x, y+2, z);
CheckBehind(phiYb, flagY, y-1, x, y-1, z, x, y-2, z);
CheckFront(phiZf, flagZ, z+1, size[2], x, y, z+1, x, y, z+2);
CheckBehind(phiZb, flagZ, z-1, x, y, z-1, x, y, z-2);
FindQuadCoef(quadCoef, phiXf, phiXb, flagX);
FindQuadCoef(quadCoef, phiYf, phiYb, flagY);
FindQuadCoef(quadCoef, phiZf, phiZb, flagZ);

double phi = quadCoef[1] +
std::sqrt(square(quadCoef[1] - 4 * quadCoef[0] * (quadCoef[2] - 1)));
phi /= (2*quadCoef[0]);
int index = getindex(x,y,z);
grid[index].value = phi;
AddClose(index);
}

void LevelSetFMGrid::CheckFront(double& phi, bool flag [], int check, int size,
int x1, int y1, int z1, int x2, int y2, int z2)
{
if(check < size) {
int index = getindex(x1,y1,z1);
if(grid[index].DoneFlag == 1) {
phi = grid[index].value;
flag[0] = 1;
if((check+1) < size) {
index = getindex(x2,y2,z2);
if(grid[index].DoneFlag == 1) {
phi = 4*phi - grid[index].value;
flag[1] = 1;
}
}
}
}
}

void LevelSetFMGrid::CheckBehind(double& phi, bool flag [], int check,
int x1, int y1, int z1, int x2, int y2, int z2)
{
if(check >= 0) {
int index = getindex(x1,y1,z1);
if(grid[index].DoneFlag == 1) {
phi = grid[index].value;
flag[2] = 1;
if((check-1) >= 0) {
index = getindex(x2,y2,z2);
if(grid[index].DoneFlag == 1) {
phi = 4*phi - grid[index].value;
flag[3] = 1;
}
}
}
}
}

void LevelSetFMGrid::FindQuadCoef(double quadCoef [], double phif, double phib, bool flag []) {
if(flag[1] & flag[3])
SOQuadCoef(quadCoef, min(phif, phib));
else if(flag[1] & !flag[2])
SOQuadCoef(quadCoef, phif);
else if(!flag[0] & flag[3])
SOQuadCoef(quadCoef, phib);
else if(flag[1] & flag[2]) {
if((phif-1) < (3*phib))
SOQuadCoef(quadCoef, phif);
else
FOQuadCoef(quadCoef, phib);
}
else if(flag[0] & flag[3]) {
if((phib-1) < (3*phif))
SOQuadCoef(quadCoef, phib);
else
FOQuadCoef(quadCoef, phif);
}
else if(flag[0] & flag[2])
FOQuadCoef(quadCoef, min(phif, phib));
else if(flag[0] & !flag[2])
FOQuadCoef(quadCoef, phif);
else if(!flag[0] & flag[2])
FOQuadCoef(quadCoef, phib);
}

void LevelSetFMGrid::SOQuadCoef(double quadCoef [], double phi) {
quadCoef[0] += 2.25;
quadCoef[1] += 1.5*phi;
quadCoef[2] += 0.25* square(phi);
}

void LevelSetFMGrid::FOQuadCoef(double quadCoef [], double phi) {
quadCoef[0]++;
quadCoef[1] += 2*phi;
quadCoef[2] += square(phi);
}*/

/*double min(double a, double b) {
if (a < b)
return a;
return b;
}*/
