#pragma once

double square(double x);

class LevelSetFMGrid3D
{
public:
	LevelSetFMGrid3D(int xResolution, int yResolution, int zResolution);
	~LevelSetFMGrid3D();
	void Set(int x, int y, int z, double value);
	double Get(int x, int y, int z);
	void Reinitialize();
private:
	void Initialize();
	void InitHeap();
	void CloseBehind(int index, int x, int y, int z);
	void CloseFront(int size, int index, int x, int y, int z);
	void FastMarch();
	void FindPhi(int x, int y, int z);
	void CheckPosExtent(int index, int size, int x, int y, int z);
	void CheckNegExtent(int index, int x, int y, int z);
	void CheckFront(double& phi, int& a, bool& flag,
		int check, int size, int x, int y, int z);
	void CheckBehind(double& phi, int& a, bool& flag,
		int check, int x, int y, int z);
	void CheckMax2(int& a, double& phi1, double phi2);
	void CheckMax3(int& a, bool& flag, double& phi1, double phi2, double phi3);
	void AddClose(int index);
	void AddToHeap(int index);
	void UpdateHeap(int index);
	int PopHeap();
	int getindex(int x, int y, int z);
	void getxyz(int index, int& x, int& y, int& z);
	struct FMContainer {
		int DoneFlag;
		int HeapPosition;
		double value;
	};
	int size[3];
	int gridLength;
	int heapSize;
	int closeSize;
	FMContainer* grid;
	std::vector<int> FMHeap;
	std::vector<int> ClosePoints;

	//void CheckFront(double& phi, bool flag [], int check, int size, 
	//	int x1, int y1, int z1, int x2, int y2, int z2);
	//void CheckBehind(double& phi, bool flag [], int check, 
	//	int x1, int y1, int z1, int x2, int y2, int z2);
	//void FindQuadCoef(double quadCoef [], double phif, double phib, bool flag []);
	//void SOQuadCoef(double quadCoef [], double phi);
	//void FOQuadCoef(double quadCoef [], double phi);
};
