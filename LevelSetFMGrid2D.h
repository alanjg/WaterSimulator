#pragma once

double square(double x);

class LevelSetFMGrid2D
{
	std::vector<int> pushed, popped;
public:
	LevelSetFMGrid2D(int xResolution, int yResolution);
	~LevelSetFMGrid2D();
	void Set(int x, int y, double value);
	double Get(int x, int y);
	void Reinitialize();
private:
	void Initialize();
	void InitHeap();
	void CloseBehind(int index, int x, int y);
	void CloseFront(int size, int index, int x, int y);
	void FastMarch();
	void FindPhi(int x, int y);
	void CheckPosExtent(int index, int size, int x, int y);
	void CheckNegExtent(int index, int x, int y);
	void CheckFront(double& phi, int& a, bool& flag, int check, int size, int x, int y);
	void CheckBehind(double& phi, int& a, bool& flag, int check, int x, int y);
	void CheckMax2(int& a, double& phi1, double phi2);
	void AddClose(int index);
	void AddToHeap(int index);
	void UpdateHeap(int index);
	int PopHeap();
	int getindex(int x, int y);
	void getxyz(int index, int& x, int& y);
	struct FMContainer 
	{
		int DoneFlag;
		int HeapPosition;
		double value;
	};
	int size[2];
	int gridLength;
	int heapSize;
	int closeSize;
	FMContainer* grid;
	std::vector<int> FMHeap;
	std::vector<int> ClosePoints;
};
