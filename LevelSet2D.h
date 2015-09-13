#pragma once

class LevelSetGrid2D;
class LevelSetFMGrid2D;
class StableFluid2D;
class ParticleSet2D;


class LevelSet2D
{
	//this should stay at 1 unless the index functions are scaled
	const int GRID_SIZE = 1;
public:
	LevelSet2D(int xResolution, int yResolution);
	~LevelSet2D();

	void MakeCircle(double x, double y, double radius);

	void Update(const StableFluid2D& f0, const StableFluid2D& f1, const StableFluid2D& f2, double timestep);
	void UpdateThreaded(const StableFluid2D& f0, const StableFluid2D& f1, const StableFluid2D& f2, double timestep);

	void Fix(const ParticleSet2D& particleSet);

	double LinearSample(double x, double y) const;
	double CubicSample(double x, double y) const;

	double Curvature(double x, double y) const;
	void normal(double x, double y, double& nx, double& ny) const;
	void gradient(double x, double y, double& gx, double& gy) const;

	void Reinitialize();

	//virtual double	eval(const Point3d& location);

	void Write(std::ostream& out);
	void CheckGrid();
private:

	// for threading
	struct WorkerData
	{
		const StableFluid2D* field;
		int iBegin;
		int iEnd;
		double timestep;
	};
	void UpdateWorker(WorkerData data);

	//the non-generic versions may be faster
	double D0(int x, int y);

	double D1x(int x, int y);
	double D1y(int x, int y);
	double D1(int x, int y, int dim);

	double D2x(int x, int y);
	double D2y(int x, int y);
	double D2(int x, int y, int dim);

	double D3x(int x, int y);
	double D3y(int x, int y);
	double D3(int x, int y, int dim);

	void gradient(double x, double y, double& gx, double& gy, double ux, double uy);

	void EulerStep(int x, int y, const StableFluid2D& grid, double timestep);
	void ENOStep(int x, int y, const StableFluid2D& grid, double timestep);
	void WENOStep(int x, int y, const StableFluid2D& grid, double timestep);

	//grid size in each dimension
	int size[2];

	LevelSetGrid2D* grid1;
	LevelSetGrid2D* grid2;
	LevelSetGrid2D* grid3;

	LevelSetFMGrid2D* gridFM;

	void GetVelocityFieldValue(int x, int y, const StableFluid2D& grid, double& ux, double& uy);
};

