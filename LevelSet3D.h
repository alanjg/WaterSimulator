#pragma once

class StableFluid;
class LevelSetGrid3D;
class LevelSetFMGrid3D;
class StableFluid3D;
class ParticleSet3D;

class LevelSet3D
{
	//this should stay at 1 unless the index functions are scaled
	const int GRID_SIZE = 1;
public:
	LevelSet3D(int xResolution, int yResolution, int zResolution);
	~LevelSet3D();

	void MakeSphere(double x, double y, double z, double radius);

	void Update(const StableFluid3D& f0, const StableFluid3D& f1, const StableFluid3D& f2, double timestep);
	void UpdateThreaded(const StableFluid3D& f0, const StableFluid3D& f1, const StableFluid3D& f2, double timestep);

	void Fix(const ParticleSet3D& particleSet);

	double LinearSample(double x, double y, double z) const;
	double CubicSample(double x, double y, double z) const;

	double Curvature(double x, double y, double z) const;
	void normal(double x, double y, double z, double& nx, double& ny, double& nz) const;
	void gradient(double x, double y, double z, double& gx, double& gy, double& gz) const;

	void Reinitialize();

	//virtual double	eval(const Point3d& location);

	void Write(std::ostream& out);
	void CheckGrid();
private:

	// for threading
	struct WorkerData
	{
		const StableFluid3D* field;
		int iBegin;
		int iEnd;
		double timestep;
	};
	void UpdateWorker(WorkerData data);

	//the non-generic versions may be faster
	double D0(int x, int y, int z);

	double D1x(int x, int y, int z);
	double D1y(int x, int y, int z);
	double D1z(int x, int y, int z);
	double D1(int x, int y, int z, int dim);

	double D2x(int x, int y, int z);
	double D2y(int x, int y, int z);
	double D2z(int x, int y, int z);
	double D2(int x, int y, int z, int dim);

	double D3x(int x, int y, int z);
	double D3y(int x, int y, int z);
	double D3z(int x, int y, int z);
	double D3(int x, int y, int z, int dim);

	void gradient(double x, double y, double z, double& gx, double& gy, double& gz, double ux, double uy, double uz);

	void EulerStep(int x, int y, int z, const StableFluid3D& grid, double timestep);
	void ENOStep(int x, int y, int z, const StableFluid3D& grid, double timestep);
	void WENOStep(int x, int y, int z, const StableFluid3D& grid, double timestep);

	//grid size in each dimension
	int size[3];

	LevelSetGrid3D* grid1;
	LevelSetGrid3D* grid2;
	LevelSetGrid3D* grid3;

	LevelSetFMGrid3D* gridFM;

	void GetVelocityFieldValue(int x, int y, int z, const StableFluid3D& grid, double& ux, double& uy, double& uz);
};

