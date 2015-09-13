#pragma once

class LevelSet3D;
class Particle3D;
class StableFluid3D;

class ParticleSet3D
{
	const double RESEED_THRESHOLD = 3.0;
	const int PARTICLES_PER_NODE = 64;

public:
	ParticleSet3D(int xResolution, int yResolution, int zResolution);
	~ParticleSet3D();

	void Update(const StableFluid3D& g0, const StableFluid3D& g1, const StableFluid3D& g2, double timestep);
	void UpdateThreaded(const StableFluid3D& g0, const StableFluid3D& g1, const StableFluid3D& g2, double timestep);
	void Reseed(const LevelSet3D& levelSet);

	typedef std::list<Particle3D*>::iterator Iterator;
	Iterator begin(int i, int j, int k) const;
	Iterator end(int i, int j, int k) const;


private:
	struct WorkerData
	{
		int iBegin, iEnd;
		const StableFluid3D *g0, *g1, *g2;
		double timestep;
	};
	void UpdateWorker(WorkerData data);
	int sizeX, sizeY, sizeZ;

	int gridLength;
	std::list<Particle3D*>* grid1;
	std::list<Particle3D*>* grid2;

	//turns a grid index into an array index
	int getindex(int i, int j, int k) const;

	//turns an array index into a grid index
	void getposition(int index, int& x, int& y, int& z) const;
};
