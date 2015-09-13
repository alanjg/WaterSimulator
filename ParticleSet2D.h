#pragma once

class LevelSet2D;
class Particle2D;
class StableFluid2D;

class ParticleSet2D
{

	const double RESEED_THRESHOLD = 3.0;
	const int PARTICLES_PER_NODE = 64;
public:
	ParticleSet2D(int xResolution, int yResolution);
	~ParticleSet2D();

	void Update(const StableFluid2D& g0, const StableFluid2D& g1, const StableFluid2D& g2, double timestep);
	void UpdateThreaded(const StableFluid2D& g0, const StableFluid2D& g1, const StableFluid2D& g2, double timestep);
	void Reseed(const LevelSet2D& levelSet);

	typedef std::list<Particle2D*>::iterator Iterator;
	Iterator begin(int i, int j) const;
	Iterator end(int i, int j) const;


private:
	struct WorkerData
	{
		int iBegin, iEnd;
		const StableFluid2D *g0, *g1, *g2;
		double timestep;
	};
	void UpdateWorker(WorkerData data);
	int sizeX, sizeY;

	int gridLength;
	std::list<Particle2D*>* grid1;
	std::list<Particle2D*>* grid2;

	//turns a grid index into an array index
	int getindex(int i, int j) const;

	//turns an array index into a grid index
	void getposition(int index, int& x, int& y) const;
};
