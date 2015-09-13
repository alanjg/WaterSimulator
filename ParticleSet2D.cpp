#include "stdafx.h"
#include "ParticleSet2D.h"
#include "LevelSet2D.h"
#include "Particle2D.h"
#include "StableFluid2D.h"

ParticleSet2D::ParticleSet2D(int xResolution, int yResolution)
{
	sizeX = xResolution;
	sizeY = yResolution;

	gridLength = sizeX*sizeY;
	grid1 = new std::list<Particle2D*>[gridLength];
	grid2 = new std::list<Particle2D*>[gridLength];
}

ParticleSet2D::~ParticleSet2D()
{
	for (int i = 0; i < gridLength; i++)
	{
		for (std::list<Particle2D*>::iterator it = grid1[i].begin(); it != grid1[i].end(); it++)
		{
			delete *it;
		}
		for (std::list<Particle2D*>::iterator it = grid2[i].begin(); it != grid2[i].end(); it++)
		{
			delete *it;
		}
	}
	delete[] grid1;
	delete[] grid2;
}

void ParticleSet2D::Update(const StableFluid2D& g0, const StableFluid2D& g1, const StableFluid2D& g2, double timestep)
{
	for (int i = 0; i < gridLength; i++)
	{
		int xi, yi;
		getposition(i, xi, yi);
		for (Iterator it = grid1[i].begin(); it != grid1[i].end();)
		{
			Particle2D* p = *it;
			grid1[i].erase(it++);

			double x, y;
			p->GetPosition(x, y);

			//get from velocity grid
			double ux0, uy0, ux1, uy1, ux2, uy2;
			g0.getvalue(x, y, ux0, uy0);
			g1.getvalue(x, y, ux1, uy1);
			g2.getvalue(x, y, ux2, uy2);
			p->Update(ux0, uy0, ux1, uy1, ux2, uy2, timestep);
			p->GetPosition(x, y);

			if (x < 0 || x >= sizeX || y < 0 || y >= sizeY)
				continue;
			grid2[getindex(int(x), int(y))].push_back(p);
		}
	}
	std::swap(grid1, grid2);
}

void ParticleSet2D::UpdateWorker(WorkerData data)
{
	for (int i = data.iBegin; i < data.iEnd; i++)
	{
		int xi, yi;
		getposition(i, xi, yi);
		for (Iterator it = grid1[i].begin(); it != grid1[i].end(); it++)
		{
			Particle2D* p = *it;

			double x, y;
			p->GetPosition(x, y);

			//get from velocity grid
			double ux0, uy0, ux1, uy1, ux2, uy2;
			(*data.g0).getvalue(x, y, ux0, uy0);
			(*data.g1).getvalue(x, y, ux1, uy1);
			(*data.g2).getvalue(x, y, ux2, uy2);
			p->Update(ux0, uy0, ux1, uy1, ux2, uy2, data.timestep);
		}
	}
}

void ParticleSet2D::UpdateThreaded(const StableFluid2D& g0, const StableFluid2D& g1, const StableFluid2D& g2, double timestep)
{
	int total = std::thread::hardware_concurrency();
	int each = gridLength / total;
	if (gridLength % total != 0)
	{
		each++;
	}
	std::vector<std::thread> threads;
	for (int i = 0; i < total; i++)
	{
		WorkerData d;
		d.g0 = &g0;
		d.g1 = &g1;
		d.g2 = &g2;
		d.timestep = timestep;
		d.iBegin = i * each;
		d.iEnd = i == total - 1 ? gridLength : (i + 1)*each;
		d.timestep = timestep;
		threads.push_back(std::thread(&ParticleSet2D::UpdateWorker, this, d));
	}

	for (unsigned int i = 0; i < threads.size(); i++)
	{
		threads[i].join();
	}
	threads.clear();
	for (int i = 0; i < gridLength; i++)
	{
		for (Iterator it = grid1[i].begin(); it != grid1[i].end();)
		{
			Particle2D* p = *it;
			grid1[i].erase(it++);

			double x, y;
			p->GetPosition(x, y);

			if (x < 0 || x >= sizeX || y < 0 || y >= sizeY)
				continue;
			grid2[getindex(int(x), int(y))].push_back(p);
		}
	}
	std::swap(grid1, grid2);
}

void ParticleSet2D::Reseed(const LevelSet2D& levelSet)
{
	for (int i = 0; i<gridLength; i++)
	{
		for (std::list<Particle2D*>::iterator it = grid1[i].begin(); it != grid1[i].end(); it++)
		{
			delete *it;
		}
		grid1[i].clear();
		int x, y;
		getposition(i, x, y);

		bool reseed = false;
		for (int dx = 0; dx<2; dx++)
		{
			for (int dy = 0; dy < 2; dy++)
			{
				if (std::fabs(levelSet.LinearSample(x + dx, y + dy)) < RESEED_THRESHOLD)
				{
					reseed = true;
				}
			}
		}
		if (reseed)
		{
			for (int j = 0; j < PARTICLES_PER_NODE; j++)
			{
				double dx = rand() / double(RAND_MAX);
				double dy = rand() / double(RAND_MAX);
				double phi = levelSet.LinearSample(x + dx, y + dy);

				grid2[i].push_back(new Particle2D(x + dx, y + dy, phi));
			}
		}
	}
	std::swap(grid1, grid2);
}

ParticleSet2D::Iterator ParticleSet2D::begin(int i, int j) const
{
	return grid1[getindex(i, j)].begin();
}

ParticleSet2D::Iterator ParticleSet2D::end(int i, int j) const
{
	return grid1[getindex(i, j)].end();
}

int ParticleSet2D::getindex(int i, int j) const
{
	return i*sizeY + j;
}

void ParticleSet2D::getposition(int index, int& x, int& y) const
{
	y = index % sizeY;
	index = index / sizeY;
	x = index;
}
