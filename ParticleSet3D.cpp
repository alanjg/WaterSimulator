#include "stdafx.h"
#include "ParticleSet3D.h"
#include "LevelSet3D.h"
#include "Particle3D.h"
#include "StableFluid3D.h"

ParticleSet3D::ParticleSet3D(int xResolution, int yResolution, int zResolution)
{
	sizeX = xResolution;
	sizeY = yResolution;
	sizeZ = zResolution;

	gridLength = sizeX*sizeY*sizeZ;
	grid1 = new std::list<Particle3D*>[gridLength];
	grid2 = new std::list<Particle3D*>[gridLength];
}

ParticleSet3D::~ParticleSet3D()
{
	for (int i = 0; i < gridLength; i++)
	{
		for (std::list<Particle3D*>::iterator it = grid1[i].begin(); it != grid1[i].end(); it++)
		{
			delete *it;
		}
		for (std::list<Particle3D*>::iterator it = grid2[i].begin(); it != grid2[i].end(); it++)
		{
			delete *it;
		}
	}
	delete[] grid1;
	delete[] grid2;
}

void ParticleSet3D::Update(const StableFluid3D& g0, const StableFluid3D& g1, const StableFluid3D& g2, double timestep)
{
	for (int i = 0; i < gridLength; i++)
	{
		int xi, yi, zi;
		getposition(i, xi, yi, zi);
		for (Iterator it = grid1[i].begin(); it != grid1[i].end();)
		{
			Particle3D* p = *it;
			grid1[i].erase(it++);

			double x, y, z;
			p->GetPosition(x, y, z);

			//get from velocity grid
			double ux0, uy0, uz0, ux1, uy1, uz1, ux2, uy2, uz2;
			g0.getvalue(x, y, z, ux0, uy0, uz0);
			g1.getvalue(x, y, z, ux1, uy1, uz1);
			g2.getvalue(x, y, z, ux2, uy2, uz2);
			p->Update(ux0, uy0, uz0, ux1, uy1, uz1, ux2, uy2, uz2, timestep);
			p->GetPosition(x, y, z);

			if (x < 0 || x >= sizeX || y < 0 || y >= sizeY || z < 0 || z >= sizeZ)
				continue;
			grid2[getindex(int(x), int(y), int(z))].push_back(p);
		}
	}
	std::swap(grid1, grid2);
}

void ParticleSet3D::UpdateWorker(WorkerData data)
{
	for (int i = data.iBegin; i < data.iEnd; i++)
	{
		int xi, yi, zi;
		getposition(i, xi, yi, zi);
		for (Iterator it = grid1[i].begin(); it != grid1[i].end(); it++)
		{
			Particle3D* p = *it;

			double x, y, z;
			p->GetPosition(x, y, z);

			//get from velocity grid
			double ux0, uy0, uz0, ux1, uy1, uz1, ux2, uy2, uz2;
			(*data.g0).getvalue(x, y, z, ux0, uy0, uz0);
			(*data.g1).getvalue(x, y, z, ux1, uy1, uz1);
			(*data.g2).getvalue(x, y, z, ux2, uy2, uz2);
			p->Update(ux0, uy0, uz0, ux1, uy1, uz1, ux2, uy2, uz2, data.timestep);
		}
	}
}

void ParticleSet3D::UpdateThreaded(const StableFluid3D& g0, const StableFluid3D& g1, const StableFluid3D& g2, double timestep)
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
		threads.push_back(std::thread(&ParticleSet3D::UpdateWorker, this, d));
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
			Particle3D* p = *it;
			grid1[i].erase(it++);

			double x, y, z;
			p->GetPosition(x, y, z);

			if (x < 0 || x >= sizeX || y < 0 || y >= sizeY || z < 0 || z >= sizeZ)
				continue;
			grid2[getindex(int(x), int(y), int(z))].push_back(p);
		}
	}
	std::swap(grid1, grid2);
}

void ParticleSet3D::Reseed(const LevelSet3D& levelSet)
{
	for (int i = 0; i<gridLength; i++)
	{
		for (std::list<Particle3D*>::iterator it = grid1[i].begin(); it != grid1[i].end(); it++)
		{
			delete *it;
		}
		grid1[i].clear();
		int x, y, z;
		getposition(i, x, y, z);

		bool reseed = false;
		for (int dx = 0; dx<2; dx++)
		{
			for (int dy = 0; dy<2; dy++)
			{
				for (int dz = 0; dz<2; dz++)
				{
					if (std::fabs(levelSet.LinearSample(x + dx, y + dy, z + dz)) < RESEED_THRESHOLD)
					{
						reseed = true;
					}
				}
			}
		}
		if (reseed)
		{
			for (int j = 0; j < PARTICLES_PER_NODE; j++)
			{
				double dx = rand() / double(RAND_MAX);
				double dy = rand() / double(RAND_MAX);
				double dz = rand() / double(RAND_MAX);
				double phi = levelSet.LinearSample(x + dx, y + dy, z + dz);

				grid2[i].push_back(new Particle3D(x + dx, y + dy, z + dz, phi));
			}
		}
	}
	std::swap(grid1, grid2);
}

ParticleSet3D::Iterator ParticleSet3D::begin(int i, int j, int k) const
{
	return grid1[getindex(i, j, k)].begin();
}

ParticleSet3D::Iterator ParticleSet3D::end(int i, int j, int k) const
{
	return grid1[getindex(i, j, k)].end();
}

int ParticleSet3D::getindex(int i, int j, int k) const
{
	return i*sizeY*sizeZ + j*sizeZ + k;
}

void ParticleSet3D::getposition(int index, int& x, int& y, int& z) const
{
	z = index % sizeZ;
	index = index / sizeZ;
	y = index % sizeY;
	index = index / sizeY;
	x = index;
}
