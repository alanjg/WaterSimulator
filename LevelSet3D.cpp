#include "stdafx.h"
#include "LevelSet3D.h"

#include "ParticleSet3D.h"
#include "Particle3D.h"
#include "StableFluid3D.h"
#include "LevelSetGrid3D.h"
#include "LevelSetFMGrid3D.h"
#include "math.h"

using std::cout;
using std::endl;
using std::min;

LevelSet3D::LevelSet3D(int xResolution, int yResolution, int zResolution)
{
	//level set rep is on the corners, extend the grid by 1
	size[0] = xResolution + 1;
	size[1] = yResolution + 1;
	size[2] = zResolution + 1;

	grid1 = new LevelSetGrid3D(size[0], size[1], size[2]);
	grid2 = new LevelSetGrid3D(size[0], size[1], size[2]);
	grid3 = new LevelSetGrid3D(size[0], size[1], size[2]);
	gridFM = new LevelSetFMGrid3D(size[0], size[1], size[2]);
}

LevelSet3D::~LevelSet3D()
{
	delete grid1;
	delete grid2;
	delete grid3;
	delete gridFM;
}

/*
double LevelSet3D::eval(const Point3d& location)
{
	double x = (location[0] + 1) * (size[0] - 1) / 2.0;
	double y = (location[1] + 1) * (size[1] - 1) / 2.0;
	double z = (location[2] + 1) * (size[2] - 1) / 2.0;
	double d = LinearSample(x, y, z);

	return d;
}
*/
void LevelSet3D::Write(std::ostream& out)
{
	grid1->Write(out);
}

void LevelSet3D::CheckGrid()
{
	for (int i = 0; i <= size[0]; i++)
	{
		for (int j = 0; j <= size[1]; j++)
		{
			for (int k = 0; k <= size[2]; k++)
			{
				double me = grid1->get(i, j, k);
				if (i > 0)
				{
					double other = grid1->get(i - 1, j, k);
					if (abs(me - other) > 1.5)
					{
						cout << "Bad grid spot at " << i << " " << j << " " << k << " " << me << " " << other << endl;
					}
				}

				if (j > 0)
				{
					double other = grid1->get(i, j - 1, k);
					if (abs(me - other) > 1.5)
					{
						cout << "Bad grid spot at " << i << " " << j << " " << k << " " << me << " " << other << endl;
					}

				}
				if (k > 0)
				{
					double other = grid1->get(i, j, k - 1);
					if (abs(me - other) > 1.5)
					{
						cout << "Bad grid spot at " << i << " " << j << " " << k << " " << me << " " << other << endl;
					}
				}
			}
		}
	}
}


void LevelSet3D::MakeSphere(double x, double y, double z, double radius)
{
	for (int i = 0; i <= size[0]; i++)
	{
		for (int j = 0; j <= size[1]; j++)
		{
			for (int k = 0; k <= size[2]; k++)
			{
				double val1 = std::sqrt(square(x - i) + square(y - j) + square(z - k)) - radius;

				double leftWall = x - radius / 3;
				double rightWall = x + radius / 3;
				double topWall = y + radius / 3;
				double bottomWall = y - radius;

				double val2 = 999999;

				if (bottomWall <= j && j <= topWall && leftWall <= i && i <= rightWall)
				{
					double top = topWall - j;
					double bottom = j - bottomWall;
					double left = i - leftWall;
					double right = rightWall - i;
					if (fabs(val2) > fabs(top))
						val2 = top;
					if (fabs(val2) > fabs(bottom))
						val2 = bottom;
					if (fabs(val2) > fabs(left))
						val2 = left;
					if (fabs(val2) > fabs(right))
						val2 = right;
					//				cout<<"In all 4"<<endl;
				}
				else if (leftWall <= i && i <= rightWall)
				{
					if (fabs(val2) > fabs(j - bottomWall))
						val2 = j - bottomWall;
					if (fabs(val2) > fabs(j - topWall))
						val2 = topWall - j;
					//			cout<<"In vertically"<<endl;
				}
				else if (bottomWall <= j && j <= topWall)
				{
					if (fabs(val2) > fabs(i - leftWall))
						val2 = i - leftWall;
					if (fabs(val2) > fabs(i - rightWall))
						val2 = rightWall - i;
					//			cout<<"In horizontally"<<endl;
				}
				else
				{
					double ul = sqrt((i - leftWall)*(i - leftWall) + (j - topWall)*(j - topWall));
					double ur = sqrt((i - rightWall)*(i - rightWall) + (j - topWall)*(j - topWall));
					double bl = sqrt((i - leftWall)*(i - leftWall) + (j - bottomWall)*(j - bottomWall));
					double br = sqrt((i - rightWall)*(i - rightWall) + (j - bottomWall)*(j - bottomWall));
					val2 = -min(min(min(ul, ur), bl), br);
					//		cout<<"Not in"<<endl;
				}

				//			val = 10*cos(i/double(size[0])*10) + 10*sin(j/double(size[1])*10) ;
				//			std::cout << val2 << std::endl; 

				grid1->get(i, j, k) = val1;
				//			grid1->get(i,j,k) = std::max(val1,val2);
				//				grid1->get(i,j,k) = j - 15;
			}
		}
	}
	grid1->UpdateBorders();
}

void LevelSet3D::UpdateWorker(WorkerData data)
{
	//Find phi(n+1)
	for (int i = data.iBegin; i < data.iEnd; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{
			for (int k = 0; k < size[2]; k++)
			{
				//				EulerStep(i,j,k,grid,timestep);
				//				ENOStep(i,j,k,grid,timestep);
				WENOStep(i, j, k, *data.field, data.timestep);
			}
		}

	}
}

void LevelSet3D::UpdateThreaded(const StableFluid3D& f0, const StableFluid3D& f1, const StableFluid3D& f2, double timestep)
{
	int total = std::thread::hardware_concurrency();
	int each = size[0] / total;
	int extra = size[0] % total;

	std::vector<std::thread> threads;
	int at = 0;
	//TVD RK THREE - see page 38 of osher & fedkiw
	for (int i = 0; i < total; i++)
	{
		WorkerData d;
		d.field = &f0;
		d.iBegin = at;
		at += each;
		if (extra > 0)
		{
			extra--;
			at++;
		}
		d.iEnd = at;
		d.timestep = timestep;
		threads.push_back(std::thread(&LevelSet3D::UpdateWorker, this, d));
	}

	for (unsigned int i = 0; i < threads.size(); i++)
	{
		threads[i].join();
	}
	threads.clear();

	grid2->UpdateBorders();
	std::swap(grid1, grid3); //grid3 holds phi(n)
	std::swap(grid1, grid2); //gird1 holds phi(n+1)

	for (int i = 0; i < total; i++)
	{
		WorkerData d;
		d.field = &f2;
		d.iBegin = i * each;
		d.iEnd = i == total - 1 ? size[0] : (i + 1)*each;
		d.timestep = timestep;
		threads.push_back(std::thread(&LevelSet3D::UpdateWorker, this, d));
	}

	for (unsigned int i = 0; i < threads.size(); i++)
	{
		threads[i].join();
	}
	threads.clear();


	grid2->UpdateBorders();

	//Set grid1 to phi(n+1/2) = .75*phi(n) + .25*phi(n+2)
	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{
			for (int k = 0; k < size[2]; k++)
			{
				grid1->get(i, j, k) = grid3->get(i, j, k)*0.75 + grid2->get(i, j, k)*0.25;
			}
		}
	}

	grid1->UpdateBorders();

	//Find phi(n+3/2)
	for (int i = 0; i < total; i++)
	{
		WorkerData d;
		d.field = &f1;
		d.iBegin = i * each;
		d.iEnd = i == total - 1 ? size[0] : (i + 1)*each;
		d.timestep = timestep;
		threads.push_back(std::thread(&LevelSet3D::UpdateWorker, this, d));
	}

	for (unsigned int i = 0; i < threads.size(); i++)
	{
		threads[i].join();
	}
	threads.clear();

	grid2->UpdateBorders();

	//Set grid1 to final phi(n+1) = phi(n) * 1/3 + phi(n+3/2) * 2/3
	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{
			for (int k = 0; k < size[2]; k++)
			{
				grid1->get(i, j, k) = grid3->get(i, j, k)*0.333333333 +
					grid2->get(i, j, k)*0.666666667;
			}
		}
	}

	grid1->UpdateBorders();
}

void LevelSet3D::Update(const StableFluid3D& f0, const StableFluid3D& f1, const StableFluid3D& f2, double timestep)
{
	//TVD RK THREE - see page 38 of osher & fedkiw

	//Find phi(n+1)
	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{
			for (int k = 0; k < size[2]; k++)
			{
				//				EulerStep(i,j,k,grid,timestep);
				//				ENOStep(i,j,k,grid,timestep);
				WENOStep(i, j, k, f0, timestep);
			}
		}
	}

	grid2->UpdateBorders();
	std::swap(grid1, grid3); //grid3 holds phi(n)
	std::swap(grid1, grid2); //grid1 holds phi(n+1)

								//Find phi(n+2)
	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{
			for (int k = 0; k < size[2]; k++)
			{
				//				EulerStep(i,j,k,grid,timestep);
				//				ENOStep(i,j,k,grid,timestep);
				WENOStep(i, j, k, f2, timestep);
			}
		}
	}

	grid2->UpdateBorders();

	//Set grid1 to phi(n+1/2) = .75*phi(n) + .25*phi(n+2)
	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{
			for (int k = 0; k < size[2]; k++)
			{
				grid1->get(i, j, k) = grid3->get(i, j, k)*0.75 + grid2->get(i, j, k)*0.25;
			}
		}
	}

	grid1->UpdateBorders();

	//Find phi(n+3/2)
	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{
			for (int k = 0; k < size[2]; k++)
			{
				//				EulerStep(i,j,k,grid,timestep);
				//				ENOStep(i,j,k,grid,timestep);
				WENOStep(i, j, k, f1, timestep);
			}
		}
	}

	grid2->UpdateBorders();

	//Set grid1 to final phi(n+1) = phi(n) * 1/3 + phi(n+3/2) * 2/3
	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{
			for (int k = 0; k < size[2]; k++)
			{
				grid1->get(i, j, k) = grid3->get(i, j, k)*0.333333333 +
					grid2->get(i, j, k)*0.666666667;
			}
		}
	}

	grid1->UpdateBorders();
}

//WE NEED ADVECTION ALONG NORMAL DIRECTION. SAME AS WENO, WITHOUGH VELOCITY GRID.
//USE NORMAL FOR VELOCITY. USE EQUATION 7.4 INSTEAD OF 3.2

void LevelSet3D::Reinitialize() {
	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{
			for (int k = 0; k < size[2]; k++) {
				gridFM->Set(i, j, k, grid1->get(i, j, k));
				//if(grid1->get(i,j,k) > -2 && grid1->get(i,j,k) < 2 && j < 10)
				//	cout << i << "   " << j << "   " << k << "   " <<
				//		grid1->get(i,j,k) << endl;
			}
		}
	}
	//cout << endl << endl;
	gridFM->Reinitialize();
	//cout << endl << endl;
	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{
			for (int k = 0; k < size[2]; k++) {
				grid1->get(i, j, k) = gridFM->Get(i, j, k);
				//if(grid1->get(i,j,k) > 0 && grid1->get(i,j,k) < 1)
				//	cout << i << "   " << j << "   " << k << "   " << 
				//		grid1->get(i,j,k) << endl;
			}
		}
	}

	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{
			for (int k = 0; k < size[2]; k++) {
				gridFM->Set(i, j, k, grid1->get(i, j, k)*(-1));
				//if(grid1->get(i,j,k)*(-1) > 0)
				//	cout << i << "   " << j << "   " << k << "   " << 
				//		grid1->get(i,j,k) << endl;
			}
		}
	}
	//cout << endl << endl;
	gridFM->Reinitialize();
	//cout << endl << endl;
	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{
			for (int k = 0; k < size[2]; k++) {
				grid1->get(i, j, k) = gridFM->Get(i, j, k)*(-1);
				//if(grid1->get(i,j,k)*(-1) > 0)
				//	cout << i << "   " << j << "   " << k << "   " << 
				//		grid1->get(i,j,k) << endl;
			}
		}
	}
	grid1->UpdateBorders();
}

void LevelSet3D::EulerStep(int x, int y, int z, const StableFluid3D& grid, double timestep)
{
	double ux, uy, uz; //obtain from velocity grid

	GetVelocityFieldValue(x, y, z, grid, ux, uy, uz);

	double gx, gy, gz;
	gradient(x, y, z, gx, gy, gz, ux, uy, uz);

	double udotgrad = gx*ux + gy*uy + gz*uz;
	grid2->get(x, y, z) = grid1->get(x, y, z) - udotgrad * timestep;
}

double LevelSet3D::D0(int x, int y, int z)
{
	return grid1->get(x, y, z);
}

double LevelSet3D::D1x(int x, int y, int z)
{
	return (D0(x + 1, y, z) - D0(x, y, z)) / GRID_SIZE;
}

double LevelSet3D::D1y(int x, int y, int z)
{
	return (D0(x, y + 1, z) - D0(x, y, z)) / GRID_SIZE;
}

double LevelSet3D::D1z(int x, int y, int z)
{
	return (D0(x, y, z + 1) - D0(x, y, z)) / GRID_SIZE;
}

double LevelSet3D::D1(int x, int y, int z, int dim)
{
	int dx[] = { 1,0,0 };
	int dy[] = { 0,1,0 };
	int dz[] = { 0,0,1 };
	return (D0(x + dx[dim], y + dy[dim], z + dz[dim]) - D0(x, y, z)) / GRID_SIZE;
}

double LevelSet3D::D2x(int x, int y, int z)
{
	return (D1x(x, y, z) - D1x(x - 1, y, z)) / (2 * GRID_SIZE);
}

double LevelSet3D::D2y(int x, int y, int z)
{
	return (D1x(x, y, z) - D1x(x, y - 1, z)) / (2 * GRID_SIZE);
}

double LevelSet3D::D2z(int x, int y, int z)
{
	return (D1x(x, y, z) - D1x(x, y, z - 1)) / (2 * GRID_SIZE);
}

double LevelSet3D::D2(int x, int y, int z, int dim)
{
	int dx[] = { 1,0,0 };
	int dy[] = { 0,1,0 };
	int dz[] = { 0,0,1 };
	return (D1(x, y, z, dim) - D1(x - dx[dim], y - dy[dim], z - dz[dim], dim)) / (2 * GRID_SIZE);
}

double LevelSet3D::D3x(int x, int y, int z)
{
	return (D2x(x + 1, y, z) - D2x(x, y, z)) / (3 * GRID_SIZE);
}

double LevelSet3D::D3y(int x, int y, int z)
{
	return (D2x(x, y + 1, z) - D2x(x, y, z)) / (3 * GRID_SIZE);
}

double LevelSet3D::D3z(int x, int y, int z)
{
	return (D2x(x, y, z + 1) - D2x(x, y, z)) / (3 * GRID_SIZE);
}

double LevelSet3D::D3(int x, int y, int z, int dim)
{
	int dx[] = { 1,0,0 };
	int dy[] = { 0,1,0 };
	int dz[] = { 0,0,1 };
	return (D2(x + dx[dim], y + dy[dim], z + dz[dim], dim) - D2(x, y, z, dim)) / (3 * GRID_SIZE);
}

void LevelSet3D::ENOStep(int x, int y, int z, const StableFluid3D& grid, double timestep)
{
	double ux, uy, uz; //obtain from velocity grid

	GetVelocityFieldValue(x, y, z, grid, ux, uy, uz);

	//we solve for these, see page 32 in Osher & Fedkiw
	double gx, gy, gz;

	int pos[3] = { x,y,z };

	int dx[] = { 1,0,0 };
	int dy[] = { 0,1,0 };
	int dz[] = { 0,0,1 };

	double u[3] = { ux,uy,uz };
	double g[3];
	for (int i = 0; i < 3; i++)
	{
		int index[3] = { x,y,z };
		int k;
		if (u[i] > 0)
		{
			index[i]--;
		}
		k = index[i];
		double q1 = D1(index[0], index[1], index[2], i);

		double dk0 = D2(index[0], index[1], index[2], i);
		double dk1 = D2(index[0] + dx[i], index[1] + dy[i], index[2] + dz[i], i);

		double c;
		int kstar;
		if (fabs(dk0) <= fabs(dk1))
		{
			c = dk0;
			index[i]--;
		}
		else
		{
			c = dk1;
		}
		kstar = index[i];
		double q2 = c*(2 * (pos[i] - k) - 1)*GRID_SIZE;

		dk0 = D3(index[0], index[1], index[2], i);
		dk1 = D3(index[0] + dx[i], index[1] + dy[i], index[2] + dz[i], i);
		double cstar;
		if (fabs(dk0) <= fabs(dk1))
		{
			cstar = dk0;
		}
		else
		{
			cstar = dk1;
		}

		double q3 = cstar*(3 * square(pos[i] - kstar) - 6 * (pos[i] - kstar) + 2)*(GRID_SIZE * GRID_SIZE);
		g[i] = q1 + q2 + q3;
	}
	gx = g[0];
	gy = g[1];
	gz = g[2];
	double udotgrad = gx*ux + gy*uy + gz*uz;
	grid2->get(x, y, z) = grid1->get(x, y, z) - udotgrad * timestep;
}

void LevelSet3D::WENOStep(int x, int y, int z, const StableFluid3D& grid, double timestep)
{
	double ux, uy, uz; //obtain from velocity grid

	GetVelocityFieldValue(x, y, z, grid, ux, uy, uz);

	//we solve for these, see page 34 in Osher & Fedkiw
	double gx, gy, gz;

	double v1, v2, v3, v4, v5;
	if (ux < 0)
	{
		v1 = D1x(x + 2, y, z);
		v2 = D1x(x + 1, y, z);
		v3 = D1x(x, y, z);
		v4 = D1x(x - 1, y, z);
		v5 = D1x(x - 2, y, z);
	}
	else
	{
		v1 = D1x(x - 3, y, z);
		v2 = D1x(x - 2, y, z);
		v3 = D1x(x - 1, y, z);
		v4 = D1x(x, y, z);
		v5 = D1x(x + 1, y, z);
	}
	gx = CalcWENO(v1, v2, v3, v4, v5);

	//now do y dimension
	if (uy < 0)
	{
		v1 = D1y(x, y + 2, z);
		v2 = D1y(x, y + 1, z);
		v3 = D1y(x, y, z);
		v4 = D1y(x, y - 1, z);
		v5 = D1y(x, y - 2, z);
	}
	else
	{
		v1 = D1y(x, y - 3, z);
		v2 = D1y(x, y - 2, z);
		v3 = D1y(x, y - 1, z);
		v4 = D1y(x, y, z);
		v5 = D1y(x, y + 1, z);
	}
	gy = CalcWENO(v1, v2, v3, v4, v5);

	//now do z dimension
	if (uz < 0)
	{
		v1 = D1z(x, y, z + 2);
		v2 = D1z(x, y, z + 1);
		v3 = D1z(x, y, z);
		v4 = D1z(x, y, z - 1);
		v5 = D1z(x, y, z - 2);
	}
	else
	{
		v1 = D1z(x, y, z - 3);
		v2 = D1z(x, y, z - 2);
		v3 = D1z(x, y, z - 1);
		v4 = D1z(x, y, z);
		v5 = D1z(x, y, z + 1);
	}
	gz = CalcWENO(v1, v2, v3, v4, v5);

	double udotgrad = gx*ux + gy*uy + gz*uz;
	grid2->get(x, y, z) = grid1->get(x, y, z) - udotgrad * timestep;
}

//fix the level set where the particles disagree with it
void LevelSet3D::Fix(const ParticleSet3D& particleSet)
{
	for (int i = 0; i < size[0] - 1; i++)
	{
		for (int j = 0; j < size[1] - 1; j++)
		{
			for (int k = 0; k < size[2] - 1; k++)
			{
				ParticleSet3D::Iterator pit, end = particleSet.end(i, j, k);
				for (pit = particleSet.begin(i, j, k); pit != end; ++pit)
				{
					Particle3D& particle = *(*pit);
					double x, y, z;
					particle.GetPosition(x, y, z);
					double phi = LinearSample(x, y, z);
					double sign = particle.Sign();

					if (phi * sign < 0)
					{
						particle.Radius(phi);
						cout << "fixing" << endl;
						//particle has crossed the boundary
						for (int dx = 0; dx < 2; dx++)
						{
							for (int dy = 0; dy < 2; dy++)
							{
								for (int dz = 0; dz < 2; dz++)
								{
									double particlePhi = particle.phi(i + dx, j + dy, k + dz);
									double& levelSetPhi = grid1->get(i + dx, j + dy, k + dz);
									if (fabs(particlePhi) < fabs(levelSetPhi))
									{
										levelSetPhi = particlePhi;
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

double LevelSet3D::LinearSample(double x, double y, double z) const
{
	double xlerp = x - int(x);
	double ylerp = y - int(y);
	double zlerp = z - int(z);

	int d[] = { 0,1 };
	double samples[2][2][2];

	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			for (int k = 0; k < 2; k++)
				samples[i][j][k] = grid1->get(int(x + i), int(y + j), int(z + k));

	//lerp in x
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			samples[0][i][j] = (1 - xlerp)*samples[0][i][j] + xlerp*samples[1][i][j];

	//lerp in y
	for (int i = 0; i < 2; i++)
		samples[0][0][i] = (1 - ylerp)*samples[0][0][i] + ylerp*samples[0][1][i];

	//lerp in z
	double final = (1 - zlerp)*samples[0][0][0] + zlerp*samples[0][0][1];

	return final;
}

double LevelSet3D::CubicSample(double x, double y, double z) const
{
	//need to do
	return 0;
}

double LevelSet3D::Curvature(double x, double y, double z) const
{
	//need to do
	return 0;
}

void LevelSet3D::normal(double x, double y, double z, double& nx, double& ny, double& nz) const
{
	gradient(x, y, z, nx, ny, nz);
	double mag = std::sqrt(nx*nx + ny*ny + nz*nz);
	if (mag > 0)
	{
		nx /= mag;
		ny /= mag;
		nz /= mag;
	}
}

void LevelSet3D::gradient(double x, double y, double z, double& gx, double& gy, double& gz) const
{
	double dx = 1;
	double dy = 1;
	double dz = 1;
	double center = LinearSample(x, y, z);
	gx = (LinearSample(x + dx, y, z) - center) / dx;
	gy = (LinearSample(x, y + dy, z) - center) / dy;
	gz = (LinearSample(x, y, z + dz) - center) / dz;
}

void LevelSet3D::gradient(double x, double y, double z, double& gx, double& gy, double& gz, double ux, double uy, double uz)
{
	double dx = 1;
	double dy = 1;
	double dz = 1;
	double center = LinearSample(x, y, z);

	if (ux < 0)
	{
		gx = (LinearSample(x + dx, y, z) - center) / dx;
	}
	else
	{
		gx = (center - LinearSample(x - dx, y, z)) / dx;
	}

	if (uy < 0)
	{
		gy = (LinearSample(x, y + dy, z) - center) / dy;
	}
	else
	{
		gy = (center - LinearSample(x, y - dy, z)) / dy;
	}

	if (uz < 0)
	{
		gz = (LinearSample(x, y, z + dz) - center) / dz;
	}
	else
	{
		gz = (center - LinearSample(x, y, z - dz)) / dz;
	}
}

void LevelSet3D::GetVelocityFieldValue(int x, int y, int z, const StableFluid3D& grid, double& ux, double& uy, double& uz)
{
	grid.getvalue(x, y, z, ux, uy, uz);
}

