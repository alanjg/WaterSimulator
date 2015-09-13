#include "stdafx.h"
#include "LevelSet2D.h"

#include "ParticleSet2D.h"
#include "Particle2D.h"
#include "StableFluid2D.h"
#include "LevelSetGrid2D.h"
#include "LevelSetFMGrid2D.h"
#include "math.h"

using std::cout;
using std::endl;
using std::min;

LevelSet2D::LevelSet2D(int xResolution, int yResolution)
{
	//level set rep is on the corners, extend the grid by 1
	size[0] = xResolution + 1;
	size[1] = yResolution + 1;

	grid1 = new LevelSetGrid2D(size[0], size[1]);
	grid2 = new LevelSetGrid2D(size[0], size[1]);
	grid3 = new LevelSetGrid2D(size[0], size[1]);
	gridFM = new LevelSetFMGrid2D(size[0], size[1]);
}

LevelSet2D::~LevelSet2D()
{
	delete grid1;
	delete grid2;
	delete grid3;
	delete gridFM;
}

/*
double LevelSet2D::eval(const Point3d& location)
{
double x = (location[0] + 1) * (size[0] - 1) / 2.0;
double y = (location[1] + 1) * (size[1] - 1) / 2.0;
double z = (location[2] + 1) * (size[2] - 1) / 2.0;
double d = LinearSample(x, y, z);

return d;
}
*/
void LevelSet2D::Write(std::ostream& out)
{
	grid1->Write(out);
}

void LevelSet2D::CheckGrid()
{
	for (int i = 0; i <= size[0]; i++)
	{
		for (int j = 0; j <= size[1]; j++)
		{
			
			double me = grid1->get(i, j);
			if (i > 0)
			{
				double other = grid1->get(i - 1, j);
				if (abs(me - other) > 1.5)
				{
					cout << "Bad grid spot at " << i << " " << j << " " << me << " " << other << endl;
				}
			}

			if (j > 0)
			{
				double other = grid1->get(i, j - 1);
				if (abs(me - other) > 1.5)
				{
					cout << "Bad grid spot at " << i << " " << j << " " << me << " " << other << endl;
				}
			}			
		}
	}
}


void LevelSet2D::MakeCircle(double x, double y, double radius)
{
	for (int i = 0; i <= size[0]; i++)
	{
		for (int j = 0; j <= size[1]; j++)
		{			
			double val1 = std::sqrt(square(x - i) + square(y - j)) - radius;

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

			grid1->get(i, j) = val1;
			//			grid1->get(i,j,k) = std::max(val1,val2);
			//				grid1->get(i,j,k) = j - 15;
			
		}
	}
	grid1->UpdateBorders();
}

void LevelSet2D::UpdateWorker(WorkerData data)
{
	//Find phi(n+1)
	for (int i = data.iBegin; i < data.iEnd; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{			
			//				EulerStep(i,j,grid,timestep);
			//				ENOStep(i,j,grid,timestep);
			WENOStep(i, j, *data.field, data.timestep);
		}
	}
}

void LevelSet2D::UpdateThreaded(const StableFluid2D& f0, const StableFluid2D& f1, const StableFluid2D& f2, double timestep)
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
		threads.push_back(std::thread(&LevelSet2D::UpdateWorker, this, d));
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
		threads.push_back(std::thread(&LevelSet2D::UpdateWorker, this, d));
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
			grid1->get(i, j) = grid3->get(i, j)*0.75 + grid2->get(i, j)*0.25;
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
		threads.push_back(std::thread(&LevelSet2D::UpdateWorker, this, d));
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
			grid1->get(i, j) = grid3->get(i, j)*0.333333333 + grid2->get(i, j)*0.666666667;
		
		}
	}

	grid1->UpdateBorders();
}

void LevelSet2D::Update(const StableFluid2D& f0, const StableFluid2D& f1, const StableFluid2D& f2, double timestep)
{
	//TVD RK THREE - see page 38 of osher & fedkiw

	//Find phi(n+1)
	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{
			//				EulerStep(i,j,grid,timestep);
			//				ENOStep(i,j,grid,timestep);
			WENOStep(i, j, f0, timestep);	
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
			//				EulerStep(i,j,grid,timestep);
			//				ENOStep(i,j,grid,timestep);
			WENOStep(i, j, f2, timestep);
		}
	}

	grid2->UpdateBorders();

	//Set grid1 to phi(n+1/2) = .75*phi(n) + .25*phi(n+2)
	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{
			grid1->get(i, j) = grid3->get(i, j)*0.75 + grid2->get(i, j)*0.25;
		}
	}

	grid1->UpdateBorders();

	//Find phi(n+3/2)
	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{
			//				EulerStep(i,j,grid,timestep);
			//				ENOStep(i,j,grid,timestep);
			WENOStep(i, j, f1, timestep);
		}
	}

	grid2->UpdateBorders();

	//Set grid1 to final phi(n+1) = phi(n) * 1/3 + phi(n+3/2) * 2/3
	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{
			grid1->get(i, j) = grid3->get(i, j)*0.333333333 + grid2->get(i, j)*0.666666667;
		}
	}

	grid1->UpdateBorders();
}

//WE NEED ADVECTION ALONG NORMAL DIRECTION. SAME AS WENO, WITHOUGH VELOCITY GRID.
//USE NORMAL FOR VELOCITY. USE EQUATION 7.4 INSTEAD OF 3.2

void LevelSet2D::Reinitialize() 
{
	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{
			gridFM->Set(i, j, grid1->get(i, j));
			//if(grid1->get(i,j) > -2 && grid1->get(i,j) < 2 && j < 10)
			//	cout << i << "   " << j << "   " <<
			//		grid1->get(i,j) << endl;	
		}
	}
	//cout << endl << endl;
	gridFM->Reinitialize();
	//cout << endl << endl;
	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{
			
			grid1->get(i, j) = gridFM->Get(i, j);
			//if(grid1->get(i,j) > 0 && grid1->get(i,j) < 1)
			//	cout << i << "   " << j << "   " <<
			//		grid1->get(i,j) << endl;
		}
	}

	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{

			gridFM->Set(i, j, grid1->get(i, j)*(-1));
			//if(grid1->get(i,j)*(-1) > 0)
			//	cout << i << "   " << j << "   " <<
			//		grid1->get(i,j) << endl;
		}
	}
	//cout << endl << endl;
	gridFM->Reinitialize();
	//cout << endl << endl;
	for (int i = 0; i < size[0]; i++)
	{
		for (int j = 0; j < size[1]; j++)
		{
			grid1->get(i, j) = gridFM->Get(i, j)*(-1);
			//if(grid1->get(i,j)*(-1) > 0)
			//	cout << i << "   " << j << "   " <<
			//		grid1->get(i,j) << endl;
		}
	}
	grid1->UpdateBorders();
}

void LevelSet2D::EulerStep(int x, int y, const StableFluid2D& grid, double timestep)
{
	double ux, uy; //obtain from velocity grid

	GetVelocityFieldValue(x, y, grid, ux, uy);

	double gx, gy;
	gradient(x, y, gx, gy, ux, uy);

	double udotgrad = gx*ux + gy*uy;
	grid2->get(x, y) = grid1->get(x, y) - udotgrad * timestep;
}

double LevelSet2D::D0(int x, int y)
{
	return grid1->get(x, y);
}

double LevelSet2D::D1x(int x, int y)
{
	return (D0(x + 1, y) - D0(x, y)) / GRID_SIZE;
}

double LevelSet2D::D1y(int x, int y)
{
	return (D0(x, y + 1) - D0(x, y)) / GRID_SIZE;
}

double LevelSet2D::D1(int x, int y, int dim)
{
	int dx[] = { 1,0 };
	int dy[] = { 0,1 };
	return (D0(x + dx[dim], y + dy[dim]) - D0(x, y)) / GRID_SIZE;
}

double LevelSet2D::D2x(int x, int y)
{
	return (D1x(x, y) - D1x(x - 1, y)) / (2 * GRID_SIZE);
}

double LevelSet2D::D2y(int x, int y)
{
	return (D1x(x, y) - D1x(x, y - 1)) / (2 * GRID_SIZE);
}

double LevelSet2D::D2(int x, int y, int dim)
{
	int dx[] = { 1,0 };
	int dy[] = { 0,1 };
	return (D1(x, y, dim) - D1(x - dx[dim], y - dy[dim], dim)) / (2 * GRID_SIZE);
}

double LevelSet2D::D3x(int x, int y)
{
	return (D2x(x + 1, y) - D2x(x, y)) / (3 * GRID_SIZE);
}

double LevelSet2D::D3y(int x, int y)
{
	return (D2x(x, y + 1) - D2x(x, y)) / (3 * GRID_SIZE);
}

double LevelSet2D::D3(int x, int y, int dim)
{
	int dx[] = { 1,0 };
	int dy[] = { 0,1 };
	return (D2(x + dx[dim], y + dy[dim], dim) - D2(x, y, dim)) / (3 * GRID_SIZE);
}

void LevelSet2D::ENOStep(int x, int y, const StableFluid2D& grid, double timestep)
{
	double ux, uy; //obtain from velocity grid

	GetVelocityFieldValue(x, y, grid, ux, uy);

	//we solve for these, see page 32 in Osher & Fedkiw
	double gx, gy;

	int pos[2] = { x,y };

	int dx[] = { 1,0 };
	int dy[] = { 0,1 };

	double u[2] = { ux,uy };
	double g[2];
	for (int i = 0; i < 2; i++)
	{
		int index[2] = { x,y };
		int k;
		if (u[i] > 0)
		{
			index[i]--;
		}
		k = index[i];
		double q1 = D1(index[0], index[1], i);

		double dk0 = D2(index[0], index[1], i);
		double dk1 = D2(index[0] + dx[i], index[1] + dy[i], i);

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

		dk0 = D3(index[0], index[1], i);
		dk1 = D3(index[0] + dx[i], index[1] + dy[i], i);
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
	double udotgrad = gx*ux + gy*uy;
	grid2->get(x, y) = grid1->get(x, y) - udotgrad * timestep;
}

void LevelSet2D::WENOStep(int x, int y, const StableFluid2D& grid, double timestep)
{
	double ux, uy; //obtain from velocity grid

	GetVelocityFieldValue(x, y, grid, ux, uy);

	//we solve for these, see page 34 in Osher & Fedkiw
	double gx, gy;

	double v1, v2, v3, v4, v5;
	if (ux < 0)
	{
		v1 = D1x(x + 2, y);
		v2 = D1x(x + 1, y);
		v3 = D1x(x, y);
		v4 = D1x(x - 1, y);
		v5 = D1x(x - 2, y);
	}
	else
	{
		v1 = D1x(x - 3, y);
		v2 = D1x(x - 2, y);
		v3 = D1x(x - 1, y);
		v4 = D1x(x, y);
		v5 = D1x(x + 1, y);
	}
	gx = CalcWENO(v1, v2, v3, v4, v5);

	//now do y dimension
	if (uy < 0)
	{
		v1 = D1y(x, y + 2);
		v2 = D1y(x, y + 1);
		v3 = D1y(x, y);
		v4 = D1y(x, y - 1);
		v5 = D1y(x, y - 2);
	}
	else
	{
		v1 = D1y(x, y - 3);
		v2 = D1y(x, y - 2);
		v3 = D1y(x, y - 1);
		v4 = D1y(x, y);
		v5 = D1y(x, y + 1);
	}
	gy = CalcWENO(v1, v2, v3, v4, v5);

	double udotgrad = gx*ux + gy*uy;
	grid2->get(x, y) = grid1->get(x, y) - udotgrad * timestep;
}

//fix the level set where the particles disagree with it
void LevelSet2D::Fix(const ParticleSet2D& particleSet)
{
	for (int i = 0; i < size[0] - 1; i++)
	{
		for (int j = 0; j < size[1] - 1; j++)
		{
			ParticleSet2D::Iterator pit, end = particleSet.end(i, j);
			for (pit = particleSet.begin(i, j); pit != end; ++pit)
			{
				Particle2D& particle = *(*pit);
				double x, y;
				particle.GetPosition(x, y);
				double phi = LinearSample(x, y);
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
							double particlePhi = particle.phi(i + dx, j + dy);
							double& levelSetPhi = grid1->get(i + dx, j + dy);
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

double LevelSet2D::LinearSample(double x, double y) const
{
	double xlerp = x - int(x);
	double ylerp = y - int(y);

	int d[] = { 0,1 };
	double samples[2][2];

	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			samples[i][j] = grid1->get(int(x + i), int(y + j));

	//lerp in x
	for (int i = 0; i < 2; i++)
		samples[0][i] = (1 - xlerp)*samples[0][i] + xlerp*samples[1][i];

	//lerp in y
	double final = (1 - ylerp)*samples[0][0] + ylerp*samples[0][1];

	return final;
}

double LevelSet2D::CubicSample(double x, double y) const
{
	//need to do
	return 0;
}

double LevelSet2D::Curvature(double x, double y) const
{
	//need to do
	return 0;
}

void LevelSet2D::normal(double x, double y, double& nx, double& ny) const
{
	gradient(x, y, nx, ny);
	double mag = std::sqrt(nx*nx + ny*ny);
	if (mag > 0)
	{
		nx /= mag;
		ny /= mag;
	}
}

void LevelSet2D::gradient(double x, double y, double& gx, double& gy) const
{
	double dx = 1;
	double dy = 1;
	double center = LinearSample(x, y);
	gx = (LinearSample(x + dx, y) - center) / dx;
	gy = (LinearSample(x, y + dy) - center) / dy;
}

void LevelSet2D::gradient(double x, double y, double& gx, double& gy, double ux, double uy)
{
	double dx = 1;
	double dy = 1;
	double center = LinearSample(x, y);

	if (ux < 0)
	{
		gx = (LinearSample(x + dx, y) - center) / dx;
	}
	else
	{
		gx = (center - LinearSample(x - dx, y)) / dx;
	}

	if (uy < 0)
	{
		gy = (LinearSample(x, y + dy) - center) / dy;
	}
	else
	{
		gy = (center - LinearSample(x, y - dy)) / dy;
	}

}

void LevelSet2D::GetVelocityFieldValue(int x, int y, const StableFluid2D& grid, double& ux, double& uy)
{
	grid.getvalue(x, y, ux, uy);
}

