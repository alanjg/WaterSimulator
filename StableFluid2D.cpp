#include "stdafx.h"
#include "StableFluid2D.h"
#include "SparseMatrix.h"

StableFluid2D::StableFluid2D()
	:viscosity(0),
	density(1),
	width(0),
	height(0),
	resX(0),
	resY(0),
	velocities(NULL),
	tempVelocities(NULL)
{
	initGrids();
}

StableFluid2D::StableFluid2D(double width_, double height_, int resX_, int resY_)
	:viscosity(0),
	density(1),
	width(width_),
	height(height_),
	resX(resX_),
	resY(resY_),
	velocities(NULL),
	tempVelocities(NULL)
{
	initGrids();
	initPoisson();
}

StableFluid2D::StableFluid2D(const StableFluid2D& toCopy)
{
	*this = toCopy;
}

StableFluid2D::~StableFluid2D()
{
	clearGrids();
}

StableFluid2D& StableFluid2D::operator =(const StableFluid2D& toCopy)
{
	viscosity = toCopy.viscosity;
	density = toCopy.density;
	if ((resX != toCopy.resX) || (resY != toCopy.resY))
	{
		clearPoisson();
		resX = toCopy.resX;
		resY = toCopy.resY;
		laplacian = toCopy.laplacian;
		preconditioner = toCopy.preconditioner;

		vel = toCopy.vel;
		r = toCopy.r;
		z = toCopy.z;
	}
	ones = toCopy.ones;
	pressure = toCopy.pressure;
	tempPressure = toCopy.tempPressure;

	if ((width != toCopy.width) || (height != toCopy.height))
	{
		clearGrids();
		width = toCopy.width;
		height = toCopy.height;
		initGrids();
	}

	for (int d = 0; d < 2; ++d)
	{
		for (int i = 0; i < resX + xOffset2[d]; ++i)
		{
			for (int j = 0; j < resY + yOffset2[d]; ++j)
			{
				tempVelocities[d][i][j] = tempVelocities[d][i][j];
				velocities[d][i][j] = toCopy.velocities[d][i][j];
			
			}
		}
	}
	return *this;
}

void StableFluid2D::step(double dt)
{
	// solve velocity
	std::swap(velocities, tempVelocities);
	advectVelocity(dt);
	std::swap(velocities, tempVelocities);

	addForcesToVelocity(dt);
	projectVelocity(dt);

}

void StableFluid2D::getvalue(double i, double j, double& ux, double& uy) const
{
	vec2d sample = linearSamp(i, j, velocities);
	ux = sample[0];
	uy = sample[1];
}

vec2d StableFluid2D::linearSamp(double x, double y, double*** v) const
{
	x *= width / (double)resX;
	y *= height / (double)resY;
	x = (x < 0 ? 0 : (x > width ? width : x));
	y = (y < 0 ? 0 : (y > height ? height : y));

	int x0 = (int)(x == width ? width - 1 : x);
	double xAlpha = x - x0;

	int y0 = (int)(y == height ? height - 1 : y);
	double yAlpha = y - y0;

	return vec2d((1 - xAlpha) * v[0][x0][y0] +
		xAlpha * v[0][x0 + 1][y0],
		(1 - yAlpha) * v[1][x0][y0] +
		yAlpha * v[1][x0][y0 + 1]);
}

vec2d StableFluid2D::linearSamp(vec2d& pos, double*** v) const
{
	return linearSamp(pos[0], pos[1], v);
}

double StableFluid2D::linearSamp(double x, double y, Vector& v) const
{
	x *= width / (double)resX;
	y *= height / (double)resY;
	x = (x < 0 ? 0 : (x > width ? width : x));
	y = (y < 0 ? 0 : (y > height ? height : y));

	int x0 = (int)(x == width ? width - 1 : x);
	double xAlpha = x - x0;

	int y0 = (int)(y == height ? height - 1 : y);
	double yAlpha = y - y0;


	double dx = (1 - xAlpha) * v[makeIndex(x0, y0)] + xAlpha * v[makeIndex(x0 + 1, y0)];
	double dy = (1 - yAlpha) * v[makeIndex(x0, y0)] + yAlpha * v[makeIndex(x0, y0 + 1)];
	return (dx + dy) / 2.0;
}

double StableFluid2D::linearSamp(vec2d& pos, Vector& v) const
{
	return linearSamp(pos[0], pos[1], v);
}

void StableFluid2D::setVisc(double viscosity_)
{
	viscosity = viscosity_;
}

void StableFluid2D::setDensity(double density_)
{
	density = density_;
}

void StableFluid2D::setSize(int width_, int height_)
{
	clearGrids();
	width = width_;
	height = height_;
	initGrids();
	initPoisson();
}

void StableFluid2D::setRes(int resX_, int resY_)
{
	clearGrids();
	resX = resX_;
	resY = resY_;
	initGrids();
	initPoisson();
}

void StableFluid2D::initGrids()
{
	clearGrids();
	velocities = new double**[3];
	tempVelocities = new double**[3];
	for (int d = 0; d < 2; ++d)
	{
		velocities[d] = new double*[resX + xOffset2[d]];
		tempVelocities[d] = new double*[resX + xOffset2[d]];
		for (int i = 0; i < resX + xOffset2[d]; ++i)
		{
			velocities[d][i] = new double[resY + yOffset2[d]];
			tempVelocities[d][i] = new double[resY + yOffset2[d]];
			for (int j = 0; j < resY + yOffset2[d]; ++j)
			{
				velocities[d][i][j] = 0.0;
				tempVelocities[d][i][j] = 0.0;	
			}
		}
	}
}

void StableFluid2D::clearGrids()
{
	if (!velocities && !tempVelocities)
		return;

	for (int d = 0; d < 2; ++d)
	{
		for (int i = 0; i < resX + xOffset2[d]; ++i)
		{
			delete[] velocities[d][i];
			delete[] tempVelocities[d][i];
		}
		delete[] velocities[d];
		delete[] tempVelocities[d];
	}
	delete[] velocities;
	delete[] tempVelocities;

	velocities = NULL;
	tempVelocities = NULL;
}

void StableFluid2D::initPoisson()
{
	clearPoisson();
	int n = resX * resY;

	for (int i = 0; i < 2; ++i)
	{
		r.push_back(Vector(n));
		z.push_back(Vector(n));
	}

	p.Resize(n);

	SparseMatrix hTranspose;
	SparseMatrix h;

	laplacian = SparseMatrix(n, n);
	hTranspose = SparseMatrix(n, n);
	h = SparseMatrix(n, n);
	preconditioner = SparseMatrix(n, n);

	pressure.Resize(n);
	tempPressure.Resize(n);
	vel.Resize(n);
	ones.Resize(n);

	for (int i = 0; i < n; ++i)
	{
		ones[i] = 1;
	}

	double rdxSqr = (resX * resX) / (width * width);
	double rdySqr = (resY * resY) / (height * height);

	// initialize values for the laplacian matrix
	for (int i = 0; i < resX; ++i)
	{
		for (int j = 0; j < resY; ++j)
	{
			int row = makeIndex(i, j);

			int x1 = (i == (resX - 1) ? row : row + 1);
			int x2 = (i == 0 ? row : row - 1);

			int y1 = (j == (resY - 1) ? row : row + resX);
			int y2 = (j == 0 ? row : row - resX);

			laplacian[row][row] += 4 * (rdxSqr + rdySqr);
			laplacian[row][x1] -= rdxSqr;
			laplacian[row][x2] -= rdxSqr;
			laplacian[row][y1] -= rdySqr;
			laplacian[row][y2] -= rdySqr;
		}
	}

	// perform Cholesky preconditioning
	// decompose laplacian into H * H^T
	SparseVector::iterator iter;

	hTranspose = laplacian.Transpose();
	SparseMatrix& A = hTranspose;
	for (int k = 0; k < n; ++k)
	{
		SparseVector::iterator aik(A[k]);
		aik.seek(k);

		if (aik.getIndex() != k)
		{
			assert(0);
		}

		double divisor = std::sqrt(*aik);

		*aik = divisor;
		aik++;

		while (aik)
		{
			*aik /= divisor;
			aik++;
		}

		SparseVector::iterator ajk(A[k]);
		ajk.seek(k + 1);

		while (ajk)
		{
			int j = ajk.getIndex();
			aik.reset();
			aik.seek(j);

			SparseVector::iterator aij(A[j]);
			aij.seek(j);


			while (aij && aik)
			{
				int i = aij.getIndex();
				if (aik.getIndex() == i)
				{
					*aij -= (*aik) * (*ajk);
					aik++;
				}

				aij++;
				if (aij)
				{
					i = aij.getIndex();

					while (aik && aik.getIndex() < i)
					{
						aik++;
					}
				}
			}
			ajk++;
		}
	}

	h = hTranspose.Transpose();

	//	cerr<<"Sparsity of laplacian is " << laplacian.Sparsity() << endl;

	preconditioner = h.FastMultiplyTranspose();

	int maxIndex = resX * resY;
	//	for (int i = 0; i < maxIndex; ++i)
	//		preconditioner[i][i] = 1;

	//	cerr<<"Sparsity of preconditioner is " << preconditioner.Sparsity() << endl;
}

void StableFluid2D::clearPoisson()
{
	laplacian.clear();
	preconditioner.clear();

	vel.clear();
	pressure.clear();
	tempPressure.clear();
	r.clear();
	z.clear();
	p.clear();
}

int StableFluid2D::makeIndex(int x, int y) const
{
	return (x + resX * y);
}

void StableFluid2D::splitIndex(int index, int& x, int& y)
{
	x = index % resX;
	y = index / resX;
}

void StableFluid2D::setBoundary()
{
	// top and bottom boundaries
	for (int j = 0; j < resX; ++j)
	{
		tempVelocities[1][j][0] = 0;
		tempVelocities[1][j][resY] = 0;
	}
	

	// left and right boundaries
	for (int i = 0; i < resY; ++i)
	{
		tempVelocities[0][0][i] = 0;
		tempVelocities[0][resX][i] = 0;
	}
}

void StableFluid2D::addForcesToVelocity(double dt)
{
	for (int i = 0; i < resX; ++i)
	{
		for (int j = 0; j <= resY; ++j)
		{
			tempVelocities[1][i][j] -= (double)(GRAV_ACCEL * dt);
		}
	}

	if (viscosity != 0)
	{
		double rdx = resX / width;
		double rdy = resY / height;

		for (int i = 0; i <= resX; ++i)
		{
			for (int j = 0; j < resY; ++j)
			{
				double v1 = (i == 0 ? tempVelocities[0][i][j] : tempVelocities[0][i - 1][j]);
				double v2 = (i == resX ? tempVelocities[0][i][j] : tempVelocities[0][i + 1][j]);
				tempVelocities[0][i][j] -= viscosity * (v1 - 2 * tempVelocities[0][i][j] + v2) * rdx;
				
			}
		}
		for (int i = 0; i < resX; ++i)
		{
			for (int j = 0; j <= resY; ++j)
			{
				double v1 = (j == 0 ? tempVelocities[1][i][j] : tempVelocities[1][i][j - 1]);
				double v2 = (j == resY ? tempVelocities[1][i][j] : tempVelocities[1][i][j + 1]);
				tempVelocities[1][i][j] -= viscosity * (v1 - 2 * tempVelocities[1][i][j] + v2) * rdy;
			}
		}
	}
}

void StableFluid2D::advectVelocity(double dt)
{
	double xStep = width / (double)resX;
	double yStep = height / (double)resY;

	double fOffsetsX[2] = { 0, yStep * 0.5 };
	double fOffsetsY[2] = { xStep * 0.5, 0 };

	vec2d pos;
	vec2d v;

	for (int d = 0; d < 2; ++d)
	{
		pos[0] = fOffsetsX[d];
		for (int i = 0; i < resX + xOffset2[d]; ++i)
		{
			pos[0] += xStep;
			pos[1] = fOffsetsY[d];
			for (int j = 0; j < resY + yOffset2[d]; ++j)
			{
				pos[1] += yStep;
				v = linearSamp(pos, tempVelocities);
				velocities[d][i][j] = linearSamp(pos - v * dt, tempVelocities)[d];
			}
		}
	}
}

void StableFluid2D::projectVelocity(double dt)
{
	setBoundary();

	Vector divergence(vel.Length());
	// fill in velocities
	double dx = width / (double)resX;
	double dy = height / (double)resY;
	for (int i = 0; i < resX; ++i)
	{
		for (int j = 0; j < resY; ++j)
		{
			int index = makeIndex(i, j);

			double div = (tempVelocities[0][i + 1][j] - tempVelocities[0][i][j]) / dx + (tempVelocities[1][i][j + 1] - tempVelocities[1][i][j]) / dy;

			vel[index] = div;
			divergence[index] = div;
			
		}
	}

	//	vel *= density / dt;
	bool useOldPCG = false;
	if (useOldPCG)
	{
		int k = 0;
		Vector laplacianTimesPressure(pressure.Length());
		laplacian.Multiply(pressure, laplacianTimesPressure);
		r[0] = vel - laplacianTimesPressure;
		double error = r[0].Magnitude2();
		while ((std::sqrt(error) > PCG_EPS) && (k < PCG_MAXITER))
		{
			//cerr << "PCG iteration " << k << ", error: " << std::sqrt(error) << endl;
			ConjugateGradient(preconditioner, z[k % 2], r[k % 2], 0.0001, 500.0);

			++k;
			int i1 = (k - 1) % 2;
			int i2 = k % 2;
			if (k == 1)
			{
				p = z[0];
			}
			else
			{
				double beta = (InnerProduct(r[i1], z[i1])) /
					(InnerProduct(r[i2], z[i2]));
				p *= beta;
				p += z[i1];
			}
			Vector laplacianTimesP(p.Length());
			laplacian.Multiply(p, laplacianTimesP);
			double alpha = (InnerProduct(r[i1], z[i1])) /
				(InnerProduct(p, laplacianTimesP));
			pressure += alpha * p;
			r[i2] = r[i1] - alpha * (laplacianTimesP);
			error = r[i2].Magnitude2();
		}
		if (std::sqrt(error) > PCG_EPS)
		{
			std::cout << "end preconditioned conj grad " << "error: " << error << std::endl;
		}
	}
	else
	{
		p = Vector(vel.Length());
		r[0] = divergence;
		//ApplyPreconditioner(r[0], z[0]);
		double sigma = InnerProduct(r[0], z[0]);

		int k = 0;
		double error = r[0].Magnitude2();
		while (k < PCG_MAXITER)
		{
			++k;
			Vector s(vel.Length());
			s = z[0];
			int k = 0;
			laplacian.Multiply(s, z[0]);
			double rho = 1.0;
			double alpha = rho / InnerProduct(s, z[0]);
			p += alpha * s;
			r[0] -= alpha * z[0];
			error = r[0].Magnitude2();
			if (error < PCG_EPS * PCG_EPS)
			{
				pressure = p;
				break;
			}
			//ApplyPreconditioner(r[0], z[0]);
			double sigmaNew = InnerProduct(r[0], z[0]);
			double beta = sigmaNew / rho;
			s = z[0] + beta * s;
			sigma = sigmaNew;
		}
		Vector laplacianTimesPressure(pressure.Length());
		laplacian.Multiply(pressure, laplacianTimesPressure);
		r[0] = vel - laplacianTimesPressure;
		error = r[0].Magnitude2();
		while ((std::sqrt(error) > PCG_EPS) && (k < PCG_MAXITER))
		{
			//cerr << "PCG iteration " << k << ", error: " << std::sqrt(error) << endl;
			ConjugateGradient(preconditioner, z[k % 2], r[k % 2], 0.0001, 500.0);

			++k;
			int i1 = (k - 1) % 2;
			int i2 = k % 2;
			if (k == 1)
			{
				p = z[0];
			}
			else
			{
				double beta = (InnerProduct(r[i1], z[i1])) /
					(InnerProduct(r[i2], z[i2]));
				p *= beta;
				p += z[i1];
			}
			Vector laplacianTimesP(p.Length());
			laplacian.Multiply(p, laplacianTimesP);
			double alpha = (InnerProduct(r[i1], z[i1])) /
				(InnerProduct(p, laplacianTimesP));
			pressure += alpha * p;
			r[i2] = r[i1] - alpha * (laplacianTimesP);
			error = r[i2].Magnitude2();
		}
		if (std::sqrt(error) > PCG_EPS)
		{
			std::cout << "end preconditioned conj grad " << "error: " << error << std::endl;
		}
	}

	/*
	mathlib::Vector PCGr, PCGd, PCGq, PCGs;
	PCGr.Resize(resX * resY * resZ);
	PCGd.Resize(resX * resY * resZ);
	PCGq.Resize(resX * resY * resZ);
	PCGs.Resize(resX * resY * resZ);

	int iter = 0;
	PCGr = vel - laplacian * pressure;
	mathlib::ConjugateGradient(preconditioner, PCGd, PCGr);
	double deltaNew = mathlib::InnerProduct(PCGr, PCGd);
	double errBound = PCG_EPS * PCG_EPS * deltaNew;

	while ((iter < PCG_MAXITER)  && (deltaNew > errBound))
	{
	if (!(iter % 3))
	cerr << "iteration " << iter << ", error: " << deltaNew << endl;
	PCGq = laplacian * PCGd;
	double alpha = deltaNew / (mathlib::InnerProduct(PCGd, PCGq));
	pressure += alpha * PCGd;
	if ((iter % 10) == 0)
	PCGr = vel - laplacian * pressure;
	else
	PCGr -= alpha * PCGq;

	mathlib::ConjugateGradient(preconditioner, PCGs, PCGr);
	double deltaOld = deltaNew;
	deltaNew = mathlib::InnerProduct(PCGr, PCGs);
	double beta = deltaNew / deltaOld;
	PCGd *= beta;
	PCGd += PCGs;
	++iter;
	}
	cout << "end preconditioned conj grad "<< "error: " << deltaNew << endl;
	*/
	std::cout << "dt:" << dt << "\t density:" << density << std::endl;
	//	dt = 1.0;
	for (int i = 0; i < resX; ++i)
	{
		for (int j = 0; j < resY; ++j)
		{
			int i1 = makeIndex(i, j);
			int i2;

			double deltaPressure;
			i2 = (i == 0 ? i1 : i1 - 1);
			velocities[0][i][j] = tempVelocities[0][i][j] - dt* 0.5 * (pressure[i1] - pressure[i2]) / dx;
			deltaPressure = pressure[i1] - pressure[i2];

			i2 = (j == 0 ? i1 : i1 - resX);
			velocities[1][i][j] = tempVelocities[1][i][j] - dt *  0.5 * (pressure[i1] - pressure[i2]) / dy;
			deltaPressure = pressure[i1] - pressure[i2];			
		}
	}
}