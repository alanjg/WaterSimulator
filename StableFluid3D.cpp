#include "stdafx.h"
#include "StableFluid3D.h"
#include "SparseMatrix.h"

StableFluid3D::StableFluid3D()
	:viscosity(0),
	density(1),
	width(0),
	height(0),
	depth(0),
	resX(0),
	resY(0),
	resZ(0),
	velocities(NULL),
	tempVelocities(NULL)
{
	initGrids();
}

StableFluid3D::StableFluid3D(double width_, double height_, double depth_, int resX_, int resY_, int resZ_)
	:viscosity(0),
	density(1),
	width(width_),
	height(height_),
	depth(depth_),
	resX(resX_),
	resY(resY_),
	resZ(resZ_),
	velocities(NULL),
	tempVelocities(NULL)
{
	initGrids();
	initPoisson();
}

StableFluid3D::StableFluid3D(const StableFluid3D& toCopy)
{
	*this = toCopy;
}

StableFluid3D::~StableFluid3D()
{
	clearGrids();
}

StableFluid3D& StableFluid3D::operator =(const StableFluid3D& toCopy)
{
	viscosity = toCopy.viscosity;
	density = toCopy.density;
	if ((resX != toCopy.resX) || (resY != toCopy.resY) || (resZ != toCopy.resZ))
	{
		clearPoisson();
		resX = toCopy.resX;
		resY = toCopy.resY;
		resZ = toCopy.resZ;
		laplacian = toCopy.laplacian;
		preconditioner = toCopy.preconditioner;

		vel = toCopy.vel;
		r = toCopy.r;
		z = toCopy.z;
	}
	ones = toCopy.ones;
	pressure = toCopy.pressure;
	tempPressure = toCopy.tempPressure;

	if ((width != toCopy.width) || (height != toCopy.height) || (depth != toCopy.depth))
	{
		clearGrids();
		width = toCopy.width;
		height = toCopy.height;
		depth = toCopy.depth;
		initGrids();
	}

	for (int d = 0; d < 3; ++d)
	{
		for (int i = 0; i < resX + xOffset3[d]; ++i)
		{
			for (int j = 0; j < resY + yOffset3[d]; ++j)
			{
				for (int k = 0; k < resZ + zOffset3[d]; ++k)
				{
					tempVelocities[d][i][j][k] = tempVelocities[d][i][j][k];
					velocities[d][i][j][k] = toCopy.velocities[d][i][j][k];
				}
			}
		}
	}
	return *this;
}

void StableFluid3D::step(double dt)
{
	// solve velocity
	std::swap(velocities, tempVelocities);
	advectVelocity(dt);
	std::swap(velocities, tempVelocities);

	addForcesToVelocity(dt);
	projectVelocity(dt);

	// solve pressure

	//	std::swap(pressure, tempPressure);
	//	advectPressure(dt);
}

void StableFluid3D::getvalue(double i, double j, double k, double& ux, double& uy, double& uz) const
{
	vec3d sample = linearSamp(i, j, k, velocities);
	ux = sample[0];
	uy = sample[1];
	uz = sample[2];
}

vec3d StableFluid3D::linearSamp(double x, double y, double z, double**** v) const
{
	x *= width / (double)resX;
	y *= height / (double)resY;
	z *= depth / (double)resZ;
	x = (x < 0 ? 0 : (x > width ? width : x));
	y = (y < 0 ? 0 : (y > height ? height : y));
	z = (z < 0 ? 0 : (z > depth ? depth : z));

	int x0 = (int)(x == width ? width - 1 : x);
	double xAlpha = x - x0;

	int y0 = (int)(y == height ? height - 1 : y);
	double yAlpha = y - y0;

	int z0 = (int)(z == depth ? depth - 1 : z);
	double zAlpha = z - z0;

	return vec3d((1 - xAlpha) * v[0][x0][y0][z0] +
		xAlpha * v[0][x0 + 1][y0][z0],
		(1 - yAlpha) * v[1][x0][y0][z0] +
		yAlpha * v[1][x0][y0 + 1][z0],
		(1 - zAlpha) * v[2][x0][y0][z0] +
		zAlpha * v[2][x0][y0][z0 + 1]);
}

vec3d StableFluid3D::linearSamp(vec3d& pos, double**** v) const
{
	return linearSamp(pos[0], pos[1], pos[2], v);
}

double StableFluid3D::linearSamp(double x, double y, double z, Vector& v) const
{
	x *= width / (double)resX;
	y *= height / (double)resY;
	z *= depth / (double)resZ;
	x = (x < 0 ? 0 : (x > width ? width : x));
	y = (y < 0 ? 0 : (y > height ? height : y));
	z = (z < 0 ? 0 : (z > depth ? depth : z));

	int x0 = (int)(x == width ? width - 1 : x);
	double xAlpha = x - x0;

	int y0 = (int)(y == height ? height - 1 : y);
	double yAlpha = y - y0;

	int z0 = (int)(z == depth ? depth - 1 : z);
	double zAlpha = z - z0;

	double dx = (1 - xAlpha) * v[makeIndex(x0, y0, z0)] + xAlpha * v[makeIndex(x0 + 1, y0, z0)];
	double dy = (1 - yAlpha) * v[makeIndex(x0, y0, z0)] + yAlpha * v[makeIndex(x0, y0 + 1, z0)];
	double dz = (1 - zAlpha) * v[makeIndex(x0, y0, z0)] + zAlpha * v[makeIndex(x0, y0, z0 + 1)];
	return (dx + dy + dz) / 3.0;
}

double StableFluid3D::linearSamp(vec3d& pos, Vector& v) const
{
	return linearSamp(pos[0], pos[1], pos[2], v);
}

void StableFluid3D::setVisc(double viscosity_)
{
	viscosity = viscosity_;
}

void StableFluid3D::setDensity(double density_)
{
	density = density_;
}

void StableFluid3D::setSize(int width_, int height_, int depth_)
{
	clearGrids();
	width = width_;
	height = height_;
	depth = depth_;
	initGrids();
	initPoisson();
}

void StableFluid3D::setRes(int resX_, int resY_, int resZ_)
{
	clearGrids();
	resX = resX_;
	resY = resY_;
	resZ = resZ_;
	initGrids();
	initPoisson();
}

void StableFluid3D::initGrids()
{
	clearGrids();
	velocities = new double***[3];
	tempVelocities = new double***[3];
	for (int d = 0; d < 3; ++d)
	{
		velocities[d] = new double**[resX + xOffset3[d]];
		tempVelocities[d] = new double**[resX + xOffset3[d]];
		for (int i = 0; i < resX + xOffset3[d]; ++i)
		{
			velocities[d][i] = new double*[resY + yOffset3[d]];
			tempVelocities[d][i] = new double*[resY + yOffset3[d]];
			for (int j = 0; j < resY + yOffset3[d]; ++j)
			{
				velocities[d][i][j] = new double[resZ + zOffset3[d]];
				tempVelocities[d][i][j] = new double[resZ + zOffset3[d]];
				for (int k = 0; k<resZ + zOffset3[d]; k++)
				{
					velocities[d][i][j][k] = 0.0;
					tempVelocities[d][i][j][k] = 0.0;
				}
			}
		}
	}
}

void StableFluid3D::clearGrids()
{

	if (!velocities && !tempVelocities)
		return;

	for (int d = 0; d < 3; ++d)
	{
		for (int i = 0; i < resX + xOffset3[d]; ++i)
		{
			for (int j = 0; j < resY + yOffset3[d]; ++j)
			{
				delete[] velocities[d][i][j];
				delete[] tempVelocities[d][i][j];
			}
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

void StableFluid3D::initPoisson()
{
	clearPoisson();
	int n = resX * resY * resZ;

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
	double rdzSqr = (resZ * resZ) / (depth * depth);

	// initialize values for the laplacian matrix
	for (int i = 0; i < resX; ++i)
	{
		for (int j = 0; j < resY; ++j)
		{
			for (int k = 0; k < resZ; ++k)
			{
				int row = makeIndex(i, j, k);

				int x1 = (i == (resX - 1) ? row : row + 1);
				int x2 = (i == 0 ? row : row - 1);

				int y1 = (j == (resY - 1) ? row : row + resX);
				int y2 = (j == 0 ? row : row - resX);

				int z1 = (k == (resZ - 1) ? row : row + resY * resX);
				int z2 = (k == 0 ? row : row - resY * resX);

				laplacian[row][row] += 6 * (rdxSqr + rdySqr + rdzSqr);
				laplacian[row][x1] -= rdxSqr;
				laplacian[row][x2] -= rdxSqr;
				laplacian[row][y1] -= rdySqr;
				laplacian[row][y2] -= rdySqr;
				laplacian[row][z1] -= rdzSqr;
				laplacian[row][z2] -= rdzSqr;
			}
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

	int maxIndex = resX * resY * resZ;
	//	for (int i = 0; i < maxIndex; ++i)
	//		preconditioner[i][i] = 1;

	//	cerr<<"Sparsity of preconditioner is " << preconditioner.Sparsity() << endl;
}

void StableFluid3D::clearPoisson()
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

int StableFluid3D::makeIndex(int x, int y, int z) const
{
	return (x + resX * (y + resY * z));
}

void StableFluid3D::splitIndex(int index, int& x, int& y, int& z)
{
	int t = resX * resY;
	x = (index % t) % resY;
	y = (index % t) / resY;
	z = index / t;
}

void StableFluid3D::setBoundary()
{
	// top and bottom boundaries
	for (int i = 0; i < resZ; ++i)
	{
		for (int j = 0; j < resX; ++j)
		{
			tempVelocities[1][j][0][i] = 0;
			tempVelocities[1][j][resY][i] = 0;
		}
	}

	// left and right boundaries
	for (int i = 0; i < resY; ++i)
	{
		for (int j = 0; j < resZ; ++j)
		{
			tempVelocities[0][0][i][j] = 0;
			tempVelocities[0][resX][i][j] = 0;
		}
	}

	// back and front boundaries
	for (int i = 0; i < resY; ++i)
	{
		for (int j = 0; j < resX; ++j)
		{
			tempVelocities[2][j][i][0] = 0;
			tempVelocities[2][j][i][resZ] = 0;
		}
	}

}

void StableFluid3D::addForcesToVelocity(double dt)
{
	for (int i = 0; i < resX; ++i)
	{
		for (int j = 0; j <= resY; ++j)
		{
			for (int k = 0; k < resZ; ++k)
			{
				tempVelocities[1][i][j][k] -= (double)(GRAV_ACCEL * dt);
			}
		}
	}

	if (viscosity != 0)
	{
		double rdx = resX / width;
		double rdy = resY / height;
		double rdz = resZ / depth;

		for (int i = 0; i <= resX; ++i)
		{
			for (int j = 0; j < resY; ++j)
			{
				for (int k = 0; k < resZ; ++k)
				{
					double v1 = (i == 0 ? tempVelocities[0][i][j][k] : tempVelocities[0][i - 1][j][k]);
					double v2 = (i == resX ? tempVelocities[0][i][j][k] : tempVelocities[0][i + 1][j][k]);
					tempVelocities[0][i][j][k] -= viscosity * (v1 - 2 * tempVelocities[0][i][j][k] + v2) * rdx;
				}
			}
		}
		for (int i = 0; i < resX; ++i)
		{
			for (int j = 0; j <= resY; ++j)
			{
				for (int k = 0; k < resZ; ++k)
				{
					double v1 = (j == 0 ? tempVelocities[1][i][j][k] : tempVelocities[1][i][j - 1][k]);
					double v2 = (j == resY ? tempVelocities[1][i][j][k] : tempVelocities[1][i][j + 1][k]);
					tempVelocities[1][i][j][k] -= viscosity * (v1 - 2 * tempVelocities[1][i][j][k] + v2) * rdy;
				}
			}
		}
		for (int i = 0; i < resX; ++i)
		{
			for (int j = 0; j < resY; ++j)
			{
				for (int k = 0; k <= resZ; ++k)
				{
					double v1 = (k == 0 ? tempVelocities[2][i][j][k] : tempVelocities[2][i][j][k - 1]);
					double v2 = (k == resZ ? tempVelocities[2][i][j][k] : tempVelocities[2][i][j][k + 1]);
					tempVelocities[2][i][j][k] -= viscosity * (v1 - 2 * tempVelocities[2][i][j][k] + v2) * rdz;
				}
			}
		}
	}
}

void StableFluid3D::advectVelocity(double dt)
{
	double xStep = width / (double)resX;
	double yStep = height / (double)resY;
	double zStep = depth / (double)resZ;

	double fOffsetsX[3] = { 0, yStep * 0.5, zStep * 0.5 };
	double fOffsetsY[3] = { xStep * 0.5, 0, zStep * 0.5 };
	double fOffsetsZ[3] = { xStep * 0.5, yStep * 0.5, 0 };

	vec3d pos;
	vec3d v;

	for (int d = 0; d < 3; ++d)
	{
		pos[0] = fOffsetsX[d];
		for (int i = 0; i < resX + xOffset3[d]; ++i)
		{
			pos[0] += xStep;
			pos[1] = fOffsetsY[d];
			for (int j = 0; j < resY + yOffset3[d]; ++j)
			{
				pos[1] += yStep;
				pos[2] = fOffsetsZ[d];
				for (int k = 0; k < resZ + zOffset3[d]; ++k)
				{
					pos[2] += zStep;

					v = linearSamp(pos, tempVelocities);
					velocities[d][i][j][k] =
						linearSamp(pos - v * dt, tempVelocities)[d];
				}
			}
		}
	}
}

void StableFluid3D::advectPressure(double dt)
{
	double xStep = width / (double)resX;
	double yStep = height / (double)resY;
	double zStep = depth / (double)resZ;

	double fOffsetsX[3] = { 0, yStep * 0.5, zStep * 0.5 };
	double fOffsetsY[3] = { xStep * 0.5, 0, zStep * 0.5 };
	double fOffsetsZ[3] = { xStep * 0.5, yStep * 0.5, 0 };

	for (int i = 0; i < resX; ++i)
	{
		for (int j = 0; j < resY; ++j)
		{
			for (int k = 0; k < resZ; ++k)
			{
				vec3d pos;
				pos[0] = i*xStep;
				pos[1] = j*yStep;
				pos[2] = k*zStep;
				vec3d v = linearSamp(pos, velocities);
				pressure[makeIndex(i, j, k)] = linearSamp(pos - v * dt, tempPressure);
			}
		}
	}
}

void StableFluid3D::projectVelocity(double dt)
{
	setBoundary();

	Vector divergence(vel.Length());
	// fill in velocities
	double dx = width / (double)resX;
	double dy = height / (double)resY;
	double dz = depth / (double)resZ;
	for (int i = 0; i < resX; ++i)
	{
		for (int j = 0; j < resY; ++j)
		{
			for (int k = 0; k < resZ; ++k)
			{
				int index = makeIndex(i, j, k);

				double div = (tempVelocities[0][i + 1][j][k] - tempVelocities[0][i][j][k]) / dx + (tempVelocities[1][i][j + 1][k] - tempVelocities[1][i][j][k]) / dy + (tempVelocities[2][i][j][k + 1] - tempVelocities[2][i][j][k]) / dz;

				vel[index] = div;
				divergence[index] = div;
			}
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
			for (int k = 0; k < resZ; ++k)
			{
				int i1 = makeIndex(i, j, k);
				int i2;

				double deltaPressure;
				i2 = (i == 0 ? i1 : i1 - 1);
				velocities[0][i][j][k] = tempVelocities[0][i][j][k] - dt* 0.5 * (pressure[i1] - pressure[i2]) / dx;
				deltaPressure = pressure[i1] - pressure[i2];

				i2 = (j == 0 ? i1 : i1 - resX);
				velocities[1][i][j][k] = tempVelocities[1][i][j][k] - dt *  0.5 * (pressure[i1] - pressure[i2]) / dy;
				deltaPressure = pressure[i1] - pressure[i2];

				i2 = (k == 0 ? i1 : i1 - resY * resX);
				velocities[2][i][j][k] = tempVelocities[2][i][j][k] - dt * 0.5 * (pressure[i1] - pressure[i2]) / dz;
				deltaPressure = pressure[i1] - pressure[i2];
			}
		}
	}
}