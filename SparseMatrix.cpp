#include "stdafx.h"
#include "SparseMatrix.h"
#include "Vector.h"
#include "SparseVector.h"

SparseMatrix::SparseMatrix()
{
	numRows_ = 0;
	numCols_ = 0;
}

SparseMatrix::SparseMatrix(int rows, int cols)
	:rows_(rows, SparseVector(cols)), numRows_(rows), numCols_(cols)
{
}

double SparseMatrix::det(int i, std::vector<bool>& use) const
{
	if (i == numCols_)
		return 1;
	double sum = 0;
	int neg = 1;
	for (int k = 0; k<numCols_; k++)
	{
		if (!use[k])
			continue;
		use[k] = false;
		sum += neg * rows_[i][k] * det(i + 1, use);
		use[k] = true;
		neg *= -1;
	}
	return sum;
}

//this is really, really slow
double SparseMatrix::Determinant() const
{
	std::vector<bool> use(numCols_, true);
	return det(0, use);
}

SparseMatrix SparseMatrix::Transpose() const
{
	SparseVector::const_iterator iter;
	SparseMatrix res(numCols_, numRows_);
	for (int i = 0; i<numRows_; i++)
	{
		iter.reset((*this)[i]);
		while (iter)
		{
			double value = *iter;
			if (value != 0.0)
			{
				int j = iter.getIndex();

				//all insertions occur at the back
				res[j].nodes_.push_back(SparseVector::SparseNode(i, value));
			}
			iter++;
		}
	}
	return res;
}

//use the fact that A * Atrans is symmetric, and only involves row-row dot products
SparseMatrix SparseMatrix::MultiplyTranspose() const
{
	SparseMatrix ret(numRows_, numRows_);
	//first phase computes A[i][j] for j >= i
	//uses n * n/2 * m operations
	for (int i = 0; i < numRows_; i++)
	{
		for (int j = i; j < numRows_; j++)
		{
			double val = InnerProduct(rows_[i], rows_[j]);
			if (val != 0)
			{
				//we're always going to be inserting at the back
				ret.rows_[i].nodes_.push_back(SparseVector::SparseNode(j, val));
			}
		}
	}
	//second phase sets A[j][i] = A[i][j], j >i
	//uses n * m operations
	for (int i = numRows_ - 1; i >= 0; i--)
	{
		SparseVector::const_iterator iter(ret[i]);
		iter.seek(i + 1);
		while (iter)
		{
			double val = *iter;
			if (val != 0)
			{
				int j = iter.getIndex();
				//this works because we build the lower triangular part backwards
				//thus, ret[j][i] will always be the first value in ret[j]
				ret[j].nodes_.push_front(SparseVector::SparseNode(i, val));
			}
			iter++;
		}
	}
	return ret;
}

//use the fact that A * Atrans is symmetric, and only involves row-row dot products
SparseMatrix SparseMatrix::FastMultiplyTranspose() const
{
	SparseMatrix ret(numRows_, numRows_);
	std::vector<std::list<int> > coverage(numCols_);
	for (int i = 0; i < numRows_; i++)
	{
		SparseVector::const_iterator it(rows_[i]);
		while (it)
		{
			coverage[it.getIndex()].push_back(i);
			it++;
		}
	}
	std::set<std::pair<int, int> > computed;
	//do pairwise collisions

	//first phase computes A[i][j] for j >= i
	//uses n * m * m operations(assuming n * m non-sparse elements of ret)
	for (unsigned int i = 0; i < coverage.size(); i++)
	{
		for (std::list<int>::iterator it = coverage[i].begin(); it != coverage[i].end(); ++it)
		{
			for (std::list<int>::iterator it2 = it; it2 != coverage[i].end(); ++it2)
			{
				if (computed.find(std::make_pair(*it, *it2)) == computed.end())
				{
					double product = InnerProduct(rows_[*it], rows_[*it2]);
					if (product != 0)
					{
						ret[*it][*it2] = product;
					}

					computed.insert(std::make_pair(*it, *it2));
				}
			}
		}
	}

	//second phase sets A[j][i] = A[i][j], j >i
	//uses n * m operations
	for (int i = numRows_ - 1; i >= 0; i--)
	{
		SparseVector::const_iterator iter(ret[i]);
		iter.seek(i + 1);
		while (iter)
		{
			double val = *iter;
			if (val != 0)
			{
				int j = iter.getIndex();
				//this works because we build the lower triangular part backwards
				//thus, ret[j][i] will always be the first value in ret[j]
				ret[j].nodes_.push_front(SparseVector::SparseNode(i, val));
			}
			iter++;
		}
	}
	return ret;
}

SparseMatrix& SparseMatrix::operator += (const SparseMatrix& rhs)
{
	for (int i = 0; i<numRows_; i++)
	{
		rows_[i] += rhs.rows_[i];
	}
	return *this;
}

SparseMatrix& SparseMatrix::operator -= (const SparseMatrix& rhs)
{
	for (int i = 0; i<numRows_; i++)
	{
		rows_[i] -= rhs.rows_[i];
	}
	return *this;
}

SparseMatrix& SparseMatrix::operator *= (const SparseMatrix& rhs)
{
	SparseMatrix temp(*this);
	SparseMatrix rhsTrans = rhs.Transpose();
	rows_ = std::vector<SparseVector >(numRows_, SparseVector(rhs.numCols_));
	numCols_ = rhs.numCols_;
	for (int i = 0; i<numRows_; i++)
	{
		for (int j = 0; j<rhs.numCols_; j++)
		{
			double sum = InnerProduct(temp[i], rhsTrans[j]);
			if (sum != 0)
				rows_[i][j] = sum;
		}
	}
	return *this;
}

SparseMatrix& SparseMatrix::operator*=(double rhs)
{
	for (int i = 0; i<numRows_; i++)
	{
		rows_[i] *= rhs;
	}
	return *this;
}

SparseMatrix& SparseMatrix::operator/=(double rhs)
{
	for (int i = 0; i<numRows_; i++)
	{
		rows_[i] /= rhs;
	}
	return *this;
}

bool SparseMatrix::operator==(const SparseMatrix& rhs) const
{
	for (int i = 0; i<numRows_; i++)
	{
		if (rows_[i] != rhs.rows_[i])
			return false;
	}
	return true;
}
bool SparseMatrix::operator!=(const SparseMatrix& rhs) const
{
	return !(operator==(rhs));
}

void SparseMatrix::MakeSparse()
{
	for (unsigned int i = 0; i < rows_.size(); i++)
	{
		rows_[i].MakeSparse();
	}
}

double SparseMatrix::Sparsity() const
{
	double total = 0.0;
	for (unsigned int i = 0; i < rows_.size(); i++)
	{
		total += rows_[i].Sparsity();
	}
	return total / numRows_;
}

void SparseMatrix::clear()
{
	for (std::vector<SparseVector >::iterator i = rows_.begin(); i != rows_.end(); ++i)
		i->clear();
}

SparseMatrix operator+(const SparseMatrix& a, const SparseMatrix& b)
{
	return SparseMatrix(a) += b;
}


SparseMatrix operator-(const SparseMatrix& a, const SparseMatrix& b)
{
	return SparseMatrix(a) -= b;
}


SparseMatrix operator*(const SparseMatrix& a, const SparseMatrix& b)
{
	return SparseMatrix(a) *= b;
}


SparseMatrix operator*(const SparseMatrix& a, double b)
{
	return SparseMatrix(a) *= b;
}


SparseMatrix operator*(double a, const SparseMatrix& b)
{
	return SparseMatrix(b) *= a;
}


SparseVector operator*(const SparseMatrix& a, const SparseVector& b)
{
	SparseVector result(a.Rows());
	for (int i = 0; i<a.numRows_; i++)
	{
		double product = InnerProduct(a[i], b);
		if (product != 0)
		{
			//always add to back
			result.nodes_.push_back(SparseVector::SparseNode(i, product));
		}
	}
	return result;
}


Vector operator*(const SparseMatrix& a, const Vector& b)
{
	Vector result(a.Rows());
	a.Multiply(b, result);
	return result;
}

void SparseMatrix::Multiply(const Vector& rhs, Vector& result) const
{
	if (rhs.Length() != result.Length()) throw new WrongDimensionException();
	for (int i = 0; i<numRows_; i++)
	{
		result[i] = InnerProduct(rows_[i], rhs);
	}
}


SparseMatrix operator/(const SparseMatrix& a, double b)
{
	return SparseMatrix(a) /= b;
}


void ConjugateGradient(const SparseMatrix& A, SparseVector& x, const SparseVector& b, double tolerance, double maxIter)
{

	int k = 0;
	SparseVector r(b - A * x);
	SparseVector p(r.Length()), w(r.Length());

	double rho1 = InnerProduct(r, r);
	double rho2 = rho1;
	double bScale = std::sqrt(InnerProduct(b, b));
	while ((std::sqrt(rho2) > (bScale * tolerance)) && (k < maxIter))
	{
		++k;
		if (k == 1)
			p = r;
		else
		{
			double beta = rho2 / rho1;
			p *= beta;
			p += r;
		}
		w = A * p;
		double alpha = rho2 / InnerProduct(p, w);
		x += alpha * p;
		r -= alpha * w;
		rho1 = rho2;
		rho2 = InnerProduct(r, r);
	}
}


void ConjugateGradient(const SparseMatrix& A, Vector& x, const Vector& b, double tolerance, double maxIter)
{
	int k = 0;
	Vector r(b - A * x);
	Vector p(r.Length()), w(r.Length());

	double rho1 = r.Magnitude2();
	double rho2 = rho1;
	double bScale = b.Magnitude();
	double errBound = bScale * tolerance;
	while ((std::sqrt(rho2) > errBound) && (k < maxIter))
	{
		++k;
		if (k == 1)
		{
			p = r;
		}
		else
		{
			double beta = rho2 / rho1;
			p *= beta;
			p += r;
		}
		w = A * p;

		double alpha = rho2 / InnerProduct(p, w);
		x += alpha * p;
		r -= alpha * w;
		rho1 = rho2;
		rho2 = r.Magnitude2();
	}
}


void JacobiIteration(const SparseMatrix& A, Vector& x, const Vector& b, int maxIter)
{
	for (int k = 0; k<maxIter; k++)
	{
		Vector xNext(x);
		for (int i = 0; i<A.Rows(); i++)
		{
			double divisor = A[i][i];
			xNext[i] = (b[i] - InnerProduct(A[i], x) + divisor*x[i]) / divisor;
		}
		x = xNext;
		Vector estimate = A * x;
		estimate -= b;
		double error = InnerProduct(estimate, estimate);
	}
}


void GaussSeidelIteration(const SparseMatrix& A, Vector& x, const Vector& b, int maxIter)
{
	for (int k = 0; k<maxIter; k++)
	{
		for (int i = 0; i<A.Rows(); i++)
		{
			double divisor = A[i][i];
			x[i] = (b[i] - InnerProduct(A[i], x) + divisor*x[i]) / divisor;
		}
		Vector estimate = A * x;
		estimate -= b;
		double error = InnerProduct(estimate, estimate);
	}
}