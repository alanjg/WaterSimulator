#pragma once

class Vector;
class SparseVector;
class SparseMatrix
{
	double det(int i, std::vector<bool>& use) const;
public:

	SparseMatrix();
	explicit SparseMatrix(int rows, int cols);
	//use default op=, copy ctr, dtr

	double Determinant() const;
	SparseMatrix Transpose() const;
	SparseMatrix MultiplyTranspose() const;
	SparseMatrix FastMultiplyTranspose() const;

	int Rows() const { return numRows_; }
	int Cols() const { return numCols_; }

	SparseVector& operator[](int index)
	{
		return rows_[index];
	}

	const SparseVector& operator[](int index) const
	{
		return rows_[index];
	}

	friend SparseMatrix operator+(const SparseMatrix& a, const SparseMatrix& b);
	friend SparseMatrix operator-(const SparseMatrix& a, const SparseMatrix& b);
	friend SparseMatrix operator*(const SparseMatrix& a, const SparseMatrix& b);
	friend SparseMatrix operator*(const SparseMatrix& a, double b);
	friend SparseMatrix operator*(double a, const SparseMatrix& b);
	friend SparseVector operator*(const SparseMatrix& a, const SparseVector& b);
	//friend Vector operator*(const SparseMatrix& a,const Vector& b);	
	friend SparseMatrix operator/(const SparseMatrix& a, double b);
	SparseMatrix& operator += (const SparseMatrix& rhs);
	SparseMatrix& operator -= (const SparseMatrix& rhs);
	SparseMatrix& operator *= (const SparseMatrix& rhs);
	SparseMatrix& operator*=(double rhs);
	SparseMatrix& operator/=(double rhs);
	void Multiply(const Vector& rhs, Vector& result) const;
	bool operator==(const SparseMatrix& rhs) const;
	bool operator!=(const SparseMatrix& rhs) const;
	void MakeSparse();
	double Sparsity() const;
	void clear();
private:
	int numRows_;
	int numCols_;
	std::vector<SparseVector > rows_;
};

SparseMatrix operator+(const SparseMatrix& a, const SparseMatrix& b);
SparseMatrix operator-(const SparseMatrix& a, const SparseMatrix& b);
SparseMatrix operator*(const SparseMatrix& a, const SparseMatrix& b);
SparseMatrix operator*(const SparseMatrix& a, double b);
SparseMatrix operator*(double a, const SparseMatrix& b);
SparseVector operator*(const SparseMatrix& a, const SparseVector& b);
//Vector operator*(const SparseMatrix& a,const Vector& b);
SparseMatrix operator/(const SparseMatrix& a, double b);
void ConjugateGradient(const SparseMatrix& A, SparseVector& x, const SparseVector& b, double tolerance = 0.0001, double maxIter = 50);
void ConjugateGradient(const SparseMatrix& A, Vector& x, const Vector& b, double tolerance = 0.0001, double maxIter = 50);
void JacobiIteration(const SparseMatrix& A, Vector& x, const Vector& b, int maxIter = 50);
void GaussSeidelIteration(const SparseMatrix& A, Vector& x, const Vector& b, int maxIter = 50);
