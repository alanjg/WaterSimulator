#pragma once

struct WrongDimensionException : public std::exception { };

class SparseMatrix;

class Vector
{
public:
	Vector() :
		Cols(0)
	{
	}

	explicit Vector(int dim) :
		rep(dim, 0), Cols(dim)
	{
	}

	double& operator[](int index) { return rep[index]; }
	const double& operator[](int index) const { return rep[index]; }
	Vector& operator+=(const Vector& rhs);
	double Normalize();
	double Magnitude() const;
	double Magnitude2() const;
	Vector& operator-=(const Vector& rhs);
	Vector& operator*=(const double& rhs);
	Vector& operator/=(const double& rhs);
	int Length() const;
	void Resize(int size);
	void clear();
private:
	std::vector<double> rep;
	int Cols;
};

Vector operator+(const Vector& lhs, const Vector& rhs);
Vector operator-(const Vector& lhs, const Vector& rhs);
Vector operator*(const Vector& lhs, const double& rhs);
Vector operator*(const double& lhs, const Vector& rhs);
Vector operator/(const Vector& lhs, const double& rhs);
Vector operator-(const Vector& lhs);
double InnerProduct(const Vector& lhs, const Vector& rhs);

