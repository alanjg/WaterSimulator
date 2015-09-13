#include "stdafx.h"
#include "Vector.h"

Vector operator+(const Vector& lhs, const Vector& rhs)
{
	return Vector(lhs) += rhs;
}

Vector operator-(const Vector& lhs, const Vector& rhs)
{
	return Vector(lhs) -= rhs;
}

Vector operator*(const Vector& lhs, const double& rhs)
{
	return Vector(lhs) *= rhs;
}

Vector operator*(const double& lhs, const Vector& rhs)
{
	return Vector(rhs) *= lhs;
}

Vector operator/(const Vector& lhs, const double& rhs)
{
	return Vector(lhs) /= rhs;
}

Vector operator-(const Vector& lhs)
{
	return Vector(lhs) *= -1;
}

double InnerProduct(const Vector& lhs, const Vector& rhs)
{
	if (lhs.Length() != rhs.Length()) throw WrongDimensionException();

	double accum(0);
	int len = lhs.Length();
	for (int i = 0; i<len; i++)
	{
		accum += lhs[i] * rhs[i];
	}
	return accum;
}

Vector& Vector::operator+=(const Vector& rhs)
{
	if (Cols != rhs.Cols) throw WrongDimensionException();

	for (int i = 0; i<Cols; i++)
	{
		rep[i] += rhs[i];
	}
	return *this;
}

double Vector::Normalize()
{
	double mag = Magnitude();
	if (mag <= 0) return 0;
	//			assert(mag > 0);
	*this /= mag;
	return mag;
}

double Vector::Magnitude() const
{
	double mag = sqrt(Magnitude2());
	return mag;
}

double Vector::Magnitude2() const
{
	double mag2 = 0;
	for (int i = 0; i<Cols; i++)
		mag2 += rep[i] * rep[i];
	return mag2;
}

Vector& Vector::operator-=(const Vector& rhs)
{
	if (Cols != rhs.Cols) throw WrongDimensionException();

	for (int i = 0; i<Cols; i++)
	{
		rep[i] -= rhs[i];
	}
	return *this;
}

Vector& Vector::operator*=(const double& rhs)
{
	for (int i = 0; i<Cols; i++)
	{
		rep[i] *= rhs;
	}
	return *this;
}

Vector& Vector::operator/=(const double& rhs)
{
	for (int i = 0; i<Cols; i++)
	{
		rep[i] /= rhs;
	}
	return *this;
}

int Vector::Length() const
{
	return Cols;
}

void Vector::Resize(int size)
{
	rep.resize(size, 0);
	fill(rep.begin(), rep.end(), 0);
	Cols = size;
}

void Vector::clear()
{
	std::fill(rep.begin(), rep.end(), 0);
}

