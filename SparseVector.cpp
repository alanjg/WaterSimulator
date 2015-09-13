#include "stdafx.h"
#include "SparseVector.h"
#include "Vector.h"

int SparseVector::Length() const { return length_; }
double& SparseVector::operator[](int index)
{
	std::list<SparseNode>::iterator it = nodes_.begin();
	while (it != nodes_.end() && it->index < index)
	{
		it++;
	}
	if (it == nodes_.end() || it->index != index)
	{
		it = nodes_.insert(it, SparseNode(index, 0));
		return it->value;
	}

	return it->value;
}
double SparseVector::operator[](int index) const
{
	std::list<SparseNode>::const_iterator it = nodes_.begin();
	while (it != nodes_.end() && it->index < index)
	{
		it++;
	}
	if (it == nodes_.end() || it->index != index)
		return 0;

	return it->value;

}

//todo - make this faster
bool SparseVector::operator==(const SparseVector& rhs) const
{
	for (int i = 0; i<rhs.length_; i++)
	{
		if (rhs[i] != (*this)[i])
			return false;
	}
	return true;
}

bool SparseVector::operator != (const SparseVector& rhs) const
{
	return !((*this) == rhs);
}


double SparseVector::at(int index) const
{
	std::list<SparseNode>::const_iterator it = nodes_.begin();
	while (it != nodes_.end() && it->index < index)
	{
		it++;
	}
	if (it == nodes_.end() || it->index != index)
		return 0;

	return it->value;
}

SparseVector SparseVector::operator *(const SparseVector& rhs) const
{
	// note: assuming two vectors have same length
	SparseVector product(length_);

	std::list<SparseVector::SparseNode>::const_iterator ita = nodes_.begin();
	std::list<SparseVector::SparseNode>::const_iterator itb = rhs.nodes_.begin();

	while (ita != nodes_.end() && itb != rhs.nodes_.end())
	{
		if (itb->index == ita->index)
		{
			product[itb->index] = ita->value * itb->value;
			++ita;
			++itb;
		}
		else if (itb->index < ita->index)
		{
			++itb;
		}
		else
		{
			++ita;
		}
	}
	return product;
}

SparseVector SparseVector::operator /(const SparseVector& rhs) const
{
	// note: assuming two vectors have same length
	SparseVector quotient(length_);

	std::list<SparseVector::SparseNode>::const_iterator ita = nodes_.begin();
	std::list<SparseVector::SparseNode>::const_iterator itb = rhs.nodes_.begin();

	while (ita != nodes_.end() && itb != rhs.nodes_.end())
	{
		if (itb->index == ita->index)
		{
			quotient[itb->index] = ita->value / itb->value;
			++ita;
			++itb;
		}
		else if (itb->index < ita->index)
		{
			++itb;
		}
		else
		{
			++ita;
		}
	}
	return quotient;
}

SparseVector& SparseVector::operator+=(const SparseVector& rhs)
{
	std::list<SparseVector::SparseNode>::iterator ita = nodes_.begin();
	std::list<SparseVector::SparseNode>::const_iterator itb = rhs.nodes_.begin();

	while (ita != nodes_.end() && itb != rhs.nodes_.end())
	{
		if (itb->index == ita->index)
		{
			ita->value += itb->value;
			++ita;
			++itb;
		}
		else if (itb->index < ita->index)
		{
			nodes_.insert(ita, SparseVector::SparseNode(itb->index, itb->value));
			++itb;
		}
		else
		{
			++ita;
		}
	}
	std::copy(itb, rhs.nodes_.end(), std::inserter(nodes_, ita));
	return *this;
}
SparseVector& SparseVector::operator-=(const SparseVector& rhs)
{
	std::list<SparseVector::SparseNode>::iterator ita = nodes_.begin();
	std::list<SparseVector::SparseNode>::const_iterator itb = rhs.nodes_.begin();

	while (ita != nodes_.end() && itb != rhs.nodes_.end())
	{
		if (itb->index == ita->index)
		{
			ita->value -= itb->value;
			++ita;
			++itb;
		}
		else if (itb->index < ita->index)
		{
			nodes_.insert(ita, SparseVector::SparseNode(itb->index, -itb->value));
			++itb;
		}
		else
		{
			++ita;
		}
	}
	std::transform(itb, rhs.nodes_.end(), inserter(nodes_, ita),
		SparseNode::negate());
	return *this;
}
SparseVector& SparseVector::operator*=(double rhs)
{
	std::transform(nodes_.begin(), nodes_.end(), nodes_.begin(),
		SparseNode::multiply(rhs));
	return *this;
}
SparseVector& SparseVector::operator/=(double rhs)
{
	transform(nodes_.begin(), nodes_.end(), nodes_.begin(),
		SparseNode::divide(rhs));
	return *this;
}

//this zeroes out this vector
void SparseVector::clear()
{
	nodes_.clear();
}

void SparseVector::MakeSparse()
{
	nodes_.erase(std::remove_if(nodes_.begin(), nodes_.end(), NodeIsZero()), nodes_.end());
}

double SparseVector::Sparsity() const
{
	return double(nodes_.size()) / double(length_);
}

void SparseVector::Resize(int size)
{
	if (size < length_)
	{
		nodes_.erase(std::remove_if(nodes_.begin(), nodes_.end(), IndexGreaterThan(size - 1)), nodes_.end());
	}
	length_ = size;
}

double InnerProduct(const SparseVector& lhs, const SparseVector& rhs)
{
	double product = 0;
	std::list<SparseVector::SparseNode>::const_iterator ita = lhs.nodes_.begin();
	std::list<SparseVector::SparseNode>::const_iterator itb = rhs.nodes_.begin();

	while (ita != lhs.nodes_.end() && itb != rhs.nodes_.end())
	{
		if (itb->index == ita->index)
		{
			product += ita->value * itb->value;
			++ita;
			++itb;
		}
		else if (itb->index < ita->index)
		{
			++itb;
		}
		else
		{
			++ita;
		}
	}
	return product;
}


double InnerProduct(const SparseVector& lhs, const Vector& rhs)
{
	if (lhs.Length() != rhs.Length()) throw WrongDimensionException();

	double product = 0;
	SparseVector::const_iterator it(lhs);

	while (it)
	{
		product += it.getVal() * rhs[it.getIndex()];

		++it;
	}

	return product;
}


double InnerProduct(const Vector& lhs, const SparseVector& rhs)
{
	if (lhs.Length() != rhs.Length()) throw WrongDimensionException();

	double product = 0;
	SparseVector::const_iterator it(rhs);

	while (it)
	{
		product += it.getVal() * lhs[it.getIndex()];

		++it;
	}

	return product;
}

SparseVector operator+(const SparseVector& a, const SparseVector& b)
{
	return SparseVector(a) += b;
}

SparseVector operator-(const SparseVector& a, const SparseVector& b)
{
	return SparseVector(a) -= b;
}

SparseVector operator*(const SparseVector& a, double b)
{
	return SparseVector(a) *= b;
}

SparseVector operator*(double a, const SparseVector& b)
{
	return SparseVector(b) *= a;
}

SparseVector operator/(const SparseVector& a, double b)
{
	return SparseVector(a) /= b;
}
