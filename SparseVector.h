#pragma once

class Vector;
class SparseMatrix;
class SparseVector
{
	friend class SparseMatrix;
private:
	struct SparseNode
	{
		SparseNode(int i, double v)
		{
			index = i;
			value = v;
		}
		int index;
		double value;
		struct multiply
		{
			double value;
			multiply(double v) { value = v; }
			SparseNode operator()(const SparseNode& node) const
			{
				return SparseNode(node.index, node.value * value);
			}
		};
		struct divide
		{
			double value;
			divide(double v) { value = v; }
			SparseNode operator()(const SparseNode& node) const
			{
				return SparseNode(node.index, node.value / value);
			}
		};
		struct negate
		{
			SparseNode operator()(const SparseNode& node) const
			{
				return SparseNode(node.index, -node.value);
			}
		};
	};
	struct NodeIsZero
	{
		bool operator()(const SparseNode& node)
		{
			return node.value == 0;
		}
	};
	struct IndexGreaterThan
	{
		int n;
		IndexGreaterThan(int i)
		{
			n = i;
		}
		bool operator()(const SparseNode& node)
		{
			return node.index > n;
		}
	};
public:
	SparseVector()
		:length_(0)
	{
	}

	explicit SparseVector(int length)
		:length_(length)
	{
	}
	//use default op=, copy ctr, and dtr

	int Length() const;
	double& operator[](int index);
	double operator[](int index) const;
	bool operator==(const SparseVector& rhs) const;
	bool operator != (const SparseVector& rhs) const;
	double at(int index) const;
	SparseVector operator *(const SparseVector& rhs) const;
	SparseVector operator /(const SparseVector& rhs) const;
	friend SparseVector operator*(const SparseMatrix& a, const SparseVector& b);
	friend double InnerProduct(const SparseVector& lhs, const SparseVector& rhs);
	friend double InnerProduct(const SparseVector& lhs, const Vector& rhs);
	friend double InnerProduct(const Vector& lhs, const SparseVector& rhs);
	SparseVector& operator+=(const SparseVector& rhs);
	SparseVector& operator-=(const SparseVector& rhs);
	SparseVector& operator*=(double rhs);
	SparseVector& operator/=(double rhs);

	//this zeroes out this vector
	void clear();

	void MakeSparse();
	double Sparsity() const;
	void Resize(int size);

	class iterator {
	public:
		iterator() {}
		iterator(SparseVector& toIter) { target = &toIter.nodes_; iter = target->begin(); }
		~iterator() {}

		double&		getVal() { return iter->value; }
		int			getIndex() { return iter->index; }

		operator bool() const { return (iter != target->end()); }
		bool		valid() const { return (iter != target->end()); }

		void		reset() { iter = target->begin(); }
		void		reset(SparseVector& toIter) { target = &toIter.nodes_; iter = target->begin(); }

		bool		seek(int index) { while (valid() && getIndex() < index) { (*this)++; } return (valid() && getIndex() == index); }

		double&		operator *	() { return iter->value; }
		iterator	operator ++	(int) { iterator i = *this; iter++; return i; }
		iterator&	operator ++	() { iter++; return *this; }

	private:
		std::list<SparseNode>::iterator	iter;
		std::list<SparseNode>* target;
	};

	class const_iterator {
	public:
		const_iterator() {}
		const_iterator(const SparseVector& toIter) { target = &toIter.nodes_; iter = target->begin(); }
		~const_iterator() {}

		double		getVal() const { return iter->value; }
		int			getIndex() const { return iter->index; }

		operator bool() const { return (iter != target->end()); }
		bool		valid() const { return (iter != target->end()); }

		void		reset() { iter = target->begin(); }
		void		reset(const SparseVector& toIter) { target = &toIter.nodes_; iter = target->begin(); }

		bool		seek(int index) { while (valid() && getIndex() < index) { (*this)++; } return (valid() && getIndex() == index); }

		double			operator *	() { return iter->value; }
		const_iterator	operator ++	(int) { const_iterator i = *this; iter++; return i; }
		const_iterator&	operator ++	() { iter++; return *this; }

	private:
		std::list<SparseNode>::const_iterator	iter;
		const std::list<SparseNode>*			target;
	};
private:
	int length_;
	std::list<SparseNode> nodes_;
};

double InnerProduct(const SparseVector& lhs, const SparseVector& rhs);
double InnerProduct(const SparseVector& lhs, const Vector& rhs);
double InnerProduct(const Vector& lhs, const SparseVector& rhs);
SparseVector operator+(const SparseVector& a, const SparseVector& b);
SparseVector operator-(const SparseVector& a, const SparseVector& b);
SparseVector operator*(const SparseVector& a, double b);
SparseVector operator*(double a, const SparseVector& b);
SparseVector operator/(const SparseVector& a, double b);