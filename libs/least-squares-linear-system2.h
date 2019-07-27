#ifndef _MATH_OPTIMIZATION_LEAST_SQUARES_LINEAR_SYSTEM_H_
#define _MATH_OPTIMIZATION_LEAST_SQUARES_LINEAR_SYSTEM_H_

#include "Eigen/Dense"
#include <map>
#include <vector>
#include <iostream>

namespace Optimization {

template<typename real>
class SparseEquation
{
public:
	std::map<unsigned int, real> a;
	real b;
	SparseEquation() : b(0) {}
};

template<typename real>
class DenseEquation
{
public:
	std::vector<real> a;
	real b;
	DenseEquation(unsigned int size) : a(size, real(0)), b(0) {}
};

template<typename real, int size = Eigen::Dynamic>
class LeastSquaresLinearSystem {
	Eigen::Matrix<real,size,size> A;
	Eigen::Matrix<real,size,1> B;
public:
	LeastSquaresLinearSystem() { A.fill(real(0)); B.fill(real(0)); }
	LeastSquaresLinearSystem(unsigned int s) : A(s,s), B(s) 
		{ A.fill(real(0)); B.fill(real(0)); }
	
	template<typename r2>
	void add_equation(const SparseEquation<r2>& equation)
	{
		for (auto i : equation.a) {
			B(i.first) += i.second*equation.b;
			for (auto j : equation.a)
				A(i.first, j.first) += i.second*j.second;
		}
	}

	void clear() { A.fill(real(0)); B.fill(real(0)); }

	template<typename r2>
	void add_equation(const DenseEquation<r2>& equation)
	{
		unsigned int ii = 0;
        for (auto i : equation.a) if (i!=real(0))
        {
			B(ii) += i*equation.b;
			unsigned int ij = 0;
            for (auto j : equation.a)
            {
				A(ii, ij) += i*j;
				++ij;
			}
			++ii;
		}
	}

	void normalize()
	{
		real max = std::max({fabs(A.maxCoeff()),fabs(B.maxCoeff()),fabs(A.minCoeff()),fabs(B.minCoeff())});
		A*=(real(1)/max);
		B*=(real(1)/max);
	}

	void operator*=(const real& v)
	{	A*=v; B*=v;       }

	LeastSquaresLinearSystem<real,size> operator*(const real& v) const
	{	LeastSquaresLinearSystem s = (*this); s*=v; return s;     }

	void operator+=(const LeastSquaresLinearSystem<real,size>& that)
	{       A+=that.A; B+=that.B;	}

	LeastSquaresLinearSystem<real,size> operator+(const LeastSquaresLinearSystem<real,size>& that)
	{       LeastSquaresLinearSystem s = (*this); s+=that; return s;   	}

	Eigen::Matrix<real,size,1> solve() const
	{	return A.ldlt().solve(B);       }
    
    /*void getEigenValuesA()
    {
        cout << "\EIGEN VALUES Matrix A\n";
        cout << A.eigenvalues();
        
    }*/
    
};



};

#endif
