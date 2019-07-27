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
		Eigen::Matrix<real, size, size> A;
		Eigen::Matrix<real, size, 1> B;
		unsigned int nequations;
	public:
		LeastSquaresLinearSystem() : nequations(0) { A.fill(real(0)); B.fill(real(0)); }
		LeastSquaresLinearSystem(unsigned int s) : A(s, s), B(s), nequations(0)
		{
			A.fill(real(0)); B.fill(real(0));
		}

		template<typename r2>
		void add_equation(const SparseEquation<r2>& equation)
		{
			++nequations;
            for (auto i : equation.a)
            {
				B(i.first) += i.second*equation.b;
				for (auto j : equation.a)
					A(i.first, j.first) += i.second*j.second;
			}
		}

		void clear() { A.fill(real(0)); B.fill(real(0)); }

		template<typename r2>
		void add_equation(const DenseEquation<r2>& equation)
		{
			++nequations;
			unsigned int ii = 0;
			for (auto i : equation.a) if (i != real(0)) {
				B(ii) += i*equation.b;
				unsigned int ij = 0;
				for (auto j : equation.a) {
					A(ii, ij) += i*j;
					++ij;
				}
				++ii;
			}
		}

		void normalize()
		{
			if (nequations > 0) {
				A *= (real(1) / real(nequations));
				B *= (real(1) / real(nequations));
				nequations = 1;
			}
		}

		void operator*=(const real& v)
		{
			A *= v; B *= v;
		}

		LeastSquaresLinearSystem<real, size> operator*(const real& v) const
		{
			LeastSquaresLinearSystem s = (*this); s *= v; return s;
		}

		void operator+=(const LeastSquaresLinearSystem<real, size>& that)
		{
			normalize();
			LeastSquaresLinearSystem<real, size> t2 = that;
			t2.normalize();
			A += t2.A; B += t2.B;
		}

		LeastSquaresLinearSystem<real, size> operator+(const LeastSquaresLinearSystem<real, size>& that)
		{
			LeastSquaresLinearSystem s = (*this); s += that; return s;
		}

		Eigen::Matrix<real, size, 1> solve() const
		{
			/**		Eigen::JacobiSVD<Eigen::Matrix<real, size, size>> svd(A);
			double cond = svd.singularValues()(0)
			/ svd.singularValues()(svd.singularValues().size()-1);
			std::cerr<<"Condition number = "<<cond<<std::endl; **/
			return A.ldlt().solve(B);
		}

		void show() const {
			std::cout << A << std::endl;
		}

		real conditionNumber() const {
			Eigen::JacobiSVD<Eigen::Matrix<real, size, size>> svd(A);
			real cond = svd.singularValues()(0)
				/ svd.singularValues()(svd.singularValues().size() - 1);

		}

		/*void getEigenValuesA()
		{
		cout << "\EIGEN VALUES Matrix A\n";
		cout << A.eigenvalues();

		}*/
        void toFile(unsigned int s) const
          {
              std::ofstream os;
              os.open("A.txt");
              for(int i=0;i<s;i++)
              {
                  for(int j=0;j<s;j++)
                      os<<A(i,j)<<" ";
                  os<<std::endl;
              }

              os.close();
              os.open("B.txt");
              os<<B;
              os.close();
          }
          void toFileB(std::string fileName) const
          {
              std::ofstream os;
              os.close();
              os.open(fileName);
              os<<B;
              os.close();
          }
	};



};

#endif
