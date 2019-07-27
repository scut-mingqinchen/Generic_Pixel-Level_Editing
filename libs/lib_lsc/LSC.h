#pragma once
// LSC.h: interface for the LSC class.
//This program demonstrates the LSC superpixel segmentation method proposed in the following paper :
//
//Jiansheng Chen, Zhengqin Li, Bo Huang, Linear Spectral Clustering Superpixel, IEEE Transactions on Image Processing, Vol. 26, No. 7, pp. 3317 - 3330, 2017.
//
//The program is free for non - commercial academic use.Any commercial use is strictly prohibited without the authors' consent. 
//
//run LSC_demo.m under Matlab to try!
//
//You may also compile the mex file yourself using "mex LSC_mex.cpp" under matlab.This demo is also available at http ://jschenthu.weebly.com/projects.html.
//
//Any problem regarding this demo, please feel free to contact us :
//
//Jiansheng Chen, jschenthu@tsinghua.edu.cn
//Zhengqin Li, zhl378@eng.ucsd.edu
//
//Enjoy!
//
//June 20, 2017
//===========================================================================
//////////////////////////////////////////////////////////////////////
#if !defined(_LSC_H_INCLUDED_)
#define _LSC_H_INCLUDED_
#include <vector>
#include<queue>
#include <string>
#include <algorithm>
#include<cmath>
#include<float.h>
using namespace std;
class lscpoint
{
public:
	int x, y;
	lscpoint(int X = 0, int Y = 0) :x(X), y(Y) {};
};
class LSCSuperpixel
{
public:
	int Label;
	int Size;
	vector<int> Neighbor;
	LSCSuperpixel(int L = 0, int S = 0) :Label(L), Size(S) {}
	vector<int> xLoc;
	vector<int> yLoc;
	friend bool operator==(LSCSuperpixel& S, int L) {
		return S.Label == L;
	}
	friend bool operator==(int L, LSCSuperpixel& S) {
		return S.Label == L;
	}
};

class LSC
{
public:
	LSC();
	virtual ~LSC();
	void RGB2XYZ(unsigned char sR, unsigned char sG, unsigned char sB, double&	X, double& Y, double& Z);
	//
	void RGB2LAB(const unsigned char& sR, const unsigned char& sG, const unsigned char& sB, unsigned char& lval, unsigned char& aval, unsigned char& bval);

	void myrgb2lab(unsigned char* r, unsigned char* g, unsigned char* b,unsigned char* L, unsigned char* A, unsigned char* B,int nRows, int nCols);
	int Seeds(int nRows, int nCols, int Row_num, int Col_num, int Row_step, int Col_step, int seed_num, lscpoint* point_array);
	//Initialize the seeds

	//Enforce Connectivity by merging very small superpixels with their neighbors
	void preEnforceConnectivity(unsigned short int* label, int nRows, int nCols);
	void EnforceConnectivity(
		float** L1,
		float** L2,
		float** a1,
		float** a2,
		float** b1,
		float** b2,
		float** x1,
		float** x2,
		float** y1,
		float** y2,
		double** W,
		unsigned short int* label,
		int threshold,
		int nRows,
		int nCols
	);
	//map pixels into ten dimensional feature space
	void Initialize(
		unsigned char* L,
		unsigned char* a,
		unsigned char* b,
		float** L1,
		float** L2,
		float** a1,
		float** a2,
		float** b1,
		float** b2,
		float** x1,
		float** x2,
		float** y1,
		float** y2,
		double** W,
		int nRows,
		int nCols,
		int StepX,
		int StepY,
		float Color,
		float Distance
	);
	//Perform weighted kmeans iteratively in the ten dimensional feature space.

	void DoSuperpixel(
		float** L1,
		float** L2,
		float** a1,
		float** a2,
		float** b1,
		float** b2,
		float** x1,
		float** x2,
		float** y1,
		float** y2,
		double** W,
		unsigned short int* label,
		lscpoint* seedArray,
		int seedNum,
		int nRows,
		int nCols,
		int StepX,
		int StepY,
		int iterationNum,
		int thresholdCoef
	);
	//count 
	int countSuperpixel(unsigned short int* label, int nRows, int nCols);

	//LSC superpixel segmentation algorithm
	void doLSC(unsigned char* R, unsigned char* G, unsigned char* B, int nRows, int nCols, int superpixelnum, double ratio, unsigned short* label);





};


#endif 