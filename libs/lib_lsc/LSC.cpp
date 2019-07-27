#include "LSC.h"
const double PI = 3.1415926;

LSC::LSC()
{
}
LSC::~LSC()
{

}
void LSC::RGB2XYZ(unsigned char sR, unsigned char sG, unsigned char sB, double&	X, double& Y, double& Z)
{
	double R = sR / 255.0;
	double G = sG / 255.0;
	double B = sB / 255.0;

	double r, g, b;

	if (R <= 0.04045)	r = R / 12.92;
	else				r = pow((R + 0.055) / 1.055, 2.4);
	if (G <= 0.04045)	g = G / 12.92;
	else				g = pow((G + 0.055) / 1.055, 2.4);
	if (B <= 0.04045)	b = B / 12.92;
	else				b = pow((B + 0.055) / 1.055, 2.4);

	X = r*0.412453 + g*0.357580 + b*0.180423;
	Y = r*0.212671 + g*0.715160 + b*0.072169;
	Z = r*0.019334 + g*0.119193 + b*0.950227;
}
void LSC::RGB2LAB(const unsigned char& sR, const unsigned char& sG, const unsigned char& sB, unsigned char& lval, unsigned char& aval, unsigned char& bval)
{
	double X, Y, Z;
	RGB2XYZ(sR, sG, sB, X, Y, Z);

	double epsilon = 0.008856;	//actual CIE standard
	double kappa = 903.3;		//actual CIE standard

	double Xr = 0.950456;	//reference white
	double Yr = 1.0;		//reference white
	double Zr = 1.088754;	//reference white

	double xr = X / Xr;
	double yr = Y / Yr;
	double zr = Z / Zr;

	double fx, fy, fz;
	if (xr > epsilon)	fx = pow(xr, 1.0 / 3.0);
	else				fx = (kappa*xr + 16.0) / 116.0;
	if (yr > epsilon)	fy = pow(yr, 1.0 / 3.0);
	else				fy = (kappa*yr + 16.0) / 116.0;
	if (zr > epsilon)	fz = pow(zr, 1.0 / 3.0);
	else				fz = (kappa*zr + 16.0) / 116.0;

	lval = (unsigned char)((116.0*fy - 16.0) / 100 * 255 + 0.5);
	aval = (unsigned char)(500.0*(fx - fy) + 128 + 0.5);
	bval = (unsigned char)(200.0*(fy - fz) + 128 + 0.5);
}
void LSC::myrgb2lab(unsigned char* r, unsigned char* g, unsigned char* b,
	unsigned char* L, unsigned char* A, unsigned char* B,
	int nRows, int nCols){
	for (int i = 0; i<nRows; i++)
		for (int j = 0; j<nCols; j++)
			RGB2LAB(r[i*nCols + j], g[i*nCols + j], b[i*nCols + j], L[i*nCols + j], A[i*nCols + j], B[i*nCols + j]);
}
int LSC::Seeds(int nRows, int nCols, int Row_num, int Col_num, int Row_step, int Col_step, int seed_num, lscpoint* point_array)
{
	int Row_remain = nRows - Row_step*Row_num;
	int Col_remain = nCols - Col_step*Col_num;
	int t1 = 1, t2 = 1;
	int count = 0;
	int centerx, centery;
    FILE *fpseed=fopen("seed.txt","w");
	for (int i = 0; i<Row_num; i++)
	{
		t2 = 1;
		for (int j = 0; j<Col_num; j++)
		{
			centerx = i*Row_step + 0.5*Row_step + t1;
			centery = j*Col_step + 0.5*Col_step + t2;
			centerx = (centerx >= nRows - 1) ? nRows - 1 : centerx;
			centery = (centery >= nCols - 1) ? nCols - 1 : centery;
			if (t2<Col_remain)
				t2++;
			point_array[count] = lscpoint(centerx, centery);
            fprintf(fpseed,"%d,%d\n",centerx,centery);
			count++;
		}
		if (t1<Row_remain)
			t1++;
	}
    fclose(fpseed);

	return count;
}
void LSC::preEnforceConnectivity(unsigned short int* label, int nRows, int nCols)
{
	const int dx8[8] = { -1, -1,  0,  1, 1, 1, 0, -1 };
	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1 };
	int adj = 0;
	int Bond = 20;
	bool **mask = new bool*[nRows];
	for (int i = 0; i<nRows; i++)
	{
		mask[i] = new bool[nCols];
		for (int j = 0; j<nCols; j++)
			mask[i][j] = 0;
	}
	vector<unsigned short> xLoc;
	vector<unsigned short> yLoc;
	for (int i = 0; i<nRows; i++)
		for (int j = 0; j<nCols; j++)
		{
			if (mask[i][j] == 0)
			{
				int L = label[i*nCols + j];
				for (int k = 0; k<8; k++)
				{
					int x = i + dx8[k];
					int y = j + dy8[k];
					if (x >= 0 && x <= nRows - 1 && y >= 0 && y <= nCols - 1)
					{
						if (mask[x][y] == 1 && label[x*nCols + y] != L)
							adj = label[x*nCols + y]; break;
					}
				}
				mask[i][j] = 1;
				xLoc.insert(xLoc.end(), i);
				yLoc.insert(yLoc.end(), j);
				int indexMarker = 0;
				while (indexMarker<xLoc.size())
				{
					int x = xLoc[indexMarker]; int y = yLoc[indexMarker];
					indexMarker++;
					int minX = (x - 1 <= 0) ? 0 : x - 1;
					int maxX = (x + 1 >= nRows - 1) ? nRows - 1 : x + 1;
					int minY = (y - 1 <= 0) ? 0 : y - 1;
					int maxY = (y + 1 >= nCols - 1) ? nCols - 1 : y + 1;
					for (int m = minX; m <= maxX; m++)
						for (int n = minY; n <= maxY; n++)
						{
							if (mask[m][n] == 0 && label[m*nCols + n] == L)
							{
								mask[m][n] = 1;
								xLoc.insert(xLoc.end(), m);
								yLoc.insert(yLoc.end(), n);
							}
						}
				}
				if (indexMarker<Bond)
				{
					for (int k = 0; k<xLoc.size(); k++)
					{
						int x = xLoc[k]; int y = yLoc[k];
						label[x*nCols + y] = adj;
					}
				}
				xLoc.clear();
				yLoc.clear();
			}
		}
	for (int i = 0; i<nRows; i++)
		delete[] mask[i];
	delete[] mask;
}

void LSC::EnforceConnectivity(
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
)
{
	unsigned char** mask = new unsigned char*[nRows];
	for (int i = 0; i<nRows; i++)
	{
		mask[i] = new unsigned char[nCols];
		for (int j = 0; j<nCols; j++)
			mask[i][j] = 0;
	}

	vector<unsigned short>strayX;
	vector<unsigned short>strayY;
	vector<unsigned short>Size;
	queue<unsigned short> xLoc;
	queue<unsigned short> yLoc;
	vector<double>centerL1;
	vector<double>centerL2;
	vector<double>centera1;
	vector<double>centera2;
	vector<double>centerb1;
	vector<double>centerb2;
	vector<double>centerx1;
	vector<double>centerx2;
	vector<double>centery1;
	vector<double>centery2;
	vector<double>centerW;
	int sLabel = -1;
	int L;
	for (int i = 0; i<nRows; i++)
		for (int j = 0; j<nCols; j++)
		{
			if (mask[i][j] == 0)
			{
				sLabel++;
				int Count = 1;
				centerL1.insert(centerL1.end(), 0);
				centerL2.insert(centerL2.end(), 0);
				centera1.insert(centera1.end(), 0);
				centera2.insert(centera2.end(), 0);
				centerb1.insert(centerb1.end(), 0);
				centerb2.insert(centerb2.end(), 0);
				centerx1.insert(centerx1.end(), 0);
				centerx2.insert(centerx2.end(), 0);
				centery1.insert(centery1.end(), 0);
				centery2.insert(centery2.end(), 0);
				centerW.insert(centerW.end(), 0);
				strayX.insert(strayX.end(), i);
				strayY.insert(strayY.end(), j);
				double Weight = W[i][j];
				centerL1[sLabel] += L1[i][j] * Weight;
				centerL2[sLabel] += L2[i][j] * Weight;
				centera1[sLabel] += a1[i][j] * Weight;
				centera2[sLabel] += a2[i][j] * Weight;
				centerb1[sLabel] += b1[i][j] * Weight;
				centerb2[sLabel] += b2[i][j] * Weight;
				centerx1[sLabel] += x1[i][j] * Weight;
				centerx2[sLabel] += x2[i][j] * Weight;
				centery1[sLabel] += y1[i][j] * Weight;
				centery2[sLabel] += y2[i][j] * Weight;
				centerW[sLabel] += W[i][j];
				L = label[i*nCols + j];
				label[i*nCols + j] = sLabel;
				mask[i][j] = 1;
				xLoc.push(i); yLoc.push(j);
				while (!xLoc.empty())
				{
					int x = xLoc.front(); xLoc.pop();
					int y = yLoc.front(); yLoc.pop();
					int minX = (x - 1 <= 0) ? 0 : x - 1;
					int maxX = (x + 1 >= nRows - 1) ? nRows - 1 : x + 1;
					int minY = (y - 1 <= 0) ? 0 : y - 1;
					int maxY = (y + 1 >= nCols - 1) ? nCols - 1 : y + 1;
					for (int m = minX; m <= maxX; m++)
						for (int n = minY; n <= maxY; n++)
						{
							if (mask[m][n] == 0 && label[m*nCols + n] == L)
							{
								Count++;
								xLoc.push(m);
								yLoc.push(n);
								mask[m][n] = 1;
								label[m*nCols + n] = sLabel;
								Weight = W[m][n];
								centerL1[sLabel] += L1[m][n] * Weight;
								centerL2[sLabel] += L2[m][n] * Weight;
								centera1[sLabel] += a1[m][n] * Weight;
								centera2[sLabel] += a2[m][n] * Weight;
								centerb1[sLabel] += b1[m][n] * Weight;
								centerb2[sLabel] += b2[m][n] * Weight;
								centerx1[sLabel] += x1[m][n] * Weight;
								centerx2[sLabel] += x2[m][n] * Weight;
								centery1[sLabel] += y1[m][n] * Weight;
								centery2[sLabel] += y2[m][n] * Weight;
								centerW[sLabel] += W[m][n];
							}
						}
				}
				Size.insert(Size.end(), Count);
				centerL1[sLabel] /= centerW[sLabel];
				centerL2[sLabel] /= centerW[sLabel];
				centera1[sLabel] /= centerW[sLabel];
				centera2[sLabel] /= centerW[sLabel];
				centerb1[sLabel] /= centerW[sLabel];
				centerb2[sLabel] /= centerW[sLabel];
				centerx1[sLabel] /= centerW[sLabel];
				centerx2[sLabel] /= centerW[sLabel];
				centery1[sLabel] /= centerW[sLabel];
				centery2[sLabel] /= centerW[sLabel];
			}
		}

	sLabel = sLabel + 1;
	int Count = 0;




	vector<int>::iterator Pointer;
	vector<LSCSuperpixel> Sarray;
	for (int i = 0; i<sLabel; i++)
	{
		if (Size[i]<threshold)
		{
			int x = strayX[i]; int y = strayY[i];
			L = label[x*nCols + y];
			mask[x][y] = 0;
			int indexMark = 0;
			LSCSuperpixel S(L, Size[i]);
			S.xLoc.insert(S.xLoc.end(), x);
			S.yLoc.insert(S.yLoc.end(), y);
			while (indexMark<S.xLoc.size())
			{
				x = S.xLoc[indexMark]; y = S.yLoc[indexMark];
				indexMark++;
				int minX = (x - 1 <= 0) ? 0 : x - 1;
				int maxX = (x + 1 >= nRows - 1) ? nRows - 1 : x + 1;
				int minY = (y - 1 <= 0) ? 0 : y - 1;
				int maxY = (y + 1 >= nCols - 1) ? nCols - 1 : y + 1;
				for (int m = minX; m <= maxX; m++)
					for (int n = minY; n <= maxY; n++)
					{
						if (mask[m][n] == 1 && label[m*nCols + n] == L)
						{
							mask[m][n] = 0;
							S.xLoc.insert(S.xLoc.end(), m);
							S.yLoc.insert(S.yLoc.end(), n);
						}
						else if (label[m*nCols + n] != L)
						{
							int NewLabel = label[m*nCols + n];
							Pointer = find(S.Neighbor.begin(), S.Neighbor.end(), NewLabel);
							if (Pointer == S.Neighbor.end())
							{
								S.Neighbor.insert(S.Neighbor.begin(), NewLabel);
							}
						}
					}

			}
			Sarray.insert(Sarray.end(), S);
		}
	}

	vector<LSCSuperpixel>::iterator S;
	vector<int>::iterator I;
	vector<int>::iterator I2;
	S = Sarray.begin();
	while (S != Sarray.end())
	{
		double MinDist = DBL_MAX;
		int Label1 = (*S).Label;
		int Label2 = -1;
		for (I = (*S).Neighbor.begin(); I != (*S).Neighbor.end(); I++)
		{
			double D = (centerL1[Label1] - centerL1[*I])*(centerL1[Label1] - centerL1[*I]) +
				(centerL2[Label1] - centerL2[*I])*(centerL2[Label1] - centerL2[*I]) +
				(centera1[Label1] - centera1[*I])*(centera1[Label1] - centera1[*I]) +
				(centera2[Label1] - centera2[*I])*(centera2[Label1] - centera2[*I]) +
				(centerb1[Label1] - centerb1[*I])*(centerb1[Label1] - centerb1[*I]) +
				(centerb2[Label1] - centerb2[*I])*(centerb2[Label1] - centerb2[*I]) +
				(centerx1[Label1] - centerx1[*I])*(centerx1[Label1] - centerx1[*I]) +
				(centerx2[Label1] - centerx2[*I])*(centerx2[Label1] - centerx2[*I]) +
				(centery1[Label1] - centery1[*I])*(centery1[Label1] - centery1[*I]) +
				(centery2[Label1] - centery2[*I])*(centery2[Label1] - centery2[*I]);
			if (D<MinDist)
			{
				MinDist = D;
				Label2 = (*I);
			}
		}
		double W1 = centerW[Label1];
		double W2 = centerW[Label2];
		double W = W1 + W2;
		centerL1[Label2] = (W2*centerL1[Label2] + W1*centerL1[Label1]) / W;
		centerL2[Label2] = (W2*centerL2[Label2] + W1*centerL2[Label1]) / W;
		centera1[Label2] = (W2*centera1[Label2] + W1*centera1[Label1]) / W;
		centera2[Label2] = (W2*centera2[Label2] + W1*centera2[Label1]) / W;
		centerb1[Label2] = (W2*centerb1[Label2] + W1*centerb1[Label1]) / W;
		centerb2[Label2] = (W2*centerb2[Label2] + W1*centerb2[Label1]) / W;
		centerx1[Label2] = (W2*centerx1[Label2] + W1*centerx1[Label1]) / W;
		centerx2[Label2] = (W2*centerx2[Label2] + W1*centerx2[Label1]) / W;
		centery1[Label2] = (W2*centery1[Label2] + W1*centery1[Label1]) / W;
		centery2[Label2] = (W2*centery2[Label2] + W1*centery2[Label1]) / W;
		centerW[Label2] = W;
		for (int i = 0; i<(*S).xLoc.size(); i++)
		{
			int x = (*S).xLoc[i]; int y = (*S).yLoc[i];
			label[x*nCols + y] = Label2;
		}
		vector<LSCSuperpixel>::iterator Stmp;
		Stmp = find(Sarray.begin(), Sarray.end(), Label2);
		if (Stmp != Sarray.end())
		{
			Size[Label2] = Size[Label1] + Size[Label2];
			if (Size[Label2] >= threshold)
			{
				Sarray.erase(Stmp);
				Sarray.erase(S);
			}
			else
			{
				(*Stmp).xLoc.insert((*Stmp).xLoc.end(), (*S).xLoc.begin(), (*S).xLoc.end());
				(*Stmp).yLoc.insert((*Stmp).yLoc.end(), (*S).yLoc.begin(), (*S).yLoc.end());
				(*Stmp).Neighbor.insert((*Stmp).Neighbor.end(), (*S).Neighbor.begin(), (*S).Neighbor.end());
				sort((*Stmp).Neighbor.begin(), (*Stmp).Neighbor.end());
				I = unique((*Stmp).Neighbor.begin(), (*Stmp).Neighbor.end());
				(*Stmp).Neighbor.erase(I, (*Stmp).Neighbor.end());
				I = find((*Stmp).Neighbor.begin(), (*Stmp).Neighbor.end(), Label1);
				(*Stmp).Neighbor.erase(I);
				I = find((*Stmp).Neighbor.begin(), (*Stmp).Neighbor.end(), Label2);
				(*Stmp).Neighbor.erase(I);
				Sarray.erase(S);
			}
		}
		else
		{
			Sarray.erase(S);
		}
		for (int i = 0; i<Sarray.size(); i++)
		{
			I = find(Sarray[i].Neighbor.begin(), Sarray[i].Neighbor.end(), Label1);
			I2 = find(Sarray[i].Neighbor.begin(), Sarray[i].Neighbor.end(), Label2);
			if (I != Sarray[i].Neighbor.end() && I2 != Sarray[i].Neighbor.end())
				Sarray[i].Neighbor.erase(I);
			else if (I != Sarray[i].Neighbor.end() && I2 == Sarray[i].Neighbor.end())
				(*I) = Label2;
		}
		S = Sarray.begin();
	}
	for (int i = 0; i<nRows; i++)
		delete[] mask[i];
	delete[] mask;
}
void LSC::DoSuperpixel(
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
)
{
	//Pre-treatment
	double** dist = new double*[nRows];
	for (int i = 0; i<nRows; i++)
		dist[i] = new double[nCols];
	double* centerL1 = new double[seedNum];
	double* centerL2 = new double[seedNum];
	double* centera1 = new double[seedNum];
	double* centera2 = new double[seedNum];
	double* centerb1 = new double[seedNum];
	double* centerb2 = new double[seedNum];
	double* centerx1 = new double[seedNum];
	double* centerx2 = new double[seedNum];
	double* centery1 = new double[seedNum];
	double* centery2 = new double[seedNum];
	double* WSum = new double[seedNum];
	int* clusterSize = new int[seedNum];



	//Initialization
	for (int i = 0; i<seedNum; i++)
	{
		centerL1[i] = 0;
		centerL2[i] = 0;
		centera1[i] = 0;
		centera2[i] = 0;
		centerb1[i] = 0;
		centerb2[i] = 0;
		centerx1[i] = 0;
		centerx2[i] = 0;
		centery1[i] = 0;
		centery2[i] = 0;
		int x = seedArray[i].x;
		int y = seedArray[i].y;
		int minX = (x - StepX / 4 <= 0) ? 0 : x - StepX / 4;
		int minY = (y - StepY / 4 <= 0) ? 0 : y - StepY / 4;
		int maxX = (x + StepX / 4 >= nRows - 1) ? nRows - 1 : x + StepX / 4;
		int maxY = (y + StepY / 4 >= nCols - 1) ? nCols - 1 : y + StepY / 4;
		int Count = 0;
		for (int j = minX; j <= maxX; j++)
			for (int k = minY; k <= maxY; k++)
			{
				Count++;
				centerL1[i] += L1[j][k];
				centerL2[i] += L2[j][k];
				centera1[i] += a1[j][k];
				centera2[i] += a2[j][k];
				centerb1[i] += b1[j][k];
				centerb2[i] += b2[j][k];
				centerx1[i] += x1[j][k];
				centerx2[i] += x2[j][k];
				centery1[i] += y1[j][k];
				centery2[i] += y2[j][k];
			}
		centerL1[i] /= Count;
		centerL2[i] /= Count;
		centera1[i] /= Count;
		centera2[i] /= Count;
		centerb1[i] /= Count;
		centerb2[i] /= Count;
		centerx1[i] /= Count;
		centerx2[i] /= Count;
		centery1[i] /= Count;
		centery2[i] /= Count;
	}



	//K-means
	for (int iteration = 0; iteration <= iterationNum; iteration++)
	{
		for (int i = 0; i<nRows; i++)
			for (int j = 0; j<nCols; j++)
				dist[i][j] = DBL_MAX;




		int minX, minY, maxX, maxY; double D;
		for (int i = 0; i<seedNum; i++)
		{
			int x = seedArray[i].x;
			int y = seedArray[i].y;
			minX = (x - (StepX) <= 0) ? 0 : x - StepX;
			minY = (y - (StepY) <= 0) ? 0 : y - StepY;
			maxX = (x + (StepX) >= nRows - 1) ? nRows - 1 : x + StepX;
			maxY = (y + (StepY) >= nCols - 1) ? nCols - 1 : y + StepY;
			for (int m = minX; m <= maxX; m++)
				for (int n = minY; n <= maxY; n++)
				{
					D = (L1[m][n] - centerL1[i])*(L1[m][n] - centerL1[i]) +
						(L2[m][n] - centerL2[i])*(L2[m][n] - centerL2[i]) +
						(a1[m][n] - centera1[i])*(a1[m][n] - centera1[i]) +
						(a2[m][n] - centera2[i])*(a2[m][n] - centera2[i]) +
						(b1[m][n] - centerb1[i])*(b1[m][n] - centerb1[i]) +
						(b2[m][n] - centerb2[i])*(b2[m][n] - centerb2[i]) +
						(x1[m][n] - centerx1[i])*(x1[m][n] - centerx1[i]) +
						(x2[m][n] - centerx2[i])*(x2[m][n] - centerx2[i]) +
						(y1[m][n] - centery1[i])*(y1[m][n] - centery1[i]) +
						(y2[m][n] - centery2[i])*(y2[m][n] - centery2[i]);
					if (D<dist[m][n])
					{
						label[m*nCols + n] = i;
						dist[m][n] = D;
					}
				}
		}



		for (int i = 0; i<seedNum; i++)
		{
			centerL1[i] = 0;
			centerL2[i] = 0;
			centera1[i] = 0;
			centera2[i] = 0;
			centerb1[i] = 0;
			centerb2[i] = 0;
			centerx1[i] = 0;
			centerx2[i] = 0;
			centery1[i] = 0;
			centery2[i] = 0;
			WSum[i] = 0;
			clusterSize[i] = 0;
			seedArray[i].x = 0;
			seedArray[i].y = 0;
		}



		for (int i = 0; i<nRows; i++)
		{
			for (int j = 0; j<nCols; j++)
			{
				int L = label[i*nCols + j];
				double Weight = W[i][j];
				centerL1[L] += Weight*L1[i][j];
				centerL2[L] += Weight*L2[i][j];
				centera1[L] += Weight*a1[i][j];
				centera2[L] += Weight*a2[i][j];
				centerb1[L] += Weight*b1[i][j];
				centerb2[L] += Weight*b2[i][j];
				centerx1[L] += Weight*x1[i][j];
				centerx2[L] += Weight*x2[i][j];
				centery1[L] += Weight*y1[i][j];
				centery2[L] += Weight*y2[i][j];
				clusterSize[L]++;
				WSum[L] += Weight;
				seedArray[L].x += i;
				seedArray[L].y += j;
			}
		}
		for (int i = 0; i<seedNum; i++)
		{
			WSum[i] = (WSum[i] == 0) ? 1 : WSum[i];
			clusterSize[i] = (clusterSize[i] == 0) ? 1 : clusterSize[i];
		}
		for (int i = 0; i<seedNum; i++)
		{
			centerL1[i] /= WSum[i];
			centerL2[i] /= WSum[i];
			centera1[i] /= WSum[i];
			centera2[i] /= WSum[i];
			centerb1[i] /= WSum[i];
			centerb2[i] /= WSum[i];
			centerx1[i] /= WSum[i];
			centerx2[i] /= WSum[i];
			centery1[i] /= WSum[i];
			centery2[i] /= WSum[i];
			seedArray[i].x /= clusterSize[i];
			seedArray[i].y /= clusterSize[i];
		}
	}



	//EnforceConnection
	int threshold = (nRows*nCols) / (seedNum*thresholdCoef);
	preEnforceConnectivity(label, nRows, nCols);
	EnforceConnectivity(L1, L2, a1, a2, b1, b2, x1, x2, y1, y2, W, label, threshold, nRows, nCols);



	//Clear Memory
	delete[]centerL1;
	delete[]centerL2;
	delete[]centera1;
	delete[]centera2;
	delete[]centerb1;
	delete[]centerb2;
	delete[]centerx1;
	delete[]centerx2;
	delete[]centery1;
	delete[]centery2;
	delete[]WSum;
	delete[]clusterSize;
	for (int i = 0; i<nRows; i++)
		delete[] dist[i];
	delete[]dist;
}
int LSC::countSuperpixel(unsigned short int* label, int nRows, int nCols)
{
	lscpoint P;
	queue<lscpoint> Q;
	int L;
	int labelNum = 0;
	int x, y, maxX, minX, maxY, minY;
	bool** mask = new bool*[nRows];
	for (int i = 0; i<nRows; i++)
	{
		mask[i] = new bool[nCols];
		for (int j = 0; j<nCols; j++)
			mask[i][j] = false;
	}

	for (int i = 0; i<nRows; i++)
		for (int j = 0; j<nCols; j++)
		{
			if (mask[i][j] == false)
			{
				L = label[i*nCols + j];
				labelNum++;
				label[i*nCols + j] = labelNum;
				P.x = i; P.y = j;
				Q.push(P);
				mask[i][j] = true;
				while (!Q.empty())
				{
					P = Q.front();
					Q.pop();
					x = P.x; y = P.y;
					minX = (x - 1 <= 0) ? 0 : x - 1;
					minY = (y - 1 <= 0) ? 0 : y - 1;
					maxX = (x + 1 >= nRows - 1) ? nRows - 1 : x + 1;
					maxY = (y + 1 >= nCols - 1) ? nCols - 1 : y + 1;
					for (int m = minX; m <= maxX; m++)
						for (int n = minY; n <= maxY; n++)
						{
							if (label[m*nCols + n] == L&&mask[m][n] == false)
							{
								P.x = m; P.y = n;
								mask[m][n] = true;
								label[m*nCols + n] = labelNum;
								Q.push(P);
							}
						}
				}
			}
		}
	for (int i = 0; i<nRows; i++)
		delete[] mask[i];
	delete[] mask;
	return labelNum;
}
void LSC::Initialize(
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
)
{
	float thetaL, thetaa, thetab, thetax, thetay;
	for (int i = 0; i<nRows; i++)
		for (int j = 0; j<nCols; j++)
		{
			thetaL = ((float)L[i*nCols + j] / (float)255)*PI / 2;
			thetaa = ((float)a[i*nCols + j] / (float)255)*PI / 2;
			thetab = ((float)b[i*nCols + j] / (float)255)*PI / 2;
			thetax = ((float)i / (float)StepX)*PI / 2;
			thetay = ((float)j / (float)StepY)*PI / 2;
			L1[i][j] = Color*cos(thetaL);
			L2[i][j] = Color*sin(thetaL);
			a1[i][j] = Color*cos(thetaa)*2.55;
			a2[i][j] = Color*sin(thetaa)*2.55;
			b1[i][j] = Color*cos(thetab)*2.55;
			b2[i][j] = Color*sin(thetab)*2.55;
			x1[i][j] = Distance*cos(thetax);
			x2[i][j] = Distance*sin(thetax);
			y1[i][j] = Distance*cos(thetay);
			y2[i][j] = Distance*sin(thetay);
		}
	double sigmaL1 = 0, sigmaL2 = 0, sigmaa1 = 0, sigmaa2 = 0, sigmab1 = 0, sigmab2 = 0, sigmax1 = 0, sigmax2 = 0, sigmay1 = 0, sigmay2 = 0;
	double size = nRows*nCols;
	for (int i = 0; i<nRows; i++)
		for (int j = 0; j<nCols; j++)
		{
			sigmaL1 += L1[i][j];
			sigmaL2 += L2[i][j];
			sigmaa1 += a1[i][j];
			sigmaa2 += a2[i][j];
			sigmab1 += b1[i][j];
			sigmab2 += b2[i][j];
			sigmax1 += x1[i][j];
			sigmax2 += x2[i][j];
			sigmay1 += y1[i][j];
			sigmay2 += y2[i][j];
		}
	sigmaL1 /= size;
	sigmaL2 /= size;
	sigmaa1 /= size;
	sigmaa2 /= size;
	sigmab1 /= size;
	sigmab2 /= size;
	sigmax1 /= size;
	sigmax2 /= size;
	sigmay1 /= size;
	sigmay2 /= size;
	for (int i = 0; i < nRows; i++) {
		for (int j = 0; j < nCols; j++)
		{
			W[i][j] = L1[i][j] * sigmaL1 +
				L2[i][j] * sigmaL2 +
				a1[i][j] * sigmaa1 +
				a2[i][j] * sigmaa2 +
				b1[i][j] * sigmab1 +
				b2[i][j] * sigmab2 +
				x1[i][j] * sigmax1 +
				x2[i][j] * sigmax2 +
				y1[i][j] * sigmay1 +
				y2[i][j] * sigmay2;
			L1[i][j] /= W[i][j];
			L2[i][j] /= W[i][j];
			a1[i][j] /= W[i][j];
			a2[i][j] /= W[i][j];
			b1[i][j] /= W[i][j];
			b2[i][j] /= W[i][j];
			x1[i][j] /= W[i][j];
			x2[i][j] /= W[i][j];
			y1[i][j] /= W[i][j];
			y2[i][j] /= W[i][j];
		}
	}
}

void LSC::doLSC(unsigned char* R, unsigned char* G, unsigned char* B, int nRows, int nCols, int superpixelnum, double ratio, unsigned short* label)
{
	//Setting Parameter
	float colorCoefficient = 20;
	float distCoefficient = colorCoefficient*ratio;
	int seedNum = superpixelnum;
	int iterationNum = 20;
	int thresholdCoef = 4;

	unsigned char *L, *a, *b;
	L = new unsigned char[nRows*nCols];
	a = new unsigned char[nRows*nCols];
	b = new unsigned char[nRows*nCols];

	myrgb2lab(R, G, B, L, a, b, nRows, nCols);


	//Produce Seeds
	int ColNum, RowNum, StepY, StepX;
	ColNum = sqrt(float(seedNum*nCols / nRows));
	RowNum = seedNum / ColNum;
	StepX = nRows / RowNum;
	StepY = nCols / ColNum;
	lscpoint *seedArray = new lscpoint[seedNum];
	int newSeedNum = Seeds(nRows, nCols, RowNum, ColNum, StepX, StepY, seedNum, seedArray);


	//Initialization
	float **L1, **L2, **a1, **a2, **b1, **b2, **x1, **x2, **y1, **y2;
	double **W;
	L1 = new float*[nRows];
	L2 = new float*[nRows];
	a1 = new float*[nRows];
	a2 = new float*[nRows];
	b1 = new float*[nRows];
	b2 = new float*[nRows];
	x1 = new float*[nRows];
	x2 = new float*[nRows];
	y1 = new float*[nRows];
	y2 = new float*[nRows];
	W = new double*[nRows];
	for (int i = 0; i<nRows; i++)
	{
		L1[i] = new float[nCols];
		L2[i] = new float[nCols];
		a1[i] = new float[nCols];
		a2[i] = new float[nCols];
		b1[i] = new float[nCols];
		b2[i] = new float[nCols];
		x1[i] = new float[nCols];
		x2[i] = new float[nCols];
		y1[i] = new float[nCols];
		y2[i] = new float[nCols];
		W[i] = new double[nCols];
	}
	Initialize(L, a, b, L1, L2, a1, a2, b1, b2, x1, x2, y1, y2, W, nRows, nCols, StepX, StepY, colorCoefficient, distCoefficient);
	delete[] L;
	delete[] a;
	delete[] b;


	//Produce Superpixel
	DoSuperpixel(L1, L2, a1, a2, b1, b2, x1, x2, y1, y2, W, label, seedArray, newSeedNum, nRows, nCols, StepX, StepY, iterationNum, thresholdCoef);
	delete[]seedArray;

	int NUMBER = countSuperpixel(label, nRows, nCols);

	//Clear Memory
	for (int i = 0; i<nRows; i++)
	{
		delete[] L1[i];
		delete[] L2[i];
		delete[] a1[i];
		delete[] a2[i];
		delete[] b1[i];
		delete[] b2[i];
		delete[] x1[i];
		delete[] x2[i];
		delete[] y1[i];
		delete[] y2[i];
		delete[] W[i];

	}
	delete[]L1;
	delete[]L2;
	delete[]a1;
	delete[]a2;
	delete[]b1;
	delete[]b2;
	delete[]x1;
	delete[]x2;
	delete[]y1;
	delete[]y2;
	delete[]W;
}


