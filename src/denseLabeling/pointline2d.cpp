
#include <iostream>
#include <math.h>
namespace PointLine {
	class CVecter2D
	{
	public:
		double x, y;
		CVecter2D() { x = y = 0.0; }
		CVecter2D(double fx, double fy) { x = fx; y = fy; }
		//-----------operator---------------------
	public:
		CVecter2D operator + (const CVecter2D& vecter)
		{
			return CVecter2D(x + vecter.x, y + vecter.y);
		}
		CVecter2D operator - (const CVecter2D& vecter) {
			return CVecter2D(x - vecter.x, y - vecter.y);
		}
		CVecter2D operator * (const double& fScalar) {
			return CVecter2D(x*fScalar, y*fScalar);
		}
		CVecter2D operator / (const double& fScalar) {
			if (fScalar) { return CVecter2D(x / fScalar, y / fScalar); }
			else { return CVecter2D(0, 0); }
		}
		void operator =  (const CVecter2D& vecter) {
			x = vecter.x;
			y = vecter.y;
		}
		bool operator ==  (const CVecter2D& vecter) {
			if (x == vecter.x && y == vecter.y)return true;
			return false;
		}
		void operator += (const CVecter2D& vecter) {
			x += vecter.x;
			y += vecter.y;
		}
		void operator -= (const CVecter2D& vecter) {
			x -= vecter.x;
			y -= vecter.y;
		}
		void operator *= (const double& fScalar) {
			x *= fScalar;
			y *= fScalar;
		}
		void operator /= (const double& fScalar) {
			if (fScalar) {
				x /= fScalar;
				y /= fScalar;
			}
			else { x = y = 0; }
		}
		//-----------Operation------------------------
	public:
		void	Normalized() {
			double fLength = GetLength();
			x = x / fLength;
			y = y / fLength;
		}
		double	Dot(const CVecter2D& vecter) {
			return (x*vecter.x + y*vecter.y);
		}
		double  GetAngle(const CVecter2D& vecter) {
			CVecter2D vec1 = vecter;
			CVecter2D vec2(x, y);
			if (vec1.GetLength() && vec2.GetLength())
			{
				vec1.Normalized();
				vec2.Normalized();
				return acos(vec1.Dot(vec2));
			}
			else { return 0; }
		}
		int	IsVertical(const CVecter2D& vecter)  //0--90,1-->90,-1--<90
		{
			double fDot = Dot(vecter);
			if (fDot == 0)return 0;
			if (fDot < 0)return -1;
			return 1;
		}
		bool IsParallel(const CVecter2D& vecter) {
			if ((x*vecter.y - y*vecter.x) == 0)return true;
			return false;
		}
		void Inverted()  //invert
		{
			x = -x;
			y = -y;
		}
		CVecter2D GetVertical(bool bClockwise = false)  //
		{
			if (bClockwise)
			{
				return CVecter2D(y, -x);
			}
			return CVecter2D(-y, x);
		}
		//-----------Properties------------------------
	public:
		bool	IsAffine()  //
		{
			if ((x + y) == 1)return true;
			return false;
		}
		bool	IsConvex() {
			if (IsAffine())
			{
				if (x >= 0 && y >= 0) { return true; }
			}
			return false;
		}
		double	GetLength() {
			return sqrt(x*x + y*y);
		}
		// function definition

	};
	class CPoint2D {
	public:
		CPoint2D() { x = y = 0.0; }
		CPoint2D(double fx, double fy) { x = fx; y = fy; }
		CPoint2D(const CPoint2D& pt) { x = pt.x; y = pt.y; }

	public:
		double x, y;
		//operation
	public:
		void operator =  (const CPoint2D& pt) {
			x = pt.x;
			y = pt.y;
		}
		CPoint2D operator + (const CVecter2D& vecter) {
			return CPoint2D(x + vecter.x, y + vecter.y);
		}
		/*CPoint2D operator - (const CVecter2D& vecter) {
			return CPoint2D(x - vecter.x, y - vecter.y);
		}*/
		CVecter2D operator - (const CPoint2D& pt) {
			return CVecter2D(x - pt.x, y - pt.y);
		}
	};
	class PointLine {
	public:
		double t ;
		double u ;
		double length;
		int nRet;
	
	public:
		PointLine() { 
			t = u = 0;
			length = 1;
			nRet = 1;
		}
		PointLine(CPoint2D p1, CPoint2D p2, CPoint2D p3, CPoint2D p4) {
			nRet = GetTwoLineMSG(p1, p2, p3, p4,t,u);
		}
		PointLine(CPoint2D p1, CPoint2D p2, CPoint2D p3) {
			nRet = GetPointToLineMSG(p1, p2, p3, length);
		}
		int GetTwoLineMSG(CPoint2D line1_pt1, CPoint2D line1_pt2, CPoint2D line2_pt1, CPoint2D line2_pt2, double& t, double& u)
		{
			//point A,B(line 1), point C,D (line 2)
			int nRet = -1;
			CVecter2D b = line1_pt2 - line1_pt1; // line 1
			CVecter2D d = line2_pt2 - line2_pt1; // line 2
			CVecter2D c = line2_pt1 - line1_pt1; // line 2 to line 1

			if (b.GetLength() < 1e-6 || d.GetLength() < 1e-6)
			{
				return -1; //not line 
			}

			if (d.IsParallel(b))  //para
			{
				nRet = -2; //Paralle
				if (d.IsParallel(c)) //d,c parellel
				{
					nRet = -3; //on same line 

					CVecter2D ad = line2_pt2 - line1_pt1; //D-A
					if (c.x != 0 && c.y != 0)
					{
						t = c.y / b.y;  // t>=0 t<=1 C between line 1 ; t<0 C on BA; t>1 C on 
						u = ad.y / b.y; // u>=0 u<=1 
					}

					if (c.x == 0)
					{
						t = c.y / b.y;
						u = ad.y / b.y;
					}
					if (c.y == 0)
					{
						t = c.x / b.x;
						u = ad.x / b.x;
					}
				}
				return nRet;
			}

			t = (c.x*d.y - c.y*d.x) / (b.x*d.y - b.y*d.x);  //0<=t<=1, intersection point is in AB; t<0 --- on extended AB; t>1,on extended BA
			u = (c.y*b.x - c.x*b.y) / (d.x*b.y - d.y*b.x);  // u CD

															//ptResult = line1_pt1 + b*t; //intersection point
			return 0;
		}
		int GetPointToLineMSG(CPoint2D pt, CPoint2D line_pt1, CPoint2D line_pt2, double& fLength)
		{
			CVecter2D v = line_pt2 - line_pt1;
			CVecter2D c = pt - line_pt1;
			double K = c.Dot(v) / v.Dot(v);
			CVecter2D Kv = v*K;  //

			int nRet = -1;  //pedal on reverse extended line v
			double fAngle = Kv.Dot(v);
			if (fAngle > 0)
			{
				nRet = 1;  //pedal on reverse extended line v
				if (Kv.GetLength() < v.GetLength())
				{
					nRet = 0;  //pedal is on line v
				}
			}

			CVecter2D vT = v.GetVertical();
			double M = c.Dot(vT) / vT.Dot(vT);
			CVecter2D MvT = vT*M;  //

			fLength = MvT.GetLength();  //distance from pt to line v
										//ptResult = pt - vT*M;  //pedal vector

			return nRet;//0---on the line;-1/1 --- one the exteded line
		}
	};
	
	
}