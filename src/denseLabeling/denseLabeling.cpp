#ifndef DENSE_LABELING
#define DENSE_LABELING
//
//  denseLabeling.cpp
//  
//
//  Created by Ana B. Cambra on 18/10/16.
//
//

#include <stdio.h>
#include <math.h>
#include <opencv2/calib3d.hpp>
#include "../superpixels/superpixels.cpp"
#include "least-squares-linear-system.h"
using namespace Eigen;
using namespace cv;
typedef Matrix<float, Dynamic, 1> MatrixFX;

class DenseLabeling : public SuperPixels
{

    //SYSTEM
    int numUnknows;
    Optimization::LeastSquaresLinearSystem<float> equations;
    Optimization::LeastSquaresLinearSystem<float> unary;
    Optimization::LeastSquaresLinearSystem<float> equal;
    Optimization::LeastSquaresLinearSystem<float> test;
    Optimization::LeastSquaresLinearSystem<float> test2;
    Optimization::LeastSquaresLinearSystem<float> test3;
    Optimization::LeastSquaresLinearSystem<float> binary;
    Optimization::LeastSquaresLinearSystem<float> binaryCOLOR;
    Optimization::LeastSquaresLinearSystem<float> binaryEdge;

    double w_unary;
    double w_equal;
    double w_color;
    double k = 3;//wb = exp(-k*diff)
    int MAX_BINARIAS = 1000000;
    double _MAX_DIFF_LAB;
    bool use_normal = true;
    int num_u = 0;
    std::vector<std::vector<int>> binaryConnect;

    Mat edgeMap;

public:

    DenseLabeling(string path, double w_u = 0.4, double w_e = 0.3, double w_c = 1.0, double diff_lab = 15.0) : SuperPixels(path),
        w_unary(w_u), w_equal(w_e), w_color(w_c), _MAX_DIFF_LAB(diff_lab)
    {
        numUnknows = ((this->maxID + 1) * 6);

        this->unary = Optimization::LeastSquaresLinearSystem<float>((maxID + 1) * 6);
        this->equal = Optimization::LeastSquaresLinearSystem<float>((maxID + 1) * 6);
        this->binaryEdge = Optimization::LeastSquaresLinearSystem<float>((maxID + 1) * 6);

        this->edgeMap = Mat::zeros(_image.rows, _image.cols, CV_8UC1);
        binaryConnect.clear();

    }
    bool setNormalFlag(bool flag) {
        return flag;
    }
    Mat* solve()
    {
        Optimization::LeastSquaresLinearSystem<float> unaryN = unary; unaryN.normalize();
        Optimization::LeastSquaresLinearSystem<float> equalN = equal; equalN.normalize();


        Optimization::LeastSquaresLinearSystem<float> B = binary + binaryEdge;

        B.normalize();


        equations = (unaryN * w_unary) + (B * (1.0 - w_unary - w_equal)) + (equalN * w_equal);

        equations.normalize();
        MatrixFX x((maxID + 1) * 6);

        x = equations.solve();
        //equations.toFile((maxID + 1) * 6);
        Mat  final = Mat::ones(_labels.rows, _labels.cols, CV_32FC1);

        for (int i = 0; i < final.cols; i++)
        {
            for (int j = 0; j < final.rows; j++)
            {

                int id = getIdFromPixel(i, j);
                final.at<float>(j, i) = (float)x(id * 6, 0)*getCoefficient(i, j, 0)
                        + (float)x(id * 6 + 1, 0)*getCoefficient(i, j, 1)
                        + (float)x(id * 6 + 2, 0)*getCoefficient(i, j, 2)
                        + (float)x(id * 6 + 3, 0)*getCoefficient(i, j, 3)
                        + (float)x(id * 6 + 4, 0)*getCoefficient(i, j, 4)
                        + (float)x(id * 6 + 5, 0)*getCoefficient(i, j, 5);
                _labels.at<float>(j, i) = final.at<float>(j, i);
            }
        }
        double min, max;
        cv::minMaxLoc(final, &min, &max);
        Mat depth(final.size(), CV_32FC1);
        depth = (float(max) - final) / (max-min);
        Mat* result = new Mat[2];//real depth depth(0-1)
        result[0]=final;
        result[1]=depth;
        return result;
    }

    //UNARY EQUATIONS
    void clearUnaries() { unary.clear(); }
    //EQUAL EQUATIONS
    void clearEqual() { equal.clear(); }


    void addGradient_Equal(int x1, int y1, int x2, int y2)
    {
        if (x1 < 0 || x1 > this->_image.cols-1 || y1 < 0 || y1 >this->_image.rows-1)
            return;
        if (x2 < 0 || x2 > this->_image.cols-1 || y2 < 0 || y2 >this->_image.rows-1)
            return;
        int x3=2*x2-x1;
        int y3=2*y2-y1;
        int id1 = getIdFromPixel(x1, y1);
        int id2 = getIdFromPixel(x2, y2);
        int id3 = getIdFromPixel(x3, y3);
        if (x3 < 0 || x3 > this->_image.cols-1 || y3 < 0 || y3 >this->_image.rows-1)
            return;


        Optimization::SparseEquation<float> eq4;

        eq4.a[id1*6+0]=getCoefficient(x1,y1,0);
        eq4.a[id1*6+1]=getCoefficient(x1,y1,1);
        eq4.a[id1*6+2]=getCoefficient(x1,y1,2);
        eq4.a[id1*6+3]=getCoefficient(x1,y1,3);
        eq4.a[id1*6+4]=getCoefficient(x1,y1,4);
        eq4.a[id1*6+5]=getCoefficient(x1,y1,5);

        eq4.a[id3*6+0]+=getCoefficient(x3,y3,0);
        eq4.a[id3*6+1]+=getCoefficient(x3,y3,1);
        eq4.a[id3*6+2]+=getCoefficient(x3,y3,2);
        eq4.a[id3*6+3]+=getCoefficient(x3,y3,3);
        eq4.a[id3*6+4]+=getCoefficient(x3,y3,4);
        eq4.a[id3*6+5]+=getCoefficient(x3,y3,5);

        eq4.a[id2*6+0]-=getCoefficient(x2,y2,0)*2;
        eq4.a[id2*6+1]-=getCoefficient(x2,y2,1)*2;
        eq4.a[id2*6+2]-=getCoefficient(x2,y2,2)*2;
        eq4.a[id2*6+3]-=getCoefficient(x2,y2,3)*2;
        eq4.a[id2*6+4]-=getCoefficient(x2,y2,4)*2;
        eq4.a[id2*6+5]-=getCoefficient(x2,y2,5)*2;

        eq4.b=0;
        equal.add_equation(eq4);

    }

    void addNormal_Equal(int x1, int y1, int x2, int y2)
    {
        //Does not add it if any of them are -1 or if both are equal.
        if (x1 < 0 || x1 > this->_image.cols || y1 < 0 || y1 >this->_image.rows)
            return;
        if (x2 < 0 || x2 > this->_image.cols || y2 < 0 || y2 >this->_image.rows)
            return;

        int id1=getIdFromPixel(x1, y1);
        int id2=getIdFromPixel(x2, y2);


        Optimization::SparseEquation<float> eq,eq1;

        //normal_X:
        eq.a[id1 * 6 + 0] = 2.0f * getCoefficient(x1, y1, 3);
        eq.a[id1 * 6 + 1] = getCoefficient(x1, y1, 4);
        eq.a[id1 * 6 + 2] = 0;
        eq.a[id1 * 6 + 3] = getCoefficient(x1, y1, 5);
        eq.a[id1 * 6 + 4] = 0;
        eq.a[id1 * 6 + 5] = 0;

        eq.a[id2 * 6 + 0] -= 2.0f * getCoefficient(x2, y2, 3);
        eq.a[id2 * 6 + 1] -= getCoefficient(x2, y2, 4);
        eq.a[id2 * 6 + 2] -= 0;
        eq.a[id2 * 6 + 3] -= getCoefficient(x2, y2, 5);
        eq.a[id2 * 6 + 4] -= 0;
        eq.a[id2 * 6 + 5] -= 0;
        eq.b = 0.0;
        equal.add_equation(eq);
        //normal_Y:
        eq1.a[id1 * 6 + 0] = 0;
        eq1.a[id1 * 6 + 1] = getCoefficient(x1, y1, 3);
        eq1.a[id1 * 6 + 2] = 2.0f * getCoefficient(x1, y1, 4);
        eq1.a[id1 * 6 + 3] = 0;
        eq1.a[id1 * 6 + 4] = getCoefficient(x1, y1, 5);
        eq1.a[id1 * 6 + 5] = 0;

        eq1.a[id2 * 6 + 0] -= 0;
        eq1.a[id2 * 6 + 1] -= getCoefficient(x2, y2, 3);
        eq1.a[id2 * 6 + 2] -= 2.0f * getCoefficient(x2, y2, 4);
        eq1.a[id2 * 6 + 3] -= 0;
        eq1.a[id2 * 6 + 4] -= getCoefficient(x2, y2, 5);
        eq1.a[id2 * 6 + 5] -= 0;
        eq1.b = 0.0;
        equal.add_equation(eq1);
    }

    void addEquation_Equal(int x1, int y1, int x2, int y2)
    {
        //Does not add it if any of them are -1 or if both are equal.
        if (x1 < 0 || x1 > this->_image.cols-1 || y1 < 0 || y1 >this->_image.rows-1)
            return;
        if (x2 < 0 || x2 > this->_image.cols-1 || y2 < 0 || y2 >this->_image.rows-1)
            return;


        int id1 = getIdFromPixel(x1, y1);
        int id2 = getIdFromPixel(x2, y2);

        Optimization::SparseEquation<float> eq;

        eq.a[id1 * 6 + 0] = getCoefficient(x1, y1, 0);
        eq.a[id1 * 6 + 1] = getCoefficient(x1, y1, 1);
        eq.a[id1 * 6 + 2] = getCoefficient(x1, y1, 2);
        eq.a[id1 * 6 + 3] = getCoefficient(x1, y1, 3);
        eq.a[id1 * 6 + 4] = getCoefficient(x1, y1, 4);
        eq.a[id1 * 6 + 5] = getCoefficient(x1, y1, 5);

        eq.a[id2 * 6 + 0] -= getCoefficient(x2, y2, 0);//use '-=' instead of '='
        eq.a[id2 * 6 + 1] -= getCoefficient(x2, y2, 1);
        eq.a[id2 * 6 + 2] -= getCoefficient(x2, y2, 2);
        eq.a[id2 * 6 + 3] -= getCoefficient(x2, y2, 3);
        eq.a[id2 * 6 + 4] -= getCoefficient(x2, y2, 4);
        eq.a[id2 * 6 + 5] -= getCoefficient(x2, y2, 5);
        eq.b = 0.0;
        equal.add_equation(eq);


    }

    void addEquation_Unary(int x, int y, float li, bool saveInput = false)
    {



        if (x < 0 || x > this->_image.cols || y < 0 || y >this->_image.rows)
            return;

        int id = getIdFromPixel(x, y);

        _labelsInput.at<float>(y, x) = li;



        if (li  > 0.0)
        {

            //COSTE UNARIO
            Optimization::SparseEquation<float> eq;
            eq.a[id * 6 + 0] = getCoefficient(x, y, 0);
            eq.a[id * 6 + 1] = getCoefficient(x, y, 1);
            eq.a[id * 6 + 2] = getCoefficient(x, y, 2);
            eq.a[id * 6 + 3] = getCoefficient(x, y, 3);
            eq.a[id * 6 + 4] = getCoefficient(x, y, 4);
            eq.a[id * 6 + 5] = getCoefficient(x, y, 5);
            eq.b = li;
            unary.add_equation(eq);
            num_u++;
        }


    }



    //compare color pixel boundaries
    //the new algorithm will use addEquations_BinariesBoundariesPerPixel() to add binary equation
    void addEquations_BinariesBoundariesPerPixel()
    {

        this->binary = Optimization::LeastSquaresLinearSystem<float>((maxID + 1) * 6);

        Mat nonZeroCoordinates;

        Mat sobel = _sobel;


        dilate(_sobel, sobel, Mat(), Point(-1, -1), 2, 1, 1);

        findNonZero((sobel != 0), nonZeroCoordinates);

        int num_b = 0;

        int size = maxID + 1;

        Mat _lab, noConnected;
        cv::cvtColor(this->_image, _lab, CV_BGR2Lab);
        cv::cvtColor(_lab, noConnected, CV_Lab2BGR);
        cv::cvtColor(noConnected, noConnected, CV_BGR2GRAY);
        cv::cvtColor(noConnected, noConnected, CV_GRAY2RGB);

        Mat mlab = _lab;



        for (int n = 0; n < (int)nonZeroCoordinates.total(); n++)
        {
            int x = nonZeroCoordinates.at<Point>(n).x;
            int y = nonZeroCoordinates.at<Point>(n).y;
            int id1 = getIdFromPixel(x, y);

            for (int i = x; i <= (x + 1); i++)
            {
                for (int j = y-1; j <= (y + 1); j++)
                {
                    int id_v = getIdFromPixel(i, j);

                    if ((i >= 0 && i <= _sobel.cols - 1) &&(j >= 0 && j <= _sobel.rows - 1)&&(id1 != id_v)&&(!(j==y-1&&i==x)))
                    {
                        float diff = sqrt(
                                    (((double)(mlab.at<Vec3b>(y, x)[0]) - (double)(mlab.at<Vec3b>(j, i)[0]))*
                                ((double)(mlab.at<Vec3b>(y, x)[0]) - (double)(mlab.at<Vec3b>(j, i)[0])))
                                + (((double)(mlab.at<Vec3b>(y, x)[1]) - (double)(mlab.at<Vec3b>(j, i)[1]))*
                                ((double)(mlab.at<Vec3b>(y, x)[1]) - (double)(mlab.at<Vec3b>(j, i)[1])))
                                + (((double)(mlab.at<Vec3b>(y, x)[2]) - (double)(mlab.at<Vec3b>(j, i)[2]))*
                                ((double)(mlab.at<Vec3b>(y, x)[2]) - (double)(mlab.at<Vec3b>(j, i)[2])))
                                );


                        double wb,wn;
                        {
                            wb = exp(-diff*k);
                            wn = exp(-(diff*diff) / 50);
                            std::vector<int> temp1,temp2;
                            temp1.push_back(x); temp1.push_back(y); temp1.push_back(i);temp1.push_back(j);
                            binaryConnect.push_back(temp1);
                        }

                        Optimization::SparseEquation<float> eq10,eq11,eq20,eq21,eq22,eq23;
                        eq10.a[id1 * 6 + 0] = getCoefficient(x, y, 0)*wb;
                        eq10.a[id1 * 6 + 1] = getCoefficient(x, y, 1)*wb;
                        eq10.a[id1 * 6 + 2] = getCoefficient(x, y, 2)*wb;
                        eq10.a[id1 * 6 + 3] = getCoefficient(x, y, 3)*wb;
                        eq10.a[id1 * 6 + 4] = getCoefficient(x, y, 4)*wb;
                        eq10.a[id1 * 6 + 5] = getCoefficient(x, y, 5)*wb;

                        eq10.a[id_v * 6 + 0] = -getCoefficient(x, y, 0)*wb;
                        eq10.a[id_v * 6 + 1] = -getCoefficient(x, y, 1)*wb;
                        eq10.a[id_v * 6 + 2] = -getCoefficient(x, y, 2)*wb;
                        eq10.a[id_v * 6 + 3] = -getCoefficient(x, y, 3)*wb;
                        eq10.a[id_v * 6 + 4] = -getCoefficient(x, y, 4)*wb;
                        eq10.a[id_v * 6 + 5] = -getCoefficient(x, y, 5)*wb;

                        binary.add_equation(eq10);

                        eq11.a[id1 * 6 + 0] = getCoefficient(i, j, 0)*wb;
                        eq11.a[id1 * 6 + 1] = getCoefficient(i, j, 1)*wb;
                        eq11.a[id1 * 6 + 2] = getCoefficient(i, j, 2)*wb;
                        eq11.a[id1 * 6 + 3] = getCoefficient(i, j, 3)*wb;
                        eq11.a[id1 * 6 + 4] = getCoefficient(i, j, 4)*wb;
                        eq11.a[id1 * 6 + 5] = getCoefficient(i, j, 5)*wb;

                        eq11.a[id_v * 6 + 0] = -getCoefficient(i, j, 0)*wb;
                        eq11.a[id_v * 6 + 1] = -getCoefficient(i, j, 1)*wb;
                        eq11.a[id_v * 6 + 2] = -getCoefficient(i, j, 2)*wb;
                        eq11.a[id_v * 6 + 3] = -getCoefficient(i, j, 3)*wb;
                        eq11.a[id_v * 6 + 4] = -getCoefficient(i, j, 4)*wb;
                        eq11.a[id_v * 6 + 5] = -getCoefficient(i, j, 5)*wb;

                        binary.add_equation(eq11);


                        eq20.a[id1 * 6 + 0] = wn * 2.0f * getCoefficient(x, y, 3);
                        eq20.a[id1 * 6 + 1] = wn*getCoefficient(x, y, 4);
                        eq20.a[id1 * 6 + 2] = 0;
                        eq20.a[id1 * 6 + 3] = wn*getCoefficient(x, y, 5);
                        eq20.a[id1 * 6 + 4] = 0;
                        eq20.a[id1 * 6 + 5] = 0;

                        eq20.a[id_v * 6 + 0] = -wn * 2.0f * getCoefficient(x, y, 3);
                        eq20.a[id_v * 6 + 1] = -wn*getCoefficient(x, y, 4);
                        eq20.a[id_v * 6 + 2] = 0;
                        eq20.a[id_v * 6 + 3] = -wn*getCoefficient(x, y, 5);
                        eq20.a[id_v * 6 + 4] = 0;
                        eq20.a[id_v * 6 + 5] = 0;

                        binary.add_equation(eq20);

                        eq21.a[id1 * 6 + 0] = wn * 2.0f * getCoefficient(i, j, 3);
                        eq21.a[id1 * 6 + 1] = wn*getCoefficient(i, j, 4);
                        eq21.a[id1 * 6 + 2] = 0;
                        eq21.a[id1 * 6 + 3] = wn*getCoefficient(i, j, 5);
                        eq21.a[id1 * 6 + 4] = 0;
                        eq21.a[id1 * 6 + 5] = 0;

                        eq21.a[id_v * 6 + 0] = -wn * 2.0f * getCoefficient(i, j, 3);
                        eq21.a[id_v * 6 + 1] = -wn*getCoefficient(i, j, 4);
                        eq21.a[id_v * 6 + 2] = 0;
                        eq21.a[id_v * 6 + 3] = -wn*getCoefficient(i, j, 5);
                        eq21.a[id_v * 6 + 4] = 0;
                        eq21.a[id_v * 6 + 5] = 0;

                        binary.add_equation(eq21);

                        eq22.a[id1 * 6 + 0] = 0;
                        eq22.a[id1 * 6 + 1] = wn*getCoefficient(x, y, 3);
                        eq22.a[id1 * 6 + 2] = wn * 2.0f * getCoefficient(x, y, 4);
                        eq22.a[id1 * 6 + 3] = 0;
                        eq22.a[id1 * 6 + 4] = wn*getCoefficient(x, y, 5);
                        eq22.a[id1 * 6 + 5] = 0;

                        eq22.a[id_v * 6 + 0] = 0;
                        eq22.a[id_v * 6 + 1] = -wn*getCoefficient(x, y, 3);
                        eq22.a[id_v * 6 + 2] = -wn * 2.0f *getCoefficient(x, y, 4);
                        eq22.a[id_v * 6 + 3] = 0;
                        eq22.a[id_v * 6 + 4] = -wn*getCoefficient(x, y, 5);
                        eq22.a[id_v * 6 + 5] = 0;

                        binary.add_equation(eq22);

                        eq23.a[id1 * 6 + 0] = 0;
                        eq23.a[id1 * 6 + 1] = wn*getCoefficient(i, j, 3);
                        eq23.a[id1 * 6 + 2] = wn * 2.0f * getCoefficient(i, j, 4);
                        eq23.a[id1 * 6 + 3] = 0;
                        eq23.a[id1 * 6 + 4] = wn*getCoefficient(i, j, 5);
                        eq23.a[id1 * 6 + 5] = 0;

                        eq23.a[id_v * 6 + 0] = 0;
                        eq23.a[id_v * 6 + 1] = -wn*getCoefficient(i, j, 3);
                        eq23.a[id_v * 6 + 2] = -wn * 2.0f * getCoefficient(i, j, 4);
                        eq23.a[id_v * 6 + 3] = 0;
                        eq23.a[id_v * 6 + 4] = -wn*getCoefficient(i, j, 5);
                        eq23.a[id_v * 6 + 5] = 0;

                        binary.add_equation(eq23);


                        num_b++;

                    }

                }
            }
        }

        //printf("these are %d binary equation\n",num_b);
    }


    //edge constaint
    void clearEdge()
    {
        edgeMap.setTo(0, edgeMap);
        binaryEdge.clear();
    }



    void addEdge(int x, int y,int x_prev,int y_prev,int dist=1)
    {
        {
            if (x < 0 || x > this->_image.cols || y < 0 || y >this->_image.rows)
                return;

            _labelsInput.at<float>(y, x) = 1;


            Mat _lab, noConnected;
            cv::cvtColor(this->_image, _lab, CV_BGR2Lab);
            cv::cvtColor(_lab, noConnected, CV_Lab2BGR);
            cv::cvtColor(noConnected, noConnected, CV_BGR2GRAY);
            cv::cvtColor(noConnected, noConnected, CV_GRAY2RGB);


            Mat mlab = _lab;

            edgeMap.at<uchar>(x, y) = 1;


            int id1 = getIdFromPixel(x, y);
            int i = 0, j = 0;
            for (i = x-2; i <= (x + 2); i++)
            {
                for (j = y-2; j <= (y + 2); j++)
                {
                    if ((i >= 0 && i <= edgeMap.cols - 1) && (j >= 0 && j <= edgeMap.rows - 1))
                    {
                        int id2 = getIdFromPixel(i, j);

                        float diff = sqrt(
                                    (((double)(mlab.at<Vec3b>(y, x)[0]) - (double)(mlab.at<Vec3b>(j, i)[0]))*
                                ((double)(mlab.at<Vec3b>(y, x)[0]) - (double)(mlab.at<Vec3b>(j, i)[0])))
                                + (((double)(mlab.at<Vec3b>(y, x)[1]) - (double)(mlab.at<Vec3b>(j, i)[1]))*
                                ((double)(mlab.at<Vec3b>(y, x)[1]) - (double)(mlab.at<Vec3b>(j, i)[1])))
                                + (((double)(mlab.at<Vec3b>(y, x)[2]) - (double)(mlab.at<Vec3b>(j, i)[2]))*
                                ((double)(mlab.at<Vec3b>(y, x)[2]) - (double)(mlab.at<Vec3b>(j, i)[2])))
                                );
                        double wb,wn;

                        wb = -exp(-diff*k);
                        wn = -exp(-(diff*diff) / 50);

                        Optimization::SparseEquation<float> eq,eq2,eq20,eq21,eq22,eq23;

                        eq.a[id1 * 6 + 0] = getCoefficient(x, y, 0)*wb;
                        eq.a[id1 * 6 + 1] = getCoefficient(x, y, 1)*wb;
                        eq.a[id1 * 6 + 2] = getCoefficient(x, y, 2)*wb;
                        eq.a[id1 * 6 + 3] = getCoefficient(x, y, 3)*wb;
                        eq.a[id1 * 6 + 4] = getCoefficient(x, y, 4)*wb;
                        eq.a[id1 * 6 + 5] = getCoefficient(x, y, 5)*wb;

                        eq.a[id2 * 6 + 0] -=getCoefficient(x, y, 0)*wb;
                        eq.a[id2 * 6 + 1] -=getCoefficient(x, y, 1)*wb;
                        eq.a[id2 * 6 + 2] -=getCoefficient(x, y, 2)*wb;
                        eq.a[id2 * 6 + 3] -=getCoefficient(x, y, 3)*wb;
                        eq.a[id2 * 6 + 4] -=getCoefficient(x, y, 4)*wb;
                        eq.a[id2 * 6 + 5] -=getCoefficient(x, y, 5)*wb;

                        binaryEdge.add_equation(eq);

                        eq2.a[id1 * 6 + 0] = getCoefficient(i, j, 0)*wb;
                        eq2.a[id1 * 6 + 1] = getCoefficient(i, j, 1)*wb;
                        eq2.a[id1 * 6 + 2] = getCoefficient(i, j, 2)*wb;
                        eq2.a[id1 * 6 + 3] = getCoefficient(i, j, 3)*wb;
                        eq2.a[id1 * 6 + 4] = getCoefficient(i, j, 4)*wb;
                        eq2.a[id1 * 6 + 5] = getCoefficient(i, j, 5) * wb;

                        eq2.a[id2 * 6 + 0] -=getCoefficient(i, j, 0)*wb;
                        eq2.a[id2 * 6 + 1] -=getCoefficient(i, j, 1)*wb;
                        eq2.a[id2 * 6 + 2] -=getCoefficient(i, j, 2)*wb;
                        eq2.a[id2 * 6 + 3] -=getCoefficient(i, j, 3)*wb;
                        eq2.a[id2 * 6 + 4] -=getCoefficient(i, j, 4)*wb;
                        eq2.a[id2 * 6 + 5] -=getCoefficient(i, j, 5) * wb;

                        binaryEdge.add_equation(eq2);



                        eq20.a[id1 * 6 + 0] = wn * 2.0f * getCoefficient(x, y, 3);
                        eq20.a[id1 * 6 + 1] = wn*getCoefficient(x, y, 4);
                        eq20.a[id1 * 6 + 2] = 0;
                        eq20.a[id1 * 6 + 3] = wn*getCoefficient(x, y, 5);
                        eq20.a[id1 * 6 + 4] = 0;
                        eq20.a[id1 * 6 + 5] = 0;

                        eq20.a[id2 * 6 + 0] -=wn * 2.0f * getCoefficient(x, y, 3);
                        eq20.a[id2 * 6 + 1] -=wn*getCoefficient(x, y, 4);
                        eq20.a[id2 * 6 + 2] -=0;
                        eq20.a[id2 * 6 + 3] -=wn*getCoefficient(x, y, 5);
                        eq20.a[id2 * 6 + 4] -=0;
                        eq20.a[id2 * 6 + 5] -=0;

                        binary.add_equation(eq20);

                        eq21.a[id1 * 6 + 0] = wn * 2.0f * getCoefficient(i, j, 3);
                        eq21.a[id1 * 6 + 1] = wn*getCoefficient(i, j, 4);
                        eq21.a[id1 * 6 + 2] = 0;
                        eq21.a[id1 * 6 + 3] = wn*getCoefficient(i, j, 5);
                        eq21.a[id1 * 6 + 4] = 0;
                        eq21.a[id1 * 6 + 5] = 0;

                        eq21.a[id2 * 6 + 0] -=wn * 2.0f * getCoefficient(i, j, 3);
                        eq21.a[id2 * 6 + 1] -=wn*getCoefficient(i, j, 4);
                        eq21.a[id2 * 6 + 2] -= 0;
                        eq21.a[id2 * 6 + 3] -=wn*getCoefficient(i, j, 5);
                        eq21.a[id2 * 6 + 4] -= 0;
                        eq21.a[id2 * 6 + 5] -=0;

                        binary.add_equation(eq21);

                        eq22.a[id1 * 6 + 0] = 0;
                        eq22.a[id1 * 6 + 1] = wn*getCoefficient(x, y, 3);
                        eq22.a[id1 * 6 + 2] = wn * 2.0f * getCoefficient(x, y, 4);
                        eq22.a[id1 * 6 + 3] = 0;
                        eq22.a[id1 * 6 + 4] = wn*getCoefficient(x, y, 5);
                        eq22.a[id1 * 6 + 5] = 0;
                        eq22.a[id2 * 6 + 0] -= 0;

                        eq22.a[id2 * 6 + 1] -=wn*getCoefficient(x, y, 3);
                        eq22.a[id2 * 6 + 2] -=wn * 2.0f *getCoefficient(x, y, 4);
                        eq22.a[id2 * 6 + 3] -= 0;
                        eq22.a[id2 * 6 + 4] -=wn*getCoefficient(x, y, 5);
                        eq22.a[id2 * 6 + 5] -= 0;

                        binary.add_equation(eq22);

                        eq23.a[id1 * 6 + 0] = 0;
                        eq23.a[id1 * 6 + 1] = wn*getCoefficient(i, j, 3);
                        eq23.a[id1 * 6 + 2] = wn * 2.0f * getCoefficient(i, j, 4);
                        eq23.a[id1 * 6 + 3] = 0;
                        eq23.a[id1 * 6 + 4] = wn*getCoefficient(i, j, 5);
                        eq23.a[id1 * 6 + 5] = 0;

                        eq23.a[id2 * 6 + 0] -=0;
                        eq23.a[id2 * 6 + 1] -=wn*getCoefficient(i, j, 3);
                        eq23.a[id2 * 6 + 2] -=wn * 2.0f * getCoefficient(i, j, 4);
                        eq23.a[id2 * 6 + 3] -= 0;
                        eq23.a[id2 * 6 + 4] -=wn*getCoefficient(i, j, 5);
                        eq23.a[id2 * 6 + 5] -= 0;

                        binary.add_equation(eq23);


                    }

                }
            }
        }
        {
            if(x==x_prev&&y==y_prev)
                return;
            if (x < 0 || x > this->_image.cols-1 || y < 0 || y >this->_image.rows-1)
                return;
            if (x_prev < 0 || x_prev > this->_image.cols-1 || y_prev < 0 || y_prev >this->_image.rows-1)
                return;
            pair<int,int> bigPoint,smallPoint;
            if(x==x_prev)
            {
                smallPoint.second=bigPoint.second=y;
                if(y>y_prev)
                {
                    bigPoint.first=x-1;
                    smallPoint.first=x+1;
                }
                else
                {
                    bigPoint.first=x+1;
                    smallPoint.first=x-1;
                }
            }
            else if(y==y_prev)
            {
                smallPoint.first=bigPoint.first=x;
                if(x>x_prev)
                {
                    bigPoint.second=y+1;
                    smallPoint.second=y-1;
                }
                else
                {
                    bigPoint.second=y-1;
                    smallPoint.second=y+1;
                }
            }
            else
            {
                if((x==x_prev-1&&y==y_prev+1)||(x==x_prev+1&&y==y_prev-1))
                {
                    bigPoint.first=x;
                    bigPoint.second=y_prev;
                    smallPoint.first=x_prev;
                    smallPoint.second=y;
                }
                else
                {
                    bigPoint.first=x_prev;
                    bigPoint.second=y;
                    smallPoint.first=x;
                    smallPoint.second=y_prev;
                }

            }
            int idBigPoint=getIdFromPixel(bigPoint.first,bigPoint.second);
            int idSmallPoint=getIdFromPixel(smallPoint.first,smallPoint.second);
            if(idBigPoint==idSmallPoint)
                return;

            Optimization::SparseEquation<float> eq;
            double wb=0.5;
            eq.a[idBigPoint * 6 + 0] = getCoefficient(bigPoint.first, bigPoint.second, 0)*wb;
            eq.a[idBigPoint * 6 + 1] = getCoefficient(bigPoint.first, bigPoint.second, 1)*wb;
            eq.a[idBigPoint * 6 + 2] = getCoefficient(bigPoint.first, bigPoint.second, 2)*wb;
            eq.a[idBigPoint * 6 + 3] = getCoefficient(bigPoint.first, bigPoint.second, 3)*wb;
            eq.a[idBigPoint * 6 + 4] = getCoefficient(bigPoint.first, bigPoint.second, 4)*wb;
            eq.a[idBigPoint * 6 + 5] = getCoefficient(bigPoint.first, bigPoint.second, 5)*wb;

            eq.a[idSmallPoint * 6 + 0] -=getCoefficient(smallPoint.first, smallPoint.second, 0)*wb;
            eq.a[idSmallPoint * 6 + 1] -=getCoefficient(smallPoint.first, smallPoint.second, 1)*wb;
            eq.a[idSmallPoint * 6 + 2] -=getCoefficient(smallPoint.first, smallPoint.second, 2)*wb;
            eq.a[idSmallPoint * 6 + 3] -=getCoefficient(smallPoint.first, smallPoint.second, 3)*wb;
            eq.a[idSmallPoint * 6 + 4] -=getCoefficient(smallPoint.first, smallPoint.second, 4)*wb;
            eq.a[idSmallPoint * 6 + 5] -=getCoefficient(smallPoint.first, smallPoint.second, 5)*wb;

            eq.b=dist;
            binary.add_equation(eq);
        }
    }


};

#endif //DENSE_LABELING
