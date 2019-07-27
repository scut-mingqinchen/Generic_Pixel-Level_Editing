#ifndef SUPERPIXELS_H
#define SUPERPIXELS_H

//#include <stdio.h>
//#include <sys/stat.h>
//#include <unistd.h>

//Opencv
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

using namespace cv;

#include "superpixel.cpp"

#include <time.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include <lib_slic/SLIC.h>
#include <imageHDR.h>
#include <lib_lsc/LSC.h>

#define PI 3.1416
#define SP_MASK 1
#define BB_SP_MASK 0
#define BB_NEIG_MASK 2


using namespace std;

class SuperPixels
{

protected:
    Mat _image; // original image in color BGR
    Mat _imageHDR;
    Mat _ids; // superpixels ids  CV_8UC1
    Mat _sobel; //MASK superpixel boundaries  UCHAR
    Mat _labelsInput; //input labeling

    Mat _labels; // labeling CV_32FC1
    Mat _coefficient[6];  //x^2 xy y^2 x y 1
    // superpixels params
    int _TAM_SP = 30;
    int _NUM_MAX_SP = 1000;

    int NUMLABELS=256;

    unsigned char _DEBUG = 1;

public:

    SuperPixel *_arraySP;
    int maxID;


    //time
    float timeSuperpixels = 0.0;

    SuperPixels()
    {
        maxID=0;
    }
    ~SuperPixels()
    {
        delete[] _arraySP;
        _image.release();
        _ids.release();
        _sobel.release();
        _labels.release();
        _labelsInput.release();
        for(int i=0;i<6;i++)
            _coefficient[i].release();

    }

    static Rect boundingbox(Mat image)
    {
        Mat nonZeroCoordinates;

        if (countNonZero(image) == 0)
            return Rect(0,0,0,0);


        findNonZero( image, nonZeroCoordinates);
        double minX=image.cols+1, minY=image.rows+1, maxX=-1,maxY=-1;

        for (int i = 0; i < nonZeroCoordinates.total(); i++)
        {
            if (nonZeroCoordinates.at<Point>(i).x <= minX) minX = nonZeroCoordinates.at<Point>(i).x;
            else if (nonZeroCoordinates.at<Point>(i).x >= maxX) maxX =  nonZeroCoordinates.at<Point>(i).x;

            if (nonZeroCoordinates.at<Point>(i).y <= minY) minY = nonZeroCoordinates.at<Point>(i).y;
            else if (nonZeroCoordinates.at<Point>(i).y >= maxY) maxY =  nonZeroCoordinates.at<Point>(i).y;
        }


        return Rect(minX,minY,(maxX-minX)+1, (maxY-minY)+1);
    }

    Mat getMask(int id)
    {
        return _arraySP[id].getMask();
    }

    Mat getMaskBB(int id)
    {
        Mat mask =_arraySP[id].getMask().clone();
        rectangle(mask,boundingbox(mask), 255,-1,cv::LINE_AA,0);
        return mask;
    }

    Mat getMaskNeigboursBB(int id)
    {
        set<int> neig = _arraySP[id].getFirstNeighbours();

        std::set<int>::iterator it;

        Mat maskN = Mat::zeros(_image.rows,_image.cols, CV_8UC1);
        for (it=neig.begin(); it!=neig.end(); ++it)
        {
            //concat mask
            bitwise_or(maskN, _arraySP[*it].getMask(),maskN);
        }

        //paint white BB neigbours
        rectangle(maskN,boundingbox(maskN), 255,-1,cv::LINE_AA,0);

        //paint black superpixel
        Mat maskID= getMaskBB(id);

        Mat mask;
        bitwise_xor(maskN,maskID,maskN);

        return maskN.clone();
    }

    void activeDEBUG(){ _DEBUG = 1;}
    void desactiveDEBUG(){ _DEBUG = 0;}

    int numPixels(int id)
    {
        return _arraySP[id].getNumPixels();
    }

    Mat getMaskNeigbours(int id)
    {
        set<int> neig = _arraySP[id].getFirstNeighbours();
        int ln=0;

        std::set<int>::iterator it;

        Mat maskN = Mat::zeros(_image.rows,_image.cols, CV_8UC1);
        for (it=neig.begin(); it!=neig.end(); ++it)
        {
            ln = _arraySP[id].addHistogramLabelSegmentation(_arraySP[*it].getLabelSegmentation());
            //printf("id: %d neig: %d l=%d  %d\n",id1,*it,ln,(int)neig.size());//getchar();
            //concat mask
            bitwise_or(maskN, _arraySP[*it].getMask(),maskN);
        }

        return maskN;

    }

    float getLabel(int x, int y)
    {
        return (float)_labels.at<float>(y,x);
    }
    float getCoefficient(int x,int y,int index)
    {
        return (float)_coefficient[index].at<float>(y,x);
    }
    int getIdFromPixel(int x, int y)
    {
        return (int)_ids.at<float>(y,x);
    }

    bool isNotNullImage()
    {
        return (_image.data != NULL);
    }

    /************************************************************************************/
    /* SuperPixels: obtain superpixels of an image (path)                               */
    /*  load _image (path) and obtain its superpixels
     *      if image_TAM_SP.sp exits -> loadFile
     *      else SLIC & save TAM_SP.sp
     *
     */
    SuperPixels(std::string path)
    {
        clock_t start = clock();
        maxID=0;

        //read image
        try{
            string extension = path.substr(path.find_last_of(".")+1,path.length());
            //printf("extension %s \n %s",extension.c_str(),path.c_str());getchar();
            if (extension == "hdr")
            {
                _imageHDR = load_hdr(path.c_str()).clone();
                //imshow("hdr",_imageHDR);
                _image = _imageHDR.clone();
                cvtColor(_image,_image,CV_BGR2RGB);
                _image.convertTo(_image,CV_8UC3,255.0,0);
            }
            else
                _image = imread(path,CV_LOAD_IMAGE_COLOR);

            if(_image.data == NULL)
            {
                printf("Image %s not found\n",path.c_str());

                return;
            }
            else
                if (_DEBUG == 1) printf("Mat _image CV_8UC1 rows %d cols %d\n",_image.rows,_image.cols);

            _ids= Mat::zeros(_image.rows,_image.cols,CV_32FC1);
            _sobel= Mat::zeros(_image.rows,_image.cols,CV_8UC1);
            _labels= Mat::ones(_ids.rows,_ids.cols,CV_32FC1)*-1;
            int ratio = (_sobel.cols >= _sobel.rows) ? _sobel.cols : _sobel.rows;
            for(int i=0;i<6;i++)
                _coefficient[i]=Mat::zeros(_ids.rows,_ids.cols,CV_32FC1);
            for(int i=0;i<_ids.rows;i++)
                for(int j=0;j<_ids.cols;j++)
                {
                    _coefficient[0].at<float>(i, j) = j*j / (float)(ratio*ratio);
                    _coefficient[1].at<float>(i, j) = j*i / (float)(ratio*ratio);
                    _coefficient[2].at<float>(i, j) = i*i / (float)(ratio*ratio);
                    _coefficient[3].at<float>(i, j) = j / (float)ratio;
                    _coefficient[4].at<float>(i, j) = i / (float)ratio;
                    _coefficient[5].at<float>(i, j) = 1.0f;
                }
            _labelsInput = Mat::zeros(_ids.rows,_ids.cols,CV_32FC1);
        }
        catch(int e)
        {
            printf("Image %s not found\n",path.c_str());
            return;
        }

        size_t found = path.find_last_of(".");
        std::string name = path.substr(0,found) + "_" + std::to_string(_TAM_SP)+".sp";
        FILE *f = fopen(name.c_str(),"r");

        if (f!=NULL)
        {
            fclose(f);

            start = clock();

            loadSuperPixels(name);

            timeSuperpixels = (float) (((double)(clock() - start)) / CLOCKS_PER_SEC);

            if (_DEBUG == 1) printf("**** TIME: load Superpixels: %f seconds\n",(float) (((double)(clock() - start)) / CLOCKS_PER_SEC) );
        }
        else
        {
            //            fclose(f);

            if (_DEBUG == 1) start = clock();

            //To-Do: ./slic_cli  --input test_image/ --contour --superpixels 100 --csv
            //calculateSLICSuperpixels(_image);
            calculateLSCSuperpixels(_image);

            if (_DEBUG == 1) printf("**** TIME: calculate Superpixels: %f seconds\n ",(float) (((double)(clock() - start)) / CLOCKS_PER_SEC) );

            //save superpixels in a file
            superpixels2file(name);

        }

        _arraySP = new SuperPixel[maxID+1];

        calculateBoundariesSuperpixels();
        initializeSuperpixels();

        imshow("superpixels",getImageSuperpixels());
        //waitKey(0);

        return;
    }//SuperPixels


    /*****************************************************************/
    /*****************************************************************/
    void setNUMLABELS(int n)
    {
        NUMLABELS = n;
    }

    Mat getImageLabelsInput() {return _labelsInput;}
    Mat getImageLabels() {return _labels;}
    Mat getImageHDR(){ return _imageHDR;}
    Mat getImage()
    {
        if(_image.data == NULL)
        {
            return Mat::zeros(100, 100, CV_8U);
        }
        else
            return _image;
    }

    Mat getImageSuperpixels(){

        Mat im= _image.clone();
        Scalar* color = new cv::Scalar( 0, 0, 255 );
        im.setTo(*color,_sobel);
        delete color;
        return im;
    }

    Mat paintSuperpixel(Mat image, int id, Scalar *color = new cv::Scalar(0,0,255)){

        Mat im =image.clone();
        // Scalar* color = new cv::Scalar( 0, 0, 255 );
        im.setTo(*color,_arraySP[id].getMask());
        /*imshow("mask",this->getMask(id));
        waitKey(0);
        destroyAllWindows();*/
        delete color;
        return im;
    }//paintSuperpixel

    Mat paintNeighboursSuperpixel(Mat image, int id){

        Mat im =image.clone();
        Scalar* color = new cv::Scalar( 0, 0, 255 );
        im.setTo(*color,_arraySP[id].getMask());

        // its neighbours
        color = new cv::Scalar( 0, 255, 0 );
        set<int> neig = _arraySP[id].getFirstNeighbours();

        std::set<int>::iterator it;
        for (it=neig.begin(); it!=neig.end(); ++it)
        {
            im.setTo(*color,_arraySP[*it].getMask());
        }

        color = new cv::Scalar( 255, 0, 0 );
        neig = _arraySP[id].getSecondNeighbours();

        for (it=neig.begin(); it!=neig.end(); ++it)
        {
            im.setTo(*color,_arraySP[*it].getMask());
        }
        delete color;

        return im;
    }//paintNeighboursSuperpixel

    /*************************************************************************************
     * calculateSLICSuperpixels
     *
     *      mat: BGR 8 bit channels
     *      fill : _ids, maxID
     *      limit total number of superpixels (maxID+1)
     */
    void calculateSLICSuperpixels(Mat mat)
    {
        // Convert matrix to unsigned int array.
        unsigned int* image = new unsigned int[mat.rows*mat.cols];
        unsigned int value = 0x0000;

        for (int i = 0; i < mat.rows; ++i) {
            for (int j = 0; j < mat.cols; ++j) {

                int b = mat.at<cv::Vec3b>(i,j)[0];
                int g = mat.at<cv::Vec3b>(i,j)[1];
                int r = mat.at<cv::Vec3b>(i,j)[2];

                value = 0x0000;
                value |= (0x00FF0000 & (r << 16));
                value |= (0x0000FF00 & (b << 8));
                value |= (0x000000FF & g);

                image[j + mat.cols*i] = value;
            }
        }

        SLIC slic;

        int* segmentation = new int[mat.rows*mat.cols];
        int numberOfLabels = 0;

        //timer.restart();

        int superpixels = (float(mat.rows*mat.cols)/(float)(_TAM_SP*_TAM_SP));
        if(superpixels > _NUM_MAX_SP) superpixels = _NUM_MAX_SP;

        //int superpixels=100;

        double compactness = 40;
        bool perturbseeds = false;
        int iterations = 30;//10;

        slic.DoSuperpixelSegmentation_ForGivenNumberOfSuperpixels(image, mat.cols, mat.rows, segmentation, numberOfLabels, superpixels, compactness, perturbseeds, iterations);

        // Convert labels.
        int** labels = new int*[mat.rows];
        for (int i = 0; i < mat.rows; ++i) {
            labels[i] = new int[mat.cols];

            for (int j = 0; j < mat.cols; ++j) {
                labels[i][j] = segmentation[j + i*mat.cols];
                if (labels[i][j] >= maxID) maxID=labels[i][j];
                _ids.at<float>(i,j)=(float)labels[i][j];
            }
        }

        delete[] image;
        delete[] segmentation;
        delete[] labels;


    }
    /*************************************************************************************
    * calculateLSCSuperpixels
    *
    *      mat: BGR 8 bit channels
    *      fill : _ids, maxID
    *      limit total number of superpixels (maxID+1)
    */
    void calculateLSCSuperpixels(Mat mat) {
        double ratio = 0.075;
        int nRows, nCols;
        unsigned char* R, *G, *B;
        unsigned short* label;
        unsigned short* output;
        nRows = mat.rows;
        nCols = mat.cols;
        int nPixels = nRows*nCols;
        R = new unsigned char[nPixels];
        G = new unsigned char[nPixels];
        B = new unsigned char[nPixels];
        label = new unsigned short[nPixels];
        //output = new unsigned short[nPixels];
        for (int i = 0; i < nRows; i++)
        {
            for (int j = 0; j < nCols; j++) {
                B[j + i*nCols] = mat.at<cv::Vec3b>(i, j)[0];
                G[j + i*nCols] = mat.at<cv::Vec3b>(i, j)[1];
                R[j + i*nCols] = mat.at<cv::Vec3b>(i, j)[2];
            }
        }
        int superpixels = (float(nPixels) / (float)(_TAM_SP*_TAM_SP));
        if (superpixels > _NUM_MAX_SP) superpixels = _NUM_MAX_SP;
        LSC lsc;
        lsc.doLSC(R, G, B,  nRows, nCols, superpixels, ratio, label);
        int** labels = new int*[nRows];
        for (int i = 0; i < nRows; i++) {
            labels[i] = new int[nCols];

            for (int j = 0; j < nCols; j++) {
                labels[i][j] = label[j + i*nCols];
                if (labels[i][j] >= maxID) maxID = labels[i][j];
                _ids.at<float>(i, j) = (float)labels[i][j];
            }
        }
        delete[] labels;
        delete[] label;
        delete[] R;
        delete[] G;
        delete[] B;
    }

    /*************************************************************************************
     * superpixels2file
     */
    void superpixels2file(string nameFile)
    {
        FILE *f;
        int w=0,h=0;
        try{
            f = fopen(nameFile.c_str(),"wb");
            h=_image.rows; w=_image.cols;

            for(int i=0;i<h;i++)
                for (int j=0; j<w; j++)
                {
                    int id =(int) _ids.at<float>(i,j);
                    fwrite(&id, sizeof id, 1, f);
                }

            fclose(f);

        }catch(int e)
        {
            printf("Exception! superpixels2file");
        }

    }//superpixels2file


    /*************************************************************************************
     * loadSLICSuperpixels
     *
     */
    void loadSuperPixels(string path)
    {
        FILE *f;
        int w=0,h=0;
        maxID = 0;

        //ifstream file (path);
        //string current_line;

        /* int i=0;
        int j=0;

        while(getline(file, current_line)){
            // Now inside each line we need to seperate the cols
            stringstream temp(current_line);
            string single_value;
            while(getline(temp,single_value,',')){
                int id =  atoi(single_value.c_str());
                if (id >= maxID) maxID = id;
                 _ids.at<float>(i,j)=(float)id;
                j=j+1;
                if (j == _image.cols)
                {
                    i= i + 1;
                    j = 0;
                }
            }
        }*/



        try{
            f = fopen(path.c_str(),"rb");
            h=_image.rows; w=_image.cols;
            _ids= Mat::zeros(h,w,CV_32FC1);

            for(int i=0;i<h;i++)
                for (int j=0; j<w; j++)
                    if(!feof(f))
                    {
                        int id;
                        fread(&id,sizeof(int),1,f);

                        if (id >= maxID) maxID = id;
                        _ids.at<float>(i,j)=(float)id;

                        //printf("%d %d %s %d\n",i,j,value.c_str() , atoi(value.c_str()));//(int)_ids.at<float>(i,j));
                    }

            fclose(f);

        }catch(int e){
            printf("Exception!");}
    }


    /*************************************************************************************
     * calculateBoundariesSuperpixels()
     *
     *
     */
    void calculateBoundariesSuperpixels()
    {

        //boundaries
        //SOBEL
        int scale = 1;
        int delta = 0;
        int ddepth =-1;// CV_16S;
        /// Generate grad_x and grad_y
        Mat grad_x, grad_y;
        Mat abs_grad_x, abs_grad_y,grad;

        /// Gradient X
        Sobel( _ids, grad_x, ddepth, 1, 0, 3, scale, delta, BORDER_DEFAULT );
        convertScaleAbs( grad_x, abs_grad_x );

        /// Gradient Y
        Sobel( _ids, grad_y, ddepth, 0, 1, 3, scale, delta, BORDER_DEFAULT );
        convertScaleAbs( grad_y, abs_grad_y );

        /// Total Gradient (approximate)
        addWeighted( abs_grad_x, 0.5, abs_grad_y, 0.5, 0, grad );
        // grad = (grad != 0);

        _sobel= (grad != 0);

        grad_x.release();
        grad_y.release();
        abs_grad_x.release();
        abs_grad_y.release();
        grad.release();

    }//calculateBoundariesSuperpixels


    /*************************************************************************************
     * initializeSuperpixels()
     * ALL labels = -1
     * obtain set of neighbours
     */
    void initializeSuperpixels()
    {
        //MASKS

        if (_ids.data != NULL )
        {
            for(int i=0; i<=maxID; i++ )
            {
                Mat mask_sp;
                mask_sp = (_ids == i);
                _arraySP[i].initialize(i,mask_sp, -1);
                mask_sp.release();
            }

            if (_DEBUG == 1) printf("Num superpixels %d\n",maxID);
        }
        else{
            printf("_ids superpixels NULL \n");
            return;
        }

        // NEIGHBOUGRS

        //first boundaries in _sobel
        for (int x = 0; x < _sobel.rows; x++)
            //for (int x = _sobel.rows; x >=0 ; --x)
        {
            for (int y = 0; y < _sobel.cols; y++)
                //for (int y = _sobel.cols; y >=0 ; --y)
            {
                //printf("x:%d y:%d\n",x,y);
                if ( _sobel.at<uchar>(x,y) == 255) //boundarie
                {
                    int id = (int)_ids.at<float>(x,y);

                    //8-neighbours
                    for (int i=(x-1); i<=(x+1) ; i++)
                    {
                        for (int j=(y-1); j<=(y+1); j++)
                        {

                            if (((x!=i) && (y!=j)) &&
                                    (i>=0 && i < _sobel.rows-1) &&
                                    (j>=0 && j < _sobel.cols-1)

                                    //solo 4 vecinas
                                    /* &&(
                                             ((i == (x-1)) && (j == (y)))  ||
                                             ((i == (x+1)) && (j == (y)))  ||
                                             ((i == (x)) && (j == (y-1)))  ||
                                             ((i == (x)) && (j == (y+1))) )//*/
                                    //solo 4 vecinas sigueintes
                                    /* &&(

                                             ((i == (x+1)) && (j == (y)))  ||
                                             ((i == (x)) && (j == (y+1))) )//*/
                                    ){
                                //add neighbours
                                int v = (int)_ids.at<float>(i,j);
                                _arraySP[id].addFirstNeighbour(v);


                            }//if neighbours
                        }//for j
                    }//for i
                }// if boundarie
            }//for y
        }//for x

        //second boundaries
        //for (int id=0; id < maxID+1; id++)
        for (int id=maxID; id >= 0; --id)
        {
            set<int> neig1 = _arraySP[id].getFirstNeighbours();
            set<int>::iterator it1;
            for (it1=neig1.begin(); it1!=neig1.end(); ++it1)
            {
                int v1 = (*it1);
                set<int> neig2 = _arraySP[v1].getFirstNeighbours();
                set<int>::iterator it2;
                for (it2=neig2.begin(); it2!=neig2.end(); ++it2)
                {
                    int v2 = (*it2);
                    bool is_in = neig1.find(v2) != neig1.end();
                    if (!is_in && v2 != id)
                        _arraySP[id].addSecondNeighbours(v2);
                }//for it2
            }//for it1
        }//for id


    }

    //////////

    Mat cropSuperpixel(Mat img, int id, float scale = 1, int maskON = 1)
    {
        //1 mask superpixel
        // 0 BB superpixel
        //2 BB N without BB superpixel

        if (maskON != BB_NEIG_MASK)
        {
            Mat roi = img(boundingbox(getMask(id))).clone();
            Size s = roi.size();
            if (maskON == SP_MASK)
            {
                Mat mask =_arraySP[id].getMask();
                Mat roiMask = mask(boundingbox(getMask(id))).clone(); //0...255
                cvtColor(roiMask,roiMask,CV_GRAY2BGR);

                if (roi.type() != roiMask.type())
                    roi.convertTo(roi, roiMask.type(),255.0,0);

                bitwise_and(roi,roiMask, roi);
            }

            resize(roi, roi, Size(s.width*scale,s.height*scale));

            return roi;
        }else //BB_SP_MASK
        {
            Mat roi = img.clone();
            rectangle(roi,boundingbox(getMask(id)),Scalar(0,0,0),-1,CV_AA,0);
            //imshow("roi",roi);

            Mat roiMask = roi(boundingbox(getMaskNeigboursBB(id))).clone();
            Size s = roiMask.size();
            resize(roiMask, roiMask, Size(s.width*scale,s.height*scale));

            return roiMask;
        }
    }//cropSuperpixel

    void calculateLabelingNeighbours(Mat seg,int numLabels)
    {
        //create image new
        Mat newLabels = Mat::zeros(_image.rows,_image.cols, CV_8UC1);

        for (int id1=0; id1 < maxID+1; id1++)
        {
            //get neighbour
            set<int> neig = _arraySP[id1].getFirstNeighbours();
            int ln=0;

            std::set<int>::iterator it;

            Mat maskN = Mat::zeros(_image.rows,_image.cols, CV_8UC1);
            for (it=neig.begin(); it!=neig.end(); ++it)
            {
                ln = _arraySP[id1].addHistogramLabelSegmentation(_arraySP[*it].getLabelSegmentation());
                //printf("id: %d neig: %d l=%d  %d\n",id1,*it,ln,(int)neig.size());//getchar();
                //concat mask
                bitwise_or(maskN, _arraySP[*it].getMask(),maskN);
            }

            _arraySP[id1].normalizeLabelFirstSegmentation((int)neig.size());
            newLabels.setTo(ln,_arraySP[id1].getMask());

            //create oriented semantic label segmentation
            _arraySP[id1].create_labelOriented(seg,numLabels,getMaskNeigboursBB(id1));
            _arraySP[id1].create_labelOrientedGlobal(seg,numLabels);

        }

        /*labelSet val(60);
        Mat l;// = Mat::zeros(seg.rows,seg.cols,CV_8UC3);
        val._DEBUG=1;
        Mat out = val.paintLabelRandom(newLabels, 60, &l).clone();
        imshow("NEIG SEGMENTATION",out);*/

    }//calculateLabelingNeighbour*/


};

#endif // SUPERPIXELS_H
