#pragma once

#include <QPixmap>
#include <QImage>
#include <string>
#include <opencv2/core/core.hpp>
#include <imageHDR.h>
#include <memory>
#include <filters/filtertonemappingmantiuk.h>

class Path 
{
    const char* st;
    std::string aux;
public:
    Path(const char* s) : st(s) { }

    const char* extension()
    {
        const char* p=st+strlen(st)-1;
        while ((p>st)&&(p[0]!='.')) p--;
        if (p[0]=='.') p++;
        return p;
    }

    const char* c_str() const { return st; }

    const char* file()
    {
        const char* p=st+strlen(st)-1;
        while ((p>st)&&((p[0]!='/')&&(p[0]!='\\'))) p--;
        if ((p[0]!='/')||(p[0]!='\\')) p++;
        return p;
    }

    static int strcasecmp(const char* s1, const char* s2)
    {
        while (*s1 != '\0' && tolower(*s1) == tolower(*s2))
        {	s1++; s2++;     }

        return tolower(*(unsigned char *) s1) - tolower(*(unsigned char *) s2);
    }

    bool extension_is(const char* ext)
    {
        return (strcasecmp(extension(), ext) == 0);
    }
};

inline void load_image(const char* filename, const std::shared_ptr<cv::Mat>& image)
{
    if (image)
    {
        if (Path(filename).extension_is("hdr"))
            *image=load_hdr(filename);
        else
            *image=cv::imread(filename);
    }
}

inline std::shared_ptr<cv::Mat> load_image(const char* filename)
{
    std::shared_ptr<cv::Mat> sol;
    if (Path(filename).extension_is("hdr")) sol=std::make_shared<cv::Mat>(load_hdr(filename));
    else sol=std::make_shared<cv::Mat>(cv::imread(filename));
    if (!sol->data) return std::shared_ptr<cv::Mat>();
    else return sol;
}

inline QPixmap qpixmap(const cv::Mat& im) {
    if (im.channels()>=3) {
        cv::Mat image;
        if (im.type()==CV_32FC3) { //If this is an HDR image we have to convert it to a 8 bit format
            cv::Mat dummy(im.rows, im.cols, CV_32F);
            dummy = cv::Scalar(1.0);
            double min, max;
            cv::minMaxLoc(im, &min, &max);
            image = FilterTonemappingMantiuk::tonemap(im, -dummy, dummy, 1, 2);
            image.convertTo(image,CV_8UC3);
            cv::cvtColor(image,image,CV_BGR2RGB);
        } else {
            cv::cvtColor(im, image, CV_BGR2RGB);
        }

        return QPixmap::fromImage(QImage(image.data, image.cols, image.rows, image.step, QImage::Format_RGB888));
    } else {
        double min, max;
        cv::minMaxLoc(im, &min, &max);
        cv::Mat image = im * (255.0/max);
            image.convertTo(image,CV_8UC1);
        return QPixmap::fromImage(QImage(image.data, image.cols, image.rows, image.step, QImage::Format_Grayscale8));
    }
}

inline void cvmat_to_qpixmap(const std::shared_ptr<cv::Mat>& mat, const std::shared_ptr<QPixmap>& pixmap)
{
    if (mat && pixmap) {
        *pixmap = qpixmap(*mat);
    }
}

