#pragma once

#include "labelimage.h"
#include "imageformats.h"

class ImageMouseBrush : public QObject
{
    Q_OBJECT
private:
    std::shared_ptr<cv::Mat>  image;   //original image
    std::shared_ptr<QPixmap>  pixmap;  //pixmap with painted strokes
    std::string path;
    bool canEdit;
    bool paint;
    QColor colorBlush;
    QPainter painter;
    FILE *fp;
public:
    ImageMouseBrush(const std::shared_ptr<cv::Mat>& _image = std::shared_ptr<cv::Mat>(), const std::shared_ptr<QPixmap>& _pixmap = std::shared_ptr<QPixmap>()) :
        image(_image), pixmap(_pixmap) { }

    void clear()
    {
        if (image)
        {
            if (!pixmap)
                pixmap = std::make_shared<QPixmap>(qpixmap(*image)); //Local ownership
            else cvmat_to_qpixmap(image, pixmap);
        }
    }

    void setImage(const std::shared_ptr<cv::Mat>& im)
    {
        image = im;
        clear();
    }

    void setPixmap(const std::shared_ptr<QPixmap>& p) //Just for the pointer
    {
        pixmap=p;
    }
    void setPath(std::string _path)
    {
        path=_path;
    }

    const std::shared_ptr<QPixmap>& getPixmap() const { return pixmap; }
    const std::shared_ptr<cv::Mat>& getImage()  const { return image; }


    //DRAW over image
    void setColorBrush(QColor color)
    {
        //cambiar color del pixel del raton
        colorBlush = color;
        canEdit    = true;
    }

    void setColorBlush(int r, int g, int b)
    {
        setColorBrush(QColor::fromRgb(r,g,b));
    }

    void unsetBrush()
    {
        canEdit=false;
    }

    void connectToLabelImage(const LabelImage* li)
    {
        if (li)
        {
            QObject::connect(li, &LabelImage::mousePixelDown,
                             this, &ImageMouseBrush::pixelDown);
            QObject::connect(li, &LabelImage::mousePixelUp,
                             this, &ImageMouseBrush::pixelUp);
            QObject::connect(li, &LabelImage::mousePixelChanged,
                             this, &ImageMouseBrush::pixelMoving);
        }

    }

    void pixelMoving(int x, int y)
    {
        if (canEdit && paint && image)
        {
            QPainter painter(pixmap.get());
            painter.setPen(Qt::NoPen);
            painter.setBrush(QBrush(colorBlush));
            painter.drawEllipse(x, y,10.0,10.0);
        }
    }

    void pixelUp(int x, int y,QMouseEvent *event)//end of clicking
    {
        paint = false;
    }

    void pixelDown(int x, int y,QMouseEvent *event)//begin clicking
    {
        if (event->button()== Qt::LeftButton)
        {
            if (!pixmap) clear();
            paint = true;
        }
        else if (event->button() == Qt::RightButton)
        {
            clear();
        }
    }


};


