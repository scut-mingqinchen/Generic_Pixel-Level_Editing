#pragma once

#include <QWidget>
#include <QScrollArea>
#include <QApplication>
#include <QLabel>
#include <QLineEdit>
#include <QSlider>
#include <QGroupBox>
#include <QDir>
#include <QTime>
#include <QDebug>
#include "labelimage.h"
#include "imagemousebrush.h"
#include "multiimageviewer.h"
#include "../filters/filter.h"
//#include <opencv2/contrib/contrib.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <math.h>
#include <opencv2/imgproc/imgproc.hpp>
#include <cstdio>
#include <cstdlib>
#include "../denseLabeling/denseLabeling.cpp"
#include "imageHDR.h"
const int slider_ticks = 1000;

class Brush {


public:
    int x_prev, y_prev;
    bool isCrossArea;
    Brush(int x_prev=-1,int y_prev=-1,bool isCrossArea=false):x_prev(x_prev),y_prev(y_prev),isCrossArea(isCrossArea){}
    virtual void onClicked(int x, int y) { }
    virtual void onMoved(int x, int y)   { }
    virtual bool shouldDrawStroke() const { return false; }

};

class BrushEditChannel : public Brush
{
    const std::vector<DenseLabeling*>& labels;
    std::vector<bool>& edited;
    std::list<std::tuple<float, int>> values;
public:
    BrushEditChannel(const std::vector<DenseLabeling*>& _labels,
                     std::vector<bool>& _edited, const std::list<std::tuple<float, int>>& _values) :
        labels(_labels), edited(_edited), values(_values) { }
    BrushEditChannel(const std::vector<DenseLabeling*>& _labels,
                     std::vector<bool>& _edited, float _value, int _channel) :
        BrushEditChannel(_labels, _edited, std::list<std::tuple<float,int>>{{std::make_tuple(_value, _channel)}}) { }

    void onClicked(int x, int y) override
    {
        for (auto t : values)
        {
            float value;
            int channel;
            std::tie(value, channel) = t;
            if (!edited[channel])
            { //We update the edited thing only on clicked (a bit extra efficiency)
                edited[channel] = true;
                labels[channel]->clearEqual();
                //printf("clearEqual\n");
                labels[channel]->clearUnaries(); //We remove the initial equation that sets up the entire thing.
            }
        }
    }

    void onMoved(int x, int y) override { //This event also happens when clicking.
        for (auto t : values)
        {
            float value; int channel;
            std::tie(value, channel) = t;
            labels[channel]->addEquation_Unary(x,y,value);
        }
    }

    bool shouldDrawStroke() const override { return true; }
};


class BrushUserDefined:public QWidget,public Brush
{
    Q_OBJECT
private:
    float value;
    const std::vector<DenseLabeling*>& labels;
    std::vector<bool>& edited;
    std::list<int> channels;
public:
    BrushUserDefined(const std::vector<DenseLabeling*>& _labels,std::vector<bool>& _edited,
                     const std::list<int> _channels,float _value=-1) :
        labels(_labels), edited(_edited),channels(_channels),value(_value)  { std::cout<<"new new new\n";}

    void onClicked(int x, int y) override
    {
        for (int channel : channels)
        {
            if (!edited[channel])
            { //We update the edited thing only on clicked (a bit extra efficiency)
                edited[channel] = true;
                labels[channel]->clearEqual();
                //printf("clearEqual\n");
                labels[channel]->clearUnaries(); //We remove the initial equation that sets up the entire thing.
            }
        }
    }

    void onMoved(int x, int y) override { //This event also happens when clicking.
        for (int channel : channels)
        {
            labels[channel]->addEquation_Unary(x,y,value);
        }
    }
    bool shouldDrawStroke() const override { return true; }
    void setValue(float value)
    {
        this->value=value;
    }
public slots:
    void valueChange(QString data)
    {
        setValue(data.toFloat());
    }


};

class BrushSetNormalEqualChannel : public Brush
{
    const std::vector<DenseLabeling*>& labels;
    std::list<int> channels;

public:
    BrushSetNormalEqualChannel(const std::vector<DenseLabeling*>& _labels,
                               const std::list<int> _channels,bool isCrossArea,int x_prev=-1,int y_prev=-1) :Brush(x_prev,y_prev,isCrossArea),
        labels(_labels), channels(_channels)  { }

    void onClicked(int x, int y) override
    {

        //printf("onClicked BrushSetGradientEqualChannel (%d,%d) ID %d\n",x,y,labels[0]->getIdFromPixel(x, y));
        for (int channel : channels)
        {
            labels[channel]->addNormal_Equal(x, y, x_prev, y_prev);
        }
        x_prev = x; y_prev = y;
    }

    void onMoved(int x, int y) override {
        for (int channel : channels) {
            // printf("onClicked BrushSetGradientEqualChannel (%d,%d) ID %d\n",x,y,labels[0]->getIdFromPixel(x, y));
            labels[channel]->addNormal_Equal(x, y, x_prev, y_prev); //Does not add it if any of them are -1 or if both are equal.
        }//
        x_prev = x; y_prev = y;
    }

    bool shouldDrawStroke() const override { return true; }
};


class BrushSetGradientEqualChannel : public Brush
{
    const std::vector<DenseLabeling*>& labels;
    std::list<int> channels;

public:
    BrushSetGradientEqualChannel(const std::vector<DenseLabeling*>& _labels,
                               const std::list<int> _channels,bool isCrossArea,int x_prev=-1,int y_prev=-1) :Brush(x_prev,y_prev,isCrossArea),
        labels(_labels), channels(_channels)  { }

    void onClicked(int x, int y) override
    {

        //printf("onClicked BrushSetGradientEqualChannel (%d,%d) ID %d\n",x,y,labels[0]->getIdFromPixel(x, y));
        for (int channel : channels)
        {
            labels[channel]->addGradient_Equal(x, y, x_prev, y_prev);
        }
        x_prev = x; y_prev = y;
    }

    void onMoved(int x, int y) override {
        for (int channel : channels) {
            // printf("onClicked BrushSetGradientEqualChannel (%d,%d) ID %d\n",x,y,labels[0]->getIdFromPixel(x, y));
            labels[channel]->addGradient_Equal(x, y, x_prev, y_prev); //Does not add it if any of them are -1 or if both are equal.
        }//
        x_prev = x; y_prev = y;
    }

    bool shouldDrawStroke() const override { return true; }
};



class BrushSetEqualChannel : public Brush
{
    const std::vector<DenseLabeling*>& labels;
    std::list<int> channels;

public:
    BrushSetEqualChannel(const std::vector<DenseLabeling*>& _labels,
                         const std::list<int> _channels,bool isCrossArea,int x_prev=-1,int y_prev=-1):Brush(x_prev,y_prev,isCrossArea),
        labels(_labels), channels(_channels)  { }

    void onClicked(int x, int y) override
    {

        for (int channel : channels)
        {

            labels[channel]->addEquation_Equal(x, y, x_prev, y_prev); //Does not add it if any of them are -1 or if both are equal.
        }//*/
        x_prev = x; y_prev = y;

    }

    void onMoved(int x, int y) override {
        for (int channel : channels) {

            labels[channel]->addEquation_Equal(x, y, x_prev, y_prev); //Does not add it if any of them are -1 or if both are equal.
        }//
        x_prev = x; y_prev = y;
    }

    bool shouldDrawStroke() const override { return true; }
};
class BrushEdgeChannel :public QWidget, public Brush {
    Q_OBJECT
private:
    const std::vector<DenseLabeling*>& labels;
    std::vector<bool>& edited;
    std::list<int> channels;
    float value;
public:
    BrushEdgeChannel(const std::vector<DenseLabeling*>& _labels, std::vector<bool>& _edited, const std::list<int> _channels) :
        Brush(),labels(_labels), edited(_edited), channels(_channels) {}
    void onClicked(int x, int y) override
    {
        for (int channel : channels)
        {
            if (!edited[channel]) {
                edited[channel] = true;
                labels[channel]->clearEdge();
            }
        }
        x_prev = x; y_prev = y;
    }
    void onMoved(int x, int y) override {
        for (int channel : channels) {
            labels[channel]->addEdge(x_prev,y_prev, x, y,value);
        }
        x_prev = x; y_prev = y;

    }

    bool shouldDrawStroke() const override { return true; }
    void setValue(float value)
    {
        this->value=value;
    }
public slots:
    void valueChange(QString data)
    {
        setValue(data.toFloat());
    }
};


class BrushPickValue : public Brush
{
    std::shared_ptr<cv::Mat> from;
    QSlider* slider;
    float min, max;
public:
    BrushPickValue(const std::shared_ptr<cv::Mat>& _from, QSlider* _slider, float _min, float _max) :
        from(_from), slider(_slider), min(_min), max(_max) { }

    void onMoved(int x, int y) override { //This event also happens when clicking.
        if (from && slider && (x>=0) && (x<from->cols) && (y>=0) && (y<from->rows))
        {
            float v = from->at<float>(y,x,0);
            slider->setValue(int(float(slider_ticks)*(v-min)/(max-min)));
        }
    }
};

class WidgetFilter : public QWidget
{
    Q_OBJECT

private:
    const Filter& filter;

    std::shared_ptr<QPixmap> strokes_pixmap;
    std::shared_ptr<cv::Mat> input_image;
    std::shared_ptr<cv::Mat> filtered_image;
    std::vector<std::shared_ptr<cv::Mat>> propagated_channels;
    std::vector<std::shared_ptr<cv::Mat>> real_depthmap;


    ImageMouseBrush *imageMouseBrush;
    MultiImageViewer *multiImageViewer;
    QButtonGroup *buttonOptions;
    QVBoxLayout  *buttonLayout;
    QLabel* info;

    int buttonIdStrokes, buttonIdInput, buttonIdFiltered,buttonIdNormal;
    float textValue;
    std::string filename;
    //Datos
    std::vector<DenseLabeling*> labels;

    int button_id;

    std::vector<QSlider*> sliderValues;

    std::vector<std::shared_ptr<Brush>> brushes;
    std::shared_ptr<Brush> activeBrush;

    void chooseButton(int id)
    {
        if ((id>=0) && (id<brushes.size()))
        {
            buttonOptions->checkedButton()->setChecked(false);
            multiImageViewer->getlabelImage()->setButtonId(id);
            if(id==5)
                multiImageViewer->getlabelImage()->setTextValue(textValue);
            int isEndOfOneEqual=-1;
            for(int i=8;i<=9;i++)
            {
                if(buttonOptions->buttons().size()<9)
                    break;
                string buttonName=buttonOptions->buttons().at(i)->text().toStdString();
                int poson=buttonName.find("on");
                if(poson<buttonName.length()&&poson>=0)
                {
                    buttonOptions->buttons().at(i)->setText((buttonName.substr(0,poson)+"off").c_str());
                    isEndOfOneEqual=i;
                    activeBrush->x_prev=-1;
                    activeBrush->y_prev=-1;
                    break;
                }
            }
            if((id==8||id==9)&&isEndOfOneEqual!=id)
            {
                    std::string buttonName=buttonOptions->checkedButton()->text().toStdString();
                    int posoff=buttonName.find("off");
                    if(posoff<buttonName.length()&&posoff>=0)
                        buttonOptions->checkedButton()->setText((buttonName.substr(0,posoff)+"on").c_str());
            }

            this->setCursor(Qt::PointingHandCursor); //Maybe have a look at the cursors somewhen
            activeBrush=brushes[id];
            if (activeBrush->shouldDrawStroke())
            {

                if(id<=4&&id>=0)
                {
                    imageMouseBrush->setColorBrush(QColor::fromHsv((360-30*id)%360,255,255));
                }
                else if(id ==5)
                {

                    float maxlen=6;
                    imageMouseBrush->setColorBrush(QColor::fromHsv((240+min(120,int(textValue/maxlen*120)))%360,255,255));
                }
                else if(6<=id&&id<=7)
                {
                    imageMouseBrush->setColorBrush(QColor::fromHsv(60+(id-6)*60,255,255));
                }
                else if(8<=id&&id<=9)
                {
                    srand(time(NULL));
                    int sp =(int(rand()))%60;
                    imageMouseBrush->setColorBrush(QColor::fromHsv(60+(id-9)*60+sp,255,255));
                }
                else if(id==10)
                {
                    imageMouseBrush->setColorBrush(QColor::fromHsv(180,255,255));
                }
            }
            else
                imageMouseBrush->unsetBrush();
        }
    }

    void updateSlider(int size)
    {
        if ((input_image) && (!input_image->empty())) processImage();
    }

    std::vector<bool> edited;

public:
    WidgetFilter(const Filter& f) : filter(f),
        strokes_pixmap(std::make_shared<QPixmap>()),
        input_image(std::make_shared<cv::Mat>()),
        filtered_image(std::make_shared<cv::Mat>()),
        propagated_channels(filter.propagatedValues().size()),
        real_depthmap(filter.propagatedValues().size()),
        imageMouseBrush(new ImageMouseBrush(input_image, strokes_pixmap)),
        multiImageViewer(new MultiImageViewer()),
        buttonOptions(new QButtonGroup()),
        buttonLayout(new QVBoxLayout()), info(new QLabel()),
        labels(filter.propagatedValues().size()),
        button_id(0),
        edited(filter.propagatedValues().size(),false)
    {
        for (std::shared_ptr<cv::Mat>& channel : propagated_channels)
            channel = std::make_shared<cv::Mat>();
        for (std::shared_ptr<cv::Mat>& channel3 : real_depthmap)
            channel3 = std::make_shared<cv::Mat>();
        QVBoxLayout *mainLayout = new QVBoxLayout();
        QHBoxLayout *layoutH = new QHBoxLayout();
        mainLayout->addLayout(layoutH);
        this->setMinimumSize(400,400);
        QVBoxLayout *sideBar = new QVBoxLayout();
        QGroupBox *brushBox = new QGroupBox(tr("Brushes"));
        brushBox->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        brushBox->setLayout(buttonLayout);
        sideBar->addWidget(brushBox);
        QVBoxLayout *sliderLayout = new QVBoxLayout();
        QGroupBox *sliderBox = new QGroupBox(tr("Parameters"));
        sliderBox->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);


        sliderBox->setLayout(sliderLayout);
        sideBar->addWidget(sliderBox);
        sideBar->addStretch();
        layoutH->addLayout(sideBar);
        layoutH->addWidget(multiImageViewer);
        QVBoxLayout *image = new QVBoxLayout;
        image->addWidget(multiImageViewer);
        //BUTTONS FOR EQUAL
        std::list<int> list_of_channels;
        for (int c = 0; c<filter.propagatedValues().size(); ++c) list_of_channels.push_back(c);

        //BUTTONS FOR STROKES
        for (auto stroke : filter.strokes())
        {
            QPushButton* button = new QPushButton(QString(stroke.name().c_str()));
            button->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
            button->setCheckable(true);
            buttonLayout->addWidget(button);
            buttonOptions->addButton(button, button_id);
            brushes.push_back(std::make_shared<BrushEditChannel>(labels, edited, stroke.values()));
            if (button_id == 0) {
                button->setChecked(true);
                chooseButton(0);
            }
            button_id++;
        }
        QLineEdit* lineEdit=new QLineEdit();
        {
            QPushButton* button = new QPushButton(QString("User-Defined"));
            button->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
            button->setCheckable(true);
            buttonLayout->addWidget(button);
            buttonOptions->addButton(button, button_id);
            auto pBrushUserDefined=std::make_shared<BrushUserDefined>(labels, edited, list_of_channels);
            brushes.push_back(pBrushUserDefined);
            if (button_id == 0) {
                button->setChecked(true);
                chooseButton(0);
            }
            button_id++;
            buttonLayout->addWidget(lineEdit);

            QObject::connect(lineEdit,&QLineEdit::textChanged,&(*pBrushUserDefined),&BrushUserDefined::valueChange);

            QObject::connect(lineEdit,&QLineEdit::textChanged,this,&WidgetFilter::valueChanged);
        }

        //SEPARATOR
        QFrame* line = new QFrame();
        line->setFrameShape(QFrame::HLine);
        line->setFrameShadow(QFrame::Sunken);
        buttonLayout->addWidget(line);

        for (int b = 0; b<2; ++b)
        {
            QPushButton* button;
            if(!b)
                button = new QPushButton(QString(("value similarity")));
            else
                button = new QPushButton(QString(("gradient similarity")));
            button->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
            button->setCheckable(true);
            buttonLayout->addWidget(button);
            buttonOptions->addButton(button, button_id);
            if(!b)
                brushes.push_back(std::make_shared<BrushSetEqualChannel>(labels, list_of_channels,false));
            else
                brushes.push_back(std::make_shared<BrushSetGradientEqualChannel>(labels, list_of_channels,false));
            if (button_id == 0)
            {
                button->setChecked(true);
                chooseButton(0);
            }
            button_id++;
        }

        for (int b = 0; b<2; ++b)
        {
            QPushButton* button;
            if(!b)
                button = new QPushButton(QString(("value similarity off")));
            else
                button = new QPushButton(QString(("gradient similarity off")));
            button->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
            button->setCheckable(true);
            buttonLayout->addWidget(button);
            buttonOptions->addButton(button, button_id);
            if(!b)
                brushes.push_back(std::make_shared<BrushSetEqualChannel>(labels, list_of_channels,true));
            else
                brushes.push_back(std::make_shared<BrushSetNormalEqualChannel>(labels, list_of_channels,true));
            if (button_id == 0)
            {
                button->setChecked(true);
                chooseButton(0);
            }
            button_id++;
        }

        QFrame* line2 = new QFrame();
        line2->setFrameShape(QFrame::HLine);
        line2->setFrameShadow(QFrame::Sunken);
        buttonLayout->addWidget(line2);
        //BUTTONES FOR GEOMETRY
        std::list<int> list_of_channels2;
        for (int c = 0; c < filter.propagatedValues().size(); ++c) list_of_channels2.push_back(c);
        string buttonName = "Edge";
        QPushButton* button = new QPushButton(QString(buttonName.c_str()));
        button->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
        button->setCheckable(true);
        buttonLayout->addWidget(button);

        buttonOptions->addButton(button, button_id);
        auto pBrushEdgeChannel=std::make_shared<BrushEdgeChannel>(labels, edited, list_of_channels);
        brushes.push_back(pBrushEdgeChannel);

        QObject::connect(lineEdit,&QLineEdit::textChanged,&(*pBrushEdgeChannel),&BrushEdgeChannel::valueChange);
        if (button_id == 0) {
            button->setChecked(true);
            chooseButton(0);
        }
        button_id++;
        for (auto fv : filter.floatValues())
        {
            QLabel* labelValue = new QLabel(fv.name().c_str());

            sliderLayout->addWidget(labelValue);
            QSlider* sliderValue = new QSlider(Qt::Horizontal);
            sliderValue->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
            sliderValue->setTickPosition(QSlider::TicksAbove);
            sliderValue->setTickInterval(slider_ticks/10);
            sliderValue->setRange(0,slider_ticks);
            sliderValue->setValue(slider_ticks/2);
            sliderValue->setEnabled(true);
            sliderValue->setTracking(false);
            //          sliderFocus->setFixedSize(180,20);
            sliderValues.push_back(sliderValue);
            sliderLayout->addWidget(sliderValue);

            QObject::connect(sliderValue,&QSlider::valueChanged,this,&WidgetFilter::updateSlider);

            if (fv.isPickable())
            {
                QPushButton* button = new QPushButton(QString("^ Pick ^"));
                button->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
                button->setCheckable(true);
                sliderLayout->addWidget(button);
                buttonOptions->addButton(button, button_id++);
                brushes.push_back(std::make_shared<BrushPickValue>(propagated_channels[fv.channelToPickFrom()],sliderValue, fv.min(), fv.max()));
            }
        }


        //add label de info
        info->setFrameStyle(QFrame::Panel | QFrame::Sunken);
        info->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
        mainLayout->addWidget(info);

        setLayout(mainLayout);

        imageMouseBrush->connectToLabelImage(multiImageViewer->labelImage());
        buttonIdInput    = multiImageViewer->add("Input",input_image);
        buttonIdStrokes  = multiImageViewer->add("Strokes",imageMouseBrush->getPixmap());
        buttonIdFiltered = multiImageViewer->add("Filtered",filtered_image);
        for (int i=0;i<propagated_channels.size();++i) {
            multiImageViewer->add(filter.propagatedValues()[i].name().c_str(),propagated_channels[i]);
        }
        //conectar señales
        //click sobre opcion, cambiar color pincel
        QObject::connect(buttonOptions, static_cast<void(QButtonGroup::*)(int)>(&QButtonGroup::buttonClicked),
                         this, &WidgetFilter::chooseButton);

        QObject::connect(multiImageViewer->labelImage(),&LabelImage::setwidgetButtonId,this,&WidgetFilter::chooseButton);

        QObject::connect(multiImageViewer->labelImage(), &LabelImage::mousePixelDown,
                         this, &WidgetFilter::clickMousePixel);
        QObject::connect(multiImageViewer->labelImage(), &LabelImage::mousePixelUp,//mousuPixelUp->unclickMOusePixel->processImage();
                         this, &WidgetFilter::unclickMousePixel);
        QObject::connect(multiImageViewer->labelImage(), &LabelImage::mousePixelChanged,
                         this, &WidgetFilter::updatePixel);
        resize(sizeHint());
    }

    void setInfo(QString sms)
    {
        info->setText(sms);
        QApplication::processEvents();
    }

    void resizeEvent(QResizeEvent* event)
    {
        if ((input_image) && (!input_image->empty())) processImage();
    }

private:
    void solveAll()
    {
        for (int i = 0; i<labels.size();++i)
        {
            Mat* results = new Mat[2];
            results = labels[i]->solve();
            *(real_depthmap[i])=results[0];
            *(propagated_channels[i])=results[1];
        }
    }

    std::vector<float> filterParameters() const {
        std::vector<float> floatValues(sliderValues.size());
        for (int i = 0; i<sliderValues.size(); ++i)
        {
            float min = filter.floatValues()[i].min();
            float max = filter.floatValues()[i].max();
            floatValues[i]=(float(sliderValues[i]->value())/float(slider_ticks))*(max-min) + min;
        }
        return floatValues;
    }


public:
    //////////////// SUPERPIXELS + BINARY EQUATIONS
    bool loadData(const QString& filename)
    {
        load_image(filename.toStdString().c_str(), input_image);
        imshow("loadData1",*input_image);
        this->filename=filename.toStdString();
        imageMouseBrush->clear();
        imageMouseBrush->setPath(filename.toStdString());

        for (int i = 0; i < edited.size(); ++i) edited[i] = false;
        //mostrar imagen en la interfaz
        multiImageViewer->setButton(buttonIdInput);
        multiImageViewer->getlabelImage()->setPath(filename.toStdString());

        for (int i = 0; i < labels.size(); ++i)
        {
            if (labels[i]) delete labels[i];
            labels[i] = new DenseLabeling(filename.toStdString());//,0.3,0.3,0.99,10.0);
            labels[i]->addEquations_BinariesBoundariesPerPixel();
            labels[i]->addEquation_Unary(0,0,filter.propagatedValues()[i].defaultValue());

        }
        solveAll();
        *filtered_image = filter.apply(*input_image, propagated_channels, filterParameters());
        setInfo("Binary equations created.");
        return true;
    }

    //////////////// INPUT UNARY EQUATIONS
    void saveData(QString dir)
    {
        processImage(true,dir.toStdString());
        setInfo("All files Saved.");
    }
public slots:
    void valueChanged(QString text)
    {
        textValue=text.toFloat();
    }


private:
    ////////////////

    //////// MOUSE
    void updatePixel(int x, int y)
    {
        if (input_image)
        {
            activeBrush->onMoved(x,y);
            if (activeBrush->shouldDrawStroke())
                multiImageViewer->setButton(buttonIdStrokes);
        }
    }

    void clickMousePixel(int x, int y,QMouseEvent * event)
    {
        if (event->button() == Qt::LeftButton)
        {
            if (input_image)
            {
                activeBrush->onClicked(x,y);
            }
        }
        else if (event->button() == Qt::RightButton)
        {
            if (input_image)
            {
                imageMouseBrush->clear();
                for (DenseLabeling* d : labels)
                {
                    d->clearUnaries();
                    d->addEquation_Unary(0,0,0.5f);
                }
                for (int i = 0; i < edited.size(); ++i)
                    edited[i] = false;
                multiImageViewer->setButton(buttonIdStrokes);
            }
        }
    }

    void unclickMousePixel(int x, int y,QMouseEvent * event,bool isTiggered,QString dir)
    {
        if(!activeBrush->isCrossArea)
        {
            activeBrush->x_prev=-1;
            activeBrush->y_prev=-1;
        }
        processImage();
    }

    void processImage(bool save=false, string dir = ".")
    {

        std::stringstream sstr;
        std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
        solveAll();
        //solve the new unary equation
        std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - start;
        sstr<<"Propagation: "<<std::setprecision(6)<<std::setw(8)<<std::fixed<<elapsed_seconds.count()<<" seconds";
        start = std::chrono::system_clock::now();
        *filtered_image = filter.apply(*input_image, propagated_channels, filterParameters());

        elapsed_seconds = std::chrono::system_clock::now() - start;
        sstr<<"        Image processing: "<<std::setprecision(6)<<std::setw(8)<<std::fixed<<elapsed_seconds.count()<<" seconds";
        setInfo(sstr.str().c_str());
        multiImageViewer->setButton(buttonIdFiltered);
        //compareWithGroundTrue(*read_depthmap[0]);
        if (save)
        {
            int p=filename.find(".");
            dir=(filename.substr(0,p)+"_result");
            std::cout<<dir<<std::endl;
            static int n=0;
            QDir *temp=new QDir;
            if(!(temp->exists(dir.c_str())))
                temp->mkdir(dir.c_str());
            for (DenseLabeling* d : labels)
            {
                Mat user = d->getImageLabelsInput() * 255.0;
                user.convertTo(user, CV_8UC1);
                string name = dir + "/user_input_" + to_string(n) + ".png";
                imwrite(name,user);

                for (int i=0;i<propagated_channels.size();++i)
                {
                    name = dir + "/propagated" + to_string(i) +"_depth" + to_string(n) +".png";
                    Mat sol = *propagated_channels[i] * 255.0;
                    sol.convertTo(sol, CV_8UC1);
                    imwrite(name,sol);
                    name = dir + "/propagated" + to_string(i) + "_normals" + to_string(n) + ".png";

                }

                name = dir + "/filtered" + to_string(n) + ".png";
                imwrite(name,*filtered_image);

                name = dir + "/user_strokes" + to_string(n) + ".png";
                strokes_pixmap->save(QString::fromStdString(name));

                n++;
            }
        }
    }


};

