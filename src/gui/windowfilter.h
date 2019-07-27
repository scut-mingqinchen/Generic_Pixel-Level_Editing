#pragma once

#include <QMainWindow>
#include <QAction>
#include <QMenu>
#include <QScrollBar>
#include <QMenuBar>
#include <QLabel>
#include <QMessageBox>
#include <QScrollArea>
#include <QDir>
#include <QFileDialog>
#include <QThread>
#include <string.h>

#include "widgetfilter.h"

class MyThread : public QThread
{
    Q_OBJECT

    WidgetFilter *gui;
    
    QString path;
    QString dir;

    bool statusFin;
    int mode;

public:
    MyThread(WidgetFilter *wf) : gui(wf) { }

    void setFile(QString file)
    {
        path = file;
        statusFin=false;
    }

    void setDir(QString file)
    {
        dir = file;
        statusFin=false;
    }

    void setMode(int i)
    {
        mode = i;
    }

protected:
    void run()
    {
        switch (mode)
        {
        case 0:
            this->gui->loadData(path);
            break;
        case 1: //save
            this->gui->saveData(dir);

        default:
            break;
        }
    }
public:
    WidgetFilter * getWidgetFilter()
    {
        return gui;
    }
};


class WindowFilter : public QMainWindow
{
    Q_OBJECT

private:
    WidgetFilter* widget;

    void open()
    {
        QString fileName = QFileDialog::getOpenFileName(this,
                                                        tr("Open image"), QDir::currentPath());

        //QApplication::setOverrideCursor(Qt::BusyCursor);

        if (fileName.isEmpty())
        {
            QMessageBox::information(this, tr("Image Viewer"),
                                     tr("Cannot load %1.").arg(fileName));
            return;
        }
        else
        {
            if (backThread->isRunning())
                QMessageBox::information(this, tr("Image Viewer"),
                                         tr("Busy"));
            else
            {
                backThread->setFile(fileName);
                backThread->setMode(0);
                backThread->start();
            }
        }
    }

    void save()
    {
        QString dir = QFileDialog::getExistingDirectory(this, tr("Open Directory"),
                                                        "/home",
                                                        QFileDialog::ShowDirsOnly
                                                        | QFileDialog::DontResolveSymlinks);
        if (backThread->isRunning())
            QMessageBox::information(this, tr("Image Viewer"),
                                     tr("Busy"));
        else
        {
            backThread->setDir(dir);
            backThread->setMode(1);
            backThread->start();
        }
        QMessageBox::information(this, tr("Image Viewer"),tr("DONE"));
    }


    QAction *openAct;
    QAction *saveAct;
    QAction *exitAct;

    QMenu *fileMenu;

    MyThread *backThread;

    void createActions()
    {
        openAct = new QAction(tr("&Open image..."), this);
        openAct->setShortcut(tr("Ctrl+O"));
        QObject::connect(openAct, &QAction::triggered, this, &WindowFilter::open);

        saveAct = new QAction(tr("&Save"), this);
        saveAct->setShortcut(tr("Ctrl+S"));
        QObject::connect(saveAct, &QAction::triggered, this, &WindowFilter::save);

        exitAct = new QAction(tr("&Exit"), this);
        exitAct->setShortcut(tr("Ctrl+Q"));
        QObject::connect(exitAct, &QAction::triggered, this, &WindowFilter::close);
    }

    void createMenus()
    {
        fileMenu = new QMenu(tr("&File"), this);
        fileMenu->addAction(openAct);
        fileMenu->addAction(saveAct);
        fileMenu->addAction(exitAct);
        menuBar()->addMenu(fileMenu);
    }

public:
    /**
    int addBrush(const char* name, float value, int channel=0)
    {
        return widget->addBrush(name, value, channel);
    }
    **/

    void setFromCommandline(int argc, char** argv)
    {
        for (int i = 1; i<(argc-1); ++i)
        {
            if (strcmp("-load",argv[i])==0) {
                backThread->setFile(QString(argv[++i]));
                backThread->setMode(0);
                backThread->start();
            }
        }

    }


    WindowFilter(const Filter& f)
    {	 
         widget = new WidgetFilter(f);
         setCentralWidget(widget);
         createActions();
         createMenus();
         setWindowTitle(tr("Image editor"));
         backThread = new MyThread(widget);
    }
    
    MyThread * getMyTHread()
    {
        return backThread;
    }


};

