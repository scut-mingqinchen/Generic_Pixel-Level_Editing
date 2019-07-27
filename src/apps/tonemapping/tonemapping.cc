#include <QApplication>

#include <gui/windowfilter.h>
#include <filters/filtertonemappingmantiuk.h>

 int main(int argc, char *argv[])
 {
     QApplication app(argc, argv);
     
     FilterTonemappingMantiuk filter;
     WindowFilter window(filter);
     window.setFromCommandline(argc,argv);
     window.show();

     return app.exec();
 }
