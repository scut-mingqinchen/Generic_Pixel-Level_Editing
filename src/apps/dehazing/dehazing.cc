#include <QApplication>

#include <gui/windowfilter.h>
#include <filters/filterdehazing.h>

#define NO_FOG 255.0
#define MIN_FOG 160.0
#define MEDIUM_FOG 96.0
#define MAX_FOG 30.0

 int main(int argc, char *argv[])
 {
     QApplication app(argc, argv);
     
     FilterDehazing filter;
     WindowFilter window(filter);
     window.setFromCommandline(argc, argv);
     window.show();

     return app.exec();
 }
