#include <QApplication>

#include <gui/windowfilter.h>
#include <filters/filterdepthoffield.h>

int main(int argc, char *argv[])
{
     QApplication app(argc, argv);
     
     FilterDepthOfField filter;
     WindowFilter window(filter);
     window.setFromCommandline(argc, argv);
     window.show();

     return app.exec();
}
