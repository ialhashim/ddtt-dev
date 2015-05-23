#include "Annotator.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    Annotator w;
    w.show();

    return a.exec();
}
