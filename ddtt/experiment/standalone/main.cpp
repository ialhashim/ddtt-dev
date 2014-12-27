#include "mainwindow.h"
#include <QApplication>
#include <QFileDialog>
#include "BatchProcess.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    //QString filename = QFileDialog::getOpenFileName(&w, "Load Jobs", "", "Jobs File (*.json)");
	QString filename = "C:/Users/Ibraheem/Desktop/StarlabPackage/jobs - Copy.json";

    QTimer::singleShot(0, [=] {
		BatchProcess * bp = new BatchProcess(filename);
		QObject::connect(bp, SIGNAL(finished()), bp, SLOT(deleteLater()));
        bp->start();
    });

    return a.exec();
}
