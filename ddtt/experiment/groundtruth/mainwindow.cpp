#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QDir>
#include <QDebug>
#include <QFile>

struct Record{
    QString shapeA, shapeB;
    QString p,q;
    QString user;
    Record(){}
    Record(QString shapeA, QString shapeB, QString p, QString q, QString user):shapeA(shapeA), shapeB(shapeB),p(p),q(q),user(user){}

    bool operator < (const Record& rhs) const {
        return (p+"-"+q+"-"+shapeA+"-"+shapeB) < (rhs.p+"-"+rhs.q+"-"+rhs.shapeA+"-"+rhs.shapeB);
    }
};

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    QVector<Record> records;

    QString datasetPath = "C:/Users/Ibraheem/Dropbox/DDTT-dataset/evaluation/userstudy";

    QDir datasetDir(datasetPath);
    QStringList subdirs = datasetDir.entryList(QDir::Dirs | QDir::NoSymLinks | QDir::NoDotAndDotDot);

    for(QString subdir : subdirs)
    {
        QDir d(datasetPath + "/" + subdir);

        auto matchFiles = d.entryList(QStringList() << "*.match", QDir::Files);
        if( matchFiles.isEmpty() ) continue;

        for(QString matchFile : matchFiles)
        {
            // Get matches from file
            QFile ff(d.absolutePath() + "/" + matchFile);
            ff.open(QFile::ReadOnly | QFile::Text);
            QTextStream in(&ff);
            auto matches = in.readAll().split("\n", QString::SkipEmptyParts);

            // Shape names
            QStringList shapes = matchFile.replace(".match","").split("-", QString::SkipEmptyParts);

            bool isSwap = shapes.front() > shapes.back();

            // Add to vector
            QVector< QPair<QString,QString> > pairMatches;
            for(QString matchLine : matches)
            {
                QStringList match = matchLine.split(" ", QString::SkipEmptyParts);

                QString p = match.front();
                QString q = match.back();

                if(isSwap) std::swap(p,q);

                pairMatches << qMakePair(p,q);
            }

            //qDebug() << pairMatches;

            // Add record
            if(isSwap) std::reverse(shapes.begin(), shapes.end());
            for(auto a : pairMatches) records << Record(shapes.front(), shapes.back(), a.first, a.second, subdir);
        }
    }

    // Debug records:
    qSort(records);
    for(Record & record : records)
    {
        qDebug() << record.p << " " << record.q << " " << record.shapeA << " " << record.shapeB << " " << record.user ;
    }

     qDebug() << "Num records =" << records.size();
}

MainWindow::~MainWindow()
{
    delete ui;
}
