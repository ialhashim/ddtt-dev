#pragma once

// Simplify runtime debugging
#include <QMessageBox>
template<typename DataType>
static inline void debugBox( DataType message ){
    QMessageBox msgBox;
    msgBox.setTextFormat(Qt::RichText);
    msgBox.setText( QString("%1").arg(message) );
    msgBox.setStandardButtons(QMessageBox::Ok);
    msgBox.exec();
}
static inline void debugBoxList( QStringList messages ){
    debugBox( messages.join("<br>") );
}
template<typename Container>
static inline void debugBoxVec( Container data ){
    QStringList l;
    for(auto d : data) l << QString("%1").arg(d);
    debugBoxList(l);
}
template<typename Container2D>
static inline void debugBoxVec2( Container2D data, int limit = -1 ){
    QStringList l;
    for(auto row : data){
        QStringList line;
        for(auto d : row) line << QString("%1").arg( d );
        l << QString("%1").arg( line.join(", ") );
        if(limit > 0 && l.size() == limit) break;
    }
    if(limit > 0 && data.size() - l.size() > 0) l << QString("... (%1) more").arg(data.size() - l.size());
    debugBoxList(l);
}

// Loading XML files
#include <QDir>
typedef QMap<QString, QVariantMap> DatasetMap;
static inline DatasetMap shapesInDataset(QString datasetPath)
{
    DatasetMap dataset;

    QDir datasetDir(datasetPath);
    QStringList subdirs = datasetDir.entryList(QDir::Dirs | QDir::NoSymLinks | QDir::NoDotAndDotDot);

    for(auto subdir : subdirs)
    {
        // Special folders
        if (subdir == "corr") continue;

        QDir d(datasetPath + "/" + subdir);

        // Check if no graph is in this folder
        if (d.entryList(QStringList() << "*.xml", QDir::Files).isEmpty()) continue;

        auto xml_files = d.entryList(QStringList() << "*.xml", QDir::Files);

        dataset[subdir]["Name"] = subdir;
        dataset[subdir]["graphFile"] = d.absolutePath() + "/" + (xml_files.empty() ? "" : xml_files.front());
        dataset[subdir]["thumbFile"] = d.absolutePath() + "/" + d.entryList(QStringList() << "*.png", QDir::Files).join("");
        dataset[subdir]["objFile"] = d.absolutePath() + "/" + d.entryList(QStringList() << "*.obj", QDir::Files).join("");
    }

    return dataset;
}
