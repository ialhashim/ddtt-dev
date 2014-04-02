#include "imagebrowser.h"
#include "ui_imagebrowser.h"
#include <QWidget>
#include <QPushButton>
#include <QFileDialog>
#include <QSettings>

#include "myimagearea.h"

MyImageArea * lastSelected = NULL;
QVector<MyImageArea*> activeViewers;

ImageBrowser::ImageBrowser(QWidget *parent) : QWidget(parent), ui(new Ui::ImageBrowser)
{
	// Anti-aliasing when using QGLWidget or subclasses
    //QGLFormat glf = QGLFormat::defaultFormat();
    //glf.setSamples(8);
    //QGLFormat::setDefaultFormat(glf);

    ui->setupUi(this);

    {
        QString datasetPath = QFileDialog::getExistingDirectory(0, "Dataset");

        // Get list of files
        QStringList filters;
        filters << "*.png" << "*.jpg" << "*.gif" << "*.jpeg" << "*.bmp";
        QStringList files = QDir(datasetPath).entryList(filters);

        foreach(QString mesh, files)
        {
            database << (datasetPath + "/" + mesh);
		}
	}

	// Settings
	{
		// Use local file
		QString path = QDir(qApp->applicationDirPath()).absoluteFilePath("meshbrowser.ini");
		QSettings qsettings(path, QSettings::IniFormat);

		// Defaults  
		if(!qsettings.allKeys().contains("nV")){
			qsettings.setValue("nV", 3);
			qsettings.setValue("nU", 2);
			qsettings.sync();
		}
 
		nV = qsettings.value("nV", 3).toInt();
		nU = qsettings.value("nU", 2).toInt();
	}

	// Create viewers
	{
		// Scroll bars
		int numPages = ceil(double(database.size()) / double(nU * nV)) - 1;
		ui->viewersScroll->setValue(0);
		ui->viewersScroll->setRange(0, numPages);
        connect(ui->viewersScroll, SIGNAL(valueChanged(int)), SLOT(loadImages()));

		// Load meshes
        loadImages();
	}

    // Connections
	//connect(ui->noneButton, &QPushButton::released, [=]() {curOp = MeshOperation::NONE_OP; ui->curLabel->setText(meshOps[curOp]);});

	connect(ui->addDelete, &QPushButton::released, [=](){ 
		if(lastSelected && deletedItems().contains(lastSelected->filename)){
			for(auto i : ui->deleteList->findItems(lastSelected->filename, Qt::MatchExactly))
				delete ui->deleteList->takeItem( ui->deleteList->row(i) );
			refreshViewers();
			return;
		}
		if(lastSelected) ui->deleteList->insertItem(0, lastSelected->filename); 
	});

	connect(ui->deleteButton, &QPushButton::released, [=](){
		for(auto filename : deletedItems()){
			QFile::remove(filename);
			for(int i = 0; i < ui->deleteList->count(); ++i){
				for(auto viewer: activeViewers){
					if(viewer->filename == filename){
						viewer->isDeleted = true;
						break;
					}
				}
			}
		}
		ui->deleteList->clear();
	});
}

QStringList ImageBrowser::deletedItems(){
	QStringList items;
	for(int i = 0; i < ui->deleteList->count(); ++i){
		QListWidgetItem* item = ui->deleteList->item(i);
		items << item->text();
	}
	return items;
}

void ImageBrowser::loadImages()
{
	lastSelected = NULL;

	QGridLayout * layout = (QGridLayout *)ui->viewersFrame->layout();
	QLayoutItem* item;
	while ( ( item = layout->takeAt( 0 ) ) != NULL ){
		delete item->widget(); delete item;
	}
	layout->setMargin(0);
	layout->setSpacing(0);

	int offset = ui->viewersScroll->value() * (nU * nV);
	int numPages = (database.size()-1) / (nU * nV);
	int pageID = (double(offset) / (database.size()-1)) * numPages;
	ui->curPage->setText(QString("Page: %1 / %2").arg(pageID).arg(numPages));
	
	activeViewers.clear();

	for(int i = 0; i < nU; i++){
		for(int j = 0; j < nV; j++){

			int idx = offset + (i*nV) + j;
			if(idx + 1 > database.size()) continue;
			QString filename = database[idx];

            //connect(viewer, &MyDrawArea::gotFocus, [=](MyDrawArea * caller) {
            //	lastSelected = caller;
            //	refreshViewers();
            //});
		}
	}

	refreshViewers();

	this->activateWindow();
	this->setFocus();
	this->raise();
}

void ImageBrowser::refreshViewers(){
	for(auto viewer : activeViewers)
	{
		viewer->parentWidget()->setStyleSheet("border: 2px solid #232323;");
		if(deletedItems().contains(viewer->filename))
			viewer->parentWidget()->setStyleSheet("border: 2px solid red;");
	}

	if(lastSelected)
	{
		lastSelected->parentWidget()->setStyleSheet("border: 2px solid rgb(255,100,100);");
		ui->selectedMesh->setText( QFileInfo(lastSelected->filename).baseName() );
	}
}

ImageBrowser::~ImageBrowser()
{
    delete ui;
}
