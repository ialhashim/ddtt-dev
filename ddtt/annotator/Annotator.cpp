#include "Annotator.h"
#include "ui_Annotator.h"
#include "StarlabDrawArea.h"
#include <QPushButton>
#include <QFileDialog>
#include <QJsonDocument>
#include <QJsonArray>

#include "mydrawarea.h"
ShapeOperation curOp = NONE_OP;
QString curLabel;
int curLabelIdx;

using namespace SurfaceMesh;
MyDrawArea * lastSelected = NULL;

QVector<MyDrawArea*> activeViewers;

QMap<QString,QStringList> lables;
QMap<QString, int> labelIndices;
QMap<QString, int> parentLabelIndices;

// Loading XML files
#include <QDir>
typedef QMap<QString, QVariantMap> DatasetMap;
DatasetMap shapesInDataset(QString datasetPath)
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

DatasetMap dataset;
QString default_folder;

Annotator::Annotator(QWidget *parent) : QWidget(parent), ui(new Ui::Annotator)
{
	curLabel = "null";
	curLabelIdx = -1;

	// Anti-aliasing when using QGLWidget or subclasses
	QGLFormat glf = QGLFormat::defaultFormat();
	glf.setSamples(8);
	QGLFormat::setDefaultFormat(glf);

	ui->setupUi(this);

	// Default colors
	{
		auto paired_colors = { "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", 
			"#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928",
			"#4D4D4D", "#5DA5DA", "#FAA43A", "#60BD68", "#F17CB0", "#B2912F", "#B276B2", "#DECF3F", "#F15854",
			"#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf" };

		int num_paired_colors = paired_colors.size();

		for (int i = 0; i < num_paired_colors; i++){
			QColor c;
			c.setNamedColor(*(paired_colors.begin() + i));
			LabelColors << c;
		}

		LabelColors << QColor(255, 97, 121) << QColor(255, 219, 88) << QColor(107, 255, 135) << QColor(255, 165, 107) << QColor(104, 126, 255) <<
						QColor(242, 5, 135) << QColor(138, 0, 242) << QColor(3, 166, 60) << QColor(242, 203, 5);
		
		ParentColors = LabelColors;

		// Fill rest with random colors
		for (int i = 0; i < 50; i++) LabelColors << starlab::qRandomColor3();
	}

    {
        // Default folder
        {
            const QString DEFAULT_DIR_KEY("default_dir");
            QString path = QDir(qApp->applicationDirPath()).absoluteFilePath("Annotator.ini");
            QSettings MySettings(path, QSettings::IniFormat);
            QString oldPath = MySettings.value(DEFAULT_DIR_KEY).toString();
            QString dirpath = QFileDialog::getExistingDirectory(this, tr("Open Directory"),
                oldPath, QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
            if (dirpath.isNull() || dirpath.isEmpty()) return;
            MySettings.setValue(DEFAULT_DIR_KEY, dirpath);
            MySettings.sync();
            default_folder = dirpath;
        }

        dataset = shapesInDataset(default_folder);

        for(auto name : dataset.keys())
        {
            auto shape = dataset[name];
            database << shape["graphFile"].toString();
		}

		// Load labels 
		{
			// Open JSON file
			auto labelsFilename = default_folder + "/labels.json";
			QFile file;
			file.setFileName(labelsFilename);
			if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) return;
			QJsonDocument jdoc = QJsonDocument::fromJson(file.readAll());
			auto json = jdoc.object();

			//auto corrArray = jdoc.toVariant().value<QVariantList>();

			auto labelsArray = json["labels"].toArray();

			for (auto l : labelsArray)
			{
				auto label = l.toObject();
				auto parent = label["parent"].toString();

				lables[parent] << label["title"].toString();
			}

			// Fill tree widget
			{
				int uc = 0;
                int pc = 0;

				for (auto parentName : lables.keys())
				{
					auto parent = new QTreeWidgetItem();
					parent->setText(0, parentName);
                    parent->setData(0, Qt::UserRole, pc);
                    parent->setBackgroundColor(0, ParentColors[pc]);

                    parentLabelIndices[parentName] = pc++;

					for (auto child : lables[parentName])
					{
						labelNames << child;

						labelIndices[child] = uc;

						auto c = new QTreeWidgetItem();
						c->setText(0, child);
						c->setData(0, Qt::UserRole, uc);
						c->setBackgroundColor(0, LabelColors[uc++]);
						parent->addChild(c);
					}
					
					ui->labelsTreeView->addTopLevelItem(parent);
				}
			}
		}
	}

	// Settings
	{
		// Use local file
		QString path = QDir(qApp->applicationDirPath()).absoluteFilePath("Annotator.ini");
		QSettings qsettings(path, QSettings::IniFormat);

		// Defaults  
		if(!qsettings.allKeys().contains("nV")){
			qsettings.setValue("nV", 4);
			qsettings.setValue("nU", 3);
			qsettings.sync();
		}
 
		nV = qsettings.value("nV", 4).toInt();
		nU = qsettings.value("nU", 3).toInt();
	}

	// Create viewers
	{
		// Scroll bars
		int numPages = ceil(double(database.size()) / double(nU * nV)) - 1;
		ui->viewersScroll->setValue(0);
		ui->viewersScroll->setRange(0, numPages);
		connect(ui->viewersScroll, SIGNAL(valueChanged(int)), SLOT(loadMeshes()));

		// Load meshes
		loadMeshes();
	}

    // Connections
    connect(ui->noneButton, &QPushButton::released, [&]() {curOp = ShapeOperation::NONE_OP; ui->curLabel->setText(shapeOps[curOp]);});
	connect(ui->rotateLeft, &QPushButton::released, [&]() {curOp = ShapeOperation::ROTATE_LEFT; ui->curLabel->setText(shapeOps[curOp]); });
	connect(ui->rotateRight, &QPushButton::released, [&]() {curOp = ShapeOperation::ROTATE_RIGHT; ui->curLabel->setText(shapeOps[curOp]); });

	connect(ui->beginAnnotateButton, &QPushButton::released, [&](){curOp = ShapeOperation::LABEL_PART; ui->curLabel->setText(shapeOps[curOp]);});
	connect(ui->labelsTreeView, &QTreeWidget::currentItemChanged, [&](){
		if (ui->labelsTreeView->selectedItems().empty()) return;
		curLabel = ui->labelsTreeView->selectedItems().front()->text(0);
		curLabelIdx = ui->labelsTreeView->selectedItems().front()->data(0, Qt::UserRole).toInt();
		ui->curLabel->setText(shapeOps[curOp] + QString(" (%1)").arg(curLabel));
	});

    connect(ui->coarseLabels, &QCheckBox::toggled, [=](){

        for(auto viewer: activeViewers){
            auto shape = viewer->m;
            for(auto n : shape->nodes)
            {
                if(!n->property.contains("label")) continue;

                auto label = n->property["label"].toString();
                auto parentLabel = label.split("-").front();

                if(ui->coarseLabels->isChecked())
                {
                    shape->setColorFor(n->id, ParentColors[parentLabelIndices[parentLabel]]);
                }
                else
                {
                    shape->setColorFor(n->id, LabelColors[labelIndices[label]]);
                }
            }
        }

        refreshViewers();
    });

    connect(ui->remainButton, &QPushButton::released, [=](){
        if(!activeViewers.empty()) return;
        for(auto viewer: activeViewers){
            QStringList untagged;
            for(auto n : viewer->m->nodes)
            {
                if(!n->property.contains("label"))
                    untagged << n->id;
            }
            QMessageBox::information(0,"Untagged", untagged.join(","));
        }
    });

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

	connect(ui->generateThumbnails, &QPushButton::released, [=](){
		foreach(QString filename, database)
        {
            auto g = QSharedPointer<Structure::Graph>(new Structure::Graph( filename ));

            MyDrawArea * viewer = new MyDrawArea(g, filename);

			if( !QFileInfo(filename).exists() ){
				viewer->isDeleted = true;
				continue;
			}

			viewer->setMinimumSize(256,256);
			viewer->setMaximumSize(256,256);
			
			viewer->bg = QColor(255,255,255);
			viewer->fg = QColor(160,160,200);
			viewer->isDrawWireframe = false;
			viewer->isDoubleLight = true;

			viewer->show();

			{
                Eigen::AlignedBox3d bbox = g->robustBBox();
                viewer->camera()->setSceneRadius(bbox.sizes().norm() * 2);
				viewer->camera()->setUpVector(qglviewer::Vec(0,0,1));
				viewer->camera()->setPosition(qglviewer::Vec(2,-2,1.5));
				viewer->camera()->lookAt(qglviewer::Vec());
				viewer->camera()->showEntireScene();
				double S = bbox.sizes().norm() * 0.2;
				bbox.extend(bbox.min() + (Vector3(-1,-1,-1) * S));
				bbox.extend(bbox.max() + (Vector3(1,1,1) * S));
				viewer->camera()->setRevolveAroundPoint(qglviewer::Vec(bbox.center()));
				viewer->camera()->setSceneCenter(qglviewer::Vec(bbox.center()));
				viewer->camera()->fitBoundingBox(qglviewer::Vec(bbox.min().data()),qglviewer::Vec(bbox.max().data()));
			}

			QFileInfo meshFileInfo( filename );
			QString meshname = meshFileInfo.baseName();
			QString path = meshFileInfo.absolutePath();
			QString thumbImgFilename = QString("%1.png").arg(path + "/" + meshname);

			qApp->processEvents();

			viewer->raise();
			viewer->makeCurrent();
			viewer->grabFrameBuffer().save( thumbImgFilename );

			qApp->processEvents();

			viewer->hide();
			delete viewer;

			qApp->processEvents();
		}

	});

}

QStringList Annotator::deletedItems(){
	QStringList items;
	for(int i = 0; i < ui->deleteList->count(); ++i){
		QListWidgetItem* item = ui->deleteList->item(i);
		items << item->text();
	}
	return items;
}

void Annotator::loadMeshes()
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

            auto g = QSharedPointer<Structure::Graph>(new Structure::Graph( filename ));

            MyDrawArea * viewer = new MyDrawArea(g, filename);

			viewer->setForegroundColor(QColor(255,255,255));

			if( !QFileInfo(filename).exists() ){
				viewer->isDeleted = true;
				continue;
			}

			// Save original
            /*Vector3VertexProperty origPnts = m->vertex_coordinates();
			for(auto v: m->vertices()) viewer->originalVertices.push_back(origPnts[v]);
			for(auto f: m->faces()){
				std::vector<Vertex> face;
				for(auto v: m->vertices(f)) face.push_back(v);
				viewer->originalFaces.push_back(face);
            }*/

			QFrame * frame = new QFrame;
			QHBoxLayout * ll = new QHBoxLayout;
			frame->setLayout(ll);
			ll->addWidget(viewer);
			layout->setMargin(2);
			layout->setSpacing(2);
			layout->addWidget(frame, i, j, 1, 1);

			activeViewers.push_back(viewer);

			viewer->setGridIsDrawn(ui->showGrids->isChecked());
			viewer->isDrawWireframe = ui->drawWireframe->isChecked();

            Eigen::AlignedBox3d bbox = g->robustBBox();

            //viewer->setGridIsDrawn();
            viewer->camera()->setSceneRadius(bbox.sizes().norm() * 2);
			viewer->camera()->setUpVector(qglviewer::Vec(0,0,1));
			viewer->camera()->setPosition(qglviewer::Vec(2,-2,1.5));
			viewer->camera()->lookAt(qglviewer::Vec());
			viewer->camera()->showEntireScene();

            double S = bbox.sizes().norm() * 0.1;
            bbox.extend(bbox.min() + (Vector3(-1,-1,-1) * S));
            bbox.extend(bbox.max() + (Vector3(1,1,1) * S));

            viewer->camera()->setRevolveAroundPoint(qglviewer::Vec(bbox.center()));
            viewer->camera()->setSceneCenter(qglviewer::Vec(bbox.center()));
			viewer->camera()->fitBoundingBox(qglviewer::Vec(bbox.min().data()),qglviewer::Vec(bbox.max().data()));

			connect(viewer, &MyDrawArea::gotFocus, [=](MyDrawArea * caller) {
				lastSelected = caller; 
				refreshViewers();
			});
		}
	}

	refreshViewers();

	this->activateWindow();
	this->setFocus();
	this->raise();
}

void Annotator::refreshViewers(){
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

Annotator::~Annotator()
{
    delete ui;
}
