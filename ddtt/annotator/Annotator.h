#pragma once

#include <QWidget>
#include <QAbstractButton>

namespace Ui {
class Annotator;
}

class Annotator : public QWidget
{
    Q_OBJECT

public:
    explicit Annotator(QWidget *parent = 0);
    ~Annotator();
	Ui::Annotator *ui;

	QStringList database;
	int nU, nV;

	QStringList deletedItems();

public slots:
	void loadMeshes();
	void refreshViewers();
};

