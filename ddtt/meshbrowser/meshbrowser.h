#ifndef MeshBrowser_H
#define MeshBrowser_H

#include <QWidget>
#include <QAbstractButton>

namespace Ui {
class MeshBrowser;
}

class MeshBrowser : public QWidget
{
    Q_OBJECT

public:
    explicit MeshBrowser(QWidget *parent = 0);
    ~MeshBrowser();
	Ui::MeshBrowser *ui;

	QStringList database;
	int nU, nV;

public slots:
	void loadMeshes();
};

#endif // MeshBrowser_H
