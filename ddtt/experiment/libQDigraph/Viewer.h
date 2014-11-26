#ifndef VIEWER_H
#define VIEWER_H

#include <QWidget>
#include <QWebFrame>
#include <QWebView>
#include <QMap>
#include <QVariant>

namespace Ui {
class Viewer;
}

class Viewer : public QWidget
{
    Q_OBJECT

public:
    explicit Viewer(QWidget *parent = 0);
    ~Viewer();
    Ui::Viewer *ui;
    QWebView * wv;
    QMap<int, QMap<QString,QVariant> >nodeProperties;

public slots:
    void insertLogItem(QString msg);

    int addNode(QString node_info);
    void addEdge(QString edge_info);
    void addEdge(int nid1, int nid2);

	void updateGraph();
	void addCSS(QString style_code);

    void doTestGraph();
    void doExpand();

signals:
    void addLogItem(QString msg);
    void nodeSelected(int nid);
    void goingToExpand();
};

#endif // VIEWER_H
