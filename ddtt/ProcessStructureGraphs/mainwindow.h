#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

	void outputResults(int G, int M, int C, QString typeShapes, QString resultsFolder);

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
