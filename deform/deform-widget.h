#pragma once

#include <QWidget>

namespace Ui { class DeformWidget; }

class DeformWidget : public QWidget
{
    Q_OBJECT

public:
    explicit DeformWidget(QWidget *parent = 0);
    ~DeformWidget();

//private:
    Ui::DeformWidget *ui;

signals:
	void shapesLoaded();
};
