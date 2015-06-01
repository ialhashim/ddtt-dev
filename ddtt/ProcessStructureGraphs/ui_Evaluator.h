/********************************************************************************
** Form generated from reading UI file 'Evaluator.ui'
**
** Created by: Qt User Interface Compiler version 5.4.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_EVALUATOR_H
#define UI_EVALUATOR_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_Evaluator
{
public:

    void setupUi(QWidget *Evaluator)
    {
        if (Evaluator->objectName().isEmpty())
            Evaluator->setObjectName(QStringLiteral("Evaluator"));
        Evaluator->resize(625, 476);

        retranslateUi(Evaluator);

        QMetaObject::connectSlotsByName(Evaluator);
    } // setupUi

    void retranslateUi(QWidget *Evaluator)
    {
        Evaluator->setWindowTitle(QApplication::translate("Evaluator", "Evaluator", 0));
    } // retranslateUi

};

namespace Ui {
    class Evaluator: public Ui_Evaluator {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_EVALUATOR_H
