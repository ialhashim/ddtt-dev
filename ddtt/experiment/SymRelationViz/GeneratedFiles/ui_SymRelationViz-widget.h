/********************************************************************************
** Form generated from reading UI file 'SymRelationViz-widget.ui'
**
** Created by: Qt User Interface Compiler version 5.4.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_SYMRELATIONVIZ_2D_WIDGET_H
#define UI_SYMRELATIONVIZ_2D_WIDGET_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_SymRelationVizWidget
{
public:
    QVBoxLayout *verticalLayout;
    QHBoxLayout *horizontalLayout;
    QLabel *label;
    QDoubleSpinBox *scaleSpinBox;
    QPushButton *loadButton;
    QSpacerItem *verticalSpacer;

    void setupUi(QWidget *SymRelationVizWidget)
    {
        if (SymRelationVizWidget->objectName().isEmpty())
            SymRelationVizWidget->setObjectName(QStringLiteral("SymRelationVizWidget"));
        SymRelationVizWidget->resize(257, 564);
        verticalLayout = new QVBoxLayout(SymRelationVizWidget);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        label = new QLabel(SymRelationVizWidget);
        label->setObjectName(QStringLiteral("label"));

        horizontalLayout->addWidget(label);

        scaleSpinBox = new QDoubleSpinBox(SymRelationVizWidget);
        scaleSpinBox->setObjectName(QStringLiteral("scaleSpinBox"));
        scaleSpinBox->setMaximum(1);
        scaleSpinBox->setSingleStep(0.1);
        scaleSpinBox->setValue(0.6);

        horizontalLayout->addWidget(scaleSpinBox);


        verticalLayout->addLayout(horizontalLayout);

        loadButton = new QPushButton(SymRelationVizWidget);
        loadButton->setObjectName(QStringLiteral("loadButton"));

        verticalLayout->addWidget(loadButton);

        verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);


        retranslateUi(SymRelationVizWidget);

        QMetaObject::connectSlotsByName(SymRelationVizWidget);
    } // setupUi

    void retranslateUi(QWidget *SymRelationVizWidget)
    {
        SymRelationVizWidget->setWindowTitle(QApplication::translate("SymRelationVizWidget", "Form", 0));
        label->setText(QApplication::translate("SymRelationVizWidget", "Scale", 0));
        loadButton->setText(QApplication::translate("SymRelationVizWidget", "Load..", 0));
    } // retranslateUi

};

namespace Ui {
    class SymRelationVizWidget: public Ui_SymRelationVizWidget {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_SYMRELATIONVIZ_2D_WIDGET_H
