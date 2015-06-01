/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.4.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *centralWidget;
    QGridLayout *gridLayout;
    QGroupBox *groupBox;
    QGridLayout *gridLayout_2;
    QRadioButton *checkbox2;
    QRadioButton *checkbox1;
    QRadioButton *checkbox3;
    QDoubleSpinBox *greedyParam;
    QSpacerItem *horizontalSpacer;
    QPushButton *conSegButton;
    QSpacerItem *horizontalSpacer_2;
    QPushButton *exportLabelsButton;
    QPushButton *extractPly;
    QPushButton *vkButton;
    QSpinBox *greedyParam2;
    QPushButton *greedyButton;
    QPushButton *visualizeButton;
    QSpacerItem *verticalSpacer;
    QPushButton *pairWiseButton;
    QPushButton *fuzzyButton;
    QPushButton *setWiseButton;
    QPushButton *clusterButton;
    QPushButton *baselineButton;
    QMenuBar *menuBar;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;
    QButtonGroup *buttonGroup;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QStringLiteral("MainWindow"));
        MainWindow->resize(362, 477);
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        gridLayout = new QGridLayout(centralWidget);
        gridLayout->setSpacing(6);
        gridLayout->setContentsMargins(11, 11, 11, 11);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        groupBox = new QGroupBox(centralWidget);
        groupBox->setObjectName(QStringLiteral("groupBox"));
        gridLayout_2 = new QGridLayout(groupBox);
        gridLayout_2->setSpacing(6);
        gridLayout_2->setContentsMargins(11, 11, 11, 11);
        gridLayout_2->setObjectName(QStringLiteral("gridLayout_2"));
        checkbox2 = new QRadioButton(groupBox);
        buttonGroup = new QButtonGroup(MainWindow);
        buttonGroup->setObjectName(QStringLiteral("buttonGroup"));
        buttonGroup->addButton(checkbox2);
        checkbox2->setObjectName(QStringLiteral("checkbox2"));
        checkbox2->setChecked(true);

        gridLayout_2->addWidget(checkbox2, 0, 1, 1, 1);

        checkbox1 = new QRadioButton(groupBox);
        buttonGroup->addButton(checkbox1);
        checkbox1->setObjectName(QStringLiteral("checkbox1"));
        checkbox1->setChecked(false);

        gridLayout_2->addWidget(checkbox1, 0, 0, 1, 1);

        checkbox3 = new QRadioButton(groupBox);
        checkbox3->setObjectName(QStringLiteral("checkbox3"));

        gridLayout_2->addWidget(checkbox3, 1, 0, 1, 1);


        gridLayout->addWidget(groupBox, 3, 0, 1, 3);

        greedyParam = new QDoubleSpinBox(centralWidget);
        greedyParam->setObjectName(QStringLiteral("greedyParam"));
        greedyParam->setDecimals(3);
        greedyParam->setMinimum(0);
        greedyParam->setMaximum(10000);
        greedyParam->setSingleStep(0.05);
        greedyParam->setValue(0.08);

        gridLayout->addWidget(greedyParam, 4, 2, 1, 1);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout->addItem(horizontalSpacer, 0, 0, 1, 1);

        conSegButton = new QPushButton(centralWidget);
        conSegButton->setObjectName(QStringLiteral("conSegButton"));

        gridLayout->addWidget(conSegButton, 2, 1, 1, 1);

        horizontalSpacer_2 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout->addItem(horizontalSpacer_2, 0, 2, 1, 1);

        exportLabelsButton = new QPushButton(centralWidget);
        exportLabelsButton->setObjectName(QStringLiteral("exportLabelsButton"));

        gridLayout->addWidget(exportLabelsButton, 1, 1, 1, 1);

        extractPly = new QPushButton(centralWidget);
        extractPly->setObjectName(QStringLiteral("extractPly"));

        gridLayout->addWidget(extractPly, 0, 1, 1, 1);

        vkButton = new QPushButton(centralWidget);
        vkButton->setObjectName(QStringLiteral("vkButton"));

        gridLayout->addWidget(vkButton, 5, 1, 1, 1);

        greedyParam2 = new QSpinBox(centralWidget);
        greedyParam2->setObjectName(QStringLiteral("greedyParam2"));
        greedyParam2->setMaximum(1000);
        greedyParam2->setValue(10);

        gridLayout->addWidget(greedyParam2, 4, 0, 1, 1);

        greedyButton = new QPushButton(centralWidget);
        greedyButton->setObjectName(QStringLiteral("greedyButton"));

        gridLayout->addWidget(greedyButton, 4, 1, 1, 1);

        visualizeButton = new QPushButton(centralWidget);
        visualizeButton->setObjectName(QStringLiteral("visualizeButton"));

        gridLayout->addWidget(visualizeButton, 6, 1, 1, 1);

        verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout->addItem(verticalSpacer, 10, 1, 1, 1);

        pairWiseButton = new QPushButton(centralWidget);
        pairWiseButton->setObjectName(QStringLiteral("pairWiseButton"));

        gridLayout->addWidget(pairWiseButton, 11, 1, 1, 1);

        fuzzyButton = new QPushButton(centralWidget);
        fuzzyButton->setObjectName(QStringLiteral("fuzzyButton"));

        gridLayout->addWidget(fuzzyButton, 7, 1, 1, 1);

        setWiseButton = new QPushButton(centralWidget);
        setWiseButton->setObjectName(QStringLiteral("setWiseButton"));

        gridLayout->addWidget(setWiseButton, 12, 1, 1, 1);

        clusterButton = new QPushButton(centralWidget);
        clusterButton->setObjectName(QStringLiteral("clusterButton"));

        gridLayout->addWidget(clusterButton, 8, 1, 1, 1);

        baselineButton = new QPushButton(centralWidget);
        baselineButton->setObjectName(QStringLiteral("baselineButton"));

        gridLayout->addWidget(baselineButton, 9, 1, 1, 1);

        MainWindow->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(MainWindow);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 362, 21));
        MainWindow->setMenuBar(menuBar);
        mainToolBar = new QToolBar(MainWindow);
        mainToolBar->setObjectName(QStringLiteral("mainToolBar"));
        MainWindow->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(MainWindow);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        MainWindow->setStatusBar(statusBar);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "MainWindow", 0));
        groupBox->setTitle(QApplication::translate("MainWindow", "GroupBox", 0));
        checkbox2->setText(QApplication::translate("MainWindow", "K-means", 0));
        checkbox1->setText(QApplication::translate("MainWindow", "DBSCAN", 0));
        checkbox3->setText(QApplication::translate("MainWindow", "K-means2", 0));
        conSegButton->setText(QApplication::translate("MainWindow", "Consistent Segmentation..", 0));
        exportLabelsButton->setText(QApplication::translate("MainWindow", "Export Labels...", 0));
        extractPly->setText(QApplication::translate("MainWindow", "Extract Ply...", 0));
        vkButton->setText(QApplication::translate("MainWindow", "VK template...", 0));
        greedyButton->setText(QApplication::translate("MainWindow", "Greedy correspondence...", 0));
        visualizeButton->setText(QApplication::translate("MainWindow", "Visualize clustering..", 0));
        pairWiseButton->setText(QApplication::translate("MainWindow", "Evaluate pair-wise...", 0));
        fuzzyButton->setText(QApplication::translate("MainWindow", "Fuzzy correspond..", 0));
        setWiseButton->setText(QApplication::translate("MainWindow", "Evaluate set-wise...", 0));
        clusterButton->setText(QApplication::translate("MainWindow", "Extract cluster...", 0));
        baselineButton->setText(QApplication::translate("MainWindow", "Baseline", 0));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
