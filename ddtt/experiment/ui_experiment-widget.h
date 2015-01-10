/********************************************************************************
** Form generated from reading UI file 'experiment-widget.ui'
**
** Created by: Qt User Interface Compiler version 5.4.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_EXPERIMENT_2D_WIDGET_H
#define UI_EXPERIMENT_2D_WIDGET_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QListWidget>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_ExperimentWidget
{
public:
    QGridLayout *gridLayout;
    QLabel *label;
    QSpacerItem *verticalSpacer_2;
    QCheckBox *isShowParts;
    QPushButton *doEnergyGuided;
    QPushButton *executeButton;
    QPushButton *clearLandmarks;
    QCheckBox *isUseDP;
    QPushButton *loadShapes;
    QCheckBox *isVisualize;
    QSpinBox *numIterations;
    QDoubleSpinBox *speed;
    QCheckBox *isShowLandmarks;
    QPushButton *saveLandmarks;
    QPushButton *searchBest;
    QPushButton *resetShapes;
    QPushButton *loadLandmarks;
    QListWidget *pathsList;
    QHBoxLayout *horizontalLayout;
    QPushButton *saveJob;
    QPushButton *loadJobs;
    QCheckBox *isOtherEnergy;
    QCheckBox *isLimitedSearch;
    QCheckBox *isAnisotropy;
    QPushButton *rotateButton;
    QPushButton *swapButton;
    QSpacerItem *verticalSpacer;
    QPushButton *clearShapes;
    QCheckBox *isUseSYMH;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label_2;
    QLineEdit *dpTopK;

    void setupUi(QWidget *ExperimentWidget)
    {
        if (ExperimentWidget->objectName().isEmpty())
            ExperimentWidget->setObjectName(QStringLiteral("ExperimentWidget"));
        ExperimentWidget->resize(298, 586);
        gridLayout = new QGridLayout(ExperimentWidget);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        label = new QLabel(ExperimentWidget);
        label->setObjectName(QStringLiteral("label"));

        gridLayout->addWidget(label, 15, 0, 1, 1);

        verticalSpacer_2 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout->addItem(verticalSpacer_2, 11, 0, 1, 3);

        isShowParts = new QCheckBox(ExperimentWidget);
        isShowParts->setObjectName(QStringLiteral("isShowParts"));
        isShowParts->setChecked(true);

        gridLayout->addWidget(isShowParts, 14, 0, 1, 1);

        doEnergyGuided = new QPushButton(ExperimentWidget);
        doEnergyGuided->setObjectName(QStringLiteral("doEnergyGuided"));

        gridLayout->addWidget(doEnergyGuided, 27, 2, 1, 1);

        executeButton = new QPushButton(ExperimentWidget);
        executeButton->setObjectName(QStringLiteral("executeButton"));
        QFont font;
        font.setPointSize(10);
        font.setBold(true);
        font.setWeight(75);
        executeButton->setFont(font);

        gridLayout->addWidget(executeButton, 16, 0, 1, 3);

        clearLandmarks = new QPushButton(ExperimentWidget);
        clearLandmarks->setObjectName(QStringLiteral("clearLandmarks"));

        gridLayout->addWidget(clearLandmarks, 7, 0, 1, 3);

        isUseDP = new QCheckBox(ExperimentWidget);
        isUseDP->setObjectName(QStringLiteral("isUseDP"));
        isUseDP->setChecked(true);

        gridLayout->addWidget(isUseDP, 21, 0, 1, 1);

        loadShapes = new QPushButton(ExperimentWidget);
        loadShapes->setObjectName(QStringLiteral("loadShapes"));

        gridLayout->addWidget(loadShapes, 0, 2, 1, 1);

        isVisualize = new QCheckBox(ExperimentWidget);
        isVisualize->setObjectName(QStringLiteral("isVisualize"));

        gridLayout->addWidget(isVisualize, 13, 0, 1, 3);

        numIterations = new QSpinBox(ExperimentWidget);
        numIterations->setObjectName(QStringLiteral("numIterations"));
        numIterations->setMaximum(1000);
        numIterations->setSingleStep(1);
        numIterations->setValue(0);

        gridLayout->addWidget(numIterations, 15, 2, 1, 1);

        speed = new QDoubleSpinBox(ExperimentWidget);
        speed->setObjectName(QStringLiteral("speed"));
        speed->setDecimals(3);
        speed->setMinimum(0.01);
        speed->setSingleStep(0.1);
        speed->setValue(0.05);

        gridLayout->addWidget(speed, 27, 0, 1, 1);

        isShowLandmarks = new QCheckBox(ExperimentWidget);
        isShowLandmarks->setObjectName(QStringLiteral("isShowLandmarks"));
        isShowLandmarks->setChecked(true);

        gridLayout->addWidget(isShowLandmarks, 9, 0, 1, 3);

        saveLandmarks = new QPushButton(ExperimentWidget);
        saveLandmarks->setObjectName(QStringLiteral("saveLandmarks"));

        gridLayout->addWidget(saveLandmarks, 5, 0, 1, 3);

        searchBest = new QPushButton(ExperimentWidget);
        searchBest->setObjectName(QStringLiteral("searchBest"));

        gridLayout->addWidget(searchBest, 30, 2, 1, 1);

        resetShapes = new QPushButton(ExperimentWidget);
        resetShapes->setObjectName(QStringLiteral("resetShapes"));

        gridLayout->addWidget(resetShapes, 1, 2, 1, 1);

        loadLandmarks = new QPushButton(ExperimentWidget);
        loadLandmarks->setObjectName(QStringLiteral("loadLandmarks"));

        gridLayout->addWidget(loadLandmarks, 6, 0, 1, 3);

        pathsList = new QListWidget(ExperimentWidget);
        pathsList->setObjectName(QStringLiteral("pathsList"));
        pathsList->setSortingEnabled(true);

        gridLayout->addWidget(pathsList, 12, 0, 1, 3);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        saveJob = new QPushButton(ExperimentWidget);
        saveJob->setObjectName(QStringLiteral("saveJob"));

        horizontalLayout->addWidget(saveJob);

        loadJobs = new QPushButton(ExperimentWidget);
        loadJobs->setObjectName(QStringLiteral("loadJobs"));

        horizontalLayout->addWidget(loadJobs);


        gridLayout->addLayout(horizontalLayout, 29, 0, 2, 1);

        isOtherEnergy = new QCheckBox(ExperimentWidget);
        isOtherEnergy->setObjectName(QStringLiteral("isOtherEnergy"));

        gridLayout->addWidget(isOtherEnergy, 18, 0, 1, 3);

        isLimitedSearch = new QCheckBox(ExperimentWidget);
        isLimitedSearch->setObjectName(QStringLiteral("isLimitedSearch"));

        gridLayout->addWidget(isLimitedSearch, 19, 0, 1, 1);

        isAnisotropy = new QCheckBox(ExperimentWidget);
        isAnisotropy->setObjectName(QStringLiteral("isAnisotropy"));
        isAnisotropy->setChecked(false);

        gridLayout->addWidget(isAnisotropy, 17, 0, 1, 3);

        rotateButton = new QPushButton(ExperimentWidget);
        rotateButton->setObjectName(QStringLiteral("rotateButton"));

        gridLayout->addWidget(rotateButton, 14, 2, 1, 1);

        swapButton = new QPushButton(ExperimentWidget);
        swapButton->setObjectName(QStringLiteral("swapButton"));

        gridLayout->addWidget(swapButton, 1, 0, 1, 1);

        verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout->addItem(verticalSpacer, 3, 0, 1, 3);

        clearShapes = new QPushButton(ExperimentWidget);
        clearShapes->setObjectName(QStringLiteral("clearShapes"));

        gridLayout->addWidget(clearShapes, 0, 0, 1, 1);

        isUseSYMH = new QCheckBox(ExperimentWidget);
        isUseSYMH->setObjectName(QStringLiteral("isUseSYMH"));

        gridLayout->addWidget(isUseSYMH, 20, 0, 1, 1);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        horizontalLayout_2->setSizeConstraint(QLayout::SetMinimumSize);
        label_2 = new QLabel(ExperimentWidget);
        label_2->setObjectName(QStringLiteral("label_2"));

        horizontalLayout_2->addWidget(label_2);

        dpTopK = new QLineEdit(ExperimentWidget);
        dpTopK->setObjectName(QStringLiteral("dpTopK"));
        dpTopK->setInputMethodHints(Qt::ImhDigitsOnly);

        horizontalLayout_2->addWidget(dpTopK);


        gridLayout->addLayout(horizontalLayout_2, 25, 0, 2, 1);


        retranslateUi(ExperimentWidget);

        QMetaObject::connectSlotsByName(ExperimentWidget);
    } // setupUi

    void retranslateUi(QWidget *ExperimentWidget)
    {
        ExperimentWidget->setWindowTitle(QApplication::translate("ExperimentWidget", "Form", 0));
        label->setText(QApplication::translate("ExperimentWidget", "Num. iterations", 0));
        isShowParts->setText(QApplication::translate("ExperimentWidget", "Show parts", 0));
        doEnergyGuided->setText(QApplication::translate("ExperimentWidget", "Energy", 0));
        executeButton->setText(QApplication::translate("ExperimentWidget", "Execute...", 0));
        clearLandmarks->setText(QApplication::translate("ExperimentWidget", "Clear landmarks", 0));
        isUseDP->setText(QApplication::translate("ExperimentWidget", "Use Dynamic Programming", 0));
        loadShapes->setText(QApplication::translate("ExperimentWidget", "Load shape..", 0));
        isVisualize->setText(QApplication::translate("ExperimentWidget", "Visualize", 0));
        isShowLandmarks->setText(QApplication::translate("ExperimentWidget", "Show landmarks", 0));
        saveLandmarks->setText(QApplication::translate("ExperimentWidget", "Save landmarks..", 0));
        searchBest->setText(QApplication::translate("ExperimentWidget", "Search..", 0));
        resetShapes->setText(QApplication::translate("ExperimentWidget", "Reset", 0));
        loadLandmarks->setText(QApplication::translate("ExperimentWidget", "Load landmarks..", 0));
        saveJob->setText(QApplication::translate("ExperimentWidget", "Save..", 0));
        loadJobs->setText(QApplication::translate("ExperimentWidget", "Load Jobs", 0));
        isOtherEnergy->setText(QApplication::translate("ExperimentWidget", "Use other energy", 0));
        isLimitedSearch->setText(QApplication::translate("ExperimentWidget", "Limited search", 0));
        isAnisotropy->setText(QApplication::translate("ExperimentWidget", "Anisotropy", 0));
        rotateButton->setText(QApplication::translate("ExperimentWidget", "Rotate..", 0));
        swapButton->setText(QApplication::translate("ExperimentWidget", "Swap shapes", 0));
        clearShapes->setText(QApplication::translate("ExperimentWidget", "Clear shapes", 0));
        isUseSYMH->setText(QApplication::translate("ExperimentWidget", "Use SYMH( random if not SYMH)", 0));
        label_2->setText(QApplication::translate("ExperimentWidget", "Top K=", 0));
        dpTopK->setInputMask(QString());
        dpTopK->setText(QApplication::translate("ExperimentWidget", "100", 0));
    } // retranslateUi

};

namespace Ui {
    class ExperimentWidget: public Ui_ExperimentWidget {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_EXPERIMENT_2D_WIDGET_H
