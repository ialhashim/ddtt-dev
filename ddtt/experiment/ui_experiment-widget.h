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
    QCheckBox *isShowLandmarks;
    QPushButton *loadShapes;
    QPushButton *clearShapes;
    QCheckBox *isAnisotropy;
    QPushButton *clearLandmarks;
    QCheckBox *isLimitedSearch;
    QPushButton *loadLandmarks;
    QPushButton *rotateButton;
    QCheckBox *isShowParts;
    QSpacerItem *verticalSpacer_2;
    QPushButton *doEnergyGuided;
    QSpinBox *numIterations;
    QDoubleSpinBox *speed;
    QCheckBox *isUseSYMH;
    QPushButton *saveLandmarks;
    QCheckBox *isOtherEnergy;
    QListWidget *pathsList;
    QSpacerItem *verticalSpacer;
    QPushButton *swapButton;
    QPushButton *resetShapes;
    QCheckBox *isVisualize;
    QPushButton *executeButton;
    QPushButton *searchBest;
    QLabel *label;
    QHBoxLayout *horizontalLayout;
    QPushButton *saveJob;
    QPushButton *loadJobs;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label_2;
    QLineEdit *dpTopK;
    QHBoxLayout *horizontalLayout_3;
    QLabel *label_3;
    QLineEdit *dpTopK_2;
    QCheckBox *isUseDP;

    void setupUi(QWidget *ExperimentWidget)
    {
        if (ExperimentWidget->objectName().isEmpty())
            ExperimentWidget->setObjectName(QStringLiteral("ExperimentWidget"));
        ExperimentWidget->resize(313, 590);
        gridLayout = new QGridLayout(ExperimentWidget);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        isShowLandmarks = new QCheckBox(ExperimentWidget);
        isShowLandmarks->setObjectName(QStringLiteral("isShowLandmarks"));
        isShowLandmarks->setChecked(true);

        gridLayout->addWidget(isShowLandmarks, 9, 0, 1, 3);

        loadShapes = new QPushButton(ExperimentWidget);
        loadShapes->setObjectName(QStringLiteral("loadShapes"));

        gridLayout->addWidget(loadShapes, 0, 2, 1, 1);

        clearShapes = new QPushButton(ExperimentWidget);
        clearShapes->setObjectName(QStringLiteral("clearShapes"));

        gridLayout->addWidget(clearShapes, 0, 0, 1, 1);

        isAnisotropy = new QCheckBox(ExperimentWidget);
        isAnisotropy->setObjectName(QStringLiteral("isAnisotropy"));
        isAnisotropy->setChecked(false);

        gridLayout->addWidget(isAnisotropy, 17, 0, 1, 3);

        clearLandmarks = new QPushButton(ExperimentWidget);
        clearLandmarks->setObjectName(QStringLiteral("clearLandmarks"));

        gridLayout->addWidget(clearLandmarks, 7, 0, 1, 3);

        isLimitedSearch = new QCheckBox(ExperimentWidget);
        isLimitedSearch->setObjectName(QStringLiteral("isLimitedSearch"));
        isLimitedSearch->setChecked(true);

        gridLayout->addWidget(isLimitedSearch, 19, 0, 1, 1);

        loadLandmarks = new QPushButton(ExperimentWidget);
        loadLandmarks->setObjectName(QStringLiteral("loadLandmarks"));

        gridLayout->addWidget(loadLandmarks, 6, 0, 1, 3);

        rotateButton = new QPushButton(ExperimentWidget);
        rotateButton->setObjectName(QStringLiteral("rotateButton"));

        gridLayout->addWidget(rotateButton, 14, 2, 1, 1);

        isShowParts = new QCheckBox(ExperimentWidget);
        isShowParts->setObjectName(QStringLiteral("isShowParts"));
        isShowParts->setChecked(true);

        gridLayout->addWidget(isShowParts, 14, 0, 1, 1);

        verticalSpacer_2 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout->addItem(verticalSpacer_2, 11, 0, 1, 3);

        doEnergyGuided = new QPushButton(ExperimentWidget);
        doEnergyGuided->setObjectName(QStringLiteral("doEnergyGuided"));

        gridLayout->addWidget(doEnergyGuided, 36, 2, 1, 1);

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

        gridLayout->addWidget(speed, 36, 0, 1, 1);

        isUseSYMH = new QCheckBox(ExperimentWidget);
        isUseSYMH->setObjectName(QStringLiteral("isUseSYMH"));

        gridLayout->addWidget(isUseSYMH, 20, 0, 1, 1);

        saveLandmarks = new QPushButton(ExperimentWidget);
        saveLandmarks->setObjectName(QStringLiteral("saveLandmarks"));

        gridLayout->addWidget(saveLandmarks, 5, 0, 1, 3);

        isOtherEnergy = new QCheckBox(ExperimentWidget);
        isOtherEnergy->setObjectName(QStringLiteral("isOtherEnergy"));

        gridLayout->addWidget(isOtherEnergy, 18, 0, 1, 3);

        pathsList = new QListWidget(ExperimentWidget);
        pathsList->setObjectName(QStringLiteral("pathsList"));
        pathsList->setSortingEnabled(true);

        gridLayout->addWidget(pathsList, 12, 0, 1, 3);

        verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout->addItem(verticalSpacer, 3, 0, 1, 3);

        swapButton = new QPushButton(ExperimentWidget);
        swapButton->setObjectName(QStringLiteral("swapButton"));

        gridLayout->addWidget(swapButton, 1, 0, 1, 1);

        resetShapes = new QPushButton(ExperimentWidget);
        resetShapes->setObjectName(QStringLiteral("resetShapes"));

        gridLayout->addWidget(resetShapes, 1, 2, 1, 1);

        isVisualize = new QCheckBox(ExperimentWidget);
        isVisualize->setObjectName(QStringLiteral("isVisualize"));

        gridLayout->addWidget(isVisualize, 13, 0, 1, 3);

        executeButton = new QPushButton(ExperimentWidget);
        executeButton->setObjectName(QStringLiteral("executeButton"));
        QFont font;
        font.setPointSize(10);
        font.setBold(true);
        font.setWeight(75);
        executeButton->setFont(font);

        gridLayout->addWidget(executeButton, 16, 0, 1, 3);

        searchBest = new QPushButton(ExperimentWidget);
        searchBest->setObjectName(QStringLiteral("searchBest"));

        gridLayout->addWidget(searchBest, 39, 2, 1, 1);

        label = new QLabel(ExperimentWidget);
        label->setObjectName(QStringLiteral("label"));

        gridLayout->addWidget(label, 15, 0, 1, 1);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        saveJob = new QPushButton(ExperimentWidget);
        saveJob->setObjectName(QStringLiteral("saveJob"));

        horizontalLayout->addWidget(saveJob);

        loadJobs = new QPushButton(ExperimentWidget);
        loadJobs->setObjectName(QStringLiteral("loadJobs"));

        horizontalLayout->addWidget(loadJobs);


        gridLayout->addLayout(horizontalLayout, 38, 0, 2, 1);

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


        gridLayout->addLayout(horizontalLayout_2, 22, 0, 1, 1);

        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        horizontalLayout_3->setSizeConstraint(QLayout::SetMinimumSize);
        label_3 = new QLabel(ExperimentWidget);
        label_3->setObjectName(QStringLiteral("label_3"));

        horizontalLayout_3->addWidget(label_3);

        dpTopK_2 = new QLineEdit(ExperimentWidget);
        dpTopK_2->setObjectName(QStringLiteral("dpTopK_2"));
        dpTopK_2->setInputMethodHints(Qt::ImhDigitsOnly);

        horizontalLayout_3->addWidget(dpTopK_2);


        gridLayout->addLayout(horizontalLayout_3, 26, 0, 1, 1);

        isUseDP = new QCheckBox(ExperimentWidget);
        isUseDP->setObjectName(QStringLiteral("isUseDP"));

        gridLayout->addWidget(isUseDP, 21, 0, 1, 1);


        retranslateUi(ExperimentWidget);

        QMetaObject::connectSlotsByName(ExperimentWidget);
    } // setupUi

    void retranslateUi(QWidget *ExperimentWidget)
    {
        ExperimentWidget->setWindowTitle(QApplication::translate("ExperimentWidget", "Form", 0));
        isShowLandmarks->setText(QApplication::translate("ExperimentWidget", "Show landmarks", 0));
        loadShapes->setText(QApplication::translate("ExperimentWidget", "Load shape..", 0));
        clearShapes->setText(QApplication::translate("ExperimentWidget", "Clear shapes", 0));
        isAnisotropy->setText(QApplication::translate("ExperimentWidget", "Anisotropy", 0));
        clearLandmarks->setText(QApplication::translate("ExperimentWidget", "Clear landmarks", 0));
        isLimitedSearch->setText(QApplication::translate("ExperimentWidget", "Limited search", 0));
        loadLandmarks->setText(QApplication::translate("ExperimentWidget", "Load landmarks..", 0));
        rotateButton->setText(QApplication::translate("ExperimentWidget", "Rotate..", 0));
        isShowParts->setText(QApplication::translate("ExperimentWidget", "Show parts", 0));
        doEnergyGuided->setText(QApplication::translate("ExperimentWidget", "Energy", 0));
        isUseSYMH->setText(QApplication::translate("ExperimentWidget", "Use SYMH( random if not SYMH)", 0));
        saveLandmarks->setText(QApplication::translate("ExperimentWidget", "Save landmarks..", 0));
        isOtherEnergy->setText(QApplication::translate("ExperimentWidget", "Use other energy", 0));
        swapButton->setText(QApplication::translate("ExperimentWidget", "Swap shapes", 0));
        resetShapes->setText(QApplication::translate("ExperimentWidget", "Reset", 0));
        isVisualize->setText(QApplication::translate("ExperimentWidget", "Visualize", 0));
        executeButton->setText(QApplication::translate("ExperimentWidget", "Execute...", 0));
        searchBest->setText(QApplication::translate("ExperimentWidget", "Search..", 0));
        label->setText(QApplication::translate("ExperimentWidget", "Num. iterations", 0));
        saveJob->setText(QApplication::translate("ExperimentWidget", "Save..", 0));
        loadJobs->setText(QApplication::translate("ExperimentWidget", "Load Jobs", 0));
        label_2->setText(QApplication::translate("ExperimentWidget", "Top K=", 0));
        dpTopK->setInputMask(QString());
        dpTopK->setText(QApplication::translate("ExperimentWidget", "20", 0));
        label_3->setText(QApplication::translate("ExperimentWidget", "Local K=", 0));
        dpTopK_2->setInputMask(QString());
        dpTopK_2->setText(QApplication::translate("ExperimentWidget", "1", 0));
        isUseDP->setText(QApplication::translate("ExperimentWidget", "Use Dynamic Programming", 0));
    } // retranslateUi

};

namespace Ui {
    class ExperimentWidget: public Ui_ExperimentWidget {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_EXPERIMENT_2D_WIDGET_H
