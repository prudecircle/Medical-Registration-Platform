/********************************************************************************
** Form generated from reading UI file 'MainWindow.ui'
**
** Created by: Qt User Interface Compiler version 5.15.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QProgressBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindowClass
{
public:
    QWidget *centralWidget;
    QHBoxLayout *horizontalLayout;
    QWidget *widget;
    QWidget *paraWidget;
    QGridLayout *gridLayout;
    QLineEdit *tx;
    QLabel *label_tz;
    QLabel *label_tx;
    QLabel *label_ty;
    QLabel *label_rz;
    QLineEdit *ry;
    QLabel *label_rx;
    QLineEdit *ty;
    QLineEdit *rx;
    QLabel *label_ry;
    QLineEdit *rz;
    QLineEdit *tz;
    QLabel *label_similarity;
    QLabel *similarity;
    QWidget *mainWidget;
    QPushButton *btn_registration;
    QComboBox *box_algorithm;
    QPushButton *btn_openX;
    QPushButton *btn_openCT;
    QWidget *imgWidget;
    QLabel *drr_label;
    QLabel *xray_label;
    QLabel *drr_label_2;
    QLabel *xray_label_2;
    QProgressBar *progressBar;
    QLabel *label_similarity_2;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *MainWindowClass)
    {
        if (MainWindowClass->objectName().isEmpty())
            MainWindowClass->setObjectName(QString::fromUtf8("MainWindowClass"));
        MainWindowClass->setWindowModality(Qt::ApplicationModal);
        MainWindowClass->resize(1239, 669);
        QPalette palette;
        QBrush brush(QColor(170, 255, 255, 255));
        brush.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::Light, brush);
        QBrush brush1(QColor(255, 255, 255, 255));
        brush1.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::Base, brush1);
        palette.setBrush(QPalette::Active, QPalette::Window, brush);
        palette.setBrush(QPalette::Inactive, QPalette::Light, brush);
        palette.setBrush(QPalette::Inactive, QPalette::Base, brush1);
        palette.setBrush(QPalette::Inactive, QPalette::Window, brush);
        palette.setBrush(QPalette::Disabled, QPalette::Light, brush);
        palette.setBrush(QPalette::Disabled, QPalette::Base, brush);
        palette.setBrush(QPalette::Disabled, QPalette::Window, brush);
        MainWindowClass->setPalette(palette);
        centralWidget = new QWidget(MainWindowClass);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        centralWidget->setMinimumSize(QSize(1200, 500));
        horizontalLayout = new QHBoxLayout(centralWidget);
        horizontalLayout->setSpacing(6);
        horizontalLayout->setContentsMargins(11, 11, 11, 11);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        widget = new QWidget(centralWidget);
        widget->setObjectName(QString::fromUtf8("widget"));
        widget->setEnabled(true);
        paraWidget = new QWidget(widget);
        paraWidget->setObjectName(QString::fromUtf8("paraWidget"));
        paraWidget->setGeometry(QRect(830, 0, 371, 261));
        QPalette palette1;
        palette1.setBrush(QPalette::Active, QPalette::Light, brush1);
        palette1.setBrush(QPalette::Active, QPalette::Base, brush1);
        palette1.setBrush(QPalette::Active, QPalette::Window, brush1);
        palette1.setBrush(QPalette::Inactive, QPalette::Light, brush1);
        palette1.setBrush(QPalette::Inactive, QPalette::Base, brush1);
        palette1.setBrush(QPalette::Inactive, QPalette::Window, brush1);
        palette1.setBrush(QPalette::Disabled, QPalette::Light, brush1);
        palette1.setBrush(QPalette::Disabled, QPalette::Base, brush1);
        palette1.setBrush(QPalette::Disabled, QPalette::Window, brush1);
        paraWidget->setPalette(palette1);
        gridLayout = new QGridLayout(paraWidget);
        gridLayout->setSpacing(6);
        gridLayout->setContentsMargins(11, 11, 11, 11);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        tx = new QLineEdit(paraWidget);
        tx->setObjectName(QString::fromUtf8("tx"));
        QFont font;
        font.setPointSize(20);
        font.setBold(false);
        font.setItalic(true);
        font.setWeight(50);
        tx->setFont(font);

        gridLayout->addWidget(tx, 0, 3, 1, 1);

        label_tz = new QLabel(paraWidget);
        label_tz->setObjectName(QString::fromUtf8("label_tz"));
        QFont font1;
        font1.setPointSize(20);
        font1.setBold(true);
        font1.setWeight(75);
        label_tz->setFont(font1);
        label_tz->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_tz, 2, 2, 1, 1);

        label_tx = new QLabel(paraWidget);
        label_tx->setObjectName(QString::fromUtf8("label_tx"));
        label_tx->setFont(font1);
        label_tx->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_tx, 0, 2, 1, 1);

        label_ty = new QLabel(paraWidget);
        label_ty->setObjectName(QString::fromUtf8("label_ty"));
        label_ty->setFont(font1);
        label_ty->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_ty, 1, 2, 1, 1);

        label_rz = new QLabel(paraWidget);
        label_rz->setObjectName(QString::fromUtf8("label_rz"));
        label_rz->setFont(font1);
        label_rz->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_rz, 2, 0, 1, 1);

        ry = new QLineEdit(paraWidget);
        ry->setObjectName(QString::fromUtf8("ry"));
        ry->setFont(font);

        gridLayout->addWidget(ry, 1, 1, 1, 1);

        label_rx = new QLabel(paraWidget);
        label_rx->setObjectName(QString::fromUtf8("label_rx"));
        label_rx->setFont(font1);
        label_rx->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_rx, 0, 0, 1, 1);

        ty = new QLineEdit(paraWidget);
        ty->setObjectName(QString::fromUtf8("ty"));
        ty->setFont(font);

        gridLayout->addWidget(ty, 1, 3, 1, 1);

        rx = new QLineEdit(paraWidget);
        rx->setObjectName(QString::fromUtf8("rx"));
        rx->setEnabled(true);
        rx->setFont(font);
        rx->setClearButtonEnabled(false);

        gridLayout->addWidget(rx, 0, 1, 1, 1);

        label_ry = new QLabel(paraWidget);
        label_ry->setObjectName(QString::fromUtf8("label_ry"));
        label_ry->setFont(font1);
        label_ry->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_ry, 1, 0, 1, 1);

        rz = new QLineEdit(paraWidget);
        rz->setObjectName(QString::fromUtf8("rz"));
        rz->setFont(font);

        gridLayout->addWidget(rz, 2, 1, 1, 1);

        tz = new QLineEdit(paraWidget);
        tz->setObjectName(QString::fromUtf8("tz"));
        tz->setFont(font);

        gridLayout->addWidget(tz, 2, 3, 1, 1);

        label_similarity = new QLabel(widget);
        label_similarity->setObjectName(QString::fromUtf8("label_similarity"));
        label_similarity->setGeometry(QRect(870, 340, 291, 141));
        QFont font2;
        font2.setPointSize(27);
        font2.setBold(true);
        font2.setWeight(75);
        label_similarity->setFont(font2);
        similarity = new QLabel(widget);
        similarity->setObjectName(QString::fromUtf8("similarity"));
        similarity->setGeometry(QRect(920, 280, 181, 51));
        QFont font3;
        font3.setPointSize(23);
        font3.setBold(true);
        font3.setWeight(75);
        similarity->setFont(font3);
        mainWidget = new QWidget(widget);
        mainWidget->setObjectName(QString::fromUtf8("mainWidget"));
        mainWidget->setGeometry(QRect(0, 0, 831, 121));
        btn_registration = new QPushButton(mainWidget);
        btn_registration->setObjectName(QString::fromUtf8("btn_registration"));
        btn_registration->setGeometry(QRect(557, 11, 171, 51));
        QPalette palette2;
        palette2.setBrush(QPalette::Active, QPalette::Button, brush);
        palette2.setBrush(QPalette::Active, QPalette::Midlight, brush);
        palette2.setBrush(QPalette::Active, QPalette::BrightText, brush);
        palette2.setBrush(QPalette::Active, QPalette::Base, brush);
        palette2.setBrush(QPalette::Active, QPalette::Window, brush);
        palette2.setBrush(QPalette::Active, QPalette::HighlightedText, brush);
        palette2.setBrush(QPalette::Active, QPalette::AlternateBase, brush);
        palette2.setBrush(QPalette::Inactive, QPalette::Button, brush);
        palette2.setBrush(QPalette::Inactive, QPalette::Midlight, brush);
        palette2.setBrush(QPalette::Inactive, QPalette::BrightText, brush);
        palette2.setBrush(QPalette::Inactive, QPalette::Base, brush);
        palette2.setBrush(QPalette::Inactive, QPalette::Window, brush);
        palette2.setBrush(QPalette::Inactive, QPalette::HighlightedText, brush);
        palette2.setBrush(QPalette::Inactive, QPalette::AlternateBase, brush);
        palette2.setBrush(QPalette::Disabled, QPalette::Button, brush);
        palette2.setBrush(QPalette::Disabled, QPalette::Midlight, brush);
        palette2.setBrush(QPalette::Disabled, QPalette::BrightText, brush);
        palette2.setBrush(QPalette::Disabled, QPalette::Base, brush);
        palette2.setBrush(QPalette::Disabled, QPalette::Window, brush);
        palette2.setBrush(QPalette::Disabled, QPalette::HighlightedText, brush);
        palette2.setBrush(QPalette::Disabled, QPalette::AlternateBase, brush);
        btn_registration->setPalette(palette2);
        btn_registration->setFont(font1);
        box_algorithm = new QComboBox(mainWidget);
        box_algorithm->addItem(QString());
        box_algorithm->addItem(QString());
        box_algorithm->addItem(QString());
        box_algorithm->setObjectName(QString::fromUtf8("box_algorithm"));
        box_algorithm->setGeometry(QRect(419, 11, 121, 51));
        QPalette palette3;
        palette3.setBrush(QPalette::Active, QPalette::Button, brush);
        palette3.setBrush(QPalette::Active, QPalette::Base, brush);
        palette3.setBrush(QPalette::Active, QPalette::AlternateBase, brush);
        palette3.setBrush(QPalette::Inactive, QPalette::Button, brush);
        palette3.setBrush(QPalette::Inactive, QPalette::Base, brush);
        palette3.setBrush(QPalette::Inactive, QPalette::AlternateBase, brush);
        palette3.setBrush(QPalette::Disabled, QPalette::Button, brush);
        QBrush brush2(QColor(224, 243, 241, 255));
        brush2.setStyle(Qt::SolidPattern);
        palette3.setBrush(QPalette::Disabled, QPalette::Base, brush2);
        palette3.setBrush(QPalette::Disabled, QPalette::AlternateBase, brush);
        box_algorithm->setPalette(palette3);
        QFont font4;
        font4.setPointSize(19);
        font4.setBold(true);
        font4.setWeight(75);
        box_algorithm->setFont(font4);
        btn_openX = new QPushButton(mainWidget);
        btn_openX->setObjectName(QString::fromUtf8("btn_openX"));
        btn_openX->setGeometry(QRect(214, 11, 191, 51));
        QPalette palette4;
        palette4.setBrush(QPalette::Active, QPalette::Button, brush);
        palette4.setBrush(QPalette::Active, QPalette::Base, brush);
        palette4.setBrush(QPalette::Active, QPalette::AlternateBase, brush);
        palette4.setBrush(QPalette::Inactive, QPalette::Button, brush);
        palette4.setBrush(QPalette::Inactive, QPalette::Base, brush);
        palette4.setBrush(QPalette::Inactive, QPalette::AlternateBase, brush);
        palette4.setBrush(QPalette::Disabled, QPalette::Button, brush);
        palette4.setBrush(QPalette::Disabled, QPalette::Base, brush2);
        palette4.setBrush(QPalette::Disabled, QPalette::AlternateBase, brush);
        btn_openX->setPalette(palette4);
        btn_openX->setFont(font1);
        btn_openCT = new QPushButton(mainWidget);
        btn_openCT->setObjectName(QString::fromUtf8("btn_openCT"));
        btn_openCT->setGeometry(QRect(9, 11, 191, 51));
        QPalette palette5;
        palette5.setBrush(QPalette::Active, QPalette::Button, brush);
        palette5.setBrush(QPalette::Active, QPalette::Base, brush);
        palette5.setBrush(QPalette::Active, QPalette::AlternateBase, brush);
        palette5.setBrush(QPalette::Inactive, QPalette::Button, brush);
        palette5.setBrush(QPalette::Inactive, QPalette::Base, brush);
        palette5.setBrush(QPalette::Inactive, QPalette::AlternateBase, brush);
        palette5.setBrush(QPalette::Disabled, QPalette::Button, brush);
        palette5.setBrush(QPalette::Disabled, QPalette::Base, brush2);
        palette5.setBrush(QPalette::Disabled, QPalette::AlternateBase, brush);
        btn_openCT->setPalette(palette5);
        btn_openCT->setFont(font1);
        imgWidget = new QWidget(widget);
        imgWidget->setObjectName(QString::fromUtf8("imgWidget"));
        imgWidget->setGeometry(QRect(10, 130, 821, 541));
        QPalette palette6;
        QBrush brush3(QColor(0, 0, 0, 255));
        brush3.setStyle(Qt::SolidPattern);
        palette6.setBrush(QPalette::Active, QPalette::Button, brush3);
        palette6.setBrush(QPalette::Active, QPalette::Light, brush3);
        palette6.setBrush(QPalette::Active, QPalette::Midlight, brush3);
        palette6.setBrush(QPalette::Active, QPalette::BrightText, brush3);
        QBrush brush4(QColor(50, 50, 50, 255));
        brush4.setStyle(Qt::SolidPattern);
        palette6.setBrush(QPalette::Active, QPalette::Base, brush4);
        QBrush brush5(QColor(34, 34, 34, 255));
        brush5.setStyle(Qt::SolidPattern);
        palette6.setBrush(QPalette::Active, QPalette::Window, brush5);
        palette6.setBrush(QPalette::Active, QPalette::AlternateBase, brush3);
        palette6.setBrush(QPalette::Inactive, QPalette::Button, brush3);
        palette6.setBrush(QPalette::Inactive, QPalette::Light, brush3);
        palette6.setBrush(QPalette::Inactive, QPalette::Midlight, brush3);
        palette6.setBrush(QPalette::Inactive, QPalette::BrightText, brush3);
        palette6.setBrush(QPalette::Inactive, QPalette::Base, brush4);
        palette6.setBrush(QPalette::Inactive, QPalette::Window, brush5);
        palette6.setBrush(QPalette::Inactive, QPalette::AlternateBase, brush3);
        palette6.setBrush(QPalette::Disabled, QPalette::Button, brush3);
        palette6.setBrush(QPalette::Disabled, QPalette::Light, brush3);
        palette6.setBrush(QPalette::Disabled, QPalette::Midlight, brush3);
        palette6.setBrush(QPalette::Disabled, QPalette::BrightText, brush3);
        palette6.setBrush(QPalette::Disabled, QPalette::Base, brush5);
        palette6.setBrush(QPalette::Disabled, QPalette::Window, brush5);
        palette6.setBrush(QPalette::Disabled, QPalette::AlternateBase, brush3);
        imgWidget->setPalette(palette6);
        imgWidget->setCursor(QCursor(Qt::CrossCursor));
        drr_label = new QLabel(imgWidget);
        drr_label->setObjectName(QString::fromUtf8("drr_label"));
        drr_label->setGeometry(QRect(20, 80, 350, 350));
        drr_label->setTextFormat(Qt::AutoText);
        drr_label->setScaledContents(true);
        drr_label->setAlignment(Qt::AlignCenter);
        xray_label = new QLabel(imgWidget);
        xray_label->setObjectName(QString::fromUtf8("xray_label"));
        xray_label->setGeometry(QRect(400, 80, 350, 350));
        xray_label->setScaledContents(true);
        xray_label->setAlignment(Qt::AlignCenter);
        drr_label_2 = new QLabel(imgWidget);
        drr_label_2->setObjectName(QString::fromUtf8("drr_label_2"));
        drr_label_2->setGeometry(QRect(0, 300, 250, 250));
        drr_label_2->setTextFormat(Qt::AutoText);
        drr_label_2->setScaledContents(true);
        drr_label_2->setAlignment(Qt::AlignCenter);
        xray_label_2 = new QLabel(imgWidget);
        xray_label_2->setObjectName(QString::fromUtf8("xray_label_2"));
        xray_label_2->setGeometry(QRect(440, 300, 250, 250));
        xray_label_2->setScaledContents(true);
        xray_label_2->setAlignment(Qt::AlignCenter);
        progressBar = new QProgressBar(widget);
        progressBar->setObjectName(QString::fromUtf8("progressBar"));
        progressBar->setGeometry(QRect(1000, 620, 191, 23));
        progressBar->setMaximum(0);
        progressBar->setValue(-1);
        label_similarity_2 = new QLabel(widget);
        label_similarity_2->setObjectName(QString::fromUtf8("label_similarity_2"));
        label_similarity_2->setGeometry(QRect(950, 560, 131, 51));

        horizontalLayout->addWidget(widget);

        MainWindowClass->setCentralWidget(centralWidget);
        mainToolBar = new QToolBar(MainWindowClass);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        mainToolBar->setEnabled(true);
        QPalette palette7;
        palette7.setBrush(QPalette::Active, QPalette::Light, brush1);
        palette7.setBrush(QPalette::Active, QPalette::Base, brush1);
        palette7.setBrush(QPalette::Active, QPalette::Window, brush1);
        palette7.setBrush(QPalette::Inactive, QPalette::Light, brush1);
        palette7.setBrush(QPalette::Inactive, QPalette::Base, brush1);
        palette7.setBrush(QPalette::Inactive, QPalette::Window, brush1);
        palette7.setBrush(QPalette::Disabled, QPalette::Light, brush1);
        palette7.setBrush(QPalette::Disabled, QPalette::Base, brush1);
        palette7.setBrush(QPalette::Disabled, QPalette::Window, brush1);
        mainToolBar->setPalette(palette7);
        MainWindowClass->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(MainWindowClass);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        MainWindowClass->setStatusBar(statusBar);

        retranslateUi(MainWindowClass);

        QMetaObject::connectSlotsByName(MainWindowClass);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindowClass)
    {
        MainWindowClass->setWindowTitle(QCoreApplication::translate("MainWindowClass", "\346\235\255\345\267\236\345\270\210\350\214\203\345\244\247\345\255\246\345\214\273\345\255\246\345\275\261\345\203\217\351\205\215\345\207\206\347\250\213\345\272\217", nullptr));
        tx->setText(QCoreApplication::translate("MainWindowClass", "0", nullptr));
        label_tz->setText(QCoreApplication::translate("MainWindowClass", "Tz", nullptr));
        label_tx->setText(QCoreApplication::translate("MainWindowClass", "Tx", nullptr));
        label_ty->setText(QCoreApplication::translate("MainWindowClass", "Ty", nullptr));
        label_rz->setText(QCoreApplication::translate("MainWindowClass", "Rz", nullptr));
        ry->setText(QCoreApplication::translate("MainWindowClass", "0", nullptr));
        label_rx->setText(QCoreApplication::translate("MainWindowClass", "Rx", nullptr));
        ty->setText(QCoreApplication::translate("MainWindowClass", "0", nullptr));
        rx->setText(QCoreApplication::translate("MainWindowClass", "90", nullptr));
        label_ry->setText(QCoreApplication::translate("MainWindowClass", "Ry", nullptr));
        rz->setText(QCoreApplication::translate("MainWindowClass", "0", nullptr));
        tz->setText(QCoreApplication::translate("MainWindowClass", "0", nullptr));
        label_similarity->setText(QString());
        similarity->setText(QCoreApplication::translate("MainWindowClass", "\347\233\270\344\274\274\345\272\246\357\274\232", nullptr));
        btn_registration->setText(QCoreApplication::translate("MainWindowClass", "\345\274\200\345\247\213\351\205\215\345\207\206", nullptr));
        box_algorithm->setItemText(0, QCoreApplication::translate("MainWindowClass", "PSO", nullptr));
        box_algorithm->setItemText(1, QCoreApplication::translate("MainWindowClass", "EO", nullptr));
        box_algorithm->setItemText(2, QCoreApplication::translate("MainWindowClass", "DE", nullptr));

        box_algorithm->setCurrentText(QCoreApplication::translate("MainWindowClass", "PSO", nullptr));
        btn_openX->setText(QCoreApplication::translate("MainWindowClass", "\350\257\273\345\217\226X\345\205\211\346\226\207\344\273\266", nullptr));
        btn_openCT->setText(QCoreApplication::translate("MainWindowClass", "\350\257\273\345\217\226CT\346\226\207\344\273\266", nullptr));
        drr_label->setText(QString());
        xray_label->setText(QString());
        drr_label_2->setText(QString());
        xray_label_2->setText(QString());
        label_similarity_2->setText(QString());
    } // retranslateUi

};

namespace Ui {
    class MainWindowClass: public Ui_MainWindowClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
