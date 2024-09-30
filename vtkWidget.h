#ifndef VTKWIDGET_H
#define WTKWIDGET_H

#include <QDialog>
#include <QWidget>
#include "vtkVolume.h"
#include "vtkRenderer.h"
//#include "InputDialog.h"
#include "vtkSmartPointer.h"
#include "vtkPiecewiseFunction.h"
#include "vtkVolumeProperty.h"
#include "vtkWindowToImageFilter.h"
#include "vtkDICOMImageReader.h"
#include "vtkImageShiftScale.h"
#include "vtkOpenGLGPUVolumeRayCastMapper.h"
#include "vtkColorTransferFunction.h"
#include "vtkBMPWriter.h"
#include "vtkFixedPointVolumeRayCastMapper.h"
#include "vtkTransform.h"
#include "vtkCamera.h"
//#include "vtkWidgetTrackballCamera.h"
#include "vtkMatrix4x4.h"

namespace Ui {
	class vtkWidget;
}

class vtkWidget : public QDialog
{
	Q_OBJECT

public:
	explicit vtkWidget(QWidget *parent = 0);
	~vtkWidget();
	//void setVolume(vtkVolume* volume);
	void setFilePath(const char* dicname);

public slots:
	void create2D();
	void addScalarOpacity();
	void deleteScalarOpacity();
	void addGradientOpacity();
	void deleteGradientOpacity();
	void ReceiveScalarData(QList<QString>* inputList);
	void ReceiveGradientData(QList<QString>* inputList);
	void changeVolume();
	void restorePara();
	void setPara();
	void myShow();
	void changeAngleAnytime(double,double,double);
	void refreshAngle(double, double, double);
	
signals:
	void sendQPixmap(QImage p);
	void sendRotate(double,double,double);
	void sendAngle(double,double,double);

private:
	const char* path;
	Ui::vtkWidget* ui;
	vtkSmartPointer<vtkRenderer> tdRenderer;
	int column;
	vtkSmartPointer<vtkVolume> drrVolume;
	//InputDialog* inputDialog;
	double x;
	double y;
	double z;
	vtkDICOMImageReader* reader;
	//vtkWidgetTrackballCamera *style;
	vtkMatrix4x4* volumeMatrix;
	bool needRefresh;

	void calAngleFromMatrix(double& x, double& y, double& z);

};
#endif // VTKWIDGET_H
