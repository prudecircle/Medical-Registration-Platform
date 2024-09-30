#pragma once
#include <stdlib.h>
#include<vector>
#include<opencv.hpp>
/*vtk File*/
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkImageViewer2.h>
#include <QVTKWidget.h>
#include <vtkDICOMImageReader.h>
#include <vtkJPEGReader.h>
#include <vtkPNGReader.h>
#include <vtkImageActor.h>
#include <vtkEventQtSlotConnect.h>
#include "vtkCommand.h"
#include <vtkOutputWindow.h>
#include "vtkSmartPointer.h"
#include "vtkImageBlend.h"
#include "vtkBMPReader.h"
#include "vtkImageCast.h"
#include "vtkBMPWriter.h"
#include "vtkFixedPointVolumeRayCastMapper.h"
#include "vtkColorTransferFunction.h"
#include "vtkPiecewiseFunction.h"
#include "vtkVolumeProperty.h"
#include "vtkExtractVOI.h"
#include <vtkDataObject.h>
#include <vtkWindows.h>
#include <qobject.h>

class DataProcess :
	public QObject
{
	Q_OBJECT
public:
	DataProcess();
	~DataProcess();
	void createXRayReader(const char* input);
	void clearXRayReader();
	vtkSmartPointer<vtkBMPReader> getXRayReader();
	void setXRayPath(std::string filename);
	std::string getXRayPath();
	bool isRotationMatrix(cv::Mat& R);
	cv::Vec3f rotationMatrixToEulerAngles(cv::Mat& R);
	cv::Mat eulerAnglesToRotationMatrix(cv::Vec3f& theta);
private:
	vtkSmartPointer<vtkBMPReader> XRayReader;
	std::string XRayFile;
	std::string DRRFile;


};