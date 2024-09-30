#pragma once
#include <cstdarg>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCenteredEuler3DTransform.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkNormalizeImageFilter.h"
#include "itkImageSeriesReader.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkRayCastInterpolateImageFunction.h"
#include "itkBMPImageIOFactory.h"
#include "itkGiplImageIOFactory.h"
#include "itkImageToVTKImageFilter.h"
#include "itkMetaDataObject.h"
#include "itkIntensityWindowingImageFilter.h"
#include "vtkImageIterator.h"
//#include "ImageAlgorithm.h"
#include "device_functions.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "dataProcess.h"
#include "CUDAdrr.cuh"
#include  <winsock2.h> 
#pragma comment(lib, "WS2_32.lib")
#pragma comment(lib, "Rpcrt4.lib")


#include <Qmessagebox>
class itkDrr //: public ImageAlgorithm
{
public:
	itkDrr(std::string input);
	~itkDrr();
	
	
	using InputImageType = typename itk::Image<float, 3>;
	using PointType = typename InputImageType::PointType;

	//void imageProcess(int count, ...);
	//void setInputName(std::string input);
	vtkDataObject* getOutput();
	itk::ImageToVTKImageFilter<itk::Image<unsigned char, 3>>::Pointer ivfilter;
	float* cpp_object3D;
	double* cudaUCharDRR();
	double* cudaUCharDRRLateral();
	void setParameters(std::vector<double> parameters);
	void setParameters(std::vector<double> parameters, int dx, int dy, float threshold);
	cv::Mat Pos2Matrix(std::vector<double> para);
private:
	bool isRead;
	itk::ImageSeriesReader< itk::Image< signed short, 3 > >::Pointer reader;
	itk::Image<float, 3>::Pointer image;
	char* input_name;

	bool m_GPUEnabled;

	//std::vector<float> parameters;
	std::vector<float> rotation;
	std::vector<float> translation;
	int dx = 1024;
	int dy = 1024;
	float threshold = 0;
};

