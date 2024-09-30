#include "dataProcess.h"

DataProcess::DataProcess()
{
	DRRFile = "tmp.bmp";
}

DataProcess::~DataProcess()
{
}

void DataProcess::createXRayReader(const char* input)
{
	XRayReader = vtkSmartPointer<vtkBMPReader>::New();
	XRayReader->SetFileName(input);
	XRayReader->Update();
	XRayFile = input;
}
void DataProcess::clearXRayReader()
{
	XRayReader = vtkSmartPointer<vtkBMPReader>::New();
}

vtkSmartPointer<vtkBMPReader> DataProcess::getXRayReader()
{
	return XRayReader;
}

void DataProcess::setXRayPath(std::string filename)
{
	XRayFile = filename;
}

std::string DataProcess::getXRayPath()
{
	return XRayFile;
}
bool DataProcess::isRotationMatrix(cv::Mat& R)
{
	cv::Mat Rt;
	transpose(R, Rt);
	cv::Mat shouldBeIdentity = Rt * R;
	cv::Mat I = cv::Mat::eye(3, 3, shouldBeIdentity.type());
	return  norm(I, shouldBeIdentity) < 1e-6;
}

// Calculates rotation matrix to euler angles
// The result is the same as MATLAB except the order
// of the euler angles ( x and z are swapped ).

cv::Vec3f DataProcess::rotationMatrixToEulerAngles(cv::Mat& R)
{
	/*if (isRotationMatrix(R)) {
		cout << " is not a rotation matrix!" << std::endl;
		abort();
	}*/
	//assert(isRotationMatrix(R));

	
	float sy = sqrt(R.at<double>(0, 0) * R.at<double>(0, 0) + R.at<double>(1, 0) * R.at<double>(1, 0));

	bool singular = sy < 1e-6; // If

	float x, y, z;
	if (!singular)
	{
		x = atan2(R.at<double>(2, 1), R.at<double>(2, 2));
		y = atan2(-R.at<double>(2, 0), sy);
		z = atan2(R.at<double>(1, 0), R.at<double>(0, 0));
	}
	else
	{
		x = atan2(-R.at<double>(1, 2), R.at<double>(1, 1));
		y = atan2(-R.at<double>(2, 0), sy);
		z = 0;
	}
	return cv::Vec3f(x, y, z);
}
cv::Mat DataProcess::eulerAnglesToRotationMatrix(cv::Vec3f& theta)
{
	// Calculate rotation about x axis
	cv::Mat R_x = (cv::Mat_<double>(3, 3) <<
		1, 0, 0,
		0, cos(theta[0]), -sin(theta[0]),
		0, sin(theta[0]), cos(theta[0])
		);
	// Calculate rotation about y axis
	cv::Mat R_y = (cv::Mat_<double>(3, 3) <<
		cos(theta[1]), 0, sin(theta[1]),
		0, 1, 0,
		-sin(theta[1]), 0, cos(theta[1])
		);
	// Calculate rotation about z axis
	cv::Mat R_z = (cv::Mat_<double>(3, 3) <<
		cos(theta[2]), -sin(theta[2]), 0,
		sin(theta[2]), cos(theta[2]), 0,
		0, 0, 1);
	// Combined rotation matrix
	cv::Mat R = R_z * R_y * R_x;
	return R;
}