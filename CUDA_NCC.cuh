#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include<iostream>
#include<opencv.hpp>
#include<vector>


using namespace cv;
using namespace std;
//#define uchar unsigned char


double* CUDA_NCC(Mat c_srcImg, uchar* f_srcImg, int PARAS_NUMS);