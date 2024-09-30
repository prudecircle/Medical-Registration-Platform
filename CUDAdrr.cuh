
#pragma once
#include <vector>
#include <iostream>
#include <opencv.hpp>
#include "CUDA_NCC.cuh"
//#include "MainWindow.h"



extern float* d_object3D;

// This variable contains the 2D output from CUDA
extern float* d_object2D;



struct Image3D
{
	float* image;
	float PixelSpacingCT[3];
	float isoCenter[3];
	int SizeCT[3];
};

struct Image2D
{
	float* image;
	int size[2]; // [rows, cols]
	//int rows, cols;
};

struct CUDAParamerters
{
	int numThreads;
	int numBlocks;
};

struct DRRParameters
{
	float stepInX[12];
	float stepInY[12];
	float corner00[12];
	float SourceWorld[12];
	float size[2];
	float threshold;
};

void loadDICOMInGPUMemory(float* cpp_object3D, int* sizeCT, float* pixelSpacing);
void loadOuputVariablesInGPUMemory(int dimX, int dimZ, int PARAS_NUMS, int flag = 0);
void loadXRayImg(std::string filename);
void freeDICOMFromGPUMemory();
void freeAuxiliaryVariablesInGPUMemory(int flag = 0);

//float* calFloatDRRwithCUDA(CUDAParamerters CUDA_Parameters, DRRParameters DRR_Parameters);
void calUCharDRRwithCUDA(CUDAParamerters CUDA_Parameters, DRRParameters DRR_Parameters, int PARAS_NUMS,bool isFrontal);



