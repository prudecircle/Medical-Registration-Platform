#include "CUDAdrr.cuh"
//#include"CUDA_NCC.cuh"
#include <cuda.h>
#include <cuda_runtime_api.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "texture_fetch_functions.h"
#include <algorithm> 
#include <iostream>
#include <vector>
#include <string>

#define KERNEL                      __global__
#define HOST                        __host__
#define DEVICE                      __device__
#define HOST_AND_DEVICE             __host__ __device__
#define DEVICE_CONST                __device__ __constant__




// This variable contains the DICOM set
float* d_object3D;

// This variable contains the 2D output from CUDA
float* d_object2D;

//X光图片
cv::Mat m_srcImg;

unsigned char* c_object2D;



// Constants depending on the DICOM
DEVICE_CONST int d_sizeCT[3];
DEVICE_CONST  float ctPixelSpacing[3];

// Constant depending on image output,单张DRR的尺寸
DEVICE_CONST int DRRImageSize[2];

// Constants dependion on the specific DRR
DEVICE_CONST  float d_DRR_Parameters[4 * 12 + 3];


// This variable contains the DICOM loaded as a Texture ( read-only, fast-cached memory)
cudaTextureObject_t tex_object3D = 0;

cudaStream_t stream1;

__global__ void cal_hist(int* hist, int* range) {
	range[0] = 0;
	range[1] = 0;
	for (int i = 0; i < 256; i++) {
		if (hist[i] != 0) {
			range[0] = i;
			break;
		}
	}
	for (int i = 255; i >= 0; i--) {
		if (hist[i] != 0) {
			range[1] = i;
			break;
		}
	}
}
__global__ void hist_uc(unsigned char* object2D, int* range, int PARAS_NUMS) {
	int total_dx = DRRImageSize[0];
	int total_dz = DRRImageSize[1];

	//Every thread calculates its own id number
	long int idx = (blockIdx.x * blockDim.x) + threadIdx.x;

	// This checks if the thread number is bigger than the amount of pixels
	if (idx >= total_dx * total_dz * PARAS_NUMS)
		return;

	float tmp = object2D[idx];

	object2D[idx] = (int)((tmp - range[0]) / (range[1] - range[0]) * 255);

}
__device__ float getPixval(cudaTextureObject_t tex_object3D, long int idx, int PARAS_NUMS) {

	//-------------kfq7.14
	int total_dx = DRRImageSize[0];
	int total_dz = DRRImageSize[1];

	int idxdz = idx / total_dx;
	int imgIdx = idx / (total_dx * total_dz);
	int dz, dx;

	float stepInX[3];
	float stepInY[3];
	float corner00[3];
	float SourceWorld[3];
	float threshold = d_DRR_Parameters[50];

	stepInX[0] = d_DRR_Parameters[imgIdx * 3];
	stepInX[1] = d_DRR_Parameters[imgIdx * 3 + 1];
	stepInX[2] = d_DRR_Parameters[imgIdx * 3 + 2];
	stepInY[0] = d_DRR_Parameters[12 + imgIdx * 3];
	stepInY[1] = d_DRR_Parameters[13 + imgIdx * 3];
	stepInY[2] = d_DRR_Parameters[14 + imgIdx * 3];
	corner00[0] = d_DRR_Parameters[24 + imgIdx * 3];
	corner00[1] = d_DRR_Parameters[25 + imgIdx * 3];
	corner00[2] = d_DRR_Parameters[26 + imgIdx * 3];
	SourceWorld[0] = d_DRR_Parameters[36 + imgIdx * 3];
	SourceWorld[1] = d_DRR_Parameters[37 + imgIdx * 3];
	SourceWorld[2] = d_DRR_Parameters[38 + imgIdx * 3];
	dz = idxdz - imgIdx * total_dz;
	dx = idx - idxdz * total_dx;


	//Calculate the spatial position of the pixel
	//drrPixelWorld_0[idx] = *corner00_0 + ((*stepInX_0)*(threadIdx.x)) + ((*stepInY_0)*(blockIdx.x));
	//drrPixelWorld_1[idx] = *corner00_1 + ((*stepInX_1)*(threadIdx.x)) + ((*stepInY_1)*(blockIdx.x));
	//drrPixelWorld_2[idx] = *corner00_2 + ((*stepInX_2)*(threadIdx.x)) + ((*stepInY_2)*(blockIdx.x));
	float drrPixelWorld[3] = { 0 };
	drrPixelWorld[0] = corner00[0] + ((stepInX[0]) * dx) + ((stepInY[0]) * dz);
	drrPixelWorld[1] = corner00[1] + ((stepInX[1]) * dx) + ((stepInY[1]) * dz);
	drrPixelWorld[2] = corner00[2] + ((stepInX[2]) * dx) + ((stepInY[2]) * dz);

	//Calculate the ray vector
	float rayVector[3] = { 0 };
	rayVector[0] = drrPixelWorld[0] - SourceWorld[0];
	rayVector[1] = drrPixelWorld[1] - SourceWorld[1];
	rayVector[2] = drrPixelWorld[2] - SourceWorld[2];

	float alpha1[3];
	float alphaN[3];
	float auxalphaMin[3] = { -2, -2, -2 };
	float auxalphaMax[3] = { 2 , 2 , 2 };


	//Calculate alphaMin and alphaMax 
	if (rayVector[2] != 0)
	{
		alpha1[0] = (0.0 - (SourceWorld[2])) / rayVector[2];
		alphaN[0] = ((d_sizeCT[2]) * (ctPixelSpacing[2]) - (SourceWorld[2])) / rayVector[2];
		auxalphaMin[0] = alphaN[0];
		auxalphaMax[0] = alpha1[0];

		if (alpha1[0] < alphaN[0])
		{
			auxalphaMin[0] = alpha1[0];
			auxalphaMax[0] = alphaN[0];
		}
	}

	if (rayVector[1] != 0)
	{
		alpha1[1] = (0.0 - (SourceWorld[1])) / rayVector[1];
		alphaN[1] = ((d_sizeCT[1]) * (ctPixelSpacing[1]) - (SourceWorld[1])) / rayVector[1];
		auxalphaMin[1] = alphaN[1];
		auxalphaMax[1] = alpha1[1];

		if (alpha1[1] < alphaN[1])
		{
			auxalphaMin[1] = alpha1[1];
			auxalphaMax[1] = alphaN[1];
		}
	}


	if (rayVector[0] != 0)
	{
		alpha1[2] = (0.0 - (SourceWorld[0])) / rayVector[0];
		alphaN[2] = ((d_sizeCT[0]) * (ctPixelSpacing[0]) - (SourceWorld[0])) / rayVector[0];
		auxalphaMin[2] = alphaN[2];
		auxalphaMax[2] = alpha1[2];

		if (alpha1[2] < alphaN[2])
		{
			auxalphaMin[2] = alpha1[2];
			auxalphaMax[2] = alphaN[2];
		}
	}


	float alphaMin;

	if (auxalphaMin[0] > auxalphaMin[1]) //x > y
	{
		alphaMin = auxalphaMin[2];
		if (auxalphaMin[0] > alphaMin) { //x > y, x > z
			alphaMin = auxalphaMin[0];
		}
	}
	else //y > x
	{
		alphaMin = auxalphaMin[2];
		if (auxalphaMin[1] > alphaMin)  //y > x, y > z
			alphaMin = auxalphaMin[1];
	}

	float alphaMax;

	if (auxalphaMax[0] < auxalphaMax[1])  // x < y
	{
		alphaMax = auxalphaMax[2];
		if (auxalphaMax[0] < alphaMax)  // x < y, x < z
			alphaMax = auxalphaMax[0];
	}
	else // y < x
	{
		alphaMax = auxalphaMax[2];
		if (auxalphaMax[1] < alphaMax)  // y < x, y < z
			alphaMax = auxalphaMax[1];
	}

	float firstIntersection[3], firstIntersectionIndex[3], firstIntersectionIndexUp[3], firstIntersectionIndexDown[3];

	//Calculate the first intersection of the ray with the planes (alphaX, alphaY and alphaZ)
	firstIntersection[0] = (SourceWorld[0]) + (alphaMin * rayVector[0]);
	firstIntersection[1] = (SourceWorld[1]) + (alphaMin * rayVector[1]);
	firstIntersection[2] = (SourceWorld[2]) + (alphaMin * rayVector[2]);

	firstIntersectionIndex[0] = firstIntersection[0] / (ctPixelSpacing[0]);
	firstIntersectionIndex[1] = firstIntersection[1] / (ctPixelSpacing[1]);
	firstIntersectionIndex[2] = firstIntersection[2] / (ctPixelSpacing[2]);


	firstIntersectionIndexUp[0] = (int)ceil(firstIntersectionIndex[0]);
	firstIntersectionIndexUp[1] = (int)ceil(firstIntersectionIndex[1]);
	firstIntersectionIndexUp[2] = (int)ceil(firstIntersectionIndex[2]);

	firstIntersectionIndexDown[0] = (int)floor(firstIntersectionIndex[0]);
	firstIntersectionIndexDown[1] = (int)floor(firstIntersectionIndex[1]);
	firstIntersectionIndexDown[2] = (int)floor(firstIntersectionIndex[2]);

	float alpha[3] = { 2,2,2 }, alphaIntersectionUp[3], alphaIntersectionDown[3];

	if (rayVector[2] != 0)
	{
		alphaIntersectionUp[2] = (firstIntersectionIndexUp[2] * (ctPixelSpacing[2]) - (SourceWorld[2])) / rayVector[2];
		alphaIntersectionDown[2] = (firstIntersectionIndexDown[2] * (ctPixelSpacing[2]) - (SourceWorld[2])) / rayVector[2];
		alpha[0] = alphaIntersectionDown[2];
		if (alphaIntersectionUp[2] > alpha[0])
			alpha[0] = alphaIntersectionUp[2];
	}

	if (rayVector[1] != 0)
	{
		alphaIntersectionUp[1] = (firstIntersectionIndexUp[1] * (ctPixelSpacing[1]) - (SourceWorld[1])) / rayVector[1];
		alphaIntersectionDown[1] = (firstIntersectionIndexDown[1] * (ctPixelSpacing[1]) - (SourceWorld[1])) / rayVector[1];
		alpha[1] = alphaIntersectionDown[1];
		if (alphaIntersectionUp[1] > alpha[1])
			alpha[1] = alphaIntersectionUp[1];
	}

	if (rayVector[0] != 0)
	{
		alphaIntersectionUp[0] = (firstIntersectionIndexUp[0] * (ctPixelSpacing[0]) - (SourceWorld[0])) / rayVector[0];
		alphaIntersectionDown[0] = (firstIntersectionIndexDown[0] * (ctPixelSpacing[0]) - (SourceWorld[0])) / rayVector[0];
		alpha[2] = alphaIntersectionDown[0];
		if (alphaIntersectionUp[0] > alpha[2])
			alpha[2] = alphaIntersectionUp[0];
	}

	float alphaU[3] = { 999,999,999 };
	//Calculate incremental values (alphaUx, alphaUx, alphaUz) when the ray intercepts the planes
	if (rayVector[2] != 0)
		alphaU[0] = (ctPixelSpacing[2]) / fabs(rayVector[2]);

	if (rayVector[1] != 0)
		alphaU[1] = (ctPixelSpacing[1]) / fabs(rayVector[1]);

	if (rayVector[0] != 0)
		alphaU[2] = (ctPixelSpacing[0]) / fabs(rayVector[0]);


	float U[3] = { -1,-1,-1 };
	// Calculate voxel index incremental values along the ray path
	if ((SourceWorld[2]) < drrPixelWorld[2])
		U[0] = 1;

	if ((SourceWorld[1]) < drrPixelWorld[1])
		U[1] = 1;

	if ((SourceWorld[0]) < drrPixelWorld[0])
		U[2] = 1;


	//Initialize the weighted sum to zero
	float d12 = 0.0, alphaCmin, alphaCminPrev;

	//Initialize the current ray position (alphaCmin)
	if (alpha[0] < alpha[1]) //x < y
	{
		alphaCmin = alpha[2];
		if (alpha[0] < alphaCmin)  //x < y, x < z
			alphaCmin = alpha[0];
	}
	else //y < x
	{
		alphaCmin = alpha[2];
		if (alpha[1] < alphaCmin)  //y < x, y < z
			alphaCmin = alpha[1];
	}

	// Initialize the current voxel index.
	float cIndexNumber[3] = { firstIntersectionIndexDown[0] , firstIntersectionIndexDown[1] , firstIntersectionIndexDown[2] };

	//The while loop represents when the ray is inside the volume
	while (alphaCmin < alphaMax)
	{
		// Store the current ray position 
		alphaCminPrev = alphaCmin;

		if ((alpha[0] <= alpha[1]) && (alpha[0] <= alpha[2])) //Ray front intercepts with x-plane. Update alphaX
		{
			alphaCmin = alpha[0];
			cIndexNumber[2] = cIndexNumber[2] + U[0];
			alpha[0] = alpha[0] + alphaU[0];
		}
		else if ((alpha[1] <= alpha[0]) && (alpha[1] <= alpha[2])) //Ray front intercepts with y-plane. Update alphaY
		{
			alphaCmin = alpha[1];
			cIndexNumber[1] = cIndexNumber[1] + U[1];
			alpha[1] = alpha[1] + alphaU[1];
		}
		else                                                                //Ray front intercepts with z-plane. Update alphaZ
		{
			alphaCmin = alpha[2];
			cIndexNumber[0] = cIndexNumber[0] + U[2];
			alpha[2] = alpha[2] + alphaU[2];
		}


		if ((cIndexNumber[0] >= 0) && (cIndexNumber[0] < (d_sizeCT[0])) &&
			(cIndexNumber[1] >= 0) && (cIndexNumber[1] < (d_sizeCT[1])) &&
			(cIndexNumber[2] >= 0) && (cIndexNumber[2] < (d_sizeCT[2])))
		{
			//If it is a valid index, get the voxel intensity

			int cIndexCoordinate[3] = { static_cast<int> (cIndexNumber[2]) ,static_cast<int> (cIndexNumber[1]) ,static_cast<int> (cIndexNumber[0]) };

			//Get current position from flat object
			long int currentPos3D = cIndexCoordinate[0] + (cIndexCoordinate[1] * (d_sizeCT[2])) + (cIndexCoordinate[2] * (d_sizeCT[2]) * (d_sizeCT[1]));

			//Retrieve intensity value from flat object
			float value = tex1Dfetch<float>(tex_object3D, currentPos3D);

			//Ignore voxels whose intensities are below the desired threshold
			if (value > threshold)
				d12 += value * (alphaCmin - alphaCminPrev);//weighted sum				
		}
	} //end of the while-loop

	float pixval = 255.0 - d12;
	if (pixval < 0)
		pixval = 255.;

	if (pixval > 255)
		pixval = 0.;

	return pixval;
}

__global__ void drrCUDA_uc(unsigned char* object2D, cudaTextureObject_t tex_object3D, int* hist, int PARAS_NUMS) {

	int total_dx = DRRImageSize[0];
	int total_dz = DRRImageSize[1];

	//Every thread calculates its own id number
	long int idx = (blockIdx.x * blockDim.x) + threadIdx.x;
	/*if (idx == 1) {
		for (int i = 0; i < 39; i++) {
			printf("d_drr[%d]=%d\n", i, d_DRR_Parameters[i]);
		}
	}*/
	// This checks if the thread number is bigger than the amount of pixels
	if (idx >= total_dx * total_dz * PARAS_NUMS)
		return;
	//--------------------kfq 7.12



	/*float pixval1 = getPixval(tex_object3D, idx, 0);
	int tmp1 = (int)pixval1;
	atomicAdd(hist + tmp1, 1);
	object2D1[idx] = (int)pixval1;

	float pixval2 = getPixval(tex_object3D, idx, 1);
	int tmp2 = (int)pixval2;
	atomicAdd(hist + tmp2, 1);
	object2D2[idx] = (int)pixval2;

	float pixval3 = getPixval(tex_object3D, idx, 2);
	int tmp3 = (int)pixval3;
	atomicAdd(hist + tmp3, 1);
	object2D3[idx] = (int)pixval3;*/







	float pixval = getPixval(tex_object3D, idx, PARAS_NUMS);
	int tmp = (int)pixval;
	atomicAdd(hist + tmp, 1);
	/*if (idx < total_dx * total_dz) {
		atomicAdd(hist1 + tmp, 1);
	}
	if (idx >= total_dx * total_dz && idx < 2 * total_dx * total_dz) {
		atomicAdd(hist2 + tmp, 1);
	}
	if (idx >= 2*total_dx * total_dz && idx < 3 * total_dx * total_dz) {
		atomicAdd(hist3 + tmp, 1);
	}
	if (idx >=3* total_dx * total_dz && idx < 4 * total_dx * total_dz) {
		atomicAdd(hist4 + tmp, 1);
	}*/
	//Assign the calculated value for the pixel to its corresponding position in the output array
	object2D[idx] = (int)pixval;
}

void loadDICOMInGPUMemory(float* cpp_object3D, int* sizeCT, float* pixelSpacing)
{
	long int object3Dsize = sizeCT[0] * sizeCT[1] * sizeCT[2];
	cudaStreamCreateWithFlags(&stream1, cudaStreamNonBlocking);


	cudaMalloc((void**)&d_object3D, object3Dsize * sizeof(float));
	cudaMemcpyAsync(d_object3D, cpp_object3D, object3Dsize * sizeof(float), cudaMemcpyHostToDevice, stream1);

	cudaMemcpyToSymbol(ctPixelSpacing, pixelSpacing, 3 * sizeof(float), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_sizeCT, sizeCT, 3 * sizeof(int), 0, cudaMemcpyHostToDevice);


	// create texture object
	cudaResourceDesc resDesc;
	memset(&resDesc, 0, sizeof(resDesc));
	resDesc.resType = cudaResourceTypeLinear;
	resDesc.res.linear.devPtr = d_object3D;
	resDesc.res.linear.desc.f = cudaChannelFormatKindFloat;
	resDesc.res.linear.desc.x = 32; // bits per channel
	resDesc.res.linear.sizeInBytes = object3Dsize * sizeof(float);

	cudaTextureDesc texDesc;
	memset(&texDesc, 0, sizeof(texDesc));
	texDesc.readMode = cudaReadModeElementType;

	// create CUDA texture object
	cudaDestroyTextureObject(tex_object3D);
	cudaCreateTextureObject(&tex_object3D, &resDesc, &texDesc, NULL);

	cudaStreamDestroy(stream1);

}
void loadOuputVariablesInGPUMemory(int dimX, int dimZ, int PARAS_NUMS, int flag)
{
	/*long int vectorSize = dimX * dimZ;
	int OutputImageSize[2] = { dimX, dimZ };*/

	//----------------------kfq7.9
	long int vectorSize = dimX * dimZ;
	int OutputImageSize[2] = { dimX, dimZ };

	if (flag == 0)
		cudaMalloc((void**)&d_object2D, PARAS_NUMS * vectorSize * sizeof(float));
	else {
		cudaMalloc((void**)&c_object2D, PARAS_NUMS * vectorSize * sizeof(unsigned char));


	}

	/*
	DRRImagesize:单张DRR的尺寸
	*/

	cudaMemcpyToSymbol(DRRImageSize, OutputImageSize, 2 * sizeof(int), 0, cudaMemcpyHostToDevice);
}

void loadXRayImg(std::string filename)
{
	m_srcImg = cv::imread(filename, cv::IMREAD_GRAYSCALE);
}

void freeDICOMFromGPUMemory()
{
	cudaFree(d_object3D);
}

void freeAuxiliaryVariablesInGPUMemory(int flag)
{
	if (flag == 0)
		cudaFree(d_object2D);
	else {
		cudaFree(c_object2D);


	}
}

void calUCharDRRwithCUDA(CUDAParamerters CUDA_Parameters, DRRParameters DRR_Parameters, int PARAS_NUMS,bool isFrontal)
{
	clock_t start, end, tstart, tend;
	start = clock();
	// size error
	//for (int i = 0; i < 9; i++) {
	//	//std::cout << "stepInX:[" << i << "]= " << DRR_Parameters.stepInX[i] << std::endl;
	//	//std::cout << "stepInY:" << DRR_Parameters.stepInY[i] << std::endl;
	//	std::cout << "corner00:" << DRR_Parameters.corner00[i] << std::endl;
	//	std::cout << "SourceWorld:" << DRR_Parameters.SourceWorld[i] << std::endl;

	//}
	cudaMemcpyToSymbol(d_DRR_Parameters, DRR_Parameters.stepInX, 51 * sizeof(float), 0, cudaMemcpyHostToDevice);

	//Block 6
	int num_Threads = CUDA_Parameters.numThreads;
	int num_Blocks = CUDA_Parameters.numBlocks;

	//------------------------------------------------------------
	//Launching the threads
	int cpu_hist[256] = { 0 };

	int* gpu_hist;

	cudaMalloc((int**)&gpu_hist, 64 * PARAS_NUMS * sizeof(int));

	cudaMemcpy(gpu_hist, cpu_hist, 64 * PARAS_NUMS * sizeof(int), cudaMemcpyHostToDevice);

	drrCUDA_uc << < num_Blocks, num_Threads >> > (c_object2D, tex_object3D, gpu_hist, PARAS_NUMS);

	cudaMemcpy(cpu_hist, gpu_hist, 64 * PARAS_NUMS * sizeof(int), cudaMemcpyDeviceToHost);

	int* gpu_range;

	cudaMalloc((int**)&gpu_range, 2 * sizeof(int));


	cal_hist << <1, 1 >> > (gpu_hist, gpu_range);

	hist_uc << <num_Blocks, num_Threads >> > (c_object2D, gpu_range,PARAS_NUMS);
	//------------------------------------------------------------
	cudaFree(gpu_hist);

	cudaFree(gpu_range);

	//Copying the result from the calculations from device to host

	cudaDeviceSynchronize();
	end = clock();
	//std::cout << "time of drr is:" << (double)(end - start) / CLOCKS_PER_SEC << std::endl;

	long int vectorSize = (int)DRR_Parameters.size[0] * (int)DRR_Parameters.size[1] * PARAS_NUMS;
	
	unsigned char* h_object2D = (unsigned char*)malloc(sizeof(unsigned char) * vectorSize);

	//
	cudaMemcpy(h_object2D, c_object2D, vectorSize * sizeof(unsigned char), cudaMemcpyDeviceToHost);
	cv::Mat image(1024 * PARAS_NUMS, 1024, CV_8UC1, (void*)h_object2D);
	
	
	if (isFrontal == 1) {
		cv::imwrite("tmp.bmp", image);
	}
	else
	{
		cv::imwrite("tmp_2.bmp", image);
	}
	
	
	//cv::imshow("img", image);
	//cv::waitKey(0);
	//cv::destroyAllWindows();

	//读取参考图像
	//MainWindow* main = new MainWindow();
	//cv::Mat m_srcImg = cv::imread(main->getXRayPath(), cv::IMREAD_GRAYSCALE);

	//double* result;// = CUDA_NCC(m_srcImg, c_object2D, PARAS_NUMS);
	free(h_object2D);
	cudaFree(c_object2D);
	//return result;
}