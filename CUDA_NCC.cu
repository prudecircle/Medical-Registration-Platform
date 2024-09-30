#include "CUDA_NCC.cuh"

//ģ��ͼ
uchar* d_m_src;

//����ͼ
//uchar* d_s_src;

//����ͼ���غ�
__device__ int s_sum[256];

//����ͼ���ؾ�ֵ
__device__ double s_mean[256];

//ģ��ͼ���غ�
__device__ int m_sum[64];

//ģ��ͼ���ؾ�ֵ
__device__ double m_mean[64];

//Э����
__device__ int cov[256];

//ģ��ͼ����
__device__ int m_d[64];

//����ͼ����
__device__ int s_d[256];

//�����ϵ��
__device__ double p[256];
__device__ double p_total[4];


//ͼ���С
__device__ __constant__ int d_Imgsize[2];


__global__ void calTotalPixVal(uchar* sample, uchar* model) {
	long int idx = (blockIdx.x * blockDim.x) + threadIdx.x;
	//printf("success%d\n",threadIdx.x);
	int n = 8;
	int sub_height = (int)(d_Imgsize[0] / n);
	int sub_width = (int)(d_Imgsize[1] / n);

	//group_idx:�ֲ�ͼƬ������
	int group_idx = (idx / 131072) * 8 + (idx - (idx / 1024) * 1024) / 128;

	//������ͳ�����غ�
	if (group_idx < 64) {
		atomicAdd(&s_sum[group_idx], (int)sample[idx]);
		atomicAdd(&m_sum[group_idx], (int)model[idx]);

	}
	else
	{
		atomicAdd(&s_sum[group_idx], (int)sample[idx]);

	}


}

__global__ void calMean(int PARAS_NUMS) {
	long int idx = (blockIdx.x * blockDim.x) + threadIdx.x;
	if (idx > PARAS_NUMS * 64) return;

	//�ֿ���������ֵ
	if (idx < 64) {
		s_mean[idx] = s_sum[idx] / (d_Imgsize[0] * d_Imgsize[1] / 64); if (s_mean[idx] == 255) s_mean[idx] = 254;
		m_mean[idx] = m_sum[idx] / (d_Imgsize[0] * d_Imgsize[1] / 64); if (m_mean[idx] == 255) m_mean[idx] = 254;
	}
	else {
		s_mean[idx] = s_sum[idx] / (d_Imgsize[0] * d_Imgsize[1] / 64); if (s_mean[idx] == 255) s_mean[idx] = 254;
	}

}



__global__ void calCovDev(uchar* sample, uchar* model, int PARAS_NUMS) {
	long int idx = (blockIdx.x * blockDim.x) + threadIdx.x;
	int n = 8;
	int sub_height = (int)(d_Imgsize[0] / n);
	int sub_width = (int)(d_Imgsize[1] / n);

	//group_idx:һ��ͼƬ�оֲ�ͼ�������
	int group_idx = (idx / 131072) * 8 + (idx - (idx / 1024) * 1024) / 128;

	if (idx > d_Imgsize[0] * d_Imgsize[1] * PARAS_NUMS) return;

	if (group_idx < 64) {
		atomicAdd(&cov[group_idx], (int)(sample[idx] - s_mean[group_idx]) * (model[idx] - m_mean[group_idx]));
		atomicAdd(&m_d[group_idx], (int)(model[idx] - m_mean[group_idx]) * (model[idx] - m_mean[group_idx]));
		atomicAdd(&s_d[group_idx], (int)(sample[idx] - s_mean[group_idx]) * (sample[idx] - s_mean[group_idx]));
	}
	else
	{
		atomicAdd(&cov[group_idx], (int)(sample[idx] - s_mean[group_idx]) * (model[idx - (int)(idx / 1048576) * 1048576] - m_mean[group_idx - (group_idx / 64) * 64]));
		atomicAdd(&s_d[group_idx], (int)(sample[idx] - s_mean[group_idx]) * (sample[idx] - s_mean[group_idx]));

	}
}



__global__ void warm(int i) {
	printf("****warm****\n");
	printf("s_sum=%d\n", s_sum[0]);
}

__global__ void calNCC(int PARAS_NUMS) {

	long int idx = (blockIdx.x * blockDim.x) + threadIdx.x;
	if (idx > PARAS_NUMS * 64) return;

	//ÿ���̼߳�����������ϵ��
	p[idx] = (double)cov[idx] / sqrt((double)s_d[idx] * m_d[idx - (idx / 64) * 64]);

}




__global__ void calTotalNCC(int PARAS_NUMS) {
	long int idx = (blockIdx.x * blockDim.x) + threadIdx.x;
	//printf("p1=%f\n", p1[0]);
	if (idx > PARAS_NUMS) return;
	for (int i = 0; i < 64; i++) {
		p_total[idx] = p[idx * 64 + i];
	}

}


double* CUDA_NCC(Mat m_srcImg, uchar* d_s_src, int PARAS_NUMS) {

	int height = m_srcImg.rows;
	int width = m_srcImg.cols;

	//����ͼ��ģ��ͼUCHAR����
	uchar* m_src = m_srcImg.data;


	int n = 8;

	int sub_height = (int)(height / n);
	int sub_width = (int)(width / n);


	//ͼ���С
	int Imgsize[2] = { height,width };
	int total_memsize = height * width * sizeof(uchar);
	int sub_Imgsize = sub_height * sub_width;



	//��ʱ����
	cudaEvent_t start, end;
	cudaEventCreate(&start);
	cudaEventCreate(&end);
	cudaEventRecord(start, 0);


	//�����Դ�
	/*cudaMalloc((void**)&p_total, 4*sizeof(double));*/
	cudaMalloc((void**)&d_m_src, total_memsize);
	/*cudaMalloc((void**)&s_sum, 256* sizeof(int));
	cudaMalloc((void**)&s_mean, 256 * sizeof(double));
	cudaMalloc((void**)&m_sum, 64 * sizeof(int));
	cudaMalloc((void**)&m_mean, 64 * sizeof(double));
	cudaMalloc((void**)&cov, 256 * sizeof(int));
	cudaMalloc((void**)&s_d, 256 * sizeof(int));
	cudaMalloc((void**)&m_d, 256 * sizeof(int));
	cudaMalloc((void**)&p, 256 * sizeof(double));*/

	int int256[256] = { 0 };
	int int64[64] = { 0 };
	double double256[256] = { 0 };
	double double64[64] = { 0 };
	double double4[4] = { 0 };
	cudaMemcpyToSymbol(s_sum, int256, 256 * sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(m_sum, int64, 64 * sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(s_mean, double256, 256 * sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(m_mean, double64, 64 * sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(cov, int256, 256 * sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(s_d, int256, 256 * sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(m_d, int64, 64 * sizeof(int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(p, double256, 256 * sizeof(double), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(p_total, double4, 4 * sizeof(double), 0, cudaMemcpyHostToDevice);

	cudaMemcpyToSymbol(d_Imgsize, Imgsize, 2 * sizeof(int), 0, cudaMemcpyHostToDevice);

	cudaMemcpy(d_m_src, m_src, total_memsize, cudaMemcpyHostToDevice);
	//cudaMemcpy(d_s_src, s_src, total_memsize, cudaMemcpyHostToDevice);


	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, 0);
	int numThreads = prop.maxThreadsPerBlock;
	int numBlocks = (int)(height * width * PARAS_NUMS / numThreads);
	//cout << "numBlocks=" << numBlocks << "   numThreads=" << numThreads << endl;

	cudaDeviceSynchronize();

	//ͳ�����غ�
	calTotalPixVal << <numBlocks, numThreads >> > (d_s_src, d_m_src);


	//�������ؾ�ֵ
	calMean << <1, 64 * PARAS_NUMS >> > (PARAS_NUMS);


	//���������Э�������
	calCovDev << <numBlocks, numThreads >> > (d_s_src, d_m_src, PARAS_NUMS);


	// cudaDeviceSynchronize();
	//��������㻥���ϵ��
	calNCC << <1, 64 * PARAS_NUMS >> > (PARAS_NUMS);

	//ƽ�����������ϵ��
	calTotalNCC << <1, PARAS_NUMS >> > (PARAS_NUMS);


	//cudaDeviceSynchronize();
	//int h_p[1] = { 0 };
	//cudaMemcpy(h_cov, cov,tile_nums* sizeof(int), cudaMemcpyDeviceToHost);


	static double* result = new double[PARAS_NUMS];
	cudaMemcpyFromSymbol(result, p_total, PARAS_NUMS * sizeof(double), 0, cudaMemcpyDeviceToHost);

	cudaEventRecord(end, 0);
	cudaEventSynchronize(start);
	cudaEventSynchronize(end);
	float time;
	cudaEventElapsedTime(&time, start, end);
	//cout << "time=" << time / 1000 << endl;
	//cout << "cov=" << h_cov[0] << endl;
	//cout << "p=" << result[0] << endl;


	//cout << h_s_sum[0] << endl;
	//cout << h_m_sum[0] << endl;

	//�ͷ��Դ�

	cudaFree(d_m_src);
	/*cudaFree(s_sum);
	cudaFree(m_sum);
	cudaFree(s_mean);
	cudaFree(m_mean);
	cudaFree(cov);
	cudaFree(s_d);
	cudaFree(m_d);
	cudaFree(p);
	cudaFree(p_total);*/

	return result;
}