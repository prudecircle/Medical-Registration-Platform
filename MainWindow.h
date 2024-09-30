#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_MainWindow.h"
#include<qfuture.h>
#include<QtConcurrent>
#include"QImageLabel.h"
#include"itkDrr.h"
#include"dataProcess.h"
#include "CUDAdrr.cuh"
#include"Particle.h"
#include"EOparticle.h"
#include"CMAES.h"
#include<HalconCpp.h>
#include<Halcon.h>
#include<math.h>

//#include "vtkWidget.h"
/*itk File*/
#include "itkCenteredEuler3DTransform.h"
/*vtk File*/
#include <vector>
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
#include "vtkBorderWidget.h"
#include "vtkBorderRepresentation.h"
#include "vtkProperty2D.h"
#include "vtkExtractVOI.h"
#include "vtkInteractorStyleSwitch.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkInteractorStyleTrackballActor.h"


//Ϊ���ܹ�ʹ��vtk��Ӵ���
#include "vtkAutoInit.h"
#define PI  3.1415926 

#define  e_0 1   //���ȣ�δʹ�ã�
#define CNT 50  //��������
#define UP_V 5//�ٶ����ֵ
#define DOWN_V -5//�ٶ���Сֵ
#define MIN -1//���ƶ���Сֵ

#define w_min 0.6//Ȩ����С
#define w_max 1.2//Ȩ�����
#define RANGE 20 //��������Ѱ�ҷ�Χ
#define VELO_MAX 5   //�����ٶ�
#define VELO_MIN -5
#define PARTICLE 10  //���Ӹ���
#define DIMEN 6   //����ά��

using namespace std;
using namespace cv;
using namespace Halcon;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    using InputImageType = typename itk::Image<float, 3>;
    using PointType = typename InputImageType::PointType;

    MainWindow(QWidget *parent = Q_NULLPTR);
    void setXRayPath(std::string path_1);
    void setTransMat(double matData[][4]);
    cv::Mat getTransMat();
    std::string getXRayPath();
    int getAlgorithm();
    //double* cudaUCharDRR(std::vector<float> rotation, std::vector<float>translation, int dx, int dy, float threshold, int PARAS_NUMS);
    void setParam(std::vector<double> parameters);             //����drr���еĲ���
    void drawDRR();
    void drawlateralDRR();
    void drawDRRs(std::vector<double> parameters);
    void drawDRRsLateral(std::vector<double> parameters);
    void viewImage(int flag);
    float calSimilarity(int flag);
    float calSimilarityLateral(int flag);
    int PARA_NUMS = 1;
    double* similarity;
    double getSimilarity_NCC(cv::Mat drr_image, cv::Mat xray_image);
    double NCC(cv::Mat drr, cv::Mat xray);
    double NCC_whole(cv::Mat drrImg, cv::Mat xrayImg);
    double contourNCC(vector<Mat> drr, vector<Mat>xray);
    void doPSO(int pn, int d);
    void initial(int pn, int d);
    void createNew();
    bool isExited(int n);
    void CMA_ES();
    double GetMutualInfo(cv::Mat ref, cv::Mat flt);
    double SSD(cv::Mat ref, cv::Mat flt);
    double SiameseSSMI(cv::Mat ref, cv::Mat flt);
    double FaceNet_PSNR(cv::Mat& ref, cv::Mat& flt);
    
    void EO();

    void DE();
    vector<double> generateRandomVector(int size, double lower_bound, double upper_bound);
    vector<vector<double>> initializePopulation(int population_size, int dimension, vector<double>& lower_bounds, vector<double>& upper_bounds);
    vector<int> argsort(vector<double>& fitness);
    vector<vector<double>> sortPopulationByFitness(vector<vector<double>>& population, vector<double>& fitness);
    vector<double> calculateFitness(vector<vector<double>>& population);
    vector<vector<double>> differentialEvolution(int population_size, int dimension, vector<double>& lower_bounds, vector<double>& upper_bounds, int max_iterations);

    int sign(double x);

    Mat ycCanny(cv::Mat input);
    Mat non_maximum_suppression_Inter(Mat dx, Mat dy);
    void trace(Mat edgeMag_nonMaxSup, Mat& edge, float lowerThresh, int r, int c, int rows, int cols);
    Mat hysteresisThreshold(Mat edgeMag_nonMaxSup, float lowerThresh, float upperThresh);
    bool checkInRange(int r, int c, int rows, int cols);

    //-------Santan���ƶȺ���
    float Similarity;
    float similarity_count(vector<Mat>ceil_img, vector<Mat>ceil_img2);
    float similarity_count1(vector<Mat>ceil_img, vector<Mat>ceil_img2);
    double getSimilarity(cv::Mat drr_image, cv::Mat xray_image);
    Mat binary_image(Mat image);
    cv::Mat blackegde_binary(Mat image);
    cv::Mat preprocessing(Mat image);
    vector<cv::Mat> ceilImgGet(Mat image);
    bool isInParRect(int x1, int y1, int x4, int y4, int x, int y);
    bool isExistEdgePtInRect(int xleft, int xright, int yup, int ydown, Mat image);
    void calMaxOfZaxis(vector<double> para);


    struct nccStruct                                          //NCC���ýṹ��
    {
        double totalPixel =0;
        double meanPixel =0;
        int variance=0;
    };

public slots:
    void openCT();
    void openXRay();
    void registration();
    void setAlgorithm();
    void updateParam();                                 //�����������²���
    void batch();

private:
    cv::Mat transMat;
    cv::Vec3f angle;
    cv::Vec3f translation;
    QImageLabel* drr_label;
    QImageLabel* xray_label;
    QImageLabel* drr_label_2;
    QImageLabel* xray_label_2;
    bool isDRR;											//�Ƿ��Ѿ�����DRR
    bool isXray;										//�Ƿ���Xrayͼ��
    bool isCT;											//�Ƿ���CTͼ��
    Ui::MainWindowClass ui;
    int dcmFlag = 0;
    //vtkWidget* vtk3dwidget;
    void clearLayer();
    void open_CT(std::string input);
    itkDrr* drr;
    DataProcess* dataProcess;                           //vtk������
    std::string xRayPath_1;
    std::string xRayPath_2;
    int method = 0;                                     //�Ż��㷨,Ĭ��PSO
    std::vector<double> parameters;

    //---------------------PSO-------------------
    Particle* particle;             // ����
    int pNum;                       // ��������
    int dim;                        // ά��
    double** p;                     // �ֲ����� ĳһ�����������̵�����λ��
    int bestIndex;                  // ȫ�����Ŷ�Ӧ�����ӵ�����ֵ
    double fg;                      // ȫ������
    double lastFG;                  // ǰһ�������̵�ȫ������
    double c1;						//ѧϰ����
    double c2;
    double PSO_w;					//Ȩ��
    vector<double> curve;

    EOparticle* C;
};

