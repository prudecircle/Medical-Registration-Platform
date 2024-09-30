#include "MainWindow.h"
#include<qpushbutton.h>
#include<qdebug.h>
#include<qfiledialog.h>
#include<qmessagebox.h>


MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
    ui.setupUi(this);
   // QPushButton* btn_openCT = new QPushButton;
    
    //btn_openCT->setParent(this);
    //btn_openCT->show();
	drr_label = new QImageLabel(ui.imgWidget);
	xray_label = new QImageLabel(ui.imgWidget);
	drr_label_2 = new QImageLabel(ui.imgWidget);
	xray_label_2 = new QImageLabel(ui.imgWidget);

	drr_label->setGeometry(20, 80, 350, 350);
	xray_label->setGeometry(400, 80, 350, 350);
	drr_label_2->setGeometry(0, 300, 250, 250);
	xray_label_2->setGeometry(440, 300, 250, 250);

	ui.progressBar->setHidden(true);
	ui.mainToolBar->setHidden(true);
	drr_label->setVisible(false);
	xray_label->setVisible(false);
	drr_label_2->setVisible(false);
	xray_label_2->setVisible(false);
	//ui.batch->setVisible(false);
	//ui.btn_trim->setVisible(false);
	

	
	isCT = false;
	isDRR = false;
	isXray = false;
    connect(ui.btn_openCT, SIGNAL(clicked()), this, SLOT(openCT()));                               //打开CT文件
	connect(ui.btn_openX, SIGNAL(clicked()), this, SLOT(openXRay()));                              //打开X光文件
	connect(ui.btn_registration, SIGNAL(clicked()), this, SLOT(registration()));                   //开始配准
	connect(ui.box_algorithm, SIGNAL(currentIndexChanged(QString)), this, SLOT(setAlgorithm()));   //选择优化算法
	//connect(ui.batch, SIGNAL(clicked()), this, SLOT(batch()));
	//参数更新区域
	connect(ui.rx, SIGNAL(returnPressed()), this, SLOT(updateParam()));
	connect(ui.ry, SIGNAL(returnPressed()), this, SLOT(updateParam()));
	connect(ui.rz, SIGNAL(returnPressed()), this, SLOT(updateParam()));
	connect(ui.tx, SIGNAL(returnPressed()), this, SLOT(updateParam()));
	connect(ui.ty, SIGNAL(returnPressed()), this, SLOT(updateParam()));
	connect(ui.tz, SIGNAL(returnPressed()), this, SLOT(updateParam()));
}

void MainWindow::openXRay() {
	QFileDialog* fileDialog = new QFileDialog(this);
	fileDialog->setWindowTitle(QStringLiteral("选择X-Ray图像"));
	fileDialog->setDirectory("");
	fileDialog->setNameFilter(tr("Images (*.bmp)"));
	fileDialog->setFileMode(QFileDialog::ExistingFiles);
	fileDialog->setViewMode(QFileDialog::Detail);

	//QString filename = QString("F:\\3\\X\\90,0,0,0,0,-100.bmp");
	//QString filename = QString("E:\\double-location_data\\spine-2\\17_tt\\x\\00003841.bmp");
	QString filename = QString("G:\\CT+X_ray\\1\\20200401\\X - Copy\\20200401201217_14.bmp");
	
	setXRayPath(string((const char*)filename.toLocal8Bit()));
	
	
		
	
		loadXRayImg(string((const char*)filename.toLocal8Bit()));
	//m_srcImg = cv::imread(filename, cv::IMREAD_GRAYSCALE);
		viewImage(2);
		

		xray_label->setVisible(true);
		xray_label->setPixmap(QPixmap::fromImage(QImage(filename)));
		xray_label->show();
		//xray_label_2->setVisible(true);
		//xray_label_2->setPixmap(QPixmap::fromImage(QImage(filename_2)));
		//xray_label_2->show();
		QMessageBox::information(this, "X Ray read complete", "X Ray read complete!");
		qDebug() << QString::fromStdString(getXRayPath()) << endl;
	//}
	//QGraphicsOpacityEffect* eff = new QGraphicsOpacityEffect;
	
	//eff->setOpacity(1.0 - getOpacity());
	//ui.xray_label->setGraphicsEffect(eff);
	
}
void MainWindow::openCT() {
	dcmFlag = 0;
	ui.progressBar->setHidden(false);
	QFileDialog* fileDialog = new QFileDialog(this);
	fileDialog->setWindowTitle(QStringLiteral("选择CT图像序列"));
	fileDialog->setDirectory("");
	fileDialog->setFileMode(QFileDialog::Directory);
	fileDialog->setViewMode(QFileDialog::Detail);
	QString filedic;
	if (fileDialog->exec()) {
		filedic = fileDialog->selectedFiles()[0];
	}
	//判断导入是否为CT文件夹
	//QString filepath = QFileDialog::getOpenFileName();
	QDir filedir(filedic);
	QFileInfoList filelist = filedir.entryInfoList();//获取文件信息列表
	//qDebug() << filedic << endl;
	for (int i = 0; i < (filelist.size()); i++) {
		QFileInfo fileinfo = filelist.at(i);
		QString filesuffix = fileinfo.suffix();
		QString filesuffix_upper = filesuffix.toUpper();
		//qDebug() << fileinfo << endl;
		//if (strcmp(filesuffix_upper.toStdString().c_str(), "DCM") == 0) {
		dcmFlag = 1;
		//break;
	//}
	//else {
		//dcmFlag = 0;
	}
	if (dcmFlag == 0) {
			//QMessageBox::warning(NULL, "Attention", "Please import the correct CT folder", QMessageBox::Ok, QMessageBox::Ok);
			QMessageBox::warning(this, QString::fromLocal8Bit("错误"), QString::fromLocal8Bit("请正确导入CT文件"));
			return;
		}
	else {
		
		clearLayer();		
		//open_CT(filedic.toStdString());
		open_CT("G:\\CT+X_ray\\resource\\jingu\\20201204robot\\ct\\ct0001");
		std::vector<double> parameters = {90, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
		setParam(parameters); 
		drawDRR();
		drawlateralDRR();
		
		viewImage(0);
		
		
		drr_label->setVisible(true);
		drr_label->setPixmap(QPixmap::fromImage(QImage("tmp.bmp")));
		drr_label->show();
		//drr_label_2->setVisible(false);
		//drr_label_2->setPixmap(QPixmap::fromImage(QImage("tmp_2.bmp")));
		//drr_label_2->show();
		ui.progressBar->setHidden(true);
		QMessageBox::information(this, "CT read complete", "CT read complete!");


		//qDebug() << "read file complete" << endl;
		//open_CT(filedic.toStdString().data());
		//QFuture<void> vtk3d = QtConcurrent::run(this, &MainApplication::vtkImitateDRR, filedic.toStdString());
		//while (!future.isFinished()) {
		//	QApplication::processEvents();
		//}
		//ui.qvtkWidget->setVisible(true);
		//viewImage(0);
		////qwait->close();
		//ui.waitBar->setHidden(true);
		//QApplication::restoreOverrideCursor();
	}
	
}

void MainWindow::clearLayer() {
	//layer = 0;
	//imageBlend = vtkSmartPointer<vtkImageBlend>::New();
	//dataProcess->clearXRayReader();
	isCT = false;
	isDRR = false;
	//isXray = false;
	//lock3d = false;
}

void MainWindow::open_CT(std::string input) {
	if (drr != NULL) {
		// 析构 释放GPU内存
		//free(rd);
		drr = NULL;
	}
	drr = new itkDrr(input);
	isCT = true;
}

void MainWindow::setXRayPath(std::string path_1) {
	xRayPath_1 = path_1;
}

void MainWindow::setTransMat(double matData[][4])
{
	transMat = cv::Mat(cv::Size(4, 4), CV_64F, matData);
	cv::Mat rotationMat(getTransMat()(cv::Rect(0, 0, 3, 3)));
	for (int i = 0; i < 3; i++) {
		translation[i] = getTransMat().at<double>(i, 3);
	}
	cv::Vec3f a = dataProcess->rotationMatrixToEulerAngles(rotationMat);
	angle = a * 180.0f / PI;
}

cv::Mat MainWindow::getTransMat()
{
	return transMat;
}

std::string MainWindow::getXRayPath() {
	return xRayPath_1;
}

int MainWindow::getAlgorithm()
{
	return method;
}

void MainWindow::registration() {
	int method = getAlgorithm();
	//qDebug() <<"current method:" << method << endl;
	switch (method) {
	case 0:
		qDebug() << "PSO" << endl;
		doPSO(PARTICLE, DIMEN);
		drawDRR();
		//drawlateralDRR();
		viewImage(0);
		drr_label->setPixmap(QPixmap::fromImage(QImage("tmp.bmp")));
		drr_label->show();
		/*drr_label_2->setPixmap(QPixmap::fromImage(QImage("tmp_2.bmp")));
		drr_label_2->show();*/
		break;

	case 1: {
		qDebug() << "EO" << endl;
		clock_t start = clock();
		EO();
		cout << "time of registration is" << (double)(clock() - start) / CLOCKS_PER_SEC << std::endl;
		viewImage(0);
		drr_label->setPixmap(QPixmap::fromImage(QImage("tmp.bmp")));
		drr_label->show();
		break;
		}
		
	case 2:{
		qDebug() << "DE" << endl;
		clock_t start = clock();
		DE();
		std::cout << "time of registration is" << (double)(clock() - start) / CLOCKS_PER_SEC << std::endl;
		viewImage(0);
		drr_label->setPixmap(QPixmap::fromImage(QImage("tmp.bmp")));
		drr_label->show();
		break;
		}
	}

}

void MainWindow::setAlgorithm()
{
	method = ui.box_algorithm->currentIndex();
}

void MainWindow::updateParam()
{
	if (isCT) {
		//qDebug() << "chufa" << endl;
		parameters.clear();
		parameters.push_back(ui.rx->text().toDouble());
		parameters.push_back(ui.ry->text().toDouble());
		parameters.push_back(ui.rz->text().toDouble());
		parameters.push_back(ui.tx->text().toDouble());
		parameters.push_back(ui.ty->text().toDouble());
		parameters.push_back(ui.tz->text().toDouble());
		setParam(parameters);
		/*for (auto iter = parameters.begin(); iter != parameters.end(); iter++) {
			qDebug() << *iter << "    ";
		}*/	
		drawDRR();
		//float similarity = calSimilarity(1);
		

		//drawlateralDRR();
		drr_label->setVisible(true);
		drr_label->setPixmap(QPixmap::fromImage(QImage("tmp.bmp")));
		drr_label->show();
		//drr_label_2->setPixmap(QPixmap::fromImage(QImage("tmp_2.bmp")));
		//drr_label_2->show();
		if (isDRR && isXray) {
			float similarity[1] = {calSimilarity(2)};
			ui.label_similarity->setText(QString::number(similarity[0]));
			//ui.label_similarity_2->setText(QString::number(similarity[1]));
			
		}
		
	}
}

void MainWindow::setParam(std::vector<double> parameters) {
	drr->setParameters(parameters);
}

void MainWindow::drawDRR()
{
	 drr->cudaUCharDRR();
	 isDRR = true;
}
void MainWindow::drawlateralDRR()
{
	drr->cudaUCharDRRLateral();
}

void MainWindow::viewImage(int flag)
{
	//此方法中的所有步骤都采用多线程来避免UI卡顿
	if (flag == 0) { //调用ITKDRR中的DRR算法
		// 导入CT图像之后使用DRR算法生成DRR图像
		//
		//isCT = true;
		//cur_DRR = -1;
		//cur_Xray = -1;
		//isDRR = true;
		//adjustPara();
	}
	else if (flag == 1) { // 展示DRR图像
		return;
	}
	else if (flag == 2) { // 展示XRay图像
		isXray = true;
	}
	
	
	if (isCT && isDRR && isXray) {
		float similarity[1] = { calSimilarity(2)};
		cout << "similarity=" << similarity[0] << std::endl;
		ui.label_similarity->setText(QString::number(similarity[0]));
	}

}

float MainWindow::calSimilarity(int flag)
{
	 static float similarity =0;
	
		cv::Mat xray_img = cv::imread(xRayPath_1, 1);
		cv::Mat drr_img = cv::imread("tmp.bmp", 1);
		if (flag == 1) {                                //NCC计算相似度
			//getSimilarity_NCC(drr_img, xray_img);
			similarity= NCC(drr_img, xray_img);
		}
		if (flag == 2) {
			similarity = getSimilarity(drr_img, xray_img);
		}
		if (flag == 3) {
			similarity= GetMutualInfo(drr_img, xray_img);
		}
		if (flag == 4) {
			similarity = getSimilarity_NCC(drr_img, xray_img);
		}
		if (flag == 5) {
			similarity = SSD(drr_img, xray_img);
		}
		if (flag == 6) {
			similarity = SiameseSSMI(drr_img, xray_img);
		}
		if (flag == 7) {
			similarity = FaceNet_PSNR(drr_img, xray_img);
		}
	cout  << similarity << std::endl;
	return similarity;
	
}
float MainWindow::calSimilarityLateral(int flag)
{
	static float similarity = 0;

	cv::Mat xray_img = cv::imread(xRayPath_2, 1);
	cv::Mat drr_img = cv::imread("tmp_2.bmp", 1);
	if (flag == 1) {                                //NCC计算相似度
		//getSimilarity_NCC(drr_img, xray_img);
		similarity = NCC(drr_img, xray_img);
	}
	if (flag == 2) {
		similarity= getSimilarity(drr_img, xray_img);
	}
	if (flag == 3) {
		similarity = GetMutualInfo(drr_img, xray_img);
	}
	if (flag == 4) {
		similarity = getSimilarity_NCC(drr_img, xray_img);
	}
	cout << similarity << std::endl;
	return similarity;

}

double MainWindow::NCC(cv::Mat drrImg, cv::Mat xrayImg)
{
	int n = 8;
	if (drrImg.size != xrayImg.size) drrImg.resize((xrayImg.cols, xrayImg.rows));
	int height = drrImg.rows;
	int width = drrImg.cols;
	int gridheight = ceil(height / 8);
	int gridwidth = ceil(width / 8);
	double p_total = 0.0; 
	double pr[64]={0};
	int m = 0;
	long double pSum = 0;
	for (int gridrows = 0; gridrows < n; gridrows++) {
		for (int gridcols = 0; gridcols < n; gridcols++) {
			nccStruct drr, x;
			int covariance = 0;
			cv::Mat tempDrrGrid = drrImg(cv::Rect(gridcols * gridwidth, gridrows * gridheight, gridwidth, gridheight));
			cv::Mat tempXGrid = xrayImg(cv::Rect(gridcols * gridwidth, gridrows * gridheight, gridwidth, gridheight));

			for (int i = 0; i < gridwidth; i++) {
				for (int j = 0; j < gridheight; j++) {
					drr.totalPixel += double(tempDrrGrid.at<uchar>(j, i));
					x.totalPixel += double(tempXGrid.at<uchar>(j, i));

				}
			}
			drr.meanPixel = drr.totalPixel / (gridheight * gridwidth);
			if (drr.meanPixel == 255) drr.meanPixel = (double)254;
			x.meanPixel = x.totalPixel / (gridheight * gridwidth);
			if (x.meanPixel == 255) x.meanPixel = 254;
			//cout << "drr平均值  " << drr.meanPixel << std::endl;
			//cout << "x平均值  " << x.meanPixel << std::endl;
			for (int i = 0; i < gridwidth; i++) {
				for (int j = 0; j < gridheight; j++) {
					drr.variance += (int)pow(double(tempDrrGrid.at<uchar>(j, i))-drr.meanPixel, 2);
					x.variance += (int)pow(double(tempXGrid.at<uchar>(j, i))-x.meanPixel, 2);
					covariance +=(int)(tempDrrGrid.at<uchar>(j, i)-drr.meanPixel) * (tempXGrid.at<uchar>(j, i)-x.meanPixel);
				}
			}

			//cout << "cov:" << covariance << std::endl;
			//cout << "drr variance" << drr.variance << std::endl;
			//cout << "x variance" << x.variance << std::endl;
			pr[m] = (double)covariance / sqrt((double)drr.variance * x.variance);
			//cout << "sqrt" << sqrt((double)drr.variance * x.variance) << std::endl;
			//cout << "results" << (double)covariance / sqrt((double)drr.variance * x.variance) << std::endl;
			//cv::imshow("drr", tempDrrGrid);
			//cv::waitKey(0);
			//cv::destroyAllWindows();
			m++;
			//cout <<"区域相似度:  "<< pr[m] << std::endl;
			pSum += (double)covariance / sqrt((double)drr.variance * x.variance);
		}
	}
	m = 0;

	//cout << "sum" << pSum << std::endl;
	p_total = pSum / (n * n);
	return p_total;
}
double MainWindow::NCC_whole(cv::Mat drrImg, cv::Mat xrayImg)
{
	int n = 8;
	if (drrImg.size != xrayImg.size) drrImg.resize((xrayImg.cols, xrayImg.rows));
	int height = drrImg.rows;
	int width = drrImg.cols;
	int gridheight = height;
	int gridwidth = width;
	double p_total = 0.0;
	//double pr[64] = { 0 };
	//int m = 0;
	long double pSum = 0;
	
			nccStruct drr, x;
			int covariance = 0;
			//cv::Mat tempDrrGrid = drrImg(cv::Rect(gridcols * gridwidth, gridrows * gridheight, gridwidth, gridheight));
			//cv::Mat tempXGrid = xrayImg(cv::Rect(gridcols * gridwidth, gridrows * gridheight, gridwidth, gridheight));

			for (int i = 0; i < gridwidth; i++) {
				for (int j = 0; j < gridheight; j++) {
					drr.totalPixel += double(drrImg.at<uchar>(i, j));
					x.totalPixel += double(xrayImg.at<uchar>(i, j));

				}
			}
			drr.meanPixel = drr.totalPixel / (gridheight * gridwidth);
			if (drr.meanPixel == 255) drr.meanPixel = (double)254;
			x.meanPixel = x.totalPixel / (gridheight * gridwidth);
			if (x.meanPixel == 255) x.meanPixel = 254;
			cout << "drr平均值  " << drr.meanPixel << std::endl;
			cout << "x平均值  " << x.meanPixel << std::endl;
			for (int i = 0; i < gridwidth; i++) {
				for (int j = 0; j < gridheight; j++) {
					drr.variance += (int)pow(double(drrImg.at<uchar>(i, j)) - drr.meanPixel, 2);
					x.variance += (int)pow(double(xrayImg.at<uchar>(i, j)) - x.meanPixel, 2);
					covariance += (int)((drrImg.at<uchar>(i, j) - drr.meanPixel) * (xrayImg.at<uchar>(i, j) - x.meanPixel));
				}
			}

			//cout << "cov:" << covariance << std::endl;
			//cout << "drr variance" << drr.variance << std::endl;
			//cout << "x variance" << x.variance << std::endl;
			//pr[m] = (double)covariance / sqrt((double)drr.variance * x.variance);
			//cout << "sqrt" << sqrt((double)drr.variance * x.variance) << std::endl;
			//cout << "results" << (double)covariance / sqrt((double)drr.variance * x.variance) << std::endl;
			//cv::imshow("drr", tempDrrGrid);
			//cv::waitKey(0);
			//cv::destroyAllWindows();
			//m++;
			//cout <<"区域相似度:  "<< pr[m] << std::endl;
			pSum = (double)covariance / sqrt(abs((double)drr.variance * x.variance));
		
	
	//m = 0;
			if (isnan(pSum) == true) pSum = 0;
	//cout << "sum" << pSum << std::endl;
	//p_total = pSum / (n * n);
	return pSum;
}

void MainWindow::doPSO(int pn, int d)
{
	clock_t start = clock();
	int n = 0;
	initial(pn, d);
	do
	{
		printf("\n\n-------第%d次迭代:------- \n", n+1);
		createNew();
		n++;
	} while (!isExited(n));
	vector<double> para;
	for (int i = 0; i < d; i++) {
	//	cout << p[bestIndex][i] << std::endl;
		para.push_back(p[bestIndex][i]);
	}
	cout << "time of registration is" << (double)(clock() - start) / CLOCKS_PER_SEC << std::endl;
	//calMaxOfZaxis(para);

	//free(particle);
	
	/*for (auto iter = curve.begin(); iter != curve.end(); iter++) {
		cout << *iter << std::endl;
	}*/
	curve.clear();
	//std::cout << p[bestIndex] << std::endl;
	//print(p[bestIndex]);
	//std::cout << "FG:  " << fg << std::endl;
}

void MainWindow::initial(int pn, int d) {
	c1 = 2;
	c2 = 2;
	fg = MIN;
	bestIndex = 0;
	pNum = pn;
	dim = d;
	particle = NULL;
	int i, j, k;
	srand((unsigned)time(NULL));//时间戳
	/*para.push_back(ui.xp_label->text().toDouble());
	para.push_back(ui.yp_label->text().toDouble());
	para.push_back(ui.zp_label->text().toDouble());
	para.push_back(ui.xr_label->text().toDouble());
	para.push_back(ui.yr_label->text().toDouble());
	para.push_back(ui.zr_label->text().toDouble());*/
	parameters.clear();
	parameters.push_back(ui.rx->text().toDouble());
	parameters.push_back(ui.ry->text().toDouble());
	parameters.push_back(ui.rz->text().toDouble());
	parameters.push_back(ui.tx->text().toDouble());
	parameters.push_back(ui.ty->text().toDouble());
	parameters.push_back(ui.tz->text().toDouble());


	//速度随机范围
	if (particle == NULL) {
		particle = new Particle[pn];         // 创建粒子数组实例对象
		for (i = 0; i < pn; ++i) {
			particle[i].setDim(d);           // 设置粒子的维数
			particle[i].position = new double[d];
			particle[i].velocity = new double[d];
			//随机初始粒子位置
			for (j = 0; j < d; j++) {
				particle[i].position[j] = rand() / (double)RAND_MAX * ((parameters[j] + 0.5 * RANGE) - (parameters[j] - 0.5 * RANGE)) + (parameters[j] - 0.5 * RANGE);
			}
			for (j = 0; j < d; j++) {
				printf("P-[%d].[%d]:%f  ", i, j, particle[i].position[j]);
			}
			printf("\n");
			for (j = 0; j < d; j++) {
				particle[i].velocity[j] = rand() / (double)RAND_MAX * (VELO_MAX - VELO_MIN) + VELO_MIN;
				printf("V-[%d].[%d]:%f   ", i, j, particle[i].velocity[j]);
			}
			printf("\n");

		}
	}
	printf("\n-------------------------------------\n");
	p = new double* [pn];                      // 粒子数
	for (i = 0; i < pn; ++i) {
		p[i] = new double[d];
		for (j = 0; j < dim; ++j) {
			p[i][j] = particle[i].position[j];  // 初始化局部最优位置
		}
			 
		
		vector<double> para;
		para.push_back(particle[i].position[0]);
		para.push_back(particle[i].position[1]);
		para.push_back(particle[i].position[2]);
		para.push_back(particle[i].position[3]);
		para.push_back(particle[i].position[4]);
		para.push_back(particle[i].position[5]);
		/*double xp = particle[i].position[0];
		double yp = particle[i].position[1];
		double zp = particle[i].position[2];
		double xr = particle[i].position[3];
		double yr = particle[i].position[4];
		double zr = particle[i].position[5];*/
		
		drawDRRs(para);
		
		float similarity = calSimilarity(2);
		particle[i].fit = similarity;
		printf("第%d个粒子的相似度: %f\n", i + 1, particle[i].fit);
		
		// 获取全局最大
		if (particle[i].fit > fg)
		{
			fg = particle[i].fit;
			bestIndex = i;                   // 改变全局最优的索引位置
		}
	}
	printf("初始化最优位置为%d， 相似度值为：fg: %f\n\n", bestIndex, fg);
	
}

void MainWindow::createNew() {
	int i, j, k;
	srand((unsigned)time(NULL));//时间戳
	// 更新速度
	for (int k = 0; k <= CNT; k++)
	{
		PSO_w = w_max - k * (w_max - w_min) / CNT;  //更新权重
		for (i = 0; i < pNum; i++)
		{
			for (j = 0; j < particle[i].dim; j++)
			{
				particle[i].velocity[j] = particle[i].velocity[j] * PSO_w +
					(rand() / (double)RAND_MAX) * c1 * (p[i][j] - particle[i].position[j]) +
					(rand() / (double)RAND_MAX) * c2 * (p[bestIndex][j] - particle[i].position[j]);

				if (particle[i].velocity[j] > VELO_MAX)
					particle[i].velocity[j] = VELO_MAX;
				else if (particle[i].velocity[j] < VELO_MIN)
					particle[i].velocity[j] = VELO_MIN;
			}
		}
	}
	// 更新位置
	for (i = 0; i < pNum; i++)
	{
		for (j = 0; j < particle[i].dim; j++)
		{
			particle[i].position[j] += particle[i].velocity[j];
			//范围限定检查
			if (particle[i].position[j] < (parameters[j] - 0.5 * RANGE))
				particle[i].position[j] = parameters[j] - 0.5 * RANGE;
			else if (particle[i].position[j] > (parameters[j] + 0.5 * RANGE))
				particle[i].position[j] = parameters[j] + 0.5 * RANGE;
		}
		// 计算该位置的函数值
		/*double xp = particle[i].position[0];
		double yp = particle[i].position[1];
		double zp = particle[i].position[2];
		double xr = particle[i].position[3];
		double yr = particle[i].position[4];
		double zr = particle[i].position[5];
		double tf = reg(xp, yp, zp, xr, yr, zr);*/
		vector<double> para;
		para.push_back(particle[i].position[0]);
		para.push_back(particle[i].position[1]);
		para.push_back(particle[i].position[2]);
		para.push_back(particle[i].position[3]);
		para.push_back(particle[i].position[4]);
		para.push_back(particle[i].position[5]);
		
		drawDRRs(para);
	
		float similarity = calSimilarity(2);
		double tf = similarity;
		if (i == 0) {
			curve.push_back(1. / tf);
		}
		
		cout << "similarity: " << tf << std::endl;
		// 检查是否局部最优
		if (tf > particle[i].fit)
		{
			particle[i].fit = tf;
			for (k = 0; k < particle[i].dim; k++)
				p[i][k] = particle[i].position[k];
			// 检查是否全局最优
			if (tf > fg)
			{
				lastFG = fg;
				fg = tf;
				bestIndex = i;
			}
		}
	}
}
bool MainWindow::isExited(int n)
{
	//迭代结束条件
	if (n >= CNT)
		return true;
	return false;
}

double MainWindow::getSimilarity(cv::Mat drr_image, cv::Mat xray_image)
{
	//std::stringstream ss(std::stringstream::in | std::stringstream::out);
	int Width = 1024, Height = 1024;
	cv::resize(drr_image, drr_image, Size(Width, Height));
	//imshow("drr", drr_image);
	vector<Mat>ceil_image = ceilImgGet(drr_image);
	//std::stringstream ss2(std::stringstream::in | std::stringstream::out);
	//imshow("x", xray_image);
	//std::string filename2 = ss2.str();
	//cv::Mat xray_image = imread(filename2, 1);//ss.str()
	if (xray_image.size() != drr_image.size()) {
		cv::resize(xray_image, xray_image, Size(Width, Height));
	}
	vector<Mat>ceil_image2 = ceilImgGet(xray_image);
	float a = similarity_count1(ceil_image, ceil_image2);
	//double a = contourNCC(ceil_image, ceil_image2);
	//waitKey(0);
	//destroyAllWindows();

	return a;
}
//double MainWindow::getSimilarity(cv::Mat drr_image, cv::Mat xray_image) {
//	int Width = 1024, Height = 1024;
//	cv::resize(drr_image, drr_image, Size(Width, Height));
//	cv::Mat output = ycCanny(drr_image); //边缘提取
//	cvtColor(xray_image, xray_image, COLOR_BGR2GRAY);	//转为灰度图
//	vector<Mat>ceil_image = { output };
//	if (xray_image.size() != drr_image.size()) {
//		cv::resize(xray_image, xray_image, Size(Width, Height));
//	}
//	vector<Mat>ceil_image2 = { xray_image };
//	float a = similarity_count(ceil_image, ceil_image2);
//	return a;
//}

//边缘检测算法
Mat MainWindow::ycCanny(cv::Mat input) {
	//输入图像

	//cv::imshow("原图",input);
	imwrite("ori.bmp", input);
	Mat image = imread(".//ori.bmp", IMREAD_GRAYSCALE);
	//cv::imshow("ori",image);
	Mat dx, dy, Gussimage;
	GaussianBlur(image, Gussimage, Size(5, 5), 0, 0);
	//cv::imwrite(" Gussimage.jpg", Gussimage);
	//imshow("高斯平滑后", Gussimage);
	cv::Sobel(Gussimage, dx, CV_32FC1, 1, 0, 3);//ksize=cv.FILTER_SCHARR
	cv::Sobel(Gussimage, dy, CV_32FC1, 0, 1, 3);
	//边缘强度的灰度级显示
	Mat edgeMag;
	cv::magnitude(dx, dy, edgeMag);
	Mat mag;
	edgeMag.convertTo(mag, CV_8UC1);
	//cv::imshow("边缘强度的灰度级显示", mag);
	//cv::imwrite("mag.jpg", mag);
	//第二步：非极大值抑制
	Mat edgeMag_nonMaxSup = non_maximum_suppression_Inter(dx, dy);
	//显示非极大值抑制后的边缘强度图
	Mat nonMaxSup;
	edgeMag_nonMaxSup.convertTo(nonMaxSup, CV_8UC1);//convertTo()函数负责转换数据类型不同的Mat，即可以将类似float型的Mat转换到imwrite()函数能够接受的类型。
	//cv::imshow("非极大值抑制", nonMaxSup);
	//cv::waitKey(0);
	//cv::imwrite("nonMaxSup.jpg", nonMaxSup);
	//第三步：双阈值的滞后阈值处理 //推荐的高与低阈值比在2：1到3：1之间、高阈值去除噪声，低阈值保留边缘信息 
	//24 48 针
	float lowerThresh = 20;//gugu 1:15 21  2:17 25  3:17 25  4:14 25 5:14 25   gupeng:1 : 14 30  2:17 27  3:12 25 
	float upperThresh = 30;//jinggu 1:10 20
	Mat edge = hysteresisThreshold(edgeMag_nonMaxSup, lowerThresh, upperThresh);//edgeMag_nonMaxSup
	//imshow("CannyEdge", edge);
	//cv::imwrite("CannyEdge.jpg", edge);

	waitKey(0);
	return edge;
}

Mat MainWindow::non_maximum_suppression_Inter(Mat dx, Mat dy)
{
	//使用平方和开方的方式计算边缘强度
	Mat edgeMag;
	cv::magnitude(dx, dy, edgeMag);
	//宽高
	int rows = dx.rows;
	int cols = dy.cols;
	Mat edgeMag_nonMaxSup = Mat::zeros(dx.size(), dx.type());
	for (int r = 1; r < rows - 1; r++)
	{
		for (int c = 1; c < cols - 1; c++)
		{
			float x = dx.at<float>(r, c);
			float y = dy.at<float>(r, c);
			if (x == 0 && y == 0)
				continue;
			float angle = atan2f(y, x) / CV_PI * 180;
			//领域内八个方向上的边缘强度
			float leftTop = edgeMag.at<float>(r - 1, c - 1);
			float top = edgeMag.at<float>(r - 1, c);
			float rightBottom = edgeMag.at<float>(r + 1, c + 1);
			float right = edgeMag.at<float>(r, c + 1);
			float rightTop = edgeMag.at<float>(r - 1, c + 1);
			float leftBottom = edgeMag.at<float>(r + 1, c - 1);
			float bottom = edgeMag.at<float>(r + 1, c);
			float left = edgeMag.at<float>(r, c - 1);
			float mag = edgeMag.at<float>(r, c);
			//左上方与上方的插值 右下方和下方的插值
			if ((angle > 45 && angle <= 90) || (angle > -135 && angle <= -90))
			{
				float ratio = x / y;
				float top = edgeMag.at<float>(r - 1, c);
				//插值
				float leftTop_top = ratio * leftTop + (1 - ratio) * top;
				float rightBottom_bottom = ratio * rightBottom + (1 - ratio) * bottom;
				if (mag > leftTop_top && mag > rightBottom_bottom)
					edgeMag_nonMaxSup.at<float>(r, c) = mag;



			}
			//右上方和上方的插值 左下方和下方的插值
			if ((angle > 90 && angle <= 135) || (angle > -90 && angle <= -45))
			{
				float ratio = abs(x / y);
				float rightTop_top = ratio * rightTop + (1 - ratio) * top;
				float leftBottom_bottom = ratio * leftBottom + (1 - ratio) * bottom;
				if (mag > rightTop_top && mag > leftBottom_bottom)
					edgeMag_nonMaxSup.at<float>(r, c) = mag;



			}
			//左上方和左方的插值 右下方和右方的插值
			if ((angle >= 0 && angle <= 45) || (angle > -180 && angle <= -135))
			{
				float ratio = y / x;
				float rightBottom_right = ratio * rightBottom + (1 - ratio) * right;
				float leftTop_left = ratio * leftTop + (1 - ratio) * left;
				if (mag > rightBottom_right && mag > leftTop_left)
					edgeMag_nonMaxSup.at<float>(r, c) = mag;


			}
			//右上方和右方的插值 左下方和左方的插值
			if ((angle > 135 && angle <= 180) || (angle > -45 && angle <= 0))
			{
				float ratio = abs(y / x);
				float rightTop_right = ratio * rightTop + (1 - ratio) * right;
				float leftBottom_left = ratio * leftBottom + (1 - ratio) * left;
				if (mag > rightTop_right && mag > leftBottom_left)
					edgeMag_nonMaxSup.at<float>(r, c) = mag;


			}
		}
	}

	return edgeMag_nonMaxSup;
}
bool MainWindow::checkInRange(int r, int c, int rows, int cols)
{
	if (r >= 0 && r < rows && c >= 0 && c < cols)
		return true;
	else
		return false;
}

//从确定边缘点出发，延长边缘
void MainWindow::trace(Mat edgeMag_nonMaxSup, Mat& edge, float lowerThresh, int r, int c, int rows, int cols)
{
	if (edge.at<uchar>(r, c) == 0)
	{
		edge.at<uchar>(r, c) = 255;
		for (int i = -1; i <= 1; i++)//1
		{
			for (int j = -1; j <= 1; j++)//1
			{
				float mag = edgeMag_nonMaxSup.at<float>(r + i, c + j);
				if (checkInRange(r + i, c + j, rows, cols) && mag >= lowerThresh)
					trace(edgeMag_nonMaxSup, edge, lowerThresh, r + i, c + j, rows, cols);
			}
		}
	}
}

//双阈值的滞后阈值处理
Mat MainWindow::hysteresisThreshold(Mat edgeMag_nonMaxSup, float lowerThresh, float upperThresh)
{
	//宽高
	int rows = edgeMag_nonMaxSup.rows;
	int cols = edgeMag_nonMaxSup.cols;
	//最后的边缘输出图
	Mat edge = Mat::zeros(Size(cols, rows), CV_8UC1);
	//滞后阈值处理
	for (int r = 1; r < rows - 1; r++)//1
	{
		for (int c = 1; c < cols - 1; c++)//1
		{
			float mag = edgeMag_nonMaxSup.at<float>(r, c);
			//大于高阈值，可作为确定边缘点被接受
			//并以该点位起始点延长边缘
			if (mag >= upperThresh)
				trace(edgeMag_nonMaxSup, edge, lowerThresh, r, c, rows, cols);
			//小于低阈值的直接被剔除
			if (mag < lowerThresh)
				edge.at<uchar>(r, c) = 0;
		}
	}
	return  edge;
}

double MainWindow::contourNCC(vector<Mat> drr, vector<Mat> xray) {
	double similarity = 0.;
	for (int t = 0; t < drr.size(); t++) {
		similarity += NCC_whole(drr[t], xray[t]);
	}
	return similarity / drr.size();
}

Mat MainWindow::binary_image(Mat image) {
	Mat pre_dst = preprocessing(image);
	//imshow("pre", pre_dst);
	Mat mask = blackegde_binary(pre_dst);
	//imshow("555", mask);
	Mat dstImage;
	Mat dotMask;
	// 初始化自适应阈值参数
	int blockSize = 35;
	int constValue = 7;
	const int maxVal = 255;

	int dotMaskConst = 35;
	/* 自适应阈值算法
	0：ADAPTIVE_THRESH_MEAN_C
	1: ADAPTIVE_THRESH_GAUSSIAN_C
	阈值类型
	0: THRESH_BINARY
	1: THRESH_BINARY_INV */
	int adaptiveMethod = 0;
	int thresholdType = 1;
	// 图像自适应阈值操作
	cv::adaptiveThreshold(pre_dst, dstImage,
		maxVal, adaptiveMethod,
		thresholdType, blockSize,
		constValue);
	cv::adaptiveThreshold(pre_dst, dotMask,
		maxVal, adaptiveMethod,
		thresholdType, blockSize,
		dotMaskConst);
	
	//imshow("dst", dstImage);
	//imshow("dotMask", dotMask);
	//waitKey(0);
	int rowNumber = dstImage.rows;
	int colNumber = dstImage.cols;
	for (int q = 0; q < rowNumber; q++)
	{
		uchar* data = mask.ptr<uchar>(q);
		uchar* data2 = dstImage.ptr<uchar>(q);

		for (int j = 0; j < colNumber; j++)
		{
			//开始处理
			if (data[j] == 0) {
				data2[j] = 0;
			}
		}
	}
	
	Mat kernel = getStructuringElement(MORPH_ELLIPSE, Size(3, 3), Point(-1, -1));
	erode(dstImage, dstImage, kernel);
	dilate(dstImage, dstImage, kernel);

	dilate(dotMask, dotMask, kernel);
	for (int q = 0; q < rowNumber; q++) {
		for (int j = 0; j < colNumber; j++) {
			if (dotMask.at<uchar>(q, j) == 255) dstImage.at<uchar>(q, j) = 0;
		}
	}
	return dstImage;

}
void MainWindow::drawDRRs(std::vector<double> parameters)
{
	setParam(parameters);
	drawDRR();
}
void MainWindow::drawDRRsLateral(std::vector<double> parameters)
{
	setParam(parameters);
	drawlateralDRR();
}
Mat MainWindow::blackegde_binary(Mat image) {
	Mat dstImage;
	threshold(image, dstImage, 45, 255, THRESH_BINARY);
	Mat kernel2 = getStructuringElement(MORPH_ELLIPSE, Size(3, 3), Point(-1, -1));
	erode(dstImage, dstImage, kernel2);
	return dstImage;
}
Mat MainWindow::preprocessing(Mat image) {
	Mat blur_dst;
	GaussianBlur(image, blur_dst, Size(5, 5), 0, 0);
	Mat gray;
	cvtColor(blur_dst, gray, cv::COLOR_BGR2GRAY);
	return gray;
}

vector<Mat> MainWindow::ceilImgGet(Mat image) {

	bool halcon = false;

	Mat binary_dst, canny_dst;
	if (halcon) {
		//Halcon算法轮廓提取
		cv::cvtColor(image, image, COLOR_BGR2GRAY);

		Halcon::Hobject hObj = Halcon::Hobject();

		/*uchar* data = new uchar[image.cols * image.rows];
		for (int i = 0; i < image.rows; i++) memcpy(data + image.cols * i, image.data + image.step * i, image.cols);
		gen_image1(&hObj, "byte", image.cols, image.rows, (Hlong)data);*/
		if (image.type() == CV_8UC1) {
			Halcon::gen_image1(&hObj, "byte", image.cols, image.rows, (Hlong)image.ptr(0));
		}
		else if (image.type() == CV_8UC3) {
			std::vector<cv::Mat> vec;
			cv::split(image, vec);

			cv::Mat imgB = vec[0];
			cv::Mat imgG = vec[1];
			cv::Mat imgR = vec[2];

			Halcon::Hobject himgR = Halcon::Hobject();
			Halcon::Hobject himgG = Halcon::Hobject();
			Halcon::Hobject himgB = Halcon::Hobject();
			Halcon::gen_image1(&himgR, "byte", image.cols, image.rows, (Hlong)imgR.ptr(0));
			Halcon::gen_image1(&himgG, "byte", image.cols, image.rows, (Hlong)imgG.ptr(0));
			Halcon::gen_image1(&himgB, "byte", image.cols, image.rows, (Hlong)imgB.ptr(0));

			compose3(himgR, himgG, himgB, &hObj);

		}

		//Halcon::HImage hImg(hObj);
		Hobject   ImaAmp, ImaDir, Edges, Skeleton;
		Hobject  Contours;

		//read_image(&Image, "F:\\santan\\data\\spine-1\\12\\x\\00002779.bmp");
		edges_image(hObj, &ImaAmp, &ImaDir, "lanser2", 0.3, "nms", 4, 20);
		threshold(ImaAmp, &Edges, 1, 255);
		Hobject  binImg;
		HTuple width;
		HTuple height;
		get_image_size(hObj, &width, &height);
		region_to_bin(Edges, &binImg, 255, 0, width, height);
		write_image(binImg, "bmp", -1, "halcon.bmp");
		binImg.~Hobject();
		Edges.~Hobject();
		hObj.~Hobject();
		ImaAmp.~Hobject();
		ImaDir.~Hobject();

		/*HTuple htCh = HTuple();
		char* cType =(char*) "";

		HImage hImg(binImg);

		hImg.ConvertImageType("byte");
		convert_image_type(binImg, &binImg, "byte");
		count_channels(binImg, &htCh);
		Hlong wid;
		Hlong hgt;
		int W, H;
		if (htCh[0].I() == 1) {
			Hlong ptr;
			get_image_pointer1(binImg, &ptr, cType, &wid, &hgt);
			W = wid;
			H = hgt;

			uchar* pdata = (uchar*)ptr;
			memcpy(canny_dst.data, pdata, W * H);
			cv::resize(canny_dst, canny_dst, Size(512, 512));
			imwrite("Halcon", canny_dst);
		}
		else if (htCh[0].I() == 3)
		{
			unsigned char* ptrR = NULL, * ptrG = NULL, * ptrB = NULL;
			unsigned char* data = NULL;
			char imgType[128] = { 0 };

			get_image_pointer3(binImg, (Hlong*)&ptrR, (Hlong*)&ptrG, (Hlong*)&ptrB, imgType, &wid, &hgt);
			W = wid;
			H = hgt;

			vector<Mat> vecM(3);
			vecM[2].create(H, W, CV_8UC1);
			vecM[1].create(H, W, CV_8UC1);
			vecM[0].create(H, W, CV_8UC1);
			uchar* pr = (uchar*)ptrR;
			uchar* pg = (uchar*)ptrG;
			uchar* pb = (uchar*)ptrB;
			memcpy(vecM[2].data, pr, W * H);
			memcpy(vecM[1].data, pg, W * H);
			memcpy(vecM[0].data, pb, W * H);
			merge(vecM, canny_dst);
			cv::cvtColor(canny_dst, canny_dst, COLOR_BGR2GRAY);
		}
		*/
		canny_dst = imread("halcon.bmp", 1);
	}
	else {
		//原生轮廓提取
		
		binary_dst = binary_image(image);
		
		Canny(binary_dst, canny_dst, 200, 300);
		//imshow("canny", canny_dst);
		//cvtColor(image, image, COLOR_BGR2GRAY);

		
	}
	
	//imshow("dst", canny_dst);
	//imwrite("1.jpg", canny_dst);
	//waitKey(0);
	int N = int(1024 * 1.0 / 40) + 1;
	int M = int(1024 * 1.0 / 40) + 1;
	int SUB_HEIGHT = int(1024 / N);
	int SUB_WIDTH = int(1024 / M);
	vector<Mat> ceil_img;

	Mat image_cut, roi_img, image_cut2, roi_img2;
	for (int j = 0; j < N; j++)
	{
		for (int i = 0; i < M; i++)
		{
			Rect rect(i * SUB_WIDTH, j * SUB_HEIGHT, SUB_WIDTH, SUB_HEIGHT);
			image_cut = Mat(canny_dst, rect);
			roi_img = image_cut.clone();
			ceil_img.push_back(roi_img);

		}
	}
	return ceil_img;

}

float MainWindow::similarity_count(vector<Mat>ceil_img, vector<Mat>ceil_img2) {
	vector<Point> edgePoint, edgePoint2;
	bool temp;
	int count = 0, edgePointNum = 0;
	Point pt, pt_1;
	for (int t = 0; t < ceil_img.size(); t++)
	{
		//
		edgePoint.clear();
		edgePoint2.clear();

		for (int row = 0; row < ceil_img[t].rows; row++)
		{
			uchar* point1 = ceil_img[t].ptr<uchar>(row);
			uchar* point2 = ceil_img2[t].ptr<uchar>(row);

			for (int col = 0; col < ceil_img[t].cols; col++)
			{
				/*if (point1[col] == 255) {
					pt = Point(row, col);
					edgePoint.push_back(pt);
				}
				if (point2[col] == 255) {
					pt_1 = Point(row, col);
					edgePoint2.push_back(pt_1);
				}*/
				if (point1[col] >= 250) {
					pt = Point(row, col);
					edgePoint.push_back(pt);//得到浮动图像和参考图像的所有边缘点
				}
				if (point2[col] >= 250) {
					pt_1 = Point(row, col);
					edgePoint2.push_back(pt_1);
				}
			}
		}
		//----------------------
		edgePointNum = edgePointNum + edgePoint.size();
		//cout << t << ":" << edgePoint.size() << endl;
		//cout << t << ":" << edgePoint2.size() << endl;
		Mat kernelX = getGaussianKernel(13, 3);
		Mat kernelY = getGaussianKernel(13, 3);
		Mat gaussianKernel = kernelX * kernelY.t();
		gaussianKernel = gaussianKernel / gaussianKernel.at<double>(6, 6);
		if (edgePoint.size() != 0 && edgePoint2.size() != 0) {
			temp = false;
			for (int a = 0; a < edgePoint.size(); a++) {
				double max = 0.;
				int xleft = edgePoint[a].x - 6;
				int xright = edgePoint[a].x + 6;
				int yup = edgePoint[a].y - 6;
				int ydown = edgePoint[a].y + 6;
				
				for (int b = 0; b < edgePoint2.size(); b++) {
					int x_2 = edgePoint2[b].x;
					int y_2 = edgePoint2[b].y;
					if (isInParRect(xleft, ydown, xright, yup, x_2, y_2) == true) {
						double value = gaussianKernel.at<double>(x_2 - xleft, y_2 - yup);
						if (value > max) max = value;
					}
				}
				count += max;
			}

		}

	}
	if (edgePointNum == 0) {
		std::cout << "Template image contour detection failed." << std::endl;
		Similarity = 0;
	}
	else {
		Similarity = float(count * 1.00000 / edgePointNum);
		cout << "Similarity_count:" << setiosflags(ios::fixed) << setprecision(6) << Similarity << std::endl;
		//cout << setiosflags(ios::fixed) << setprecision(6) << Similarity << endl;

	}
	return Similarity;
}

float MainWindow::similarity_count1(vector<Mat>ceil_img, vector<Mat>ceil_img2) {
	vector<Point> edgePoint, edgePoint2;

	double count = 0.0;
	int edgePointNum = 0;
	float Similarity;
	Point pt, pt_1;
	for (int t = 0; t < ceil_img.size(); t++)
	{
		//
		edgePoint.clear();
		edgePoint2.clear();

		for (int row = 0; row < ceil_img[t].rows; row++)
		{
			uchar* point1 = ceil_img[t].ptr<uchar>(row);
			uchar* point2 = ceil_img2[t].ptr<uchar>(row);

			for (int col = 0; col < ceil_img[t].cols; col++)
			{
			
				if (point1[col] >= 250) {
					pt = Point(row, col);
					edgePoint.push_back(pt);//得到浮动图像和参考图像的所有边缘点
				}
				if (point2[col] >= 250) {
					pt_1 = Point(row, col);
					edgePoint2.push_back(pt_1);
				}
			}
		}
		edgePointNum = edgePointNum + edgePoint.size();
		if (edgePoint.size() != 0) {
			for (int a = 0; a < edgePoint.size(); a++) {

				int x1left = edgePoint[a].x - 4;
				int x1right = edgePoint[a].x + 4;
				int y1up = edgePoint[a].y - 4;
				int y1down = edgePoint[a].y + 4;

				int x2left = edgePoint[a].x - 6;
				int x2right = edgePoint[a].x + 6;
				int y2up = edgePoint[a].y - 6;
				int y2down = edgePoint[a].y + 6;

				int x3left = edgePoint[a].x - 10;
				int x3right = edgePoint[a].x + 10;
				int y3up = edgePoint[a].y - 10;
				int y3down = edgePoint[a].y + 10;

				int x4left = edgePoint[a].x - 15;
				int x4right = edgePoint[a].x + 15;
				int y4up = edgePoint[a].y -15;
				int y4down = edgePoint[a].y + 15;
				if (isExistEdgePtInRect(x1left, x1right, y1up, y1down, ceil_img2[t]) == true) {
					count += 1;
					continue;
				}
				else if (isExistEdgePtInRect(x2left, x2right, y2up, y2down, ceil_img2[t]) == true) {
					count += 0.6;
					continue;
				}
				else if (isExistEdgePtInRect(x3left, x3right, y3up, y3down, ceil_img2[t]) == true) {
					count += 0.3;
					continue;
				}
				else if (isExistEdgePtInRect(x4left, x4right, y4up, y4down, ceil_img2[t]) == true) {
					count += 0.1;
					continue;
				}
			}
		}
	}
	if (edgePointNum == 0) {
		std::cout << "Template image contour detection failed." << std::endl;
		Similarity = 0;
	}
	else {
		Similarity = float(count * 1.00000 / edgePointNum);
		cout << "Similarity:" << setiosflags(ios::fixed) << setprecision(6) << Similarity << std::endl;
		//cout << setiosflags(ios::fixed) << setprecision(6) << Similarity << endl;
	}
	return Similarity;
}
bool MainWindow::isExistEdgePtInRect(int xleft, int xright, int yup, int ydown, Mat image) {
	if (xleft < 0) xleft = 0;
	if (xright >= image.cols) xright = image.cols - 1;
	if (yup < 0) yup = 0;
	if (ydown >= image.rows) ydown = image.rows - 1;
	for (int cols = xleft; cols < xright; cols++) {
		for (int rows = yup; rows < ydown; rows++) {
			if (image.at<uchar>(cols, rows) == 255) return true;
		}
	}
	return false;
}

bool MainWindow::isInParRect(int x1, int y1, int x4, int y4, int x, int y) {
	if (x <= x1) {
		return false;
	}
	if (x >= x4) {
		return false;
	}
	if (y >= y1) {
		return false;
	}
	if (y <= y4) {
		return false;
	}
	else {
		return true;
	}
}


double MainWindow::getSimilarity_NCC(cv::Mat drr_image, cv::Mat xray_image) {
	cvtColor(drr_image, drr_image, COLOR_BGR2GRAY);
	cvtColor(xray_image, xray_image, COLOR_BGR2GRAY);
	cv::equalizeHist(drr_image, drr_image);
	cv::equalizeHist(xray_image, xray_image);
	cv::GaussianBlur(drr_image, drr_image, cv::Size(5, 5), 5, 5);
	cv::GaussianBlur(xray_image, xray_image, cv::Size(5, 5), 5, 5);
	//cv::threshold(drr_image, drr_image, 100, 255, cv::THRESH_BINARY);
	//cv::threshold(xray_image, xray_image, 150, 255, cv::THRESH_BINARY);
	//cv::adaptiveThreshold(drr_image, drr_image, 255, cv::ADAPTIVE_THRESH_MEAN_C, cv::THRESH_BINARY, 35, 5);
	//cv::adaptiveThreshold(xray_image, xray_image, 255, cv::ADAPTIVE_THRESH_MEAN_C, cv::THRESH_BINARY, 35, 5);
	//cv::imshow("drr",drr_image);
	//cv::imshow("x", xray_image);
	//cv::waitKey(0);
	//cv::destroyAllWindows();

	return NCC_whole(drr_image,xray_image);
}

void MainWindow::CMA_ES() {

	CMAES<double> evo;
	double* arFunvals, * const* pop, * xfinal, * current_x;

	// Initialize everything
	const int dim = 6;
	double xstart[dim];
	int a = 0, b = 180;
	for (int i = 0; i < dim; i++) 
		xstart[i] = (rand() % (b - a + 1)) + a;

	double stddev[dim];
	for (int i = 0; i < dim; i++) 
		stddev[i] = 0.5;
	Parameters<double> parameters;
	// TODO Adjust parameters here
	parameters.init(dim, xstart, stddev);
	arFunvals = evo.init(parameters);

	std::cout << evo.sayHello() << std::endl;

	// Iterate until stop criterion holds
	while (!evo.testForTermination())
	{
		// Generate lambda new search points, sample population
		pop = evo.samplePopulation(); // Do not change content of pop


		// evaluate the new search points using fitfun from above
		for (int i = 0; i < evo.get(CMAES<double>::Lambda); ++i) {
			vector<double> para;
			para.push_back((double)pop[i][0]);
			para.push_back((double)pop[i][1]);
			para.push_back((double)pop[i][2]);
			para.push_back((double)pop[i][3]);
			para.push_back((double)pop[i][4]);
			para.push_back((double)pop[i][5]);
			setParam(para);
			drawDRR();
			arFunvals[i] = calSimilarity(2);
			//arFunvals[i] = cost(pop[i], (int)evo.get(CMAES<double>::Dimension));
						current_x = evo.getNew(CMAES<double>::XMean);
			//printf("value = %f\t", 1 / arFunvals[i]);
			cout << "value=" << calSimilarity(2) << std::endl;
			//printf("x = %f,%f,%f,%f,%f,%f\n", pop[i][0], pop[i][1], pop[i][2],pop[i][3], pop[i][4], pop[i][5]);
		}
		// update the search distribution used for sampleDistribution()
		evo.updateDistribution(arFunvals);
	}
	std::cout << "Stop:" << std::endl << evo.getStopMessage();
	evo.writeToFile(CMAES<double>::WKResume, "resumeevo1.dat"); // write resumable state of CMA-ES

	// get best estimator for the optimum, xmean
	xfinal = evo.getNew(CMAES<double>::XMean); // "XBestEver" might be used as well
	vector<double> para;
	para.push_back((double)xfinal[0]);
	para.push_back((double)xfinal[1]);
	para.push_back((double)xfinal[2]);
	para.push_back((double)xfinal[3]);
	para.push_back((double)xfinal[4]);
	para.push_back((double)xfinal[5]);
	setParam(para);
	//printf("%f,%f,%f,%f,%f,%f\n", xfinal);

	/*ui.rx->setText(QString::number((xfinal[0], 'd', 2));
	ui.ry->setText(QString::number(adjustValue(xfinal[1]), 'd', 2));
	ui.rz->setText(QString::number(adjustValue(xfinal[2]), 'd', 2));
	ui.tx->setText(QString::number(adjustValue(xfinal[3]), 'd', 2));
	ui.ty->setText(QString::number(adjustValue(xfinal[41]), 'd', 2));
	ui.tz->setText(QString::number(adjustValue(xfinal[5]), 'd', 2));*/


	// do something with final solution and finally release memory
	delete[] xfinal;
}

double MainWindow::GetMutualInfo(cv::Mat ref, cv::Mat flt)
{
	cvtColor(ref, ref, COLOR_BGR2GRAY);
	cvtColor(flt, flt, COLOR_BGR2GRAY);

	cv::equalizeHist(ref, ref);
	cv::equalizeHist(flt, flt);
	cv::Size ksize(7, 7);
	GaussianBlur(ref, ref, ksize, 5);
	GaussianBlur(flt, flt, ksize, 5);
	cv::Mat joint_histogram(256, 256, CV_64FC1, cv::Scalar(0));

	for (int i = 0; i < ref.cols; ++i) {
		for (int j = 0; j < ref.rows; ++j) {
			int ref_intensity = ref.at<uchar>(j, i);
			int flt_intensity = flt.at<uchar>(j, i);
			joint_histogram.at<double>(ref_intensity, flt_intensity) = joint_histogram.at<double>(ref_intensity, flt_intensity) + 1;
			double v = joint_histogram.at<double>(ref_intensity, flt_intensity);
		}
	}

	for (int i = 0; i < 256; ++i) {
		for (int j = 0; j < 256; ++j) {
			joint_histogram.at<double>(j, i) = joint_histogram.at<double>(j, i) / (1.0 * ref.rows * ref.cols);
			double v = joint_histogram.at<double>(j, i);
		}
	}

	double entropy = 0.0;
	for (int i = 0; i < 256; ++i) {
		for (int j = 0; j < 256; ++j) {
			double v = joint_histogram.at<double>(j, i);
			if (v > 0.000000000000001) {
				entropy += v * log(v) / log(2);
			}
		}
	}
	entropy *= -1;

	std::vector<double> hist_ref(256, 0.0);
	for (int i = 0; i < joint_histogram.rows; ++i) {
		for (int j = 0; j < joint_histogram.cols; ++j) {
			hist_ref[i] += joint_histogram.at<double>(i, j);
		}
	}

	//cv::Size ksize2(5, 0);
	//  cv::GaussianBlur(hist_ref, hist_ref, ksize2, 5);

	std::vector<double> hist_flt(256, 0.0);
	for (int i = 0; i < joint_histogram.cols; ++i) {
		for (int j = 0; j < joint_histogram.rows; ++j) {
			hist_flt[i] += joint_histogram.at<double>(j, i);
		}
	}

	//   cv::GaussianBlur(hist_flt, hist_flt, ksize2, 5);

	double entropy_ref = 0.0;
	for (int i = 0; i < 256; ++i) {
		if (hist_ref[i] > 0.000000000001) {
			entropy_ref += hist_ref[i] * log2(hist_ref[i]);

		}
	}
	entropy_ref *= -1;
	//std::cout << entropy_ref << "~~ ";

	double entropy_flt = 0.0;
	for (int i = 0; i < 256; ++i) {
		if (hist_flt[i] > 0.000000000001) {
			entropy_flt += hist_flt[i] * log2(hist_flt[i]);
		}
	}
	entropy_flt *= -1;
	// std::cout << entropy_flt << "++ ";

	double mutual_information = 2 * (entropy_flt + entropy_ref - entropy)/(entropy_flt+entropy_ref);
	return mutual_information;
}


//————————————————————————————以下为DE——————————————————————————————————————————————————————————
//函数用于生成指定大小的随机向量，其中每个元素都在给定的上下界范围内
vector<double> MainWindow::generateRandomVector(int size, double lower_bound, double upper_bound) {
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> dis(lower_bound, upper_bound);

	vector<double> result(size);
	for (int i = 0; i < size; ++i) {
		result[i] = dis(gen);
	}
	return result;
}

//用于初始化种群。它根据给定的种群大小、维度和上下界范围，生成一个二维向量，表示具有随机初始值的种群。
vector<vector<double>> MainWindow::initializePopulation(int population_size, int dimension,
	vector<double>& lower_bounds, vector<double>& upper_bounds) {
	vector<vector<double>> population(population_size, vector<double>(dimension));
	int range = 10;
	srand((unsigned)time(NULL));
	parameters.clear();
	parameters.push_back(ui.rx->text().toDouble());
	parameters.push_back(ui.ry->text().toDouble());
	parameters.push_back(ui.rz->text().toDouble());
	parameters.push_back(ui.tx->text().toDouble());
	parameters.push_back(ui.ty->text().toDouble());
	parameters.push_back(ui.tz->text().toDouble());
	vector<double> paraRange(parameters);
	for (int i = 0; i < population_size; ++i) {
		for (int j = 0; j < dimension; ++j) {
			population[i][j] = parameters[j] + (2.0 * rand() / (double)RAND_MAX - 1.0) * range;
		}
	}
	return population;
}

//计算种群中每个个体的适应度值。它调用目标函数，并将目标函数的返回值作为个体的适应度值。
vector<double> MainWindow::calculateFitness(vector<vector<double>>& population) {
	int population_size = population.size();
	vector<double> fitness(population_size);
	vector<double> para;
	//输出一下生成的种群值
	for (int i = 0; i < population_size; ++i) {
		vector<double> temp;
		temp.clear();	

		temp.push_back(population[i][0]);
		temp.push_back(population[i][1]);
		temp.push_back(population[i][2]);
		temp.push_back(population[i][3]);
		temp.push_back(population[i][4]);
		temp.push_back(population[i][5]);

		drawDRRs(temp);
		float s1 = calSimilarity(2);
		fitness[i] = s1;
	}
	return fitness;
}

//返回一个按适应度值排序后的索引向量。它使用sort函数和自定义的排序比较函数，按照适应度值从小到大对索引进行排序。
vector<int> MainWindow::argsort(vector<double>& fitness) {
	vector<int> indices(fitness.size());
	iota(indices.begin(), indices.end(), 0);
	sort(indices.begin(), indices.end(), [&fitness](int a, int b) { return fitness[a] < fitness[b]; });
	return indices;
}

//据适应度值对种群进行排序
vector<vector<double>> MainWindow::sortPopulationByFitness(vector<vector<double>>& population, vector<double>& fitness) {
	vector<int> sorted_indices = argsort(fitness);
	reverse(sorted_indices.begin(), sorted_indices.end());
	vector<vector<double>> sorted_population(population.size(), vector<double>(population[0].size()));
	for (int i = 0; i < population.size(); ++i) {
		sorted_population[i] = population[sorted_indices[i]];
	}
	return sorted_population;
}

vector<vector<double>> MainWindow::differentialEvolution(int population_size, int dimension, vector<double>& lower_bounds, vector<double>& upper_bounds, int max_iterations) {
	//“CR”参数控制亲本在后代生成中的影响。值越高，表示父代的影响越小。
	//“F”参数衡量被选择来计算突变值的一组解对的影响。
	vector<double> crossover_probability = { 0.8, 0.3 };          //交叉概率  c 0.8	0.3 
	vector<double> scaling_factor = { 0.5,0.8 };                  //缩放因子  f 0.5	0.8

	//初始化种群
	vector<vector<double>> population = initializePopulation(population_size, dimension, lower_bounds, upper_bounds);
	//计算初始种群的适应度值
	vector<double> fitness = calculateFitness(population);
	population = sortPopulationByFitness(population, fitness);

	double best_fitness = fitness[0];
	vector<double> best_position = population[0];

	//使用差分进化算子对种群进行迭代更新，生成新的种群，并计算新种群的适应度值
	for (int iteration = 0; iteration < max_iterations; ++iteration) {
		cout << "---------------DE第" << iteration + 1 << "次迭代----------------" << std::endl;
		vector<vector<double>> new_population(population_size, vector<double>(dimension));
		vector<double> new_fitness(population_size);

		for (int i = 0; i < population_size; ++i) {
			//选择其他三个个体作为参考
			vector<double> indices = generateRandomVector(5, 0, population_size - 1);
			vector<double> x1 = population[indices[0]];
			vector<double> x2 = population[indices[1]];
			vector<double> x3 = population[indices[2]];
			vector<double> x4 = population[indices[3]];
			vector<double> x5 = population[indices[4]];
			//计算差分向量mutant
			vector<double> mutant(dimension);
			for (int j = 0; j < dimension; ++j) {
				if (iteration >= 0 && iteration <= max_iterations / 2) {
					//mutant[j] = x1[j] + scaling_factor[0] * (x2[j] - x3[j]);																				//DE/rand/1	
					mutant[j] = x1[j] + scaling_factor[0] * (x2[j] - x3[j]) + scaling_factor[0] * (x4[j] - x5[j]);											//DE/rand/2	
					//mutant = x1 + f/2(x1-x2-x3)//有待参考
				}
				else {
					//mutant[j] = x1[j] + scaling_factor[1] * (x2[j] - x3[j]);														 		        		//DE/rand/1
					mutant[j] = x1[j] + scaling_factor[1] * (x2[j] - x3[j]) + scaling_factor[1] * (x4[j] - x5[j]);											//DE/rand/2		
				}
			}
			//通过交叉概率的控制，将差分向量的部分维度与当前个体的对应维度进行交叉，生成试验向量trial。
			//如果随机生成的概率值小于等于交叉概率，则使用差分向量的对应维度值，否则使用当前个体的对应维度值。
			vector<double> trial(dimension);
			for (int j = 0; j < dimension; ++j) {
				if ((iteration >= 0 && iteration <= max_iterations / 2) && generateRandomVector(1, 0, 1)[0] <= crossover_probability[0])
					/*if(generateRandomVector(1, 0, 1)[0] <= crossover_probability[0])*/
				{
					trial[j] = mutant[j];
				}
				else if (generateRandomVector(1, 0, 1)[0] <= crossover_probability[1]) {
					trial[j] = mutant[j];
				}
				else {
					trial[j] = population[i][j];
				}
			}

			//使用生成的试验向量trial替换当前个体的参数向量，并计算新个体的适应度值。
			new_population[i] = trial;
			this->drawDRRs(trial);
			float s1 = calSimilarity(2);
			new_fitness[i] = s1;
		}

		for (int i = 0; i < population_size; ++i) {
			if (new_fitness[i] > fitness[i]) {
				population[i] = new_population[i];
				fitness[i] = new_fitness[i];
			}
		}

		population = sortPopulationByFitness(population, fitness);

		if (fitness[0] > best_fitness) {
			best_fitness = fitness[0];
			best_position = population[0];
			//使用数组记录最好与次好的两组fitness和position
		}
	}

	return { best_position };
}

void MainWindow::DE() {								//差分进化算法
	int population_size = 10;
	int dimension = 6;
	vector<double> lower_bounds(dimension, -100);
	vector<double> upper_bounds(dimension, 100);
	int max_iterations = 50;

	vector<vector<double>> result = differentialEvolution(population_size, dimension, lower_bounds, upper_bounds, max_iterations);
	loadOuputVariablesInGPUMemory(1024, 1024, 1, 1);
	std::cout << "**best result**:" << std::endl;
	drawDRRs(result[0]);
	float s1 = calSimilarity(2);
	freeAuxiliaryVariablesInGPUMemory(1);
}

void MainWindow::EO() {                         //基于遗传算法的平衡优化器
	int a1 = 1;                                  //exploration系数，控制全局搜索能力
	int a2 = 1;                                  //exploitation系数，控制局部搜索能力
	int GP = 0.5;                                //Genetation probability,生成速率控制浓度更新的参与概率
	int Max_iter = 25;
	int pn = 10;                             //粒子数
	int dim = 6;
	int range =20;
	C = NULL;
	srand((unsigned)time(NULL));
	parameters.clear();
	parameters.push_back(ui.rx->text().toDouble());
	parameters.push_back(ui.ry->text().toDouble());
	parameters.push_back(ui.rz->text().toDouble());
	parameters.push_back(ui.tx->text().toDouble());
	parameters.push_back(ui.ty->text().toDouble());
	parameters.push_back(ui.tz->text().toDouble());
	vector<double> paraRange(parameters);
	//initialize
	if (C == NULL) {
		C = new EOparticle[pn];
		for (int i = 0; i < pn; i++) {
			C[i].value = new double[dim];
			for (int j = 0; j < dim; j++) {
				C[i].value[j] = parameters[j] + (2.0 * rand() / (double)RAND_MAX - 1.0) * range;
			}
		}
	}

	//初始化候选粒子
	double* Ceq1 = new double[dim];
	double* Ceq2 = new double[dim];
	double* Ceq3 = new double[dim];
	double* Ceq4 = new double[dim];

	double Ceq1_fit = 99999.;
	double Ceq2_fit = 99999.;
	double Ceq3_fit = 99999.;
	double Ceq4_fit = 99999.;

	vector<double> convergenceCurve;
	//进入迭代
	int V = 1;
	//double* fit_best;
	//EOparticle* C_best;

	double* Ceq;
	double C_best[6] = { 0. };
	double fit_best = 9999.;
	loadOuputVariablesInGPUMemory(1024, 1024, pn * Max_iter, 1);
	for (int iter = 0; iter < Max_iter; iter++) {
		cout << "---------------第" << iter + 1 << "次迭代----------------" << std::endl;

		double* fitness = new double[pn];
		clock_t part1start, part1end, part2start, part2end;
		part1start = clock();
		vector<double> para;
		for (int i = 0; i < pn; i++) {
			vector<double> temp;
			temp.clear();
			para.push_back(C[i].value[0]);
			para.push_back(C[i].value[1]);
			para.push_back(C[i].value[2]);
			para.push_back(C[i].value[3]);
			para.push_back(C[i].value[4]);
			para.push_back(C[i].value[5]);

			temp.push_back(C[i].value[0]);
			temp.push_back(C[i].value[1]);
			temp.push_back(C[i].value[2]);
			temp.push_back(C[i].value[3]);
			temp.push_back(C[i].value[4]);
			temp.push_back(C[i].value[5]);

			drawDRRs(temp);
			
			fitness[i] = calSimilarity(7);

		}

		part1end = clock();
		//fitness = drawDRRs(para);

		part2start = clock();
		for (int i = 0; i < pn; i++) {
			fitness[i] = -fitness[i];
			//cout << "s======" << fitness[i] << std::endl;
		}
		part1start = clock();
		for (int i = 0; i < pn; i++) {
			//范围限定检查
			for (int j = 0; j < dim; j++) {
				if (C[i].value[j] < paraRange[j] - range) C[i].value[j] = paraRange[j] - range;
				else if (C[i].value[j] > paraRange[j] + range) C[i].value[j] = paraRange[j] + range;
			}

			if (fitness[i] < fit_best) {
				fit_best = fitness[i];
				for (int j = 0; j < dim; j++) {
					C_best[j] = C[i].value[j];
				}
				cout << "entered" << std::endl;
			}
			cout << "best_fit=" << 1. - fit_best << std::endl;
			cout << "C_best" << std::endl;

			for (int j = 0; j < dim; j++) {
				cout << C_best[j] << std::endl;
			}

			convergenceCurve.push_back(fitness[i]);
			
			if (fitness[i] < Ceq1_fit) {
				Ceq1_fit = fitness[i];
				Ceq1 = C[i].value;
			}
			else if (fitness[i] > Ceq1_fit && fitness[i] < Ceq2_fit) {
				Ceq2_fit = fitness[i];
				Ceq2 = C[i].value;
			}
			else if (fitness[i] > Ceq2_fit && fitness[i] < Ceq3_fit) {
				Ceq3_fit = fitness[i];
				Ceq3 = C[i].value;
			}
			else if (fitness[i] > Ceq3_fit && fitness[i] < Ceq4_fit) {
				Ceq4_fit = fitness[i];
				Ceq4 = C[i].value;
			}
		}
		//memory saving

		if (iter == 0) {
			/*	fit_best = fitness[];
				C_best = new EOparticle[pn];
				C_best = C;*/
		}

		//构建Cave
		double* Cave = new double[dim];
		for (int j = 0; j < dim; j++) {
			Cave[j] = (Ceq1[j] + Ceq2[j] + Ceq3[j] + Ceq4[j]) / 4;
		}
		//构建平衡池
		vector<double*> C_pool = { Ceq1,Ceq2,Ceq3,Ceq4,Cave };

		double t = pow((1 - iter / Max_iter), (a2 * (iter / Max_iter)));      //  Eq(9)
		double* lambda = new double[dim];
		double* r = new double[dim];
		double* F = new double[dim];
		double* GCP = new double[dim];
		double* G0 = new double[dim];
		double* G = new double[dim];
		for (int i = 0; i < pn; i++) {
			for (int j = 0; j < dim; j++) {
				lambda[j] = rand() / (double)RAND_MAX;
				r[j] = rand() / (double)RAND_MAX;
				F[j] = a1 * sign(r[j] - 0.5) * (exp(-lambda[j] * t) - 1);
				double r1 = rand() / (double)RAND_MAX;
				double r2 = rand() / (double)RAND_MAX;
				GCP[j] = 0.5 * r1 * (r2 >= GP);
			}

			Ceq = C_pool[int(5 * (rand() / (double)RAND_MAX))];
			for (int j = 0; j < dim; j++) {
				G0[j] = GCP[j] * (Ceq[j] - lambda[j] * C[i].value[j]);
				G[j] = G0[j] * F[j];
				C[i].value[j] = Ceq[j] + (C[i].value[j] - Ceq[j]) * F[j] + G[j] / (lambda[j] * V) * (1 - F[j]);
			}
		}
		part2end = clock();
		//cout << "parallel time:" << (double)(part2end - part1start) / CLOCKS_PER_SEC << std::endl;
		convergenceCurve.push_back(Ceq1_fit);
	}
	for (auto iter = convergenceCurve.begin(); iter != convergenceCurve.end(); iter++) {
		cout << *iter << std::endl;
	}
	vector<double> para;
	cout << "最优参数：" << std::endl;
	cout << "best similarity:" << -fit_best << std::endl;
	for (int j = 0; j < dim; j++) {
		para.push_back(C_best[j]);
	}
	loadOuputVariablesInGPUMemory(1024, 1024, 1, 1);

	/*double temp = para[0];
	temp += 10.0;*/
	//para[0] += 10;
	//para[5] -= 10; 
	//para[4] -= 10;//手动修改的地方1
	drawDRRs(para);
	//drawDRRsLateral(para);
	freeAuxiliaryVariablesInGPUMemory(1);
}

int MainWindow::sign(double x) {
	if (x > 0) 
		return 1; 
	if (x < 0) 
		return -1; 
	return 0;
}

void MainWindow::batch() {
	int rxstart = 70;
	int rxend = 115;
	int rystart = -20;
	int ryend = 25;
	int rzstart = -20;
	int rzend = 25;
	int xstart = -20;
	int xend = 20;
	int ystart = -20;
	int yend = 40;
	int zstart = -140;
	int zend = -270;
	ui.rx->setText(QString::number(100));
	ui.ry->setText(QString::number(0));
	ui.rz->setText(QString::number(0));
	ui.tx->setText(QString::number(-5));
	ui.ty->setText(QString::number(-10));
	ui.tz->setText(QString::number(-250));

	for (int i = rxstart; i <= rxend; i += 5) {
		ui.rx->setText(QString::number(i));
		updateParam();
	}
	cout << std::endl;
	ui.rx->setText(QString::number(90));

	for (int i = rystart; i <= ryend; i += 5) {
		ui.ry->setText(QString::number(i));
		updateParam();
	}
	cout << std::endl;
	ui.ry->setText(QString::number(0));

	for (int i = rzstart; i <= rzend; i += 5) {
		ui.rz->setText(QString::number(i));
		updateParam();
	}
	cout << std::endl;
	ui.rz->setText(QString::number(0));
	
	for (int i = xstart; i <= xend; i += 10) {
		ui.tx->setText(QString::number(i));
		updateParam();
	}
	ui.tx->setText(QString::number(-5));
	cout << std::endl;

	for (int y = ystart; y >= yend; y -= 10) {
		ui.ty->setText(QString::number(y));
		updateParam();
	}
	ui.ty->setText(QString::number(-80));
	cout << std::endl;

	for (int z = zstart; z >= zend; z -= 10) {
		ui.tz->setText(QString::number(z));
		updateParam();
	}
	ui.tz->setText(QString::number(-210));
	cout << std::endl;
}

double MainWindow::SSD(cv::Mat ref, cv::Mat flt) {
	int height = ref.rows;
	int width = ref.cols;
	double ssd = 0;
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			ssd += pow((double)ref.at<uchar>(i, j) - (double)flt.at<uchar>(i, j), 2);
		}
	}
	ssd = ssd / (height * width);
	return ssd;
}

double MainWindow::FaceNet_PSNR(cv::Mat& ref, cv::Mat& flt) {
	//注意，当两幅图像一样时这个函数计算出来的psnr为0 
	Mat s1;
	absdiff(ref, flt, s1);
	s1.convertTo(s1, CV_32F);//转换为32位的float类型，8位不能计算平方  
	s1 = s1.mul(s1);
	Scalar s = sum(s1);  //计算每个通道的和  
	double sse = s.val[0] + s.val[1] + s.val[2];
	if (sse <= 1e-10) // for small values return zero  
		return 0;
	else
	{
		double mse = sse / (double)(ref.channels() * ref.total()); //  sse/(w*h*3)  
		double psnr = 10.0 * log10((255 * 255) / mse);
		return psnr;
	}
}

double MainWindow::SiameseSSMI(cv::Mat ref, cv::Mat flt) {
	const double C1 = 6.5025, C2 = 58.5225;     //用于稳定计算的参数

	int d = CV_32F;

	Mat I1, I2;
	ref.convertTo(I1, d);
	flt.convertTo(I2, d);

	//计算图像的平方
	Mat I2_2 = I2.mul(I2);
	Mat I1_2 = I1.mul(I1);
	Mat I1_I2 = I1.mul(I2);

	//对图像进行高斯模糊操作，计算图像的均值
	Mat mu1, mu2;
	GaussianBlur(I1, mu1, Size(11, 11), 1.5);
	GaussianBlur(I2, mu2, Size(11, 11), 1.5);

	Mat mu1_2 = mu1.mul(mu1);
	Mat mu2_2 = mu2.mul(mu2);
	Mat mu1_mu2 = mu1.mul(mu2);

	Mat sigma1_2, sigma2_2, sigma12;
	GaussianBlur(I1_2, sigma1_2, Size(11, 11), 1.5);
	sigma1_2 -= mu1_2;
	GaussianBlur(I2_2, sigma2_2, Size(11, 11), 1.5);
	sigma2_2 -= mu2_2;
	GaussianBlur(I1_I2, sigma12, Size(11, 11), 1.5);
	sigma12 -= mu1_mu2;

	//计算图像的方差   
	Mat t1, t2, t3;
	t1 = 2 * mu1_mu2 + C1;
	t2 = 2 * sigma12 + C2;
	t3 = t1.mul(t2);
	t1 = mu1_2 + mu2_2 + C1;
	t2 = sigma1_2 + sigma2_2 + C2;
	t1 = t1.mul(t2);

	Mat ssim_map;
	divide(t3, t1, ssim_map);
	Scalar mssim = mean(ssim_map);

	return mssim[0];
}


void MainWindow::calMaxOfZaxis(vector<double> para) {
	//cout << "calculating Z axis..." << std::endl;
	//斐波那契法求Z轴极值
	double epsilon = 0.0001, f_max;
	vector<double>fn(2);
	fn[0] = 1; fn[1] = 1;
	int i;
	double a_fib = para[5] - 20.0, b_fib = para[5] + 20.0, f_lambda = 0, f_mu = 0;
	for (i = 2; fn[i - 1] < 1.0 / epsilon; i++)
	{
		fn.push_back(fn[i - 2] + fn[i - 1]);
	}
	//n_fib：迭代次数
	int n_fib = i-1;
	
	double lambda = a_fib + (fn[n_fib - 2] / fn[n_fib]) * (b_fib - a_fib), mu = a_fib + (fn[n_fib - 1] / fn[n_fib]) * (b_fib - a_fib);
	para[5] = lambda;
	drawDRRsLateral(para);
	f_lambda = calSimilarityLateral(2);
	para[5] = mu;
	drawDRRsLateral(para);
	f_mu =1.- calSimilarityLateral(2);

	for (int i = 0; i < n_fib - 2; i++)
	{
		if (f_lambda < f_mu)
		{
			a_fib = a_fib, b_fib = mu, mu = lambda, f_mu = f_lambda;
			lambda = a_fib + (fn[n_fib - i - 3] / fn[n_fib - 1 - i]) * (b_fib - a_fib);
			para[5] = lambda;
			drawDRRsLateral(para);
			f_lambda =1.- calSimilarityLateral(2);
		}
		else
		{
			a_fib = lambda, b_fib = b_fib, lambda = mu, f_lambda = f_mu;
			mu = a_fib + (fn[n_fib - i - 2] / fn[n_fib - 1 - i]) * (b_fib - a_fib);
			para[5] = mu;
			drawDRRsLateral(para);
			f_mu =1.- calSimilarityLateral(2);
		}
		//cout << fixed << setw(12) << setprecision(5) << a_fib << fixed << setw(12) << setprecision(5) << lambda << fixed << setw(12) << setprecision(5) << mu << fixed << setw(12) << setprecision(5) << b_fib << endl;
	}
	mu = lambda + 0.001;
	para[5] = mu;
	drawDRRsLateral(para);
	f_mu = 1. - calSimilarityLateral(2);
	//f_mu = fx(mu);
	if (f_lambda < f_mu)
	{
		a_fib = a_fib, b_fib = mu;
		para[5] = a_fib + 0.5 * (b_fib - a_fib);
		drawDRRsLateral(para);
		f_max =  calSimilarityLateral(2);
	}
	else
	{
		a_fib = lambda, b_fib = b_fib; 
		para[5] = a_fib + 0.5 * (b_fib - a_fib);
		drawDRRsLateral(para);
		f_max =  calSimilarityLateral(2);
	}
	para[5] = a_fib + 0.5 * (b_fib - a_fib);
	setParam(para);
}