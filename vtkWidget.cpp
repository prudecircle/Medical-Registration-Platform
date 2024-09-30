#include "vtkWidget.h"
//#include "ui_vtkWidget.h"
#include "vtkRenderWindow.h"
#include "vtkImageViewer2.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkMatrix4x4.h"
#include "vtkMatrix3x3.h"
#include "vtkAxesActor.h"
#include <QDebug>

vtkWidget::vtkWidget(QWidget *parent)
	: QDialog(parent), ui(new Ui::vtkWidget)
{
	ui->setupUi(this);
	// 参数设置
	x = 0.0;
	y = 0.0;
	z = 0.0;
	needRefresh = true;

	//相机设置
	style = vtkWidgetTrackballCamera::New();
	style->setWidget(this);
	vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::New();
	iren->SetInteractorStyle(style);
	ui->qvtkWidget->GetRenderWindow()->SetInteractor(iren);

	// table设置
	QStringList header;
	header << QString::fromLocal8Bit("阈值") << QString::fromLocal8Bit("透明度");
	ui->scalarTable->setHorizontalHeaderLabels(header);
	ui->scalarTable->setSelectionBehavior(QAbstractItemView::SelectRows);
	ui->scalarTable->setSelectionMode(QAbstractItemView::SingleSelection); //设置选择模式，选择单行
	ui->scalarTable->setRowCount(3);
	ui->scalarTable->setItem(0, 0, new QTableWidgetItem("-16.4458"));
	ui->scalarTable->setItem(0, 1, new QTableWidgetItem("0.0"));
	ui->scalarTable->setItem(1, 0, new QTableWidgetItem("340.0"));
	ui->scalarTable->setItem(1, 1, new QTableWidgetItem("0.0"));
	ui->scalarTable->setItem(2, 0, new QTableWidgetItem("3000.0"));
	ui->scalarTable->setItem(2, 1, new QTableWidgetItem("1.0"));

	ui->gradientTable->setHorizontalHeaderLabels(header);
	ui->gradientTable->setSelectionBehavior(QAbstractItemView::SelectRows);
	ui->gradientTable->setSelectionMode(QAbstractItemView::SingleSelection); //设置选择模式，选择单行
	ui->gradientTable->setRowCount(2);
	ui->gradientTable->setItem(0, 0, new QTableWidgetItem("0"));
	ui->gradientTable->setItem(0, 1, new QTableWidgetItem("1.0"));
	ui->gradientTable->setItem(1, 0, new QTableWidgetItem("255.0"));
	ui->gradientTable->setItem(1, 1, new QTableWidgetItem("1.0"));

	connect(ui->plane_btn, SIGNAL(clicked()), this, SLOT(create2D()));
	connect(ui->add_btn1, SIGNAL(clicked()), this, SLOT(addScalarOpacity()));
	connect(ui->delete_btn1, SIGNAL(clicked()), this, SLOT(deleteScalarOpacity()));
	connect(ui->add_btn2, SIGNAL(clicked()), this, SLOT(addGradientOpacity()));
	connect(ui->delete_btn2, SIGNAL(clicked()), this, SLOT(deleteGradientOpacity()));
	connect(ui->restore_btn, SIGNAL(clicked()), this, SLOT(restorePara()));

	//connect(ui->xr_edit, SIGNAL(returnPressed()), this, SLOT(setPara()));
	//connect(ui->yr_edit, SIGNAL(returnPressed()), this, SLOT(setPara()));
	//connect(ui->zr_edit, SIGNAL(returnPressed()), this, SLOT(setPara()));
	connect(ui->dx_lineEdit_3,SIGNAL(returnPressed()),this,SLOT(setPara()));

	connect(this, SIGNAL(sendAngle(double,double,double)), this, SLOT(changeAngleAnytime(double,double,double)));
}

vtkWidget::~vtkWidget()
{
	delete ui;
}


void vtkWidget::setFilePath(const char * dicname)
{

	double time = (double)clock() / CLOCKS_PER_SEC;
	path = dicname;
	reader = vtkDICOMImageReader::New();	
	reader->SetDirectoryName(dicname);
	//reader->SetDataSpacing(0.40039101, 0.40039101, 0.6250000);
	reader->Update();

	std::cout << "vtk readfile :" << (double)clock() / CLOCKS_PER_SEC - time << std::endl;
	time = (double)clock() / CLOCKS_PER_SEC;
	vtkImageShiftScale* scale = vtkImageShiftScale::New();
	
	scale->SetInputData(reader->GetOutput());
	scale->SetOutputScalarTypeToUnsignedChar();
	scale->Update();

	vtkSmartPointer<vtkOpenGLGPUVolumeRayCastMapper> volumeMapper =
		vtkSmartPointer<vtkOpenGLGPUVolumeRayCastMapper>::New();
	volumeMapper->SetInputData(reader->GetOutput());

	vtkSmartPointer<vtkPiecewiseFunction> compositeOpacity =
		vtkSmartPointer<vtkPiecewiseFunction>::New();
	compositeOpacity->AddPoint(-16.4458, 0.0);
	compositeOpacity->AddPoint(340.0, 0.0);
	compositeOpacity->AddPoint(3000.0, 1.0);

	vtkSmartPointer<vtkPiecewiseFunction> volumeGradientOpacity =
		vtkSmartPointer<vtkPiecewiseFunction>::New();
	volumeGradientOpacity->AddPoint(0.0, 1.0);
	volumeGradientOpacity->AddPoint(255.0, 1.0);

	vtkSmartPointer<vtkColorTransferFunction> color =
		vtkSmartPointer<vtkColorTransferFunction>::New();
	color->AddRGBPoint(-1000.0, 0.8, 0.8, 0.8);
	color->AddRGBPoint(3000.0, 0.0, 0.0, 0.0);


	vtkSmartPointer<vtkVolumeProperty> volumeProperty =
		vtkSmartPointer<vtkVolumeProperty>::New();
	volumeProperty->SetScalarOpacity(compositeOpacity);
	volumeProperty->SetGradientOpacity(volumeGradientOpacity);
	volumeProperty->SetColor(color);
	volumeProperty->SetInterpolationTypeToLinear();
	volumeProperty->SetShade(false);
	volumeProperty->SetAmbient(0.1);
	volumeProperty->SetDiffuse(0.9);
	volumeProperty->SetSpecular(0.2);
	volumeProperty->SetSpecularPower(10.0);
	
	drrVolume = vtkSmartPointer<vtkVolume>::New();

	drrVolume->SetMapper(volumeMapper);
	drrVolume->SetProperty(volumeProperty);

	style->setVolume(drrVolume);

	drrVolume->SetOrigin(drrVolume->GetCenter());
	//drrVolume->SetOrientation(xr, yr, zr);
	drrVolume->RotateZ(180);
	
	needRefresh = true;
	std::cout << "during vtk:" << (double)clock() / CLOCKS_PER_SEC - time << std::endl;
}

void vtkWidget::addScalarOpacity()
{
	//初始化inputDialog
	inputDialog = new InputDialog();
	inputDialog->setModal(true); //总在最前
	connect(inputDialog, SIGNAL(sendDataList(QList<QString>*)), this, SLOT(ReceiveScalarData(QList<QString>*)));
	inputDialog->show();
	inputDialog->exec();
}

void vtkWidget::deleteScalarOpacity() {
	int i = ui->scalarTable->currentRow();
	ui->scalarTable->removeRow(i);

	changeVolume();
}

void vtkWidget::addGradientOpacity()
{
	//初始化inputDialog
	inputDialog = new InputDialog();
	inputDialog->setModal(true); //总在最前
	connect(inputDialog, SIGNAL(sendDataList(QList<QString>*)), this, SLOT(ReceiveGradientData(QList<QString>*)));
	inputDialog->show();
	inputDialog->exec();
}

void vtkWidget::deleteGradientOpacity()
{
	int i = ui->gradientTable->currentRow();
	ui->gradientTable->removeRow(i);

	changeVolume();
}

void vtkWidget::ReceiveGradientData(QList<QString>* inputList)
{
	ui->gradientTable->insertRow(ui->gradientTable->rowCount());
	QTableWidgetItem *newItem = new QTableWidgetItem();
	newItem->setText(inputList->at(0));
	ui->gradientTable->setItem(ui->gradientTable->rowCount() - 1, 0, newItem);

	QTableWidgetItem *newItem1 = new QTableWidgetItem();
	newItem1->setText(inputList->at(1));
	ui->gradientTable->setItem(ui->gradientTable->rowCount() - 1, 1, newItem1);

	changeVolume();
}

void vtkWidget::ReceiveScalarData(QList<QString>* inputList)
{
	ui->scalarTable->insertRow(ui->scalarTable->rowCount());
	QTableWidgetItem *newItem = new QTableWidgetItem();
	newItem->setText(inputList->at(0));
	ui->scalarTable->setItem(ui->scalarTable->rowCount() - 1, 0, newItem);

	QTableWidgetItem *newItem1 = new QTableWidgetItem();
	newItem1->setText(inputList->at(1));
	ui->scalarTable->setItem(ui->scalarTable->rowCount() - 1, 1, newItem1);

	changeVolume();
}


void vtkWidget::changeVolume() {
	vtkSmartPointer<vtkPiecewiseFunction> volumeScalarOpacity = vtkSmartPointer<vtkPiecewiseFunction>::New();
	for (int i = 1; i < ui->scalarTable->rowCount(); i++) {
		double x = ui->scalarTable->item(i, 0)->text().toDouble();
		double y = ui->scalarTable->item(i, 1)->text().toDouble();
		volumeScalarOpacity->AddPoint(x, y);
	}
	vtkSmartPointer<vtkPiecewiseFunction> volumeGradientOpacity =
		vtkSmartPointer<vtkPiecewiseFunction>::New();
	for (int i = 1; i < ui->gradientTable->rowCount(); i++) {
		double x = ui->gradientTable->item(i, 0)->text().toDouble();
		double y = ui->gradientTable->item(i, 1)->text().toDouble();
		volumeGradientOpacity->AddPoint(x, y);
	}
	vtkSmartPointer<vtkColorTransferFunction> color =
		vtkSmartPointer<vtkColorTransferFunction>::New();
	color->AddRGBPoint(-1000.0, 0.8, 0.8, 0.8);
	color->AddRGBPoint(3000.0, 0.0, 0.0, 0.0);
	vtkSmartPointer<vtkVolumeProperty> volumeProperty =
		vtkSmartPointer<vtkVolumeProperty>::New();
	volumeProperty->SetScalarOpacity(volumeScalarOpacity);
	volumeProperty->SetGradientOpacity(volumeGradientOpacity);
	volumeProperty->SetColor(color);
	volumeProperty->SetInterpolationTypeToLinear();
	volumeProperty->SetShade(false);
	volumeProperty->SetAmbient(0.1);
	volumeProperty->SetDiffuse(0.9);
	volumeProperty->SetSpecular(0.2);
	volumeProperty->SetSpecularPower(10.0);
	drrVolume->SetProperty(volumeProperty);
	tdRenderer->DrawOn();
	ui->qvtkWidget->GetRenderWindow()->Render();
}

void vtkWidget::restorePara()
{
	ui->scalarTable->setRowCount(3);
	ui->scalarTable->setItem(0, 0, new QTableWidgetItem("-16.4458"));
	ui->scalarTable->setItem(0, 1, new QTableWidgetItem("0.0"));
	ui->scalarTable->setItem(1, 0, new QTableWidgetItem("340.0"));
	ui->scalarTable->setItem(1, 1, new QTableWidgetItem("0.0"));
	ui->scalarTable->setItem(2, 0, new QTableWidgetItem("3000.0"));
	ui->scalarTable->setItem(2, 1, new QTableWidgetItem("1.0"));

	ui->gradientTable->setRowCount(2);
	ui->gradientTable->setItem(0, 0, new QTableWidgetItem("0"));
	ui->gradientTable->setItem(0, 1, new QTableWidgetItem("1.0"));
	ui->gradientTable->setItem(1, 0, new QTableWidgetItem("255.0"));
	ui->gradientTable->setItem(1, 1, new QTableWidgetItem("1.0"));

	//ui->xr_edit->setText("0");
	//ui->yr_edit->setText("0");
	//ui->zr_edit->setText("0");
	//drrVolume->SetOrientation(0, 0, 180);

	ui->dx_lineEdit->setText("1024");
	ui->dx_lineEdit_3->setText("170");
	ui->dy_lineEdit->setText("1024");
	
	changeVolume();
	setPara();
}

void vtkWidget::setPara()
{
	
	// 要调整的角度
	//double xr = ui->xr_edit->text().toDouble();
	//double yr = ui->yr_edit->text().toDouble();
	//double zr = ui->zr_edit->text().toDouble();
	
	/*
	// 与当前角度的差
	double xd, yd, zd;
	xd = xr - x;
	yd = yr - y;
	zd = zr - z;
	// 存储调整后的角度
	x = xr;
	y = yr;
	z = zr;
	*/

	//drrVolume->SetScale(20);

	//drrVolume->SetOrigin(drrVolume->GetCenter());
	//drrVolume->SetOrientation(xr, yr, zr);
	//drrVolume->RotateZ(zr);
	//drrVolume->RotateX(xr);
	//drrVolume->RotateY(yr);
	
	//ui->xr_edit->setText("0");
	//ui->yr_edit->setText("0");
	//ui->zr_edit->setText("0");


	//double scale = ui->dx_lineEdit_3->text().toDouble();
	//drrVolume->SetScale(scale/100);//缩放比例

	tdRenderer->DrawOn();

	ui->qvtkWidget->GetRenderWindow()->Render();
}

void vtkWidget::changeAngleAnytime(double a,double b,double c) {
	//calAngleFromMatrix(x, y, z);
	//x -= xe;
	//y -= ye;
	//z += ze;
	x = a;
	y = b;
	z = c;
	ui->xr_edit->setText(QString::number(a, 'g', 4));
	ui->yr_edit->setText(QString::number(b, 'g', 4));
	ui->zr_edit->setText(QString::number(c, 'g', 4));
}

void vtkWidget::refreshAngle(double a,double b,double c) {
	emit sendAngle(a,b,c);
}

void vtkWidget::create2D() {

	vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter =
		vtkSmartPointer<vtkWindowToImageFilter>::New();
	windowToImageFilter->SetInput(ui->qvtkWidget->GetRenderWindow());
	windowToImageFilter->SetInputBufferTypeToRGBA();
	windowToImageFilter->ReadFrontBufferOff();
	windowToImageFilter->Update();

	/*vtkSmartPointer<vtkBMPWriter> writer =
		vtkSmartPointer<vtkBMPWriter>::New();
	writer->SetFileName("vtk.bmp");
	writer->SetInputConnection(windowToImageFilter->GetOutputPort());
	writer->Write();
	*/

	vtkImageData* imgData = windowToImageFilter->GetOutput();
	

	int width = imgData->GetDimensions()[0];
	int height = imgData->GetDimensions()[1];

	QImage image(width, height, QImage::Format_RGB32);
	QRgb *rgbPtr =
		reinterpret_cast<QRgb *>(image.bits()) + width * (height - 1);
	unsigned char *colorsPtr =
		reinterpret_cast<unsigned char *>(imgData->GetScalarPointer());

	// Loop over the vtkImageData contents.
	for (int row = 0; row < height; row++)
	{
		for (int col = 0; col < width; col++)
		{
			// Swap the vtkImageData RGB values with an equivalent QColor
			*(rgbPtr++) = QColor(colorsPtr[0], colorsPtr[1], colorsPtr[2]).rgb();
			colorsPtr += imgData->GetNumberOfScalarComponents();
		}

		rgbPtr -= width * 2;
	}
	
	/*vtkSmartPointer<vtkImageViewer2> viewer =
		vtkSmartPointer<vtkImageViewer2>::New();
	viewer->SetInputConnection(reader->GetOutputPort());
	viewer->SetColorLevel(400);
	viewer->SetColorWindow(1000);

	vtkSmartPointer<vtkRenderWindowInteractor> rwi =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	viewer->SetupInteractor(rwi);
	rwi->Start();*/

	width = ui->dx_lineEdit->text().toInt();
	height = ui->dy_lineEdit->text().toInt();
	image = image.scaled(width, height, Qt::IgnoreAspectRatio, Qt::SmoothTransformation);

	emit sendQPixmap(image);
	emit sendRotate(x, y, z);
}

void vtkWidget::myShow() {
	show();

	if (needRefresh) {
		tdRenderer = vtkSmartPointer<vtkRenderer>::New();
		tdRenderer->SetBackground(1, 1, 1);

		ui->qvtkWidget->GetRenderWindow()->AddRenderer(tdRenderer);

		tdRenderer->AddActor(drrVolume);

		tdRenderer->GetRenderWindow()->Render();
		tdRenderer->DrawOn();

		needRefresh = false;
	}

	//setPara();
	vtkCamera * camera = tdRenderer->GetActiveCamera();
}

void vtkWidget::calAngleFromMatrix(double& xe, double& ye, double& ze) {
	vtkSmartPointer<vtkMatrix4x4> matrix = vtkSmartPointer<vtkMatrix4x4>::New();
	drrVolume->GetMatrix(matrix);
	double* tmpx;
	tmpx = drrVolume->GetOrientation();
	
	//std::cout << "vtk-orientation: "<< tmpx[0] << " " << tmpx[1] << " " << tmpx[2] << endl;
	//matrix->Print(std::cout);
	
	float sx = (float)sqrt(matrix->GetElement(2, 0)*matrix->GetElement(2, 0)
		+ matrix->GetElement(2, 2)*matrix->GetElement(2, 2));
	bool singular = sx < 1e-6;
	double xd, yd, zd;
	if (!singular) {
		xd = atan2(matrix->GetElement(2, 1), sx);
		yd = -atan2(matrix->GetElement(2, 0), matrix->GetElement(2, 2));
		zd = -atan2(matrix->GetElement(0, 1), matrix->GetElement(1, 1));
	}
	else {
		xd = atan2(-matrix->GetElement(2, 1), sx);
		if (xd > 0) {
			zd = atan2(matrix->GetElement(0, 2), matrix->GetElement(0, 0));
			yd = 0;
		}
		else {
			zd = atan2(matrix->GetElement(1, 0), matrix->GetElement(0, 0));
			yd = 0;
		}
		
	}
	xe = xd * (180.0 / acos(-1));
	ye = yd * (180.0 / acos(-1));
	ze = zd * (180.0 / acos(-1));
	std::cout << "vtk angle: " << xe << " " << ye << " " << ze << endl;
	//std::cout << "vtk matrix: " << endl;
	//std::cout << matrix->GetElement(0, 0) << " " << matrix->GetElement(0, 1) << " " << matrix->GetElement(0, 2) << endl;
	//std::cout << matrix->GetElement(1, 0) << " " << matrix->GetElement(1, 1) << " " << matrix->GetElement(1, 2) << endl;
	//std::cout << matrix->GetElement(2, 0) << " " << matrix->GetElement(2, 1) << " " << matrix->GetElement(2, 2) << endl;

}
