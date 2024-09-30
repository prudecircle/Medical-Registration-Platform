#include "itkDrr.h"



itkDrr::itkDrr(std::string input)
{
	isRead = false;
	reader = nullptr;

	//定义像素类型，图像类型，三维有符号数，定义指针
	typedef float PixelType;
	const unsigned int Dimension = 3;
	typedef itk::Image< PixelType, Dimension > ImageType;
	typedef itk::ImageSeriesReader< ImageType > ReaderType;

	//声明读、写 DICOM 图 像 的 itk::GDCMImageIO对象
	//itk::GDCMSeriesFileNames对象将生成并将构成所有体数据的切片的文件名进行排序
	typedef itk::GDCMImageIO ImageIOType;
	typedef itk::GDCMSeriesFileNames NamesGeneratorType;
	ImageIOType::Pointer gdcmIO = ImageIOType::New();
	NamesGeneratorType::Pointer namesGenerator = NamesGeneratorType::New();

	bool result = reader.IsNull();
	double time = (double)clock() / CLOCKS_PER_SEC;
	//设置读取路径
	//用文件名发生器生成被读的文件名和被写的文件名
	namesGenerator->SetInputDirectory(input);
	const ReaderType::FileNamesContainer& filenames = namesGenerator->GetInputFileNames();

	//设置DICOM图像IO对象和被读的文件名的列表
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetImageIO(gdcmIO);
	reader->SetFileNames(filenames);
	try
	{
		reader->Update();
	}
	catch (const itk::ExceptionObject& err)
	{
		std::string strerr = err.GetDescription();
		return;
	}

	image = reader->GetOutput();
	isRead = true;
	ImageType::SizeType sizeCT = image->GetLargestPossibleRegion().GetSize();
	InputImageType::SpacingType ctPixelSpacing = image->GetSpacing();
	cpp_object3D = (float*)malloc(sizeof(float) * sizeCT[0] * sizeCT[1] * sizeCT[2]);
	for (int i = 0; i < sizeCT[0]; i++)
		for (int j = 0; j < sizeCT[1]; j++)
			for (int k = 0; k < sizeCT[2]; k++) {
				itk::Index<3> index;
				index[0] = i;
				index[1] = j;
				index[2] = k;
				cpp_object3D[k + j * (sizeCT[2]) + i * (sizeCT[1] * sizeCT[2])] = image->GetPixel(index);
			}

	int SizeCT[3];
	SizeCT[0] = sizeCT[0];
	SizeCT[1] = sizeCT[1];
	SizeCT[2] = sizeCT[2];
	float PixelSpacingCT[3];
	PixelSpacingCT[0] = ctPixelSpacing[0];
	PixelSpacingCT[1] = ctPixelSpacing[1];
	PixelSpacingCT[2] = ctPixelSpacing[2];
	loadDICOMInGPUMemory(cpp_object3D, SizeCT, PixelSpacingCT);
	free(cpp_object3D);
	cpp_object3D = NULL;
	//std::cout << "readfile time:" << (double)clock() / CLOCKS_PER_SEC - time << "s" << std::endl;

}


itkDrr::~itkDrr()
{
}
cv::Mat itkDrr::Pos2Matrix(std::vector<double> para)
{
	// Calculate rotation about x axis
	cv::Mat R_x = (cv::Mat_<double>(3, 3) <<
		1, 0, 0,
		0, cos(para[0]), -sin(para[0]),
		0, sin(para[0]), cos(para[0])
		);
	// Calculate rotation about y axis
	cv::Mat R_y = (cv::Mat_<double>(3, 3) <<
		cos(para[1]), 0, sin(para[1]),
		0, 1, 0,
		-sin(para[1]), 0, cos(para[1])
		);
	// Calculate rotation about z axis
	cv::Mat R_z = (cv::Mat_<double>(3, 3) <<
		cos(para[2]), -sin(para[2]), 0,
		sin(para[2]), cos(para[2]), 0,
		0, 0, 1);
	// Combined rotation matrix
	cv::Mat R = R_z * R_y * R_x;
	cv::Mat translation = (cv::Mat_<double>(3, 1) << para[3], para[4], para[5]);
	
	cv::Mat A = Mat(Size(3, 4), CV_64F);
	cv::hconcat(R, translation, A);
	cv::Mat B = Mat(Size(4, 4), CV_64F);
	cv::Mat zero = (cv::Mat_<double>(1, 4) << 0, 0, 0, 1);
	
	cv::vconcat(A, zero, B);
	return B;
}

double* itkDrr::cudaUCharDRR() {
	clock_t start;
	start = clock();

	DRRParameters para;
	cout << "drr para:  ";
	for (int i = 0; i < 3; i++) {
		cout << rotation[i]<<"   ";
	}
	for (int i = 0; i < 3; i++) {
		cout << translation[i] << "   ";
	}
	cout << std::endl;
	double output_size[3];
	for (int i = 0; i < 1; i++) {
		InputImageType::PointType ori;
		ori[0] = 0.0;
		ori[1] = 0.0;
		ori[2] = 0.0;
		image->SetOrigin(ori);

		//----------------kfq 4.20
		 
		typedef itk::CenteredEuler3DTransform<double> TransformType;
		TransformType::Pointer transform = TransformType::New();
		TransformType::OutputVectorType translations;
		const double dtr = (atan(1.0) * 4.0) / 180.0;
		//transform->SetTranslation(translation);
		transform->SetRotation(dtr * rotation[0], dtr * rotation[1], dtr * rotation[2]);
		translations[0] = translation[0];
		translations[1] = translation[1];
		translations[2] = translation[2];

		transform->SetTranslation(translations);
		itk::Matrix<double> transformMat = transform->GetMatrix();
		//cout << "transform Matrix=" << transformMat << std::endl;
		itk::Vector<double> SourceWorld;

		typename InputImageType::SizeType sizeCT;
		typename InputImageType::RegionType regionCT;
		typename InputImageType::SpacingType ctPixelSpacing;
		typename InputImageType::PointType ctOrigin;

		ctPixelSpacing = image->GetSpacing();
		ctOrigin = image->GetOrigin();
		regionCT = image->GetLargestPossibleRegion();
		sizeCT = regionCT.GetSize();
		InputImageType::SizeType imSize = image->GetBufferedRegion().GetSize();

		TransformType::InputPointType isocenter;
		isocenter[0] = ctOrigin[0] + ctPixelSpacing[0] * static_cast<double>(imSize[0]) / 2.0;
		isocenter[1] = ctOrigin[1] + ctPixelSpacing[1] * static_cast<double>(imSize[1]) / 2.0;
		isocenter[2] = ctOrigin[2] + ctPixelSpacing[2] * static_cast<double>(imSize[2]) / 2.0;
		transform->SetCenter(isocenter);
		

		TransformType::Pointer m_InverseTransform = TransformType::New();
		m_InverseTransform->SetComputeZYX(true);
		TransformType::Pointer m_ComposeTransform = TransformType::New();
		m_ComposeTransform->SetComputeZYX(true);
		TransformType::Pointer m_GantryRoTransform = TransformType::New();
		m_GantryRoTransform->SetComputeZYX(true);
		m_GantryRoTransform->SetIdentity();
		TransformType::Pointer m_CamShiftTransform = TransformType::New();
		m_CamShiftTransform->SetComputeZYX(true);
		m_CamShiftTransform->SetIdentity();
		TransformType::Pointer m_CamRotTransform = TransformType::New();
		m_CamRotTransform->SetComputeZYX(true);
		m_CamRotTransform->SetIdentity();
		float dtrf = (atan(1.0) * 4.0) / 180.0;
		m_CamRotTransform->SetRotation(dtrf * (-90.0), 0., 0.);

		m_ComposeTransform->SetIdentity();
		m_ComposeTransform->Compose(transform, 0);
		float rprojection = 90.;
		float angle = -dtr * rprojection;
		m_GantryRoTransform->SetRotation(angle, 2 * angle, 0.);
		m_GantryRoTransform->SetCenter(isocenter);
		m_ComposeTransform->Compose(m_GantryRoTransform, 0);
		TransformType::OutputVectorType focalpointtranslation;
		focalpointtranslation[0] = -isocenter[0];
		//focalpointtranslation[1] = 400 - isocenter[1];
		focalpointtranslation[1] = 400 - isocenter[1];
		focalpointtranslation[2] = -isocenter[2];
		//std::cout << isocenter[0] << " " << isocenter[1] << " " << isocenter[2] << endl;
		//std::cout << focalpointtranslation[0] << " " << focalpointtranslation[1] << " " << focalpointtranslation[2] << endl;
		m_CamShiftTransform->SetTranslation(focalpointtranslation);
		m_ComposeTransform->Compose(m_CamShiftTransform, 0);
		m_ComposeTransform->Compose(m_CamRotTransform, 0);



		m_ComposeTransform->GetInverse(m_InverseTransform);

		itk::Matrix<double> m_matrix = m_InverseTransform->GetMatrix();
		itk::Vector<double, 3> m_offset = m_InverseTransform->GetOffset();
		
		std::vector<double> pose_para1 = { 0,0,0,ctOrigin[0],ctOrigin[1],ctOrigin[2] };
		std::vector<double> pose_para2 = { 0,0,0,0,0,1020};
		std::vector<double> pose_para3 = { CV_PI,0,0,dx * ctPixelSpacing[0] * 0.5f,dy * ctPixelSpacing[1] * 0.5f,0 };
		
		cv::Mat PosMat = (cv::Mat_<double>(4, 4) <<
			-0.998934, 0.033668, -0.0315872, -11.9511,
			0.03284, 0.0373304, -0.998763, -809.549,
			- 0.0324472, -0.998736, -0.0383962, -766.221,
			0, 0, 0, 1);
		cv::Mat frontal2lateralcb = (cv::Mat_<double>(4, 4) <<
			0.0121391, -0.00461972, -0.999916, 388.875,
			0.00666581, 0.999967, -0.00453903, -20.8345,
			0.999904, -0.00661015, 0.0121695, 349.604,
			0, 0, 0, 1);
		cout << "m_matrix" << m_matrix << std::endl;
		cout << "m_offset" << m_offset << std::endl;
		cv::Mat itkMat = (cv::Mat_<double>(4, 4) <<
			m_matrix[0][0], m_matrix[0][1], m_matrix[0][2], m_offset[0],
			m_matrix[1][0], m_matrix[1][1], m_matrix[1][2], m_offset[1],
			m_matrix[2][0], m_matrix[2][1], m_matrix[2][2], m_offset[2],
			0, 0, 0, 1);
		cv::Mat trans_mat = Pos2Matrix(pose_para3) * Pos2Matrix(pose_para2) * PosMat * Pos2Matrix(pose_para1);
		
		cv::invert(trans_mat, trans_mat);


		
		
		
		PointType m_SourcePoint;
		m_SourcePoint[0] = 0.;
		m_SourcePoint[1] = 0.;
		m_SourcePoint[2] = 0.;
		PointType tmp_SourceWorld;
		m_InverseTransform->TransformPoint(m_SourcePoint);
		tmp_SourceWorld = m_matrix * m_SourcePoint + m_offset;

		itk::Matrix<double> m_direction = image->GetDirection();

		//output_2d_image
		float projectRange = ctPixelSpacing[0] * sizeCT[0] > ctPixelSpacing[2] * sizeCT[2] ? ctPixelSpacing[0] * sizeCT[0] : ctPixelSpacing[2] * sizeCT[2];

		double output_spacing[3];
		output_spacing[0] = 0.3125;//sx
		output_spacing[1] = 0.3125;//sy
		//output_spacing[0] = 0.21;//sx
		//output_spacing[1] = 0.21;//sy
		//output_spacing[0] = (double)(projectRange*1.15)/1024;//sx
		//output_spacing[1] = (double)(projectRange * 1.15) / 1024;//sy
		output_spacing[2] = 1;

		output_size[0] = this->dx;//dx
		output_size[1] = this->dy;//dy
		output_size[2] = 1.0;

		itk::Vector<double> origin;
		origin[0] = -output_spacing[0] * (output_size[0] - 1.) / 2.;
		origin[1] = -output_spacing[1] * (output_size[1] - 1.) / 2.;
		origin[2] = -400;

		itk::Vector<double> corner00;
		itk::Vector<double> cornerOrigin01;
		itk::Vector<double> corner01;
		itk::Vector<double> cornerOrigin10;
		itk::Vector<double> corner10;

		cornerOrigin01[0] = origin[0] + output_spacing[0] * (output_size[0] - 1.);
		cornerOrigin01[1] = origin[1];
		cornerOrigin01[2] = origin[2];

		cornerOrigin10[0] = origin[0];
		cornerOrigin10[1] = origin[1] + output_spacing[1] * (output_size[1] - 1.);
		cornerOrigin10[2] = origin[2];

		corner00 = m_matrix * origin + m_offset;
		corner01 = m_matrix * cornerOrigin01 + m_offset;
		corner10 = m_matrix * cornerOrigin10 + m_offset;
		SourceWorld[0] = tmp_SourceWorld[0];
		SourceWorld[1] = tmp_SourceWorld[1];
		SourceWorld[2] = tmp_SourceWorld[2];

		itk::Vector<double> stepInX;
		stepInX[0] = (corner01[0] - corner00[0]) / (output_size[0]);
		stepInX[1] = (corner01[1] - corner00[1]) / (output_size[0]);
		stepInX[2] = (corner01[2] - corner00[2]) / (output_size[0]);
		itk::Vector<double> stepInY;
		stepInY[0] = (corner10[0] - corner00[0]) / (output_size[1]);
		stepInY[1] = (corner10[1] - corner00[1]) / (output_size[1]);
		stepInY[2] = (corner10[2] - corner00[2]) / (output_size[1]);


		para.stepInX[i * 3] = stepInX[0];
		para.stepInX[i * 3 + 1] = stepInX[1];
		para.stepInX[i * 3 + 2] = stepInX[2];
		para.stepInY[i * 3] = stepInY[0];
		para.stepInY[i * 3 + 1] = stepInY[1];
		para.stepInY[i * 3 + 2] = stepInY[2];
		para.corner00[i * 3] = corner00[0];
		para.corner00[i * 3 + 1] = corner00[1];
		para.corner00[i * 3 + 2] = corner00[2];
		para.SourceWorld[i * 3] = SourceWorld[0];
		para.SourceWorld[i * 3 + 1] = SourceWorld[1];
		para.SourceWorld[i * 3 + 2] = SourceWorld[2];
		para.size[0] = output_size[0];
		para.size[1] = output_size[1];
		para.threshold = this->threshold;

		Image3D image3d;
		image3d.image = cpp_object3D;
		image3d.isoCenter[0] = isocenter[0];
		image3d.isoCenter[1] = isocenter[1];
		image3d.isoCenter[2] = isocenter[2];
		image3d.PixelSpacingCT[0] = ctPixelSpacing[0];
		image3d.PixelSpacingCT[1] = ctPixelSpacing[1];
		image3d.PixelSpacingCT[2] = ctPixelSpacing[2];
		image3d.SizeCT[0] = sizeCT[0];
		image3d.SizeCT[1] = sizeCT[1];
		image3d.SizeCT[2] = sizeCT[2];
	}


	CUDAParamerters cudaPara;
	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, 0);
	cudaPara.numThreads = prop.maxThreadsPerBlock;
	cudaPara.numBlocks = (int)ceil((float)output_size[0] * output_size[1]  / cudaPara.numThreads);


	//clock_t start, end;
	//start = clock();

	loadOuputVariablesInGPUMemory((int)output_size[0], (int)output_size[1], 1, 1);
	
	calUCharDRRwithCUDA(cudaPara, para, 1,true);




	freeAuxiliaryVariablesInGPUMemory(1);
	//end = clock();
	//cout << "time of process is:" << (double)(end - start) / CLOCKS_PER_SEC << endl;


	//cv::Mat image((int)output_size[0] * PARAS_NUMS, (int)output_size[1], CV_8UC1, (void*)object2d);










	//cout << "time:" << (double)(clock() - start) / CLOCKS_PER_SEC << endl;




	return 0;
}
double* itkDrr::cudaUCharDRRLateral() {
	clock_t start;
	start = clock();

	DRRParameters para;
	cout << "drr para:  ";
	for (int i = 0; i < 3; i++) {
		cout << rotation[i] << "   ";
	}
	for (int i = 0; i < 3; i++) {
		cout << translation[i] << "   ";
	}
	cout << std::endl;
	double output_size[3];
	for (int i = 0; i < 1; i++) {
		InputImageType::PointType ori;
		ori[0] = 0.0;
		ori[1] = 0.0;
		ori[2] = 0.0;
		image->SetOrigin(ori);


		typedef itk::CenteredEuler3DTransform<double> TransformType;
		TransformType::Pointer transform = TransformType::New();
		TransformType::OutputVectorType translations;
		const double dtr = (atan(1.0) * 4.0) / 180.0;
		//transform->SetTranslation(translation);
		transform->SetRotation(dtr * rotation[0], dtr * rotation[1], dtr * rotation[2]);
		translations[0] = translation[0];
		translations[1] = translation[1];
		translations[2] = translation[2];

		transform->SetTranslation(translations);
		itk::Matrix<double> transformMat = transform->GetMatrix();
		//cout << "transform Matrix=" << transformMat << std::endl;
		itk::Vector<double> SourceWorld;

		typename InputImageType::SizeType sizeCT;
		typename InputImageType::RegionType regionCT;
		typename InputImageType::SpacingType ctPixelSpacing;
		typename InputImageType::PointType ctOrigin;

		ctPixelSpacing = image->GetSpacing();
		ctOrigin = image->GetOrigin();
		regionCT = image->GetLargestPossibleRegion();
		sizeCT = regionCT.GetSize();
		InputImageType::SizeType imSize = image->GetBufferedRegion().GetSize();

		TransformType::InputPointType isocenter;
		isocenter[0] = ctOrigin[0] + ctPixelSpacing[0] * static_cast<double>(imSize[0]) / 2.0;
		isocenter[1] = ctOrigin[1] + ctPixelSpacing[1] * static_cast<double>(imSize[1]) / 2.0;
		isocenter[2] = ctOrigin[2] + ctPixelSpacing[2] * static_cast<double>(imSize[2]) / 2.0;
		transform->SetCenter(isocenter);


		TransformType::Pointer m_InverseTransform = TransformType::New();
		m_InverseTransform->SetComputeZYX(true);
		TransformType::Pointer m_ComposeTransform = TransformType::New();
		m_ComposeTransform->SetComputeZYX(true);
		TransformType::Pointer m_GantryRoTransform = TransformType::New();
		m_GantryRoTransform->SetComputeZYX(true);
		m_GantryRoTransform->SetIdentity();
		TransformType::Pointer m_CamShiftTransform = TransformType::New();
		m_CamShiftTransform->SetComputeZYX(true);
		m_CamShiftTransform->SetIdentity();
		TransformType::Pointer m_CamRotTransform = TransformType::New();
		m_CamRotTransform->SetComputeZYX(true);
		m_CamRotTransform->SetIdentity();
		float dtrf = (atan(1.0) * 4.0) / 180.0;
		m_CamRotTransform->SetRotation(dtrf * (-90.0), 0., 0.);

		m_ComposeTransform->SetIdentity();
		m_ComposeTransform->Compose(transform, 0);
		float rprojection = 90.;
		float angle = -dtr * rprojection;
		m_GantryRoTransform->SetRotation(angle, 2 * angle, 0.);
		m_GantryRoTransform->SetCenter(isocenter);
		m_ComposeTransform->Compose(m_GantryRoTransform, 0);
		TransformType::OutputVectorType focalpointtranslation;
		focalpointtranslation[0] = -isocenter[0];
		//focalpointtranslation[1] = 400 - isocenter[1];
		focalpointtranslation[1] = 400 - isocenter[1];
		focalpointtranslation[2] = -isocenter[2];
		//std::cout << isocenter[0] << " " << isocenter[1] << " " << isocenter[2] << endl;
		//std::cout << focalpointtranslation[0] << " " << focalpointtranslation[1] << " " << focalpointtranslation[2] << endl;
		m_CamShiftTransform->SetTranslation(focalpointtranslation);
		m_ComposeTransform->Compose(m_CamShiftTransform, 0);
		m_ComposeTransform->Compose(m_CamRotTransform, 0);



		m_ComposeTransform->GetInverse(m_InverseTransform);

		itk::Matrix<double> m_matrix = m_InverseTransform->GetMatrix();
		itk::Vector<double, 3> m_offset = m_InverseTransform->GetOffset();








		cv::Mat frontalMat = (cv::Mat_<double>(4, 4) <<
			m_matrix[0][0], m_matrix[0][1], m_matrix[0][2], m_offset[0],
			m_matrix[1][0], m_matrix[1][1], m_matrix[1][2], m_offset[1],
			m_matrix[2][0], m_matrix[2][1], m_matrix[2][2], m_offset[2],
			0, 0, 0, 1);
		cv::Mat transMat = (cv::Mat_<double>(4, 4) <<
			-0.025631840169871, 0.9995580007914995, -0.01506815649542876, 49.54787381914251,
			-0.99705852603129, -0.02665099793168167, -0.0718546052672259, 186.012167498308,
			-0.07222448305540932, 0.01318214416459025, 0.9973007657790259, -11.09006937237569,
			0, 0, 0, 1);
		cv::Mat lateralMat = transMat * frontalMat;

		m_matrix[0][0] = lateralMat.at<double>(0, 0);
		m_matrix[0][1] = lateralMat.at<double>(0, 1);
		m_matrix[0][2] = lateralMat.at<double>(0, 2);
		m_matrix[1][0] = lateralMat.at<double>(1, 0);
		m_matrix[1][1] = lateralMat.at<double>(1, 1);
		m_matrix[1][2] = lateralMat.at<double>(1, 2);
		m_matrix[2][0] = lateralMat.at<double>(2, 0);
		m_matrix[2][1] = lateralMat.at<double>(2, 1);
		m_matrix[2][2] = lateralMat.at<double>(2, 2);
		m_offset[0] = lateralMat.at<double>(0, 3);
		m_offset[1] = lateralMat.at<double>(1, 3);
		m_offset[2] = lateralMat.at<double>(2, 3);
		






		PointType m_SourcePoint;
		m_SourcePoint[0] = 0.;
		m_SourcePoint[1] = 0.;
		m_SourcePoint[2] = 0.;
		PointType tmp_SourceWorld;
		m_InverseTransform->TransformPoint(m_SourcePoint);
		tmp_SourceWorld = m_matrix * m_SourcePoint + m_offset;

		itk::Matrix<double> m_direction = image->GetDirection();

		//output_2d_image
		float projectRange = ctPixelSpacing[0] * sizeCT[0] > ctPixelSpacing[2] * sizeCT[2] ? ctPixelSpacing[0] * sizeCT[0] : ctPixelSpacing[2] * sizeCT[2];

		double output_spacing[3];
		output_spacing[0] = 0.3125;//sx
		output_spacing[1] = 0.3125;//sy
		//output_spacing[0] = 0.21;//sx
		//output_spacing[1] = 0.21;//sy
		//output_spacing[0] = (double)(projectRange*1.15)/1024;//sx
		//output_spacing[1] = (double)(projectRange * 1.15) / 1024;//sy
		output_spacing[2] = 1;

		output_size[0] = this->dx;//dx
		output_size[1] = this->dy;//dy
		output_size[2] = 1.0;

		itk::Vector<double> origin;
		origin[0] = -output_spacing[0] * (output_size[0] - 1.) / 2.;
		origin[1] = -output_spacing[1] * (output_size[1] - 1.) / 2.;
		origin[2] = -400;

		itk::Vector<double> corner00;
		itk::Vector<double> cornerOrigin01;
		itk::Vector<double> corner01;
		itk::Vector<double> cornerOrigin10;
		itk::Vector<double> corner10;

		cornerOrigin01[0] = origin[0] + output_spacing[0] * (output_size[0] - 1.);
		cornerOrigin01[1] = origin[1];
		cornerOrigin01[2] = origin[2];

		cornerOrigin10[0] = origin[0];
		cornerOrigin10[1] = origin[1] + output_spacing[1] * (output_size[1] - 1.);
		cornerOrigin10[2] = origin[2];

		corner00 = m_matrix * origin + m_offset;
		corner01 = m_matrix * cornerOrigin01 + m_offset;
		corner10 = m_matrix * cornerOrigin10 + m_offset;
		SourceWorld[0] = tmp_SourceWorld[0];
		SourceWorld[1] = tmp_SourceWorld[1];
		SourceWorld[2] = tmp_SourceWorld[2];

		itk::Vector<double> stepInX;
		stepInX[0] = (corner01[0] - corner00[0]) / (output_size[0]);
		stepInX[1] = (corner01[1] - corner00[1]) / (output_size[0]);
		stepInX[2] = (corner01[2] - corner00[2]) / (output_size[0]);
		itk::Vector<double> stepInY;
		stepInY[0] = (corner10[0] - corner00[0]) / (output_size[1]);
		stepInY[1] = (corner10[1] - corner00[1]) / (output_size[1]);
		stepInY[2] = (corner10[2] - corner00[2]) / (output_size[1]);


		para.stepInX[i * 3] = stepInX[0];
		para.stepInX[i * 3 + 1] = stepInX[1];
		para.stepInX[i * 3 + 2] = stepInX[2];
		para.stepInY[i * 3] = stepInY[0];
		para.stepInY[i * 3 + 1] = stepInY[1];
		para.stepInY[i * 3 + 2] = stepInY[2];
		para.corner00[i * 3] = corner00[0];
		para.corner00[i * 3 + 1] = corner00[1];
		para.corner00[i * 3 + 2] = corner00[2];
		para.SourceWorld[i * 3] = SourceWorld[0];
		para.SourceWorld[i * 3 + 1] = SourceWorld[1];
		para.SourceWorld[i * 3 + 2] = SourceWorld[2];
		para.size[0] = output_size[0];
		para.size[1] = output_size[1];
		para.threshold = this->threshold;

		Image3D image3d;
		image3d.image = cpp_object3D;
		image3d.isoCenter[0] = isocenter[0];
		image3d.isoCenter[1] = isocenter[1];
		image3d.isoCenter[2] = isocenter[2];
		image3d.PixelSpacingCT[0] = ctPixelSpacing[0];
		image3d.PixelSpacingCT[1] = ctPixelSpacing[1];
		image3d.PixelSpacingCT[2] = ctPixelSpacing[2];
		image3d.SizeCT[0] = sizeCT[0];
		image3d.SizeCT[1] = sizeCT[1];
		image3d.SizeCT[2] = sizeCT[2];
	}


	CUDAParamerters cudaPara;
	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, 0);
	cudaPara.numThreads = prop.maxThreadsPerBlock;
	cudaPara.numBlocks = (int)ceil((float)output_size[0] * output_size[1]  / cudaPara.numThreads);


	//clock_t start, end;
	//start = clock();

	loadOuputVariablesInGPUMemory((int)output_size[0], (int)output_size[1], 1, 1);
	calUCharDRRwithCUDA(cudaPara, para, 1,false);




	freeAuxiliaryVariablesInGPUMemory(1);
	//end = clock();
	//cout << "time of process is:" << (double)(end - start) / CLOCKS_PER_SEC << endl;


	//cv::Mat image((int)output_size[0] * PARAS_NUMS, (int)output_size[1], CV_8UC1, (void*)object2d);


	return 0;
}

void itkDrr::setParameters(std::vector<double> parameters, int dx, int dy, float threshold)
{
	clock_t start = clock();
	for (int m = 0; m < 1; m++) {
		vector<double> sub_parameters(parameters.begin(), parameters.begin() + 6);

		for (int i = 0; i < parameters.size(); i++) {
			if (i % 6 < 3) {
				rotation.push_back(parameters[i]);
			}
			else
			{
				translation.push_back(parameters[i]);
			}
		}

		if (rotation.size() != translation.size()) {
			cout << "amount of rotation doesn't match translation!" << endl;
			
		}
		if (rotation.size() > 12) {
			cout << "too many parameters!" << endl;
			
		}
		if (rotation.size() % 3 != 0) {
			cout << "missing parameters!" << endl;
		}
	}
	this->dx = dx;
	this->dy = dy;
	this->threshold = threshold;
}

void itkDrr::setParameters(std::vector<double> parameters)
{
	
	rotation.clear();
	translation.clear();
	for (int m = 0; m < 1; m++) {
		//vector<double> sub_parameters(parameters.begin(), parameters.begin() + 6);

		for (int i = 0; i < parameters.size(); i++) {
			if (i % 6 < 3) {
				rotation.push_back(parameters[i]);
			}
			else
			{
				translation.push_back(parameters[i]);
			}
		}

		if (rotation.size() != translation.size()) {
			cout << "amount of rotation doesn't match translation!" << endl;

		}
		if (rotation.size() > 12) {
			cout << "too many parameters!" << endl;

		}
		if (rotation.size() % 3 != 0) {
			cout << "missing parameters!" << endl;
		}
	}
	//cout <<"translation========" << translation[0] << std::endl;
}

//void itkDrr::setInputName(std::string input) {
//	// 重新导入CT序列后需要对reader缓存进行清除
//	//reader = nullptr;
//	//isRead = false;
//
//	////定义像素类型，图像类型，三维有符号数，定义指针
//	//typedef signed short PixelType;
//	//const unsigned int Dimension = 3;
//	//typedef itk::Image< PixelType, Dimension > ImageType;
//	//typedef itk::ImageSeriesReader< ImageType > ReaderType;
//
//	////声明读、写 DICOM 图 像 的 itk::GDCMImageIO对象
//	////itk::GDCMSeriesFileNames对象将生成并将构成所有体数据的切片的文件名进行排序
//	//typedef itk::GDCMImageIO ImageIOType;
//	//typedef itk::GDCMSeriesFileNames NamesGeneratorType;
//	//ImageIOType::Pointer gdcmIO = ImageIOType::New();
//	//NamesGeneratorType::Pointer namesGenerator = NamesGeneratorType::New();
//
//	//if (!isRead) {
//	//	bool result = reader.IsNull();
//	//	double time = (double)clock() / CLOCKS_PER_SEC;
//	//	//设置读取路径
//	//	//用文件名发生器生成被读的文件名和被写的文件名
//	//	namesGenerator->SetInputDirectory(input);
//	//	const ReaderType::FileNamesContainer& filenames = namesGenerator->GetInputFileNames();
//
//	//	//设置DICOM图像IO对象和被读的文件名的列表
//	//	ReaderType::Pointer reader = ReaderType::New();
//	//	reader->SetImageIO(gdcmIO);
//	//	reader->SetFileNames(filenames);
//	//	try
//	//	{
//	//		reader->Update();
//	//	}
//	//	catch (const itk::ExceptionObject & err)
//	//	{
//	//		std::string strerr = err.GetDescription();
//	//		return;
//	//	}
//
//	//	image = reader->GetOutput();
//	//	isRead = true;
//	//	//qDebug() << "itk read complete" << endl;
//	//	//std::cout << "itk readfile :" << (double)clock() / CLOCKS_PER_SEC - time << std::endl;
//	//}
//
//	//-----------------
//
//	
//
//	//ImageType::SizeType sizeCT = image->GetLargestPossibleRegion().GetSize();
//	//InputImageType::SpacingType ctPixelSpacing = image->GetSpacing();
//	//free(cpp_object3D);
//}

vtkDataObject * itkDrr::getOutput()
{
	return ivfilter->GetOutput();
}

//void itkDrr::imageProcess(int count, ...)
//{
//
//	bool ok;
//	va_list arg_ptr;
//	va_start(arg_ptr, count);
//	//va_arg(arg_ptr, int);
//	// 轴变换
//	float rx = va_arg(arg_ptr, double);
//	float ry = va_arg(arg_ptr, double);
//	float rz = va_arg(arg_ptr, double);
//	// 相机变换参数
//	float tx = va_arg(arg_ptr, double);
//	float ty = va_arg(arg_ptr, double);
//	float tz = va_arg(arg_ptr, double);
//	//相对物体的旋转中心
//	float cx = 0;
//	float cy = 0;
//	float cz = 0;
//	//float sid = 400.;
//	float length = va_arg(arg_ptr, double); // 源像距
//	float distance = 0.0;
//	int dx = va_arg(arg_ptr, int);//int dx = 501; //像素数
//	int dy = va_arg(arg_ptr, int);//int dy = 501;
//
//	
//	float sx = va_arg(arg_ptr, double);//像素间距
//	float sy = va_arg(arg_ptr, double);
//
//
//	//2d投影位置 无需修改
//	float o2Dx = 0;
//	float o2Dy = 0;
//	//阈值
//	double threshold = va_arg(arg_ptr, double);//double threshold = 0;
//
//	int ww = va_arg(arg_ptr, int);	//窗宽
//	int wl = va_arg(arg_ptr, int);	//窗位
//
//	va_end(arg_ptr);
//
//	// Software Guide : BeginLatex
//	//
//	// Although we generate a 2D projection of the 3D volume for the
//	// purposes of the interpolator both images must be three dimensional.
//	//
//	// Software Guide : EndLatex
//
//	// Software Guide : BeginCodeSnippet
//	constexpr unsigned int Dimension = 3;
//	using InputPixelType = short;
//	using OutputPixelType = unsigned char;
//
//	using InputImageType = itk::Image<InputPixelType, Dimension>;
//	using OutputImageType = itk::Image<OutputPixelType, Dimension>;
//
//	// Software Guide : EndCodeSnippet
//
//	// Software Guide : BeginLatex
//	//
//	// For the purposes of this example we assume the input volume has
//	// been loaded into an \code{itk::Image image}.
//	//
//	// Software Guide : EndLatex
//
//	if (!isRead)
//	{
//		// 需要读文件
//		return;
//	}
//
//	// Software Guide : BeginLatex
//	//
//	// Creation of a \code{ResampleImageFilter} enables coordinates for
//	// each of the pixels in the DRR image to be generated. These
//	// coordinates are used by the \code{RayCastInterpolateImageFunction}
//	// to determine the equation of each corresponding ray which is cast
//	// through the input volume.
//	//
//	// Software Guide : EndLatex
//
//	// Software Guide : BeginCodeSnippet
//	using FilterType = itk::ResampleImageFilter<InputImageType, InputImageType>;
//
//	FilterType::Pointer filter = FilterType::New();
//
//	filter->SetInput(image);
//	filter->SetDefaultPixelValue(0);
//	// Software Guide : EndCodeSnippet
//
//	// Software Guide : BeginLatex
//	//
//	// An Euler transformation is defined to position the input volume.
//	// The \code{ResampleImageFilter} uses this transform to position the
//	// output DRR image for the desired view.
//	//
//	// Software Guide : EndLatex
//
//	// Software Guide : BeginCodeSnippet
//	using TransformType = itk::CenteredEuler3DTransform<double>;
//
//	TransformType::Pointer transform = TransformType::New();
//
//	//transform->SetComputeZYX(true);
//
//	TransformType::OutputVectorType translation;
//
//	translation[0] = tx;
//	translation[1] = ty;
//	translation[2] = tz;
//
//	// constant for converting degrees into radians
//	const double dtr = (std::atan(1.0) * 4.0) / 180.0;
//
//	
//	transform->SetTranslation(translation);
//	transform->SetRotation(dtr * rx, dtr * ry, dtr * rz);
//
//
//	InputImageType::PointType   imOrigin = image->GetOrigin();
//	InputImageType::SpacingType imRes = image->GetSpacing();
//
//	using InputImageRegionType = InputImageType::RegionType;
//	using InputImageSizeType = InputImageRegionType::SizeType;
//
//	InputImageRegionType imRegion = image->GetBufferedRegion();
//	InputImageSizeType   imSize = imRegion.GetSize();
//
//	imOrigin[0] += imRes[0] * static_cast<double>(imSize[0]) / 2.0;
//	imOrigin[1] += imRes[1] * static_cast<double>(imSize[1]) / 2.0;
//	imOrigin[2] += imRes[2] * static_cast<double>(imSize[2]) / 2.0;
//
//	TransformType::InputPointType center;
//	center[0] = cx + imOrigin[0];
//	center[1] = cy + imOrigin[1];
//	center[2] = cz + imOrigin[2];
//
//	transform->SetCenter(center);
//
//	// Software Guide : EndCodeSnippet
//
//	// Software Guide : BeginLatex
//	//
//	// The \code{RayCastInterpolateImageFunction} is instantiated and passed the
//	// transform object. The \code{RayCastInterpolateImageFunction} uses this
//	// transform to reposition the x-ray source such that the DRR image
//	// and x-ray source move as one around the input volume. This coupling
//	// mimics the rigid geometry of the x-ray gantry.
//	//
//	// Software Guide : EndLatex
//
//	// Software Guide : BeginCodeSnippet
//	using InterpolatorType =
//		itk::RayCastInterpolateImageFunction<InputImageType, double>;
//	InterpolatorType::Pointer interpolator = InterpolatorType::New();
//	interpolator->SetTransform(transform);
//	// Software Guide : EndCodeSnippet
//
//	// Software Guide : BeginLatex
//	//
//	// We can then specify a threshold above which the volume's
//	// intensities will be integrated.
//	//
//	// Software Guide : EndLatex
//
//	// Software Guide : BeginCodeSnippet
//	interpolator->SetThreshold(threshold);
//	// Software Guide : EndCodeSnippet
//
//	// Software Guide : BeginLatex
//	//
//	// The ray-cast interpolator needs to know the initial position of the
//	// ray source or focal point. In this example we place the input
//	// volume at the origin and halfway between the ray source and the
//	// screen. The distance between the ray source and the screen
//	// is the "source to image distance" \code{sid} and is specified by
//	// the user.
//	//
//	// Software Guide : EndLatex
//
//	// Software Guide : BeginCodeSnippet
//	InterpolatorType::InputPointType focalpoint;
//
//	focalpoint[0] = center[0];//imOrigin[0];
//	focalpoint[1] = center[1];// imOrigin[1];
//	focalpoint[2] = length - distance; // imOrigin[2] - sid / 2.;
//	std::cout << "focalpoint: " << focalpoint[0] << " " << focalpoint[1] << " " << focalpoint[2] << endl;
//
//	interpolator->SetFocalPoint(focalpoint);
//	// Software Guide : EndCodeSnippet
//
//	// Software Guide : BeginLatex
//	//
//	// Having initialised the interpolator we pass the object to the
//	// resample filter.
//	//
//	// Software Guide : EndLatex
//
//	// Software Guide : BeginCodeSnippet
//	//interpolator->Print(std::cout);
//	
//	filter->SetInterpolator(interpolator);
//	filter->SetTransform(transform);
//	// Software Guide : EndCodeSnippet
//
//	// Software Guide : BeginLatex
//	//
//	// The size and resolution of the output DRR image is specified via the
//	// resample filter.
//	//
//	// Software Guide : EndLatex
//
//	// Software Guide : BeginCodeSnippet
//
//	// setup the scene
//	InputImageType::SizeType size;
//
//	size[0] = dx; // number of pixels along X of the 2D DRR image
//	size[1] = dy; // number of pixels along Y of the 2D DRR image
//	size[2] = 1;  // only one slice
//
//	filter->SetSize(size);
//
//	InputImageType::SpacingType spacing;
//
//	spacing[0] = sx;  // pixel spacing along X of the 2D DRR image [mm]
//	spacing[1] = sy;  // pixel spacing along Y of the 2D DRR image [mm]
//	spacing[2] = 1.0; // slice thickness of the 2D DRR image [mm]
//	std::cout << "spacing: " << spacing[0] << " " << spacing[1] << " " << spacing[2] << endl;
//
//	filter->SetOutputSpacing(spacing);
//
//	// Software Guide : EndCodeSnippet
//
//	// Software Guide : BeginLatex
//	//
//	// In addition the position of the DRR is specified. The default
//	// position of the input volume, prior to its transformation is
//	// half-way between the ray source and screen and unless specified
//	// otherwise the normal from the "screen" to the ray source passes
//	// directly through the centre of the DRR.
//	//
//	// Software Guide : EndLatex
//
//	// Software Guide : BeginCodeSnippet
//
//	double origin[Dimension];
//
//	origin[0] = center[0] + o2Dx - sx * ((double)dx ) / 2.;
//	origin[1] = center[1] + o2Dy - sy * ((double)dy ) / 2.;
//	origin[2] = imOrigin[2] + length / 2.;
//	std::cout << "origin: " << origin[0] << " " << origin[1] << " " << origin[2] << endl;
//
//	filter->SetOutputOrigin(origin);
//	filter->Update();
//	// Software Guide : EndCodeSnippet
//
//	// create writer
//
//		// Software Guide : BeginLatex
//		//
//		// The output of the resample filter can then be passed to a writer to
//		// save the DRR image to a file.
//		//
//		// Software Guide : EndLatex
//
//		// Software Guide : BeginCodeSnippet
//	itk::GiplImageIOFactory::RegisterOneFactory();
//	using RescaleFilterType =
//		itk::RescaleIntensityImageFilter<InputImageType, OutputImageType>;
//	RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
//	rescaler->SetOutputMinimum(0);
//	rescaler->SetOutputMaximum(255);
//	rescaler->SetInput(filter->GetOutput());
//	rescaler->Update();
//
//	//窗宽窗位
//	typedef itk::IntensityWindowingImageFilter <OutputImageType, OutputImageType> IntensityWindowingImageFilterType;
//	IntensityWindowingImageFilterType::Pointer intensityFilter = IntensityWindowingImageFilterType::New();
//	intensityFilter->SetInput(rescaler->GetOutput());
//	intensityFilter->SetWindowLevel(ww, wl);	//修改窗宽窗位
//	intensityFilter->SetOutputMinimum(0);
//	intensityFilter->SetOutputMaximum(255);
//	intensityFilter->Update();
//
//	///////////////
//	double time = (double)clock() / CLOCKS_PER_SEC;
//	ivfilter = itk::ImageToVTKImageFilter<OutputImageType>::New();
//	ivfilter->SetInput(intensityFilter->GetOutput());
//	ivfilter->Update();
//	std::cout << "during itk :" << (double)clock() / CLOCKS_PER_SEC - time << std::endl;
//	///////////////
//	// Software Guide : EndCodeSnippet
//	
//
//
//	// 颜色取反
//	
//	vtkImageIterator<unsigned char> iter(ivfilter->GetOutput(), ivfilter->GetOutput()->GetExtent());
//	while (!iter.IsAtEnd()) {
//		unsigned char* inSI = iter.BeginSpan();
//		unsigned char* inSIEnd = iter.EndSpan();
//		while (inSI != inSIEnd) {
//			*inSI = 255 - *inSI;
//			++inSI;
//		}
//		iter.NextSpan();
//	}
//	
//
//
//	return;
//}
