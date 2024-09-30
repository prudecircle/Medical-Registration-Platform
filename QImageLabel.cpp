#include "QImageLabel.h"

#include "vtkWidget.h"
#include<QDebug>
QImageLabel::QImageLabel(QWidget* p) : QLabel(p) {
	QGraphicsOpacityEffect* eff = new QGraphicsOpacityEffect;
	this->setScaledContents(true);
	eff->setOpacity(0.7);
	this->setGraphicsEffect(eff);
	
}

void QImageLabel::setMoveable(bool flag) {
	moveable = flag;
}

void QImageLabel::mousePressEvent(QMouseEvent* event)
{
	if (event->button() == Qt::LeftButton)
	{
		m_move = true;
		m_startPoint = event->globalPos();
		m_windowPoint = this->frameGeometry().topLeft();
	}
}

void QImageLabel::mouseMoveEvent(QMouseEvent* event)
{
	if (event->buttons() & Qt::LeftButton)
	{
		//移动中的鼠标位置相对于初始位置的相对位置.
		QPoint relativePos = event->globalPos() - m_startPoint;
		//然后移动窗体即可.
		this->move(m_windowPoint + relativePos);
	}
}

void QImageLabel::mouseReleaseEvent(QMouseEvent* event)
{
	if (event->button() == Qt::LeftButton)
	{
		//改变移动状态.
		m_move = false;
	}
}

void QImageLabel::wheelEvent(QWheelEvent* event) {
	if (event->delta() > 0) {
		setImageSize(imgWidth * 1.1, imgHeight * 1.1);
		QImage tmp = image.scaled(imgWidth, imgHeight, Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
		this->resize(imgWidth, imgHeight);
		this->setPixmap(QPixmap::fromImage(tmp));
	}
	else if (event->delta() < 0) {
		setImageSize(imgWidth * 0.9, imgHeight * 0.9);
		QImage tmp = image.scaled(imgWidth, imgHeight, Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
		this->resize(imgWidth, imgHeight);
		this->setPixmap(QPixmap::fromImage(tmp));
	}
}

void QImageLabel::setImagePixmap(QPixmap map)
{
	pixmap = map;
	imgWidth = map.width();
	imgHeight = map.height();
	this->setPixmap(map);
}

void QImageLabel::setImage(QImage img)
{
	image = img;
	imgWidth = img.width();
	imgHeight = img.height();
	this->setPixmap(QPixmap::fromImage(img));
	this->resize(imgWidth, imgHeight);
}

void QImageLabel::reverseColor() {
	unsigned char* data = image.bits();
	int totalBytes = image.byteCount();
	for (int i = 0; i < totalBytes; i++) {
		*data = 255 - *data;
		++data;
	}
	image.save("tmp.bmp");		// 如果不缓存 修改透明度后会有异常
	image = QImage("tmp.bmp");
	QImage tmp = image.scaled(imgWidth, imgHeight, Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
	this->resize(imgWidth, imgHeight);
	this->setPixmap(QPixmap::fromImage(tmp));
}

void QImageLabel::setImageSize(int w, int h)
{
	imgWidth = w;
	imgHeight = h;
}

QImage QImageLabel::copyImage(int x, int y, int width, int height) {
	QImage tmp = image.scaled(imgWidth, imgHeight, Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
	QImage copyimg = tmp.copy(x, y, width, height);
	return copyimg;
}

void QImageLabel::setSize(const QRect& rect) {
	setGeometry(rect);
	setImageSize(rect.width(), rect.height());
	QImage tmp = image.scaled(imgWidth, imgHeight, Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
	this->resize(imgWidth, imgHeight);
	this->setPixmap(QPixmap::fromImage(tmp));
}