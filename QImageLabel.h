#pragma once
#pragma once
#ifndef QIMAGELABEL_H
#define QIMAGELABEL_H

#include <qlabel.h>
#include <QMouseEvent>
#include<qgraphicseffect.h>
class QImageLabel :
	public QLabel
{
	Q_OBJECT
public:
	QImageLabel(QWidget* parent = 0);
	void setImagePixmap(QPixmap map);
	void setImage(QImage img);
	QImage copyImage(int x, int y, int width, int height);
	void setSize(const QRect& rect);
	void reverseColor();
	void setMoveable(bool moveable);

protected:
	void mousePressEvent(QMouseEvent* event);
	void mouseMoveEvent(QMouseEvent* event);
	void mouseReleaseEvent(QMouseEvent* event);
	void wheelEvent(QWheelEvent* event);
	void setImageSize(int w, int h);

private:
	bool moveable;                     //设置图片是否可拖动
	bool m_move;
	QPoint m_startPoint;
	QPoint m_windowPoint;
	int action;
	QPixmap pixmap;
	QImage image;
	int imgWidth = 0;
	int imgHeight = 0;
};

#endif //QIMAGELABEL_H
