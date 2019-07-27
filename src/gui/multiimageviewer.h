#pragma once

#include <QPushButton>
#include <QButtonGroup>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include "labelimage.h"
#include <vector>

class MultiImageViewer : public QWidget {
    Q_OBJECT

private:
    	QHBoxLayout  button_layout;
    	QButtonGroup button_group;
    	LabelImage   label_image;
    	int          button_id;

	std::vector<std::shared_ptr<cv::Mat>> images;
	std::vector<std::shared_ptr<QPixmap>> pixmaps;	

	void buttonClicked(int id) 
	{
        if ((id>=0) && (id<pixmaps.size()))
        {
            if (images[id])
            {
                if (!pixmaps[id])
                    pixmaps[id] = std::make_shared<QPixmap>(qpixmap(*images[id]));
                else
                    cvmat_to_qpixmap(images[id], pixmaps[id]);
			}
			label_image.set(*(pixmaps[id]));
		}
	}

public:
	MultiImageViewer(QWidget *parent = 0) :
       		button_group(parent), button_id(0)
	{
		QWidget* button_widget = new QWidget();
		button_widget->setLayout(&button_layout);

		QVBoxLayout *layout = new QVBoxLayout();
		layout->addWidget(button_widget);
		layout->addWidget(&label_image);
		this->setLayout(layout);

		QObject::connect(&button_group, static_cast<void(QButtonGroup::*)(int)>(&QButtonGroup::buttonClicked),
				this, &MultiImageViewer::buttonClicked);
	}

	void setButton(int id) {
		QAbstractButton* button = button_group.button(id);
		if (button) button->setChecked(true);
		buttonClicked(id);
	}

	int add(const char* label, const std::shared_ptr<cv::Mat>& image)
	{
		QPushButton* button = new QPushButton(QString(label));
		button->setCheckable(true);
		button_layout.addWidget(button);
		button_group.addButton(button, button_id);
		images.push_back(image);
		pixmaps.push_back(std::shared_ptr<QPixmap>());
		if (button_id == 0) {
			button->setChecked(true);
		}
		return button_id++;
	}

	int add(const char* label, const std::shared_ptr<QPixmap>& pixmap)
	{
		QPushButton* button = new QPushButton(QString(label));
		button->setCheckable(true);
		button_layout.addWidget(button);
		button_group.addButton(button, button_id);
		images.push_back(std::shared_ptr<cv::Mat>());
		pixmaps.push_back(pixmap);
		if (button_id == 0) {
			button->setChecked(true);
		}
		return button_id++;
	}

	const LabelImage* labelImage() const { return &label_image; }
    	LabelImage* getlabelImage() {return &label_image;}
};

