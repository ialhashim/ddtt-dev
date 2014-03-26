#pragma once

#include <QGraphicsProxyWidget>

struct DeformationPath;

class DeformPathItemWidget : public QGraphicsProxyWidget
{
public:
	DeformPathItemWidget(DeformationPath * usedPath = NULL);

	DeformationPath * path;
};
