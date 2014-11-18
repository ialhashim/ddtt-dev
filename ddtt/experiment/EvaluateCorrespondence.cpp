#include "EvaluateCorrespondence.h"

Q_DECLARE_METATYPE(Vector3);
Q_DECLARE_METATYPE(Array2D_Vector4d);

void EvaluateCorrespondence::prepare(Structure::ShapeGraph * shape)
{
	double smallest_length = DBL_MAX;
	for (auto n : shape->nodes) smallest_length = std::min(smallest_length, n->area());

	int density_count = 4;
	double resolution = smallest_length / density_count;

	shape->property["sampling_resolution"].setValue(resolution);

	for (auto n : shape->nodes)
	{
		Array2D_Vector4d samples_coords = n->discretizedPoints(resolution);
		n->property["samples_coords"].setValue(samples_coords);

		n->property["orig_diagonal"].setValue(n->diagonal());
		n->property["orig_start"].setValue(n->startPoint());
	}

	shape->property["corrEvalReady"].setValue(true);
}

void EvaluateCorrespondence::evaluate(QStringList fixedNodes, Structure::ShapeGraph *shape)
{
	for (auto n : shape->nodes)
	{
		Array2D_Vector4d samples_coords = n->property["samples_coords"].value<Array2D_Vector4d>();

		
	}
}
