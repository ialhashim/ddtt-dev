#include "EvaluateCorrespondence.h"

Q_DECLARE_METATYPE(Vector3);
Q_DECLARE_METATYPE(Array2D_Vector4d);

Array1D_Vector3 EvaluateCorrespondence::spokesFromLink(Structure::Link * l)
{
	auto samples1 = l->n1->property["samples_coords"].value<Array2D_Vector4d>();
	auto samples2 = l->n2->property["samples_coords"].value<Array2D_Vector4d>();

	Array1D_Vector3 spokes;
	for (auto rowi : samples1) for (auto ci : rowi)
		for (auto rowj : samples2) for (auto cj : rowj)
			spokes.push_back(l->n1->position(ci) - l->n2->position(cj));

	return spokes;
}

Array1D_Vector3 EvaluateCorrespondence::spokesFromShape(Structure::ShapeGraph * shape)
{
	Array1D_Vector3 allspokes;
	QMap<QString, bool> seen;
	for (auto l : shape->edges)
	{
		// Ignore duplicate edges
		QString key = (l->n1->id < l->n2->id) ? (l->n1->id + l->n2->id) : (l->n2->id + l->n1->id);
		if (seen[key]) continue;
		seen[key] = true;

		for (auto s : EvaluateCorrespondence::spokesFromLink(l)) allspokes.push_back(s);
	}
	return allspokes;
}

void EvaluateCorrespondence::prepare(Structure::ShapeGraph * shape)
{
	double smallest_length = DBL_MAX;
	for (auto n : shape->nodes) smallest_length = std::min(smallest_length, n->area());

	double sum_length = 0;
	for (auto n : shape->nodes) sum_length += n->area();
	double avg_length = sum_length / shape->nodes.size();

	int density_count = 3;
	double resolution = avg_length / density_count;

	shape->property["sampling_resolution"].setValue(resolution);

	// Unary properties:
	for (auto n : shape->nodes)
	{
		Array2D_Vector4d samples_coords = n->discretizedPoints(resolution);
		n->property["samples_coords"].setValue(samples_coords);

		n->property["orig_diagonal"].setValue(n->diagonal());
		n->property["orig_start"].setValue(n->startPoint());
	}

	// Ignore duplicate edges
	QMap<QString, bool> seen;

	// Binary properties:
	for (auto l : shape->edges)
	{
		QString key = (l->n1->id < l->n2->id) ? (l->n1->id + l->n2->id) : (l->n2->id + l->n1->id);
		if (seen[key]) continue;
		seen[key] = true;
		l->property["orig_spokes"].setValue(EvaluateCorrespondence::spokesFromLink(l));
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
