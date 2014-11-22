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

		// Special case: degenerate sheet
		if (n->type() == Structure::SHEET && samples_coords.empty()){
			double step = 0.25;
			for (double u = 0; u <= 1.0; u += step){
				Array1D_Vector4d row;
				for (double v = 0; v <= 1.0; v += step)
					row.push_back(Eigen::Vector4d(u,v,0,0));
				samples_coords.push_back(row);
			}
		}

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

double EvaluateCorrespondence::evaluate(Structure::ShapeGraph *shape)
{
	QVector<double> feature_vector;

	// Ignore duplicate edges
	QMap<QString, bool> seen;

	// Binary properties:
	for (auto l : shape->edges)
	{
		QString key = (l->n1->id < l->n2->id) ? (l->n1->id + l->n2->id) : (l->n2->id + l->n1->id);
		if (seen[key]) continue;
		seen[key] = true;

		auto original = l->property["orig_spokes"].value<Array1D_Vector3>();
		auto current = EvaluateCorrespondence::spokesFromLink(l);

		for (size_t i = 0; i < original.size(); i++)
		{
			double v = original[i].normalized().dot(current[i].normalized());
			feature_vector << v;
		}
	}

	Eigen::VectorXd original_v = Eigen::VectorXd::Ones(feature_vector.size());

	Eigen::VectorXd v(feature_vector.size());
	for (size_t d = 0; d < feature_vector.size(); d++) v[d] = feature_vector[d];

	return v.norm() / original_v.norm();
}
