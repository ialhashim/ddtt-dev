#include "EvaluateCorrespondence.h"

#include "hausdorff.h"

Array1D_Vector3 EvaluateCorrespondence::spokesFromLink(Structure::Link * l)
{
	auto samples1 = l->n1->property["samples_coords"].value<Array2D_Vector4d>();
	auto samples2 = l->n2->property["samples_coords"].value<Array2D_Vector4d>();

	if (samples1.empty()) samples1 = EvaluateCorrespondence::sampleNode(l->n1, 0);
	if (samples2.empty()) samples1 = EvaluateCorrespondence::sampleNode(l->n2, 0);

	Array1D_Vector3 spokes;
	for (auto rowi : samples1) for (auto ci : rowi)
		for (auto rowj : samples2) for (auto cj : rowj)
			spokes.push_back(l->n1->position(ci) - l->n2->position(cj));

	return spokes;
}

Array2D_Vector4d EvaluateCorrespondence::sampleNode(Structure::Node * n, double resolution)
{
	if (resolution == 0) resolution = n->area() / 3;

	Array2D_Vector4d samples_coords = n->discretizedPoints(resolution);

	// Special case: degenerate sheet
	if (n->type() == Structure::SHEET && samples_coords.empty()){
		double step = 0.25;
		for (double u = 0; u <= 1.0; u += step){
			Array1D_Vector4d row;
			for (double v = 0; v <= 1.0; v += step)
				row.push_back(Eigen::Vector4d(u, v, 0, 0));
			samples_coords.push_back(row);
		}
	}

	n->property["samples_coords"].setValue(samples_coords);

	n->property["orig_diagonal"].setValue(n->diagonal());
	n->property["orig_start"].setValue(n->startPoint());

	return samples_coords;
}

void EvaluateCorrespondence::prepare(Structure::ShapeGraph * shape)
{
	double sum_length = 0;
	for (auto n : shape->nodes) sum_length += n->area();
	double avg_length = sum_length / shape->nodes.size();

	int density_count = 3;
	double resolution = avg_length / density_count;

	shape->property["sampling_resolution"].setValue(resolution);

	for (auto n : shape->nodes)
		EvaluateCorrespondence::sampleNode(n, resolution);

	// Sample curve/sheet
	for (auto l : shape->edges)
		l->property["orig_spokes"].setValue(EvaluateCorrespondence::spokesFromLink(l));
}

double EvaluateCorrespondence::evaluate(Structure::ShapeGraph *shape)
{
	QVector<double> feature_vector;

	// Binary properties:
	for (auto l : shape->edges)
	{
		auto original = l->property["orig_spokes"].value<Array1D_Vector3>();
		auto current = EvaluateCorrespondence::spokesFromLink(l);

		// Scale by difference in closest distance of before / after
		auto minDist = [](Array1D_Vector3& pts){double m = DBL_MAX; for (auto& p : pts) m = std::min(m, p.norm()); return m; };
		double s_orig = minDist(original);
		double s_curr = minDist(current);
		double scale = std::min(s_orig,s_curr) / std::max(s_orig, s_curr);

		bool isAssignedNull = l->n1->property["isAssignedNull"].toBool() ||
							  l->n2->property["isAssignedNull"].toBool();

		for (size_t i = 0; i < original.size(); i++)
		{
			double v = original[i].normalized().dot(current[i].normalized());

			if (isAssignedNull) v = 0;

			// Penalize bad values e.g. nan
			if (!std::isfinite(v)) v = 0;

			// Scaling
			v *= scale;

			feature_vector << v;
		}
	}

	Eigen::VectorXd original_v = Eigen::VectorXd::Ones(feature_vector.size());

	Eigen::VectorXd v(feature_vector.size());
	for (size_t d = 0; d < feature_vector.size(); d++) v[d] = feature_vector[d];

	double v_norm = v.norm(), original_norm = original_v.norm();

	// Value of 1 means perfect match, anything lower is bad
	double similarity = std::min(v_norm, original_norm) / std::max(v_norm, original_norm);
	double cost = 1.0 - similarity;

	return cost;
}

QMap<QString, NanoKdTree*> EvaluateCorrespondence::kdTreesNodes(Structure::ShapeGraph * shape)
{
	QMap < QString, NanoKdTree* > result;

	auto buildKdTree = [](Structure::Node * n, const Array2D_Vector4d & coords){
		auto t = new NanoKdTree;
		for (auto row : coords) for(auto coord : row) t->addPoint(n->position(coord));
		t->build();
		return t;
	};

	for (auto n : shape->nodes)
	{
		auto coords = n->property["samples_coords"].value<Array2D_Vector4d>();
		if (coords.empty()) coords = EvaluateCorrespondence::sampleNode(n, 0);
		result[n->id] = buildKdTree(n, coords);
	}

	return result;
}

QMap<QString, QMap<QString, double> > EvaluateCorrespondence::hausdroffDistance(Structure::ShapeGraph * shapeA, Structure::ShapeGraph * shapeB)
{
	QMap<QString, QMap<QString, double> > result;

	// Initialize search acceleration
	auto trees_a = EvaluateCorrespondence::kdTreesNodes(shapeA);
	auto trees_b = EvaluateCorrespondence::kdTreesNodes(shapeB);

	for (auto nA : shapeA->nodes)
	{
		for (auto nB : shapeB->nodes)
		{
			double distance = hausdroff::distance(trees_a[nA->id], trees_b[nB->id]);
			result[nA->id][nB->id] = distance;
		}
	}

	// Clean up
	for (auto t : trees_a) delete t;
	for (auto t : trees_b) delete t;

	return result;
}
