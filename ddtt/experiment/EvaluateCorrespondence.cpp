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

double EvaluateCorrespondence::evaluate(Energy::SearchNode * searchNode)
{
	QVector<double> feature_vector;

	auto shape = searchNode->shapeA.data();
	auto targetShape = searchNode->shapeB.data();

	// Binary properties:
	for (auto l : shape->edges)
	{
		auto original = l->property["orig_spokes"].value<Array1D_Vector3>();
		auto current = EvaluateCorrespondence::spokesFromLink(l);

		/// Connection: scale by difference in closeness distance
		double connection_weight = 1.0;
		{
			auto minMaxDist = [](Array1D_Vector3& pts){
				double mn = DBL_MAX, mx = -DBL_MAX;
				for (auto& p : pts) { double d = p.norm(); mn = std::min(mn, d); mx = std::max(mx, d); }
				return std::make_pair(mn, mx);
			};

			auto bounds_orig = minMaxDist(original);
			auto bounds_curr = minMaxDist(current);

			if (bounds_curr.first > bounds_orig.first){
				double ratio = std::min(1.0, bounds_curr.first / bounds_orig.second);
				connection_weight = 1.0 - ratio;
			}
		}

		/// Geometric: scale by geometric similarity
		double geom_weight = 1.0;
		{
			auto fixedNodeChangedType = [&](Structure::Node * n){
				bool isFixedAlready = searchNode->mapping.contains(n->id);
				if (!isFixedAlready) return false;
				if (n->type() == targetShape->getNode(searchNode->mapping[n->id])->type()) return false;
				return true;
			};
			
			auto ratio = [&](Structure::Node * n){
				if (fixedNodeChangedType(n)){
					double rs = 1.0 / shape->relationOf(n->id).parts.size();
					double rt = 1.0 / targetShape->relationOf(searchNode->mapping[n->id]).parts.size();
					return rs == rt ? 0 : std::min(rs, rt);
				}
				return 0.0;
			};

			double r = std::max( ratio(l->n1), ratio(l->n2) );
			geom_weight = 1.0 - r;
		}

		/// Structure: favor diversity in assignments
		double diversity_weight = 1.0;
		{
			if (searchNode->mapping.contains(l->n1->id) && searchNode->mapping.contains(l->n2->id))
				if (searchNode->mapping[l->n1->id] == searchNode->mapping[l->n2->id])
					diversity_weight = 0.5;
		}

		for (size_t i = 0; i < original.size(); i++)
		{
			double v = original[i].normalized().dot(current[i].normalized());

			// Bad values e.g. nan are broken edge values
			if (!std::isfinite(v)) v = 0;

			// Bound: beyond 90 degrees the original edge is broken 
			if (v < 0) v = 0;

			// Scaling
			v *= connection_weight * geom_weight * diversity_weight;

			feature_vector << v;
		}
	}

	Eigen::VectorXd original_v = Eigen::VectorXd::Ones(feature_vector.size());

	Eigen::VectorXd v(feature_vector.size());
	for (size_t d = 0; d < feature_vector.size(); d++) v[d] = feature_vector[d];

	double v_norm = v.norm(), original_norm = original_v.norm();

	// Value of 1 means perfect match, anything lower is bad
	double similarity = v_norm / original_norm;
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

double EvaluateCorrespondence::RMSD(Structure::ShapeGraph * shapeA, Structure::ShapeGraph * shapeB)
{
	// Collect point sets X, Y
	auto sampleShape = [&](Structure::ShapeGraph * shape){
		std::vector<Vector3> samples;
		for (auto n : shape->nodes){
			auto coords = n->property["samples_coords"].value<Array2D_Vector4d>();
			if (coords.empty()) coords = EvaluateCorrespondence::sampleNode(n, 0);
			for (auto row : coords) for (auto c : row) samples.push_back( n->position(c) );
		}
		return samples;
	};

	auto X = sampleShape(shapeA);
	auto Y = sampleShape(shapeB);

	// Point query acceleration
	auto buildKdTree = []( std::vector<Vector3> samples ){
		auto t = new NanoKdTree;
		for (auto p : samples) t->addPoint(p);
		t->build();
		return QSharedPointer<NanoKdTree>(t);
	};
	QSharedPointer<NanoKdTree> X_tree = buildKdTree(X);
	QSharedPointer<NanoKdTree> Y_tree = buildKdTree(Y);

	auto dist_x_Y = [](QSharedPointer<NanoKdTree> X_tree, QSharedPointer<NanoKdTree> Y_tree){
		auto & X = X_tree->cloud.pts, Y = Y_tree->cloud.pts;
		int n = X.size();
		double sum = 0;
		for (int i = 0; i < n; i++)
		{
			auto & p = X[i];
			KDResults matches;
			Y_tree->k_closest(p, 1, matches); // min_j
			sum += matches.front().second; // (Euclidean distance)^2
		}
		return sum;
	};

	int n = (int) X.size(); // or Y.size()??

	double a = (dist_x_Y(X_tree, Y_tree) + dist_x_Y(Y_tree, X_tree));
	double b = 2.0 * n;

	return sqrt( a / b );
}
