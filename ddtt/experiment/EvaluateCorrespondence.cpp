#include "EvaluateCorrespondence.h"

#include "hausdorff.h"

#include "helpers/PhysicsHelper.h"

int EvaluateCorrespondence::numSamples = 3;

Array1D_Vector3 coupledNormals(Structure::Link * l){
	Array1D_Vector3 normals;
	Eigen::Vector4d c1(0,0,0,0), c2(1,1,0,0), midc(0.5,0.5,0,0);
	
	Vector3 v1 = (l->n1->position(c2) - l->n1->position(c1)).normalized();
	Vector3 v2 = (l->n2->position(midc) - l->n1->position(c1)).normalized();

	Vector3 v3 = (l->n2->position(c2) - l->n2->position(c1)).normalized();
	Vector3 v4 = (l->n1->position(midc) - l->n2->position(c1)).normalized();

	normals.push_back(v1.cross(v2).normalized());
	normals.push_back(v3.cross(v4).normalized());

	return normals;
}

Array1D_Vector3 EvaluateCorrespondence::spokesFromLink(Structure::ShapeGraph * shape, Structure::Link * l)
{
	auto samples1 = l->n1->property["samples_coords"].value<Array2D_Vector4d>();
	auto samples2 = l->n2->property["samples_coords"].value<Array2D_Vector4d>();

	if (samples1.empty()) samples1 = EvaluateCorrespondence::sampleNode(shape, l->n1, 0);
	if (samples2.empty()) samples1 = EvaluateCorrespondence::sampleNode(shape, l->n2, 0);

	Array1D_Vector3 spokes;
	for (auto rowi : samples1) for (auto ci : rowi)
		for (auto rowj : samples2) for (auto cj : rowj)
			spokes.push_back(l->n1->position(ci) - l->n2->position(cj));

	return spokes;
}

Array2D_Vector4d EvaluateCorrespondence::sampleNode(Structure::ShapeGraph * shape, Structure::Node * n, double resolution)
{
	// Auto-set resolution when asked
	if (resolution == 0) 
		resolution = ((n->type() == Structure::SHEET) ? n->length() / 4.0 : n->length()) * (1.0 / numSamples);

	Array2D_Vector4d samples_coords = n->discretizedPoints(resolution);

	// Special case: degenerate sheet
	if (n->type() == Structure::SHEET && samples_coords.empty()){
		double step = 1.0 / numSamples;
		for (double u = 0; u <= 1.0; u += step){
			Array1D_Vector4d row;
			for (double v = 0; v <= 1.0; v += step)
				row.push_back(Eigen::Vector4d(u, v, 0, 0));
			samples_coords.push_back(row);
		}
	}

	n->property["samples_coords"].setValue(samples_coords);

	// Unary node properties
	{
		n->property["orig_diagonal"].setValue(n->diagonal());
		n->property["orig_start"].setValue(n->startPoint());
		n->property["orig_length"].setValue(n->length());

		// Volume
		double volume = 0.0;
		auto mesh = n->property["mesh"].value< QSharedPointer<SurfaceMeshModel> >();
		if (mesh->property("volume").isValid())
		{
			volume = mesh->property("volume").toDouble();
		}
		else
		{
			volume = PhysicsHelper(mesh.data()).volume();
			mesh->setProperty("volume", volume);
		}
			
		n->property["orig_volume"].setValue(volume);
	}

	return samples_coords;
}

void EvaluateCorrespondence::prepare(Structure::ShapeGraph * shape)
{
	double sum_length = 0;
	for (auto n : shape->nodes) sum_length += n->area();
	double avg_length = sum_length / shape->nodes.size();

	double resolution = avg_length / numSamples;

	shape->property["sampling_resolution"].setValue(resolution);

	for (auto n : shape->nodes)
		EvaluateCorrespondence::sampleNode(shape, n, resolution);

	// Sample curve/sheet
	for (auto l : shape->edges)
	{
		l->property["orig_spokes"].setValue(EvaluateCorrespondence::spokesFromLink(shape, l));
		l->property["orig_centroid_dir"].setValue(Vector3((l->n1->center() - l->n2->center()).normalized()));
	}
}

QMap<QString, NanoKdTree*> EvaluateCorrespondence::kdTreesNodes(Structure::ShapeGraph * shape)
{
	QMap < QString, NanoKdTree* > result;

	auto buildKdTree = [](Structure::Node * n, const Array2D_Vector4d & coords){
		auto t = new NanoKdTree;
		for (auto row : coords) for (auto coord : row) t->addPoint(n->position(coord));
		t->build();
		return t;
	};

	for (auto n : shape->nodes)
	{
		auto coords = n->property["samples_coords"].value<Array2D_Vector4d>();
		if (coords.empty()) coords = EvaluateCorrespondence::sampleNode(shape, n, 0);
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
			if (coords.empty()) coords = EvaluateCorrespondence::sampleNode(shape, n, 0);
			for (auto row : coords) for (auto c : row) samples.push_back(n->position(c));
		}
		return samples;
	};

	auto X = sampleShape(shapeA);
	auto Y = sampleShape(shapeB);

	// Point query acceleration
	auto buildKdTree = [](std::vector<Vector3> samples){
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

	int n = (int)X.size(); // or Y.size()??

	double a = (dist_x_Y(X_tree, Y_tree) + dist_x_Y(Y_tree, X_tree));
	double b = 2.0 * n;

	return sqrt(a / b);
}

double EvaluateCorrespondence::evaluate2(Energy::SearchNode * searchNode)
{
	auto shape = searchNode->shapeA.data();
	auto targetShape = searchNode->shapeB.data();

	QVector<double> feature_vector;

	feature_vector << 0;

	Eigen::Map<Eigen::VectorXd> v(&feature_vector[0], feature_vector.size());
	Eigen::VectorXd original_v = Eigen::VectorXd::Ones(v.size());
	double v_norm = v.norm(), original_norm = original_v.norm();
	double similarity = v_norm / original_norm;
	double cost = 1.0 - similarity;
	return cost;

	return 0;
}

double EvaluateCorrespondence::evaluate(Energy::SearchNode * searchNode)
{
	QVector<double> feature_vector;

	auto shape = searchNode->shapeA.data();
	auto targetShape = searchNode->shapeB.data();

	// Logging
	QVariantMap costMap;

	// Binary properties:
	for (auto l : shape->edges)
	{
		auto original_spokes = l->property["orig_spokes"].value<Array1D_Vector3>();
		auto current_spokes = EvaluateCorrespondence::spokesFromLink(shape, l);

		/// (a) Connection: scale by difference in closeness distance
		double connection_weight = 1.0;
		{
			auto minMaxDist = [](Array1D_Vector3& vectors){
				double mn = DBL_MAX, mx = -DBL_MAX;
				for (auto& p : vectors) { double d = p.norm(); mn = std::min(mn, d); mx = std::max(mx, d); }
				return std::make_pair(mn, mx);
			};

			auto bounds_orig = minMaxDist(original_spokes);
			auto bounds_curr = minMaxDist(current_spokes);

			double connection_weight_near = 1.0, connection_weight_far = 1.0;

			// Near point got further
			if (bounds_curr.first > bounds_orig.first){
				double ratio = std::min(1.0, bounds_curr.first / bounds_orig.second);
				connection_weight_near = 1.0 - ratio;
			}

			// Far point got closer
			//if (bounds_curr.second < bounds_orig.second * 1.25){
			//	connection_weight_far = bounds_curr.second / bounds_orig.second;
			//}

			connection_weight = std::min(connection_weight_near, connection_weight_far);
		}

		/// (b) Geometric: scale by geometric (topological?) similarity
		double split_weight = 1.0;
		{
			auto volumeRatio = [&](Structure::Node * n){
				double ratio = 1.0;

				if (!n->property["isManyMany"].toBool())
				{
					if (n->property["isSplit"].toBool())
					{
						double partVolume = n->length();

						QStringList targetIDs = targetShape->relationOf(searchNode->mapping[n->id]).parts;
						double targetVolume = 0.0;
						for (auto tpartID : targetIDs) targetVolume += targetShape->getNode(tpartID)->length();

						ratio = targetVolume / std::max(partVolume, targetVolume);
					}
					if (n->property["isMerged"].toBool())
					{
						auto targetNode = targetShape->getNode(searchNode->mapping[n->id]);
						double targetVolume = targetNode->length();

						QStringList sourceIDs = n->property["groupParts"].toStringList();
						double sourceVolume = 0.0;
						for (auto partID : sourceIDs)
						{
							auto nj = shape->getNode(partID);
							sourceVolume += nj->length();
						}

						ratio = sourceVolume / std::max(targetVolume, sourceVolume);
					}
				}

				return ratio;
			};

			split_weight = std::min(volumeRatio(l->n1), volumeRatio(l->n2));
		}

		/// (b) Geometric: scale by geometric (topological?) similarity
		/*double split_weight = 1.0;
		{
			auto ratio = [&](Structure::Node * n){
				if (n->property["isSplit"].toBool() || n->property["isMerged"].toBool()){
					double rs = 1.0 / shape->relationOf(n->id).parts.size();
					double rt = 1.0 / targetShape->relationOf(searchNode->mapping[n->id]).parts.size();
					return rs == rt ? 0 : std::min(rs, rt);
				}
				return 0.0;
			};

			double r = std::max(ratio(l->n1), ratio(l->n2));
			split_weight = 1.0 - r;
		}*/

		QVector < double > link_vector;

		for (size_t i = 0; i < original_spokes.size(); i++)
		{
			double v = original_spokes[i].normalized().dot(current_spokes[i].normalized());

			// Bad values e.g. nan are broken feature values
			if (!std::isfinite(v)) v = 0;

			// Bound: the original edge is broken beyond a certain point
			v = std::max(0.0, v);

			// Scaling
			v *= connection_weight * split_weight;

			link_vector << v;
		}

		for (auto v : link_vector) feature_vector << v;

		double sum_link_vector = 0; for (auto v : link_vector) sum_link_vector += v;
		double avg_link_vector = sum_link_vector / link_vector.size();
		//feature_vector << avg_link_vector;

		// Logging:
		qSort(link_vector);
		costMap[l->id].setValue(QString("[%1, %2, avg = %6] w_conn = %3 / w_geom = %4 / w_diverse = %5 / w_context = %7")
			.arg(link_vector.front()).arg(link_vector.back()).arg(connection_weight).arg(split_weight).arg(1).arg(avg_link_vector).arg(1));
	}

	// Logging:
	shape->property["costs"].setValue(costMap);

	Eigen::Map<Eigen::VectorXd> v(&feature_vector[0], feature_vector.size());
	Eigen::VectorXd original_v = Eigen::VectorXd::Ones(v.size());

	double v_norm = v.norm(), original_norm = original_v.norm();

	// Value of 1 means perfect matching to original shape, anything lower is bad
	double similarity = v_norm / original_norm;
	double cost = 1.0 - similarity;

	return cost;
}
