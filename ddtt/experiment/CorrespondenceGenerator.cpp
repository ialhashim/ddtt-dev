#include "CorrespondenceGenerator.h"

typedef QVector< QVector<size_t> > Assignments;

CorrespondenceGenerator::CorrespondenceGenerator(Structure::ShapeGraph *shapeA, Structure::ShapeGraph *shapeB) : a(shapeA), b(shapeB)
{

}

Paths CorrespondenceGenerator::generate()
{
	Paths result;

	auto boxA = a->bbox(), boxB = b->bbox();

	QStringList nodesA, nodesB;
	for (auto n : a->nodes) nodesA << n->id;
	for (auto n : b->nodes) nodesB << n->id;
	for (auto g : a->groups) for (auto nid : g) nodesA.removeAll(nid);
	for (auto g : b->groups) for (auto nid : g) nodesB.removeAll(nid);
	Structure::NodeGroups tmpA, tmpB;
	for (auto nid : nodesA) tmpA.push_back( QVector<QString>() << nid );
	for (auto nid : nodesB) tmpB.push_back( QVector<QString>() << nid );
	for (auto g : a->groups) tmpA.push_back(g);
	for (auto g : b->groups) tmpB.push_back(g);
	a->groups = tmpA;
	b->groups = tmpB;

	/// Build similarity matrix of groups:
	Eigen::MatrixXd similiarity = Eigen::MatrixXd::Ones(a->groups.size(), b->groups.size() + 1);

	for (size_t i = 0; i < similiarity.rows(); i++){
		double minVal = DBL_MAX;
		for (size_t j = 0; j < similiarity.cols(); j++){
			double val = 0;

			if (j < b->groups.size())
			{
				Eigen::AlignedBox3d groupBoxA, groupBoxB;
				for (auto sid : a->groups[i]) groupBoxA.extend(a->getNode(sid)->bbox(1.01));
				for (auto sid : b->groups[j]) groupBoxB.extend(b->getNode(sid)->bbox(1.01));

				Vector3 relativeCenterA = (groupBoxA.center() - boxA.min()).array() / boxA.sizes().array();
				Vector3 relativeCenterB = (groupBoxB.center() - boxB.min()).array() / boxB.sizes().array();

				val = (relativeCenterA - relativeCenterB).norm();

				minVal = std::min(minVal, val);
			}

			similiarity(i, j) = val;
		}
		double maxVal = similiarity.row(i).maxCoeff();
		if (maxVal == 0) maxVal = 1.0;

		for (size_t j = 0; j < b->groups.size(); j++)
			similiarity(i, j) = (similiarity(i, j) - minVal) / maxVal;
	}

	std::vector< std::vector<float> > data;
	for (int i = 0; i < similiarity.rows(); i++){
		std::vector<float> dataRow;
		for (int j = 0; j < similiarity.cols(); j++)
			dataRow.push_back(similiarity(i, j));
		data.push_back(dataRow);
	}
	//if (false) showTableColorMap(data, true); // DEBUG

	// Parameters:
	double similarity_threshold = 0.4;

	if (similiarity.rows() < 4) similarity_threshold = 1.0;

	// Collect good candidates
	Assignments candidates;
	for (size_t i = 0; i < similiarity.rows(); i++){
		QVector<size_t> candidate;
		for (size_t j = 0; j < similiarity.cols(); j++){
			if (similiarity(i, j) < similarity_threshold)
				candidate << j;
		}
		candidates << candidate;
	}

	int count_threshold = 1;

	Assignments assignments;
	if (candidates.size() > 7)
	{
		debugBox(QString("%1 candidates. Too complex (lower similairty?)").arg(candidates.size()));
		cart_product(assignments, candidates, 10000);
		count_threshold = 12;
	}
	else
	{
		cart_product(assignments, candidates);
	}

	// Filter assignments
	Assignments filtered;
	for (auto & a : assignments)
	{
		QMap<size_t, size_t> counts;
		bool accept = true;
		auto NOTHING_SEGMENT = similiarity.cols() - 1;

		for (auto i : a) counts[i]++;
		for (auto k : counts.keys())
		{
			if (k != NOTHING_SEGMENT && counts[k] > count_threshold){
				accept = false;
				break;
			}
		}

		if (accept) filtered << a;
	}

	for (auto & assignment : filtered)
	{
		QVector < QStringList > landmarksA, landmarksB;

		for (size_t i = 0; i < assignment.size(); i++)
		{
			if (assignment[i] == b->groups.size()) continue;

			auto grpA = a->groups[i];
			auto grpB = b->groups[assignment[i]];

			// Resolve many-to-many
			if (grpA.size() > 1 && grpB.size() > 1)
			{
				QVector<QString> tmpA, tmpB;

				Eigen::AlignedBox3d groupBoxA, groupBoxB;
				for (auto sid : grpA) groupBoxA.extend(a->getNode(sid)->bbox(1.01));
				for (auto sid : grpB) groupBoxB.extend(b->getNode(sid)->bbox(1.01));

				for (auto nodeid_i : grpA){
					auto nodeA = a->getNode(nodeid_i);

					QMap < QString, double > dists;
					for (auto nodeid_j : grpB){
						auto nodeB = b->getNode(nodeid_j);

						Vector3 posA = (nodeA->position(Eigen::Vector4d(0.5, 0.5, 0, 0)) - boxA.min()).array() / boxA.sizes().array();
						Vector3 posB = (nodeB->position(Eigen::Vector4d(0.5, 0.5, 0, 0)) - boxB.min()).array() / boxB.sizes().array();

						dists[nodeid_j] = (posA - posB).norm();
					}
					auto nodeid_j = sortQMapByValue(dists).first().second;

					landmarksA << (QStringList() << nodeid_i);
					landmarksB << (QStringList() << nodeid_j);
				}
			}
			else
			{
				landmarksA << QStringList::fromVector(grpA);
				landmarksB << QStringList::fromVector(grpB);
			}
		}

		result.push_back( qMakePair(landmarksA, landmarksB) );
	}

    return result;
}
