#include "DeformEnergy2.h"

#define _USE_MATH_DEFINES 
#include <math.h>

Array2D_Vector4d DeformEnergy2::sideCoordinates = Structure::ShapeGraph::computeSideCoordinates();

DeformEnergy2::DeformEnergy2(Structure::ShapeGraph *shapeA, Structure::ShapeGraph *shapeB,
                             const QVector<QStringList> &a_landmarks, const QVector<QStringList> &b_landmarks,
							 bool debugging) : a(shapeA), b(shapeB)
{
	// Perform topological operations needed for sheet-curves correspondence case:
	mappedPartsA = a_landmarks; mappedPartsB = b_landmarks;
	for (size_t i = 0; i < mappedPartsA.size(); i++)
	{
		auto & la = mappedPartsA[i];
		auto & lb = mappedPartsB[i];

		bool isOneA = la.size() == 1, isOneB = lb.size() == 1;
		bool isSheetA = a->getNode(la.front())->type() == Structure::SHEET;
		bool isSheetB = b->getNode(lb.front())->type() == Structure::SHEET;
		bool isOneToMany = (isOneA != isOneB), isOneIsSheet = (isSheetA || isSheetB), isOtherNotSheet = !isSheetA || !isSheetB;

		// Map between parts from shape A to shape B
		for (size_t j = 0; j < la.size(); j++)
		{
			mappingAB[la[j]] = lb;

			// Align parts in their parameter space
			for (auto tid : lb) Structure::ShapeGraph::correspondTwoNodes(la[j], a, tid, b);
		}

		if (isOneToMany && isOneIsSheet && isOtherNotSheet)
		{
			Structure::ShapeGraph * graph = (isOneA) ? shapeB : shapeA;
			QStringList * landmarks = (isOneA) ? &lb : &la;

			QString newnode = Structure::ShapeGraph::convertCurvesToSheet(graph, *landmarks, sideCoordinates);
			Structure::Node * newNode = graph->getNode(newnode);

			landmarks->clear();
			landmarks->push_back(newnode);
		}
	}

	// Collect set A of all parts corresponded:
	QStringList A;
	for (auto l : mappedPartsA) for (auto nid : l) A << nid;

	// Collect set B of parts sharing edges in both source and target:
	ShapeEdges B, B_missed;
	for (size_t i = 0; i < mappedPartsA.size(); i++){
		for (size_t j = i+1; j < mappedPartsA.size(); j++){
			auto partsShareEdge = [&](Structure::ShapeGraph * graph, const QStringList& partsA, const QStringList& partsB)
			{
				QStringList allParts;
				for (auto partID : partsA) for (auto nid : partID.split(",")) allParts << nid;
				for (auto partID : partsB) for (auto nid : partID.split(",")) allParts << nid;
				for (size_t u = 0; u < allParts.size(); u++)
					for (size_t v = u + 1; v < allParts.size(); v++)
						if (graph->shareEdge(allParts[u], allParts[v])) return true;
				return false;
			};

			if ( partsShareEdge(a, mappedPartsA[i], mappedPartsA[j]) )
			{
				bool isExistOnBoth = partsShareEdge(b, mappedPartsB[i], mappedPartsB[j]);

				for (size_t n = 0; n < mappedPartsA[i].size(); n++)
				{
					for (size_t m = 0; m < mappedPartsA[j].size(); m++)
					{
						if (!a->shareEdge(mappedPartsA[i][n], mappedPartsA[j][m])) continue;

						if (isExistOnBoth)
							B << qMakePair(mappedPartsA[i][n], mappedPartsA[j][m]);
						else
							B_missed << qMakePair(mappedPartsA[i][n], mappedPartsA[j][m]);
					}
				}
			}
		}
	}

	// Evaluate energy terms
	double theta = 0, omega = 0, kappa = 0, rho = 0;
	{
		theta = computeAngles(B, B_missed);
		omega = computeEdges(B, B_missed);
		kappa = computeContext(A);
		rho = computeSymmetry(A);
	}

	// Structural energy:
	double w_s = 1.0;
	double E_s = theta + omega + kappa + rho;
	this->energyTerms["w_s"].setValue(w_s);
	this->energyTerms["E_s"].setValue(E_s);

	this->energyTerms["_t"].setValue(theta);
	this->energyTerms["_o"].setValue(omega);
	this->energyTerms["_k"].setValue(kappa);
	this->energyTerms["_r"].setValue(rho);

	// Geometric energy:
	double w_g = 1.0;
	double E_g = 1.0;
	this->energyTerms["w_g"].setValue(w_g);
	this->energyTerms["E_g"].setValue(E_g);

	// Total Energy:
	double E_total = (w_s * E_s) + (w_g * E_g);

	// Regularization energy:
	double w_reg = 1.0;
	double E_reg = E_regularizer(a_landmarks);
	this->energyTerms["w_reg"].setValue( w_reg );
	this->energyTerms["E_reg"].setValue( E_reg );

	this->total_energy = E_total + (w_reg * E_reg);

	// Clean up after any previously applied topological operations:
	{
		auto cleanUpMerges = [&](Structure::ShapeGraph * graph, QVector<QStringList> mappedParts){
			for (auto landmark : mappedParts){
				bool isMerged = false;
				for (auto & nid : landmark) if (nid.contains(",")){ isMerged = true; break; }
				if (!isMerged) continue;
				graph->removeNode(landmark.front());
			}
		};

		cleanUpMerges(a, mappedPartsA);
		cleanUpMerges(b, mappedPartsB);
	}
}

/* ===================	*/
/*  Energy Terms:		*/
/* ===================	*/
double DeformEnergy2::computeAngles(const ShapeEdges & B, const ShapeEdges & B_missed)
{
	if (B.size() == 0) return 1.0;

	QMap<QString,double> angleDissimilarity;

	// Compare angle using normals at edges of parametric parts:
	auto computeAngle = [&](Structure::ShapeGraph * shape, const PairParts & parts){
		auto edge = shape->getEdge(parts.first, parts.second);

		Vector3 pos;
		std::vector<Eigen::Vector3d> frame1, frame2;
		edge->n1->get(edge->getCoord(edge->n1->id).front(), pos, frame1);
		edge->n2->get(edge->getCoord(edge->n2->id).front(), pos, frame2);

		Vector3 v1 = edge->n1->type() == Structure::SHEET ? frame1.back() : frame1.front();
		Vector3 v2 = edge->n2->type() == Structure::SHEET ? frame2.back() : frame2.front();

		double angle = acos(abs(v1.dot(v2)));

		return angle;
	};

	// Go over all edge relations:
	for (size_t i = 0; i < B.size(); i++)
	{
		auto corrEdges = correspondedEdges(B[i]);

		QVector<double> anglesA, anglesB;
		for (auto edge : corrEdges.first) anglesA << computeAngle(a, edge);
		for (auto edge : corrEdges.second) anglesB << computeAngle(b, edge);

		// Compare the angles, and use the minimum difference
		double minDiff = DBL_MAX;
		for (size_t u = 0; u < anglesA.size(); u++){
			for (size_t v = 0; v < anglesB.size(); v++){
				minDiff = std::min(minDiff, abs(anglesA[u] - anglesB[u]));
			}
		}

		double angle_difference = minDiff / (M_PI * 0.5);

		angleDissimilarity[B[i].first + "-" + B[i].second] = angle_difference;
	}

	// Penalties of missed relations:
	for (size_t i = 0; i < B_missed.size(); i++)
	{
		angleDissimilarity[B_missed[i].first + "-" + B_missed[i].second] = 1.0;
	}

	// Sum up angle dissimilarities
	double sum = 0;
	for (auto dissimilarity : angleDissimilarity) sum += dissimilarity;

	double theta = sum / (B.size() + B_missed.size());
	return theta;
}

double DeformEnergy2::computeEdges(const ShapeEdges & B, const ShapeEdges & B_missed)
{
	if (B.size() == 0) return 1.0;

	QMap<QString, double> edgeDissimilarity;

	auto getCoordinates = [&](Structure::ShapeGraph * shape, const PairParts & parts, bool isDropV1, bool isDropV2){
		auto edge = shape->getEdge(parts.first, parts.second);
		auto c1 = edge->getCoord(parts.first).front();
		auto c2 = edge->getCoord(parts.second).front();

		if (isDropV1) c1[1] = 0;
		if (isDropV2) c2[1] = 0;

		return qMakePair(c1, c2);
	};

	// Go over all edge relations:
	for (size_t i = 0; i < B.size(); i++)
	{
		auto corrEdges = correspondedEdges(B[i]);

		// We need to resolve this for sheet parts (default: disable 'v' coordinate)
		bool isDropV1 = true, isDropV2 = true;

		QVector< QPair<Eigen::Vector4d, Eigen::Vector4d> > coordsA, coordsB;
		for (auto edge : corrEdges.first) coordsA << getCoordinates(a, edge, isDropV1, isDropV2);
		for (auto edge : corrEdges.second) coordsB << getCoordinates(b, edge, isDropV1, isDropV2);

		QPair<double, double> minDiff(DBL_MAX, DBL_MAX);
		for (size_t u = 0; u < coordsA.size(); u++){
			for (size_t v = 0; v < coordsB.size(); v++){
				minDiff.first = std::min(minDiff.first, (coordsA[u].first - coordsB[u].first).norm());
				minDiff.second = std::min(minDiff.second, (coordsA[u].second - coordsB[u].second).norm());
			}
		}

		QPair<double, double> edge_difference (minDiff.first, minDiff.second);

		// Average of edge distance
		edgeDissimilarity[B[i].first + "-" + B[i].second] = (edge_difference.first + edge_difference.second) / 2.0;
	}

	// Penalties of missed relations:
	for (size_t i = 0; i < B_missed.size(); i++)
	{
		edgeDissimilarity[B_missed[i].first + "-" + B_missed[i].second] = 1.0;
	}

	// Sum up edge dissimilarities
	double sum = 0;
	for (auto dissimilarity : edgeDissimilarity) sum += dissimilarity;

	double omega = sum / (B.size() + B_missed.size());
	return omega;
}

double DeformEnergy2::computeContext(const QStringList & A)
{
	if (A.size() == 0) return 1.0;

	// Compute spatial range of a part along the upright axis (conventionally z-axis):
	auto partInterval = [&](Structure::ShapeGraph * shape, const QString & part){
		std::vector<double> z_values;
		auto node = shape->getNode(part);
		z_values.push_back(node->position(Eigen::Vector4d(0, 0, 0, 0)).z());
		z_values.push_back(node->position(Eigen::Vector4d(0.5, 0.5, 0, 0)).z());
		z_values.push_back(node->position(Eigen::Vector4d(1, 1, 0, 0)).z());
		auto minmax = std::minmax_element(z_values.begin(), z_values.end());
		return qMakePair(*minmax.first, *minmax.second);
	};

	// Context computation:
	typedef QVector< QPair<double, double> > VecPairDoubles;
	auto context = [&](const VecPairDoubles & intervals, size_t i, size_t j){
		bool isBelow = intervals[i].second < intervals[j].first;
		bool isAbove = intervals[i].first > intervals[j].second;
		if (isBelow) return -1;
		else if (isAbove) return 1;
		else return 0;
	};

	// Pre-compute parts intervals on source and target shapes:
	VecPairDoubles partsIntervalA, partsIntervalB;
	for (size_t i = 0; i < A.size(); i++) partsIntervalA << partInterval(a, A[i]);
	for (size_t i = 0; i < A.size(); i++) partsIntervalB << partInterval(b, mappingAB[A[i].split(",").front()].front());

	// Compute part context dissimilarity
	QMap<QString, double> contextDissimilarity;

	for (size_t i = 0; i < A.size(); i++)
	{
		QVector<double> partContextA, partContextB;
		for (size_t j = 0; j < A.size(); j++)
		{
			partContextA << context(partsIntervalA, i, j);
			partContextB << context(partsIntervalB, i, j);
		}

		// Compute dissimilarity ratio
		int numDissimilar = 0;
		for (size_t j = 0; j < A.size(); j++)
			if (partContextA[j] != partContextB[j])
				numDissimilar++;

		contextDissimilarity[A[i]] = double(numDissimilar) / A.size();
	}

	// Sum up context dissimilarities
	double sum = 0;
	for (auto dissimilarity : contextDissimilarity) sum += dissimilarity;

	double kappa = sum / (A.size());
	return kappa;
}

double DeformEnergy2::computeSymmetry(const QStringList & A)
{
	if (A.size() == 0) return 1.0;

	// Groups from source with members in A 
	QMap < QString, QVector<QString> > foundGroups;
	for (auto group : a->groups){
		for (auto nid : group){
			if (A.contains(nid)){
				QStringList key;
				for (auto nid : group) key << nid;
				foundGroups[key.join("-")] = group;
				break;
			}
		}
	}

	// Compute group dissimilarity:
	QMap<QString, double> contextDissimilarity;
	for (auto groupKey : foundGroups.keys())
	{
		auto & group = foundGroups[groupKey];

		int beforeCount = group.size();

		// Find its equivalent
		QString targetPart = mappingAB[group.front()].front();
		int afterCount = b->groupsOf(targetPart).front().size();

		double ratio = double(std::min(beforeCount, afterCount)) / std::max(beforeCount, afterCount);

		contextDissimilarity[groupKey] = 1.0 - ratio;
	}

	// Sum up group dissimilarities
	double sum = 0;
	for (auto dissimilarity : contextDissimilarity) sum += dissimilarity;

	double rho = sum / (A.size());
	return rho;
}

double DeformEnergy2::E_regularizer(const QVector<QStringList> & correspondedSet)
{
	QStringList sourceNodes;
	for (auto l : correspondedSet) for (auto nid : l) sourceNodes << nid;

	// Get shape A nodes as groups
	auto grps = a->groups;
	for (auto n : a->nodes){
		bool isFound = false;
		for (auto g : grps) if (g.contains(n->id)) isFound = true;
		if (!isFound && !n->id.contains(",")) grps.push_back(QVector<QString>() << n->id);
	}

	// Count corresponded groups
	int correspondedGroups = 0;
	for (auto g : grps){
		bool isFound = false;
		for (auto nid : sourceNodes){
			if (g.contains(nid)){
				isFound = true;
				break;
			}
		}
		if (isFound) correspondedGroups++;
	}

	double ratio = double(correspondedGroups) / grps.size();
	if (ratio == 0) ratio = 0.001; // nothing was corresponded

	return 1.0 / ratio;
}

// Utility:
QPair< QVector<PairParts>, QVector<PairParts> > DeformEnergy2::correspondedEdges(const PairParts & partsA)
{
	QVector < PairParts > pairEdgesA, pairEdgesB;

	// From source:
	for (auto ni : partsA.first.split(","))
		for (auto nj : partsA.second.split(","))
			if (a->shareEdge(ni, nj)) pairEdgesA << PairParts(ni, nj);

	// From target:
	for (auto edge : pairEdgesA){
		for (auto ti : mappingAB[edge.first])
			for (auto tj : mappingAB[edge.second])
				if (b->shareEdge(ti, tj)) pairEdgesB << PairParts(ti, tj);
	}

	assert(!pairEdgesA.isEmpty() && !pairEdgesB.isEmpty());

	// Contains duplicates in one case
	return qMakePair(pairEdgesA, pairEdgesB);
}
