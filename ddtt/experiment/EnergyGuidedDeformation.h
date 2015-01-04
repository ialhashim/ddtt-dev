#pragma once
#undef SearchPath

#include "ShapeGraph.h"
#include "RenderObjectExt.h"

#include <QStack>
#include "tree.hh"

namespace Energy
{
	typedef QVector< QPair<QStringList, QStringList> > Assignments;
	static QString null_part = "NULL_PART";
	typedef QSet<QString> QSetString;

	struct SearchNode
	{
		double cost, energy;
		QSetString fixed, current, unassigned;
		Assignments assignments;
		QMap<QString, QString> mapping;
		QSharedPointer<Structure::ShapeGraph> shapeA, shapeB;
		QMap<QString, QVariant> property;
		int num_children;

		SearchNode(QSharedPointer<Structure::ShapeGraph> shapeA = QSharedPointer<Structure::ShapeGraph>(), 
			QSharedPointer<Structure::ShapeGraph> shapeB = QSharedPointer<Structure::ShapeGraph>(),
			const QSetString & fixed = QSetString(), const Assignments & assignments = Assignments(),
			const QSetString & unassigned = QSetString(), const QMap<QString, QString> & mapping = QMap<QString, QString>(),
			double cost = std::numeric_limits<double>::max(), double energy = 0.0)
			: shapeA(shapeA), shapeB(shapeB), fixed(fixed), assignments(assignments), unassigned(unassigned), mapping(mapping), 
			cost(cost), energy(energy), num_children(0){}

		QSetString fixedOnTarget(){ QSetString result; for (auto a : assignments) for (auto p : a.second) result << p; return result; }
		QSetString unassignedList(){
			QSetString result;
			for (auto n : shapeA->nodes){
				bool isUnassigned = true;
				if (fixed.contains(n->id)){ isUnassigned = false; }
				else { for (auto a : assignments) if (a.first.contains(n->id)){ isUnassigned = false; break; } }
				if(isUnassigned) result << n->id;
			}
			return result; 
		}

		bool operator<(const SearchNode & path) const { return energy < path.energy; }
	};

	typedef tree<SearchNode> SearchTree;

	struct GuidedDeformation
	{
		QSharedPointer<Structure::ShapeGraph> origShapeA, origShapeB;
		QVector< SearchTree > searchTrees;

		static void preprocess(Structure::ShapeGraph * shapeA, Structure::ShapeGraph * shapeB);
		void searchAll(Structure::ShapeGraph * shapeA, Structure::ShapeGraph * shapeB, QVector<Energy::SearchNode> & roots);

		static void topologicalOpeartions(Structure::ShapeGraph *shapeA, Structure::ShapeGraph *shapeB, QStringList & la, QStringList & lb);
		static void applyDeformation(Structure::ShapeGraph *shapeA, Structure::ShapeGraph *shapeB, const QStringList & la, const QStringList & lb, const QSetString & fixed, bool isSaveKeyframes = false);
		static void postDeformation(Structure::ShapeGraph * shape, const QSet<QString> & fixed);

		static void applyAssignment(Energy::SearchNode * path, bool isSaveKeyframes);
		static QVector<Energy::SearchNode> suggestChildren(Energy::SearchNode & path);

		QVector<Energy::SearchNode*> solutions();
		QVector<Energy::SearchNode*> parents();
		QVector<Energy::SearchNode*> childrenOf(Energy::SearchNode * path);
		QVector<Energy::SearchNode*> getEntirePath(Energy::SearchNode * path);

		void applySearchPath(QVector<Energy::SearchNode*> path);
	};
}
