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

	struct SearchNode
	{
		double cost;
		QStringList fixed, current, unassigned;
		Assignments assignments;
		QMap<QString, QString> mapping;
		Structure::ShapeGraph *shapeA, *shapeB;
		QMap<QString, QVariant> property;
		int num_children;

		SearchNode(Structure::ShapeGraph * shapeA = NULL, Structure::ShapeGraph * shapeB = NULL, 
			const QStringList & fixed = QStringList(), const Assignments & assignments = Assignments(), 
			const QStringList & unassigned = QStringList(), const QMap<QString, QString> & mapping = QMap<QString, QString>(),
			double cost = std::numeric_limits<double>::max()) : shapeA(shapeA), shapeB(shapeB), fixed(fixed), 
			assignments(assignments), unassigned(unassigned), mapping(mapping), cost(cost), num_children(0){}

		QStringList fixedOnTarget(){ QStringList result; for (auto a : assignments) for (auto p : a.second) result << p; return result; }
		QStringList unassignedList(){ 
			QStringList result;  
			for (auto n : shapeA->nodes){
				bool isUnassigned = true;
				if (fixed.contains(n->id)){ isUnassigned = false; }
				else { for (auto a : assignments) if (a.first.contains(n->id)){ isUnassigned = false; break; } }
				if(isUnassigned) result.push_back(n->id);
			}
			return result; 
		}

		bool operator<(const SearchNode & path) const { return cost < path.cost; }
	};

	typedef tree<SearchNode> SearchTree;

	struct GuidedDeformation
	{
		Structure::ShapeGraph *origShapeA, *origShapeB;
		QVector< SearchTree > searchTrees;

		static void preprocess(Structure::ShapeGraph * shapeA, Structure::ShapeGraph * shapeB);
		void searchAll(Structure::ShapeGraph * shapeA, Structure::ShapeGraph * shapeB, QVector<Energy::SearchNode> & roots);

		static void topologicalOpeartions(Structure::ShapeGraph *shapeA, Structure::ShapeGraph *shapeB, QStringList & la, QStringList & lb);
		static void applyDeformation(Structure::ShapeGraph *shapeA, Structure::ShapeGraph *shapeB, const QStringList & la, const QStringList & lb, const QStringList & fixed, bool isSaveKeyframes = false);
		static void postDeformation(Structure::ShapeGraph * shape, const QStringList & fixed);

		static void applyAssignment(Energy::SearchNode * path, bool isSaveKeyframes);
		static QVector<Energy::SearchNode> suggestChildren(const Energy::SearchNode & path);

		QVector<Energy::SearchNode*> solutions();
		QVector<Energy::SearchNode*> parents();
		QVector<Energy::SearchNode*> childrenOf(Energy::SearchNode * path);
		QVector<Energy::SearchNode*> getEntirePath(Energy::SearchNode * path);

		void applySearchPath(const QVector<Energy::SearchNode*> & path);
	};
}
