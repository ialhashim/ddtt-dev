#pragma once
#include <QStack>

#include "ShapeGraph.h"
#include "RenderObjectExt.h"

namespace Energy
{
	static QString null_part = "NULL_PART";
	typedef QStack< QPair<QStringList, QStringList> > AssignmentsStack;

	struct SearchPath{
		Structure::ShapeGraph *shapeA, *shapeB;
		QStringList fixed;
		AssignmentsStack assignments;
		QStringList unassigned;
		QVector<SearchPath> children;
		double cost;

		SearchPath(Structure::ShapeGraph * shapeA = NULL, Structure::ShapeGraph * shapeB = NULL, const QStringList & fixed = QStringList(),
			const AssignmentsStack & assignments = AssignmentsStack(), const QStringList & unassigned = QStringList(),
			double cost = 0) : shapeA(shapeA), shapeB(shapeB), fixed(fixed), assignments(assignments), unassigned(unassigned), cost(cost){}

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

		bool operator<(const SearchPath & path){ return cost < path.cost; }
	};

	struct GuidedDeformation{
		QVector<SearchPath> search_paths;

		void searchAll();		
		void explore(SearchPath & path);

		// Utility:
		static void topologicalOpeartions(Structure::ShapeGraph *shapeA, Structure::ShapeGraph *shapeB,
			QStringList & la, QStringList & lb);
		static void applyDeformation(Structure::ShapeGraph *shapeA, Structure::ShapeGraph *shapeB, 
			const QStringList & la, const QStringList & lb, const QStringList & fixed);
	};
}
