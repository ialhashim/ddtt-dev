#pragma once
#include <QStack>

#include "ShapeGraph.h"
#include "RenderObjectExt.h"

#undef SearchPath

#include "tree.hh"

namespace Energy
{
	static QString null_part = "NULL_PART";
	typedef QVector< QPair<QStringList, QStringList> > Assignments;

	struct SearchPath{
		QVector<SearchPath> children;
		double cost;
		QStringList fixed, current, unassigned;
		Assignments assignments;
		QMap<QString, QString> mapping;
		Structure::ShapeGraph *shapeA, *shapeB;
		QMap<QString, QVariant> property;

		SearchPath(Structure::ShapeGraph * shapeA = NULL, Structure::ShapeGraph * shapeB = NULL, const QStringList & fixed = QStringList(),
			const Assignments & assignments = Assignments(), const QStringList & unassigned = QStringList(), const QMap<QString, QString> & mapping = QMap<QString, QString>(),
			double cost = 0) : shapeA(shapeA), shapeB(shapeB), fixed(fixed), assignments(assignments), unassigned(unassigned), mapping(mapping), cost(cost){}

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

		// Extract as a tree structure
		typedef SearchPath* SearchNode;
		static void exploreAsTree(tree<SearchNode> & t, tree<SearchNode>::iterator & parent, Energy::SearchPath * current){
			auto cur_itr = t.append_child(parent, current);
			for (auto & p : current->children)
				SearchPath::exploreAsTree(t, cur_itr, &p);
		}
		static tree<Energy::SearchPath*> exploreAsTree(QVector<Energy::SearchPath> & paths){
			tree<Energy::SearchPath*> t;
			for (auto & current : paths){
				auto root = t.insert(t.begin(), &current);
				for (auto & child : current.children)
					SearchPath::exploreAsTree(t, root, &child);
			}
			return t;
		}
	};

	struct GuidedDeformation{
		QVector<SearchPath> search_paths;

		void searchAll();		
		void explore(SearchPath & path);

		QVector<Energy::SearchPath*> solutions();
		QVector<Energy::SearchPath*> parents();

		// Utility:
		static void topologicalOpeartions(Structure::ShapeGraph *shapeA, Structure::ShapeGraph *shapeB,
			QStringList & la, QStringList & lb);
		static void applyDeformation(Structure::ShapeGraph *shapeA, Structure::ShapeGraph *shapeB, 
			const QStringList & la, const QStringList & lb, const QStringList & fixed);
	};
}
