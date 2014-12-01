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
		static QVector <Energy::SearchPath*> getEntirePath(SearchPath * selected_path, QVector<SearchPath> & paths){
			auto t = SearchPath::exploreAsTree(paths);

			// Search for selected solution:
			auto itr = t.begin();
			for (; itr != t.end(); itr++) if (*itr == selected_path) break;

			// Find ancestors
			tree < Energy::SearchPath::SearchNode >::iterator current = itr;
			QVector <Energy::SearchPath*> entire_path;
			while (current != 0){
				entire_path.push_front(*current);
				current = t.parent(current);
			}

			return entire_path;
		}

		// Smaller memory footprint
		typedef byte PartIndex;
		QVector< QPair<QVector<PartIndex>, QVector<PartIndex> > > _assignments;
		QVector< PartIndex > _fixed, _current, _unassigned;
		QMap < PartIndex, PartIndex > _mapping;
		void compress(const QMap<QString, PartIndex> & mapA, const QMap<QString, PartIndex> & mapB);
		void decompress(const QMap<PartIndex, QString> & mapA, const QMap<PartIndex, QString> & mapB);
	};

	struct GuidedDeformation{
		QVector<SearchPath> search_paths;
		Structure::ShapeGraph *origShapeA, *origShapeB;

		void searchAll();		
		void explore(SearchPath & path);

		QVector<Energy::SearchPath*> solutions();
		QVector<Energy::SearchPath*> parents();

		void applySearchPath(const QVector<Energy::SearchPath*> & path);

		// Utility:
		void topologicalOpeartions(Structure::ShapeGraph *shapeA, Structure::ShapeGraph *shapeB,
			QStringList & la, QStringList & lb);
		void applyDeformation(Structure::ShapeGraph *shapeA, Structure::ShapeGraph *shapeB, 
			const QStringList & la, const QStringList & lb, const QStringList & fixed, bool isSaveKeyframes = false);

		// Compression:
		QMap<SearchPath::PartIndex, QString> idxPartMapA, idxPartMapB;
		QMap<QString, SearchPath::PartIndex> partIdxMapA, partIdxMapB;
	};
}
