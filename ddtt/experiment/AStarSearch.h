#pragma once

#include "EnergyGuidedDeformation.h"
#include "EvaluateCorrespondence.h"

#include "stlastar.h"

namespace AStar
{
	class PathSearchNode : public Energy::SearchNode{
	public:
		PathSearchNode(const Energy::SearchNode & n = Energy::SearchNode()) : Energy::SearchNode(n){}

		/* Heuristic */
		double GoalDistanceEstimate( PathSearchNode & )
		{
			auto allDists = EvaluateCorrespondence::hausdroffDistance(shapeA.data(), shapeB.data());

			double sum = 0;

			for (auto partID : allDists.keys())
			{
				auto dists = allDists[partID];

				double curSum = 0;
				for (auto d : dists) curSum += d;
				sum += curSum / dists.size();
			}
			
			return sum / allDists.size();
		}

		/* Actual cost (arc cost) */
		double GetCost( PathSearchNode &successor )
		{
			Energy::GuidedDeformation::applyAssignment(&successor, false);
			return successor.cost;
		}

		bool GetSuccessors( AStarSearch<PathSearchNode> *astarsearch )
		{
			auto suggestions = Energy::GuidedDeformation::suggestChildren(*this);

			for (auto suggestion : suggestions)
				astarsearch->AddSuccessor(PathSearchNode(suggestion));

			return suggestions.size();
		}

		bool IsGoal()
		{
			return unassigned.empty();
		}

		bool IsSameState( PathSearchNode &rhs )
		{
			if (this->mapping.empty()) return false;
			return this->mapping == rhs.mapping;
		}
	};

	std::vector< std::vector<Energy::SearchNode> > search(Energy::SearchNode & root)
	{
		// Prepare shapes
		Energy::GuidedDeformation::preprocess(root.shapeA.data(), root.shapeB.data());

		Energy::SearchNode start;
		start.shapeA = QSharedPointer<Structure::ShapeGraph>(new Structure::ShapeGraph(*root.shapeA.data()));
		start.shapeB = QSharedPointer<Structure::ShapeGraph>(new Structure::ShapeGraph(*root.shapeB.data()));
		start.assignments = root.assignments;
		start.unassigned = root.unassignedList();
		Energy::GuidedDeformation::applyAssignment(&start, true);

		auto startCopy = new Energy::SearchNode(start);

		PathSearchNode startNode(*startCopy);

		std::vector< std::vector<Energy::SearchNode> > result;

		/// Perform search:
		{
			int num_solutions = 1000;
			int max_open_set = 10000;

			AStarSearch<PathSearchNode> astarsearch(num_solutions);

			astarsearch.SetStartAndGoalStates(startNode, startNode);

			unsigned int SearchState;
			unsigned int SearchSteps = 0;

			do
			{
				SearchState = astarsearch.SearchStep(num_solutions, max_open_set);
				SearchSteps++;
			} while (SearchState == AStarSearch<PathSearchNode>::SEARCH_STATE_SEARCHING);

			// Convert to original type
			for (auto & row : astarsearch.solutions){
				std::vector<Energy::SearchNode> r;
				for (auto & element : row) r.push_back(element);
				result.push_back(r);
			}
		}

		return result;
	}
}
