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

		bool GetSuccessors( AStarSearch<PathSearchNode> *astarsearch, PathSearchNode *parent_node )
		{
			auto suggestions = Energy::GuidedDeformation::suggestChildren(*this);

			for (auto suggestion : suggestions)
			{
				auto NewNode = PathSearchNode(suggestion);
				astarsearch->AddSuccessor(NewNode);
			}

			return true;
		}

		bool IsGoal( PathSearchNode & nodeGoal )
		{
			return unassigned.empty();
		}

		bool IsSameState( PathSearchNode &rhs )
		{
			return this->mapping == rhs.mapping;
		}
	};

	Energy::SearchNode * search(Energy::SearchNode & start)
	{
		// Prepare shapes
		Energy::GuidedDeformation::preprocess(start.shapeA.data(), start.shapeB.data());
		Energy::GuidedDeformation::applyAssignment(&start, true);

		auto startCopy = new Energy::SearchNode(start);

		PathSearchNode startNode(*startCopy);

		AStarSearch<PathSearchNode> astarsearch;

		astarsearch.SetStartAndGoalStates(startNode, startNode);

		unsigned int SearchState;
		unsigned int SearchSteps = 0;

		do
		{
			SearchState = astarsearch.SearchStep();
			SearchSteps++;
		} 
		while (SearchState == AStarSearch<PathSearchNode>::SEARCH_STATE_SEARCHING);

		auto * node = astarsearch.GetSolutionNode();

		PathSearchNode final_node = node->m_UserState;
		while (node = node->child)
		{
			final_node = node->m_UserState;
		}

		astarsearch.FreeSolutionNodes();

		return new Energy::SearchNode(new Structure::ShapeGraph(*final_node.shapeA.data()), new Structure::ShapeGraph(*final_node.shapeB.data()), final_node.fixed,
			final_node.assignments, final_node.unassigned, final_node.mapping, final_node.cost, final_node.energy);
	}
}
