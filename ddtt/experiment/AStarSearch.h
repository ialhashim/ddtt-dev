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
			auto allDists = EvaluateCorrespondence::hausdroffDistance(shapeA, shapeB);

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
			Energy::GuidedDeformation::applyAssignment(&successor, true);
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

			// Clean up:
			if (suggestions.size())
			{
				delete shapeA;
				delete shapeB;
				shapeA = shapeB = NULL;
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
		Energy::GuidedDeformation::preprocess(start.shapeA, start.shapeB);
		Energy::GuidedDeformation::applyAssignment(&start, true);
		PathSearchNode startNode(start);

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

		return new Energy::SearchNode(final_node.shapeA, final_node.shapeB, final_node.fixed, 
			final_node.assignments, final_node.unassigned, final_node.mapping, final_node.cost);
	}
}
