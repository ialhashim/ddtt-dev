#include "ParticleDeformer.h"
#include "ParticleCorresponder.h"
#include "myglobals.h"

#include "BasicTable.h"

ParticleDeformer::ParticleDeformer(ParticleMesh *pmeshA, ParticleMesh *pmeshB): sA(pmeshA), sB(pmeshB)
{
	generateAssignments();
}

void ParticleDeformer::generateAssignments()
{
	auto segsA = sA->property["segments"].value<Segments>();
	auto segsB = sB->property["segments"].value<Segments>();

	//auto idA = segsA.keys(), idB = segsB.keys();
	auto boxA = sA->bbox(), boxB = sB->bbox();

	auto groupsA = sA->property["groups"].value< std::vector< std::vector<size_t> > >();
	auto groupsB = sB->property["groups"].value< std::vector< std::vector<size_t> > >();

	/// Build similarity matrix of groups:
	Eigen::MatrixXd similiarity = Eigen::MatrixXd::Ones( groupsA.size(), std::max(groupsA.size() + 1, groupsB.size() + 1) );

	for(size_t i = 0; i < similiarity.rows(); i++){
		double minVal = DBL_MAX;
		for(size_t j = 0; j < similiarity.cols(); j++){
			double val = 0;

			if( j < groupsB.size() )
			{
				Eigen::AlignedBox3d groupBoxA, groupBoxB;

				for(auto sid : groupsA[i]) groupBoxA.extend( segsA[sid].property["bbox"].value<Eigen::AlignedBox3d>() );
				for(auto sid : groupsB[j]) groupBoxB.extend( segsB[sid].property["bbox"].value<Eigen::AlignedBox3d>() );

				Vector3 relativeCenterA = (groupBoxA.center()-boxA.min()).array() / boxA.sizes().array();
				Vector3 relativeCenterB = (groupBoxB.center()-boxB.min()).array() / boxB.sizes().array();

				val = ( relativeCenterA - relativeCenterB ).norm();

				minVal = std::min(minVal, val);
			}

			similiarity(i,j) = val;
		}
		double maxVal = similiarity.row(i).maxCoeff();
		if(maxVal == 0) maxVal = 1.0;

		for(size_t j = 0; j < groupsB.size(); j++)
			similiarity(i,j) = (similiarity(i,j) - minVal) / maxVal;
	}

	similiarity.array() /= similiarity.maxCoeff();

	std::vector< std::vector<float> > data;
	for(int i = 0; i < similiarity.rows(); i++){
		std::vector<float> dataRow;
		for(int j = 0; j < similiarity.cols(); j++) 
			dataRow.push_back( similiarity(i,j) );
		data.push_back(dataRow);
	}
	showTableColorMap(data, true);

	// Find Acceptable Permutations
	{
		double threshold = 0.5;

		bool isTranspose = false;
		if( similiarity.rows() > similiarity.cols() ) isTranspose = true;
		if( isTranspose ) similiarity = similiarity.transpose();

		std::vector< std::vector<size_t> > possibleSets;

		// Fill set
		std::vector< size_t > idxSet;
		for(size_t i = 0; i < similiarity.cols(); i++) idxSet.push_back(i);

		// Go over permutations
		while ( std::next_permutation(idxSet.begin(), idxSet.begin() + similiarity.rows()) )
		{
			// Evaluate set
			bool accept = true;
			for(size_t i = 0; i < idxSet.size(); i++){
				if( similiarity(i, idxSet[i]) > threshold ){
					accept = false;
					break;
				}
			}

			if( accept ) possibleSets.push_back( idxSet );
		}

		debugBox(possibleSets.size());
		debugBoxVec2( possibleSets, 20 );
	}
}
