#include "CorrespondenceGenerator.h"
#include "CorrespondencePrepare.h"
#include "PartCorresponder.h"
#include "myglobals.h"
#include "BasicTable.h"

CorrespondenceGenerator::CorrespondenceGenerator(ParticleMesh * pmeshA, ParticleMesh * pmeshB) : sA(pmeshA), sB(pmeshB)
{
	computedAssignments = segmentAssignFromGroupAssign( generateGroupAssignments() );
}

Assignments CorrespondenceGenerator::generateGroupAssignments()
{		
	Assignments generatedAssignments;

	auto & segsA = sA->property["segments"].value<Segments>();
	auto & segsB = sB->property["segments"].value<Segments>();

	auto boxA = sA->bbox(), boxB = sB->bbox();

	auto & groupsA = sA->property["groups"].value< std::vector< std::vector<size_t> > >();
	auto & groupsB = sB->property["groups"].value< std::vector< std::vector<size_t> > >();

	/// Build similarity matrix of groups:
	Eigen::MatrixXd similiarity = Eigen::MatrixXd::Ones( groupsA.size(), groupsB.size() + 1 );

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

	std::vector< std::vector<float> > data;
	for(int i = 0; i < similiarity.rows(); i++){
		std::vector<float> dataRow;
		for(int j = 0; j < similiarity.cols(); j++) 
			dataRow.push_back( similiarity(i,j) );
		data.push_back(dataRow);
	}
	if(false) showTableColorMap(data, true); // DEBUG

	// Parameters:
	double similarity_threshold = 0.25;
	int count_threshold = 1;

	// Collect good candidates
	Assignments candidates;
	for(size_t i = 0; i < similiarity.rows(); i++){
		QVector<size_t> candidate;
		for(size_t j = 0; j < similiarity.cols(); j++){
			if(similiarity(i,j) < similarity_threshold)
				candidate << j;
		}
		candidates << candidate;
	}

	Assignments assignments;
	cart_product(assignments, candidates);

	// Filter assignments
	Assignments filtered;
	for(auto & a : assignments)
	{
		QMap<size_t,size_t> counts;
		bool accept = true;
		auto NOTHING_SEGMENT = similiarity.cols()-1;

		for(auto i : a) counts[i]++;
		for(auto k : counts.keys()){
			if(k != NOTHING_SEGMENT && counts[k] > count_threshold){
				accept = false;
				break;
			}
		}

		if( accept ) filtered << a;
	}

	// DEBUG:
	if(true) saveToTextFile("assignments.txt", vecToString2(filtered));

	generatedAssignments = filtered;
	
	return generatedAssignments;
}

QVector<Pairings> CorrespondenceGenerator::segmentAssignFromGroupAssign( Assignments groupAssignments )
{
	QVector<Pairings> segmentAssignments;

	auto & groupsA = sA->property["groups"].value< std::vector< std::vector<size_t> > >();
	auto & groupsB = sB->property["groups"].value< std::vector< std::vector<size_t> > >();

	auto & segsA = sA->property["segments"].value<Segments>();
	auto & segsB = sB->property["segments"].value<Segments>();

	// Pre-compute bounding boxes for groups
	QVector<Eigen::AlignedBox3d> groupBoxA(groupsA.size()), groupBoxB(groupsB.size());
	for(size_t i = 0; i < groupsA.size(); i++ ) 
		for(auto sid : groupsA[i]) groupBoxA[i].extend( segsA[sid].property["bbox"].value<Eigen::AlignedBox3d>() );
	for(size_t j = 0; j < groupsB.size(); j++ ) 
		for(auto sid : groupsB[j]) groupBoxB[j].extend( segsB[sid].property["bbox"].value<Eigen::AlignedBox3d>() );

	// Transfer group assignments into member assignments
	for( auto groupAssign : groupAssignments )
	{
		Pairings segAssign;

		for(size_t i = 0; i < groupAssign.size(); i++)
		{
			size_t j = groupAssign[i];
			if(j > groupsB.size()-1) continue;

			auto & gA = groupsA[i];
			auto & gB = groupsB[j];

			QSet<size_t> seenB;

			// Forward:
			for(auto sid : gA){
				auto & boxI = segsA[sid].property["bbox"].value<Eigen::AlignedBox3d>();
				Vector3 centerI = (boxI.center() - groupBoxA[i].min()).array() / groupBoxA[i].sizes().array();

				QMap<double, size_t> dists;

				for(auto sjd : gB){
					auto & boxJ = segsB[sjd].property["bbox"].value<Eigen::AlignedBox3d>();
					Vector3 centerJ = (boxJ.center() - groupBoxB[j].min()).array() / groupBoxB[j].sizes().array();
					double d = (centerI - centerJ).norm();
					dists[d] = sjd;
				}

				auto sjd = dists.values().front();
				seenB << sjd;

				segAssign << qMakePair(sid, sjd);
			}

			// Backward:
			for(auto sjd : gB){
				if(seenB.contains(sjd)) continue;

				auto & boxJ = segsB[sjd].property["bbox"].value<Eigen::AlignedBox3d>();
				Vector3 centerJ = (boxJ.center() - groupBoxB[j].min()).array() / groupBoxB[j].sizes().array();

				QMap<double, size_t> dists;

				for(auto sid : gA){
					auto & boxI = segsA[sid].property["bbox"].value<Eigen::AlignedBox3d>();
					Vector3 centerI = (boxI.center() - groupBoxA[i].min()).array() / groupBoxA[i].sizes().array();
					double d = (centerI - centerJ).norm();
					dists[d] = sid;
				}

				auto sid = dists.values().front();
				segAssign << qMakePair(sid, sjd);
			}
		}

		segmentAssignments << segAssign;
	}

	return segmentAssignments;
}
