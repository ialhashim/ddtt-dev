#include "DeformationPath.h"

DeformationPath::DeformationPath(){
    gcorr = NULL;
    weight = 0.0;
    idx = i = si = 0;
}

void DeformationPath::normalizeErrors()
{
	double minval = DBL_MAX;
	double maxval = -DBL_MAX;

	for(auto error : errors){
		minval = std::min(minval, error);
		maxval = std::max(maxval, error);
	}

	for(auto & error : errors) error = (error - minval) / (maxval - minval);
}

double DeformationPath::minWeight(const std::vector<DeformationPath> &paths)
{
    double minval = DBL_MAX;
    for(auto & p : paths) minval = std::min(minval, p.weight);
    return minval;
}

double DeformationPath::maxWeight(const std::vector<DeformationPath> &paths)
{
    double maxval = -DBL_MAX;
    for(auto & p : paths) maxval = std::max(maxval, p.weight);
    return maxval;
}

void DeformationPath::normalize(std::vector<DeformationPath> & paths)
{
	double minval = minWeight(paths);
	double maxval = maxWeight(paths);
	double range = maxval - minval;

	for(auto & p : paths) p.weight = (p.weight - minval) / range;
}

void DeformationPath::execute()
{
	scheduler = QSharedPointer<Scheduler>( new Scheduler );
	blender = QSharedPointer<TopoBlender>( new TopoBlender( gcorr, scheduler.data() ) );
	scheduler->executeAll();

	property["isReady"].setValue( true );
}
