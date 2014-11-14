#pragma once

#include "PropagateProximity.h"

class Propagate
{
public:
    static void propagateSymmetry(const QStringList & fixedNodes, Structure::Graph * graph);
};
