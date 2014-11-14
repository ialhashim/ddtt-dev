#pragma once

#include <QStringList>
namespace Structure{ struct Graph; }

struct PropagateProximity
{
    static void prepareForProximity(Structure::Graph * graph);
    static void propagateProximity(const QStringList & fixedNodes, Structure::Graph * graph);
};
