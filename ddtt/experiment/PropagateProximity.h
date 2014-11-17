#pragma once

#include <QStringList>
namespace Structure{ struct Graph; struct ShapeGraph; }

struct PropagateProximity
{
    static void prepareForProximity(Structure::Graph * graph);
	static void propagate(const QStringList & fixedNodes, Structure::ShapeGraph * graph);
};
