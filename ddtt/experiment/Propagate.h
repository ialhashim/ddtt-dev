#pragma once
#include <QStringList>
namespace Structure{ class Graph; }

class Propagate
{
public:
    static void propagateProximity(const QStringList & fixedNodes, Structure::Graph * graph);
    static void propagateSymmetry(const QStringList & fixedNodes, Structure::Graph * graph);
};
