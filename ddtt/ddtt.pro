TEMPLATE = subdirs

# Libraries
SUBDIRS += NURBS
SUBDIRS += GlSplatRendererLib
SUBDIRS += Reconstruction
SUBDIRS += StructureGraphLib
SUBDIRS += experiment/libQDigraph
SUBDIRS += ProcessStructureGraphs

# Main plugin
SUBDIRS += experiment

# Standalone tool
SUBDIRS += experiment/standalone
