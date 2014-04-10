TEMPLATE = subdirs

# Libraries
SUBDIRS += NURBS
SUBDIRS += GlSplatRendererLib
SUBDIRS += Reconstruction
SUBDIRS += StructureGraphLib
SUBDIRS += AuctionLIB

# Main plugin
SUBDIRS += ddtt-plugin functional

# Aux. plugins
SUBDIRS += empty-mesh
