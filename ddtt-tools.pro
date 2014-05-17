TEMPLATE = subdirs
CONFIG += ordered

# Libraries
#SUBDIRS += NURBS

# Plugins
#SUBDIRS += voxel_resampler

SUBDIRS  = agd \         	# Average geodesic distance
    symmetry
SUBDIRS += bdf   		# Biharmonic distance function
SUBDIRS += sdf   		# Shape Diameter Function
SUBDIRS += curvature	# Curvature
SUBDIRS += repair		# Basic mesh repair
SUBDIRS += segmentation	# Segmentation via skeletons
SUBDIRS += symmetry		# Basic symmetry analysis

#SUBDIRS += test        # Performance test

# Dependency map
