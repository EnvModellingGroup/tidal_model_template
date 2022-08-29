import numpy as np
import qgis.core
import sys
import qmesh3

# EDIT ME #
shp_file = "../mesh/simulation_boundaries.shp"

#############################################

output_mask = "../mesh/mask.shp"
#Reading in the shapefile describing the domain boundaries, and creating a gmsh file.
boundaries = qmesh.vector.shapefileTools.Shapes()
boundaries.fromFile(shp_file)
loopShapes = qmesh.vector.shapefileTools.identifyLoops(boundaries,
          isGlobal=False, defaultPhysID=1000,
          fixOpenLoops=True)
polygonShapes = qmesh.vector.shapefileTools.identifyPolygons(loopShapes, smallestMeshedArea=50000000, 
                                                                 meshedAreaPhysID = 1)
polygonShapes.writeFile(output_mask)
