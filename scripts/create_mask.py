import numpy as np
import qgis.core
import sys
import qmesh3
import argparse

def main():

    parser = argparse.ArgumentParser(
         prog="create_mask",
         description="""Create a polygon from a simulation boundary shapefile"""
    )
    parser.add_argument(
        '-v', 
        '--verbose', 
        action='store_true', 
        help="Verbose output: mainly progress reports.",
        default=False
    )
    parser.add_argument(
        'shapefile',
        help='The shapefile used in qmesh when making your mesh'
    )
    parser.add_argument(
        'output',
        help='The output mask shapefile'
    )


    args = parser.parse_args()
    verbose = args.verbose    
    shp_file = args.shapefile
    output_mask = args.output

    #Reading in the shapefile describing the domain boundaries, and creating a gmsh file.
    boundaries = qmesh3.vector.shapefileTools.Shapes()
    boundaries.fromFile(shp_file)
    loopShapes = qmesh3.vector.shapefileTools.identifyLoops(boundaries,
              isGlobal=False, defaultPhysID=1000,
              fixOpenLoops=True)
    polygonShapes = qmesh3.vector.shapefileTools.identifyPolygons(loopShapes, smallestMeshedArea=50000000, 
                                                                     meshedAreaPhysID = 1)
    polygonShapes.writeFile(output_mask)

if __name__ == "__main__":
    main()
