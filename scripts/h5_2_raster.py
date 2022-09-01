import shutil
from thetis import *
import os.path
import sys
import math
import numpy as np
import argparse


mapping = {
        'tidal_range': "TidalRange",
        'min_fs' : 'MinFS',
        'max_fs' : 'MaxFS',
        }

#create function
def main():

    parser = argparse.ArgumentParser(
         prog="create_raster_from_h5",
         description="""Will create a raster files for a h5 file"""
    )
    parser.add_argument(
            '--resolution',
            help="Set the resolution of the XY file. Default is 5000m",
            type=float,
            default=5000
            )
    parser.add_argument(
            'input_file',
            metavar='input_file',
            help='The h5 file to rasterise. Without the .h5'
            )
    parser.add_argument(
            'output_file',
            metavar='output_file',
            help='The output xyz file. We add the xyz extension'
            )
    parser.add_argument(
            'mesh',
            metavar='mesh',
            help='The mesh file.'
            )
    parser.add_argument(
            '--min_max',
            type=float,
            nargs=4,
            help='The min and max of x and y in order: minx maxx miny maxy. Default is whole domain'
            )
    parser.add_argument(
            '--velocity',
            action="store_true",
            help='Your h5 file is a velocity, so a U and V raster will be produced',
            )
    parser.add_argument(
            '-v',
            '--verbose',
            action='store_true',
            help="Verbose output: mainly progress reports.",
            default=False
            )
    args = parser.parse_args()
    verbose = args.verbose
    resolution = args.resolution
    input_file = args.input_file
    meshfile = args.mesh
    output_file = args.output_file
    velocity = args.velocity
    
    #load Mesh
    mesh2d = Mesh(meshfile)

    min_max = args.min_max
    if min_max == None:
        # produces coords per core. Damn.
        min_coords = np.amin(mesh2d.coordinates.dat.data_ro, axis=0)
        max_coords = np.amax(mesh2d.coordinates.dat.data_ro, axis=0)
        x_min = min_coords[0]
        y_min = min_coords[1]
        x_max = max_coords[0]
        y_max = max_coords[1]
        comm = COMM_WORLD
        # so we need to do an all reduce
        import mpi4py.MPI
        global_x_max = np.zeros(1)
        global_x_min = np.zeros(1)
        global_y_max = np.zeros(1)
        global_y_min = np.zeros(1)
        comm.Barrier()
        comm.Allreduce(x_max, global_x_max, op=MPI.MAX)
        comm.Allreduce(x_min, global_x_min, op=MPI.MIN)
        comm.Allreduce(y_max, global_y_max, op=MPI.MAX)
        comm.Allreduce(y_min, global_y_min, op=MPI.MIN)
        comm.Barrier()
        x_min = global_x_min[0]
        x_max = global_x_max[0]
        y_min = global_y_min[0]
        y_max = global_y_max[0]
    else:
        x_min = min_max[0]
        x_max = min_max[1]
        y_min = min_max[2]
        y_max = min_max[3]

    if verbose:
        PETSc.Sys.Print("Coords: ", x_min, y_min, x_max, y_max)

    # we create a list of 2D points. These are the centre points of the raster we want
    # hence the resolution/2
    # arange won't go past the max, so we don't need to account for resolution at this end
    x_coords = np.arange(x_min+(resolution/2.), x_max, resolution)
    y_coords = np.arange(y_min+(resolution/2), y_max, resolution)
    raster_coords = []
    # loop y, then x so we get the xyz format right
    for y in y_coords:
        for x in x_coords:
            raster_coords.append([x,y])

    # create our function, note this is assumed to be DG
    P1_2d = FunctionSpace(mesh2d, 'DG', 1)
    function = Function(P1_2d)

    # create empty 2D arrays - same size as the raster.
    # One for u, one for v
    # these arrays should match the size of the *output* raster file
    if (velocity):
        u_data_set = np.empty(len(raster_coords))
        v_data_set = np.empty(len(raster_coords))
    else:
        data_set = np.empty(len(raster_coords))

    # loop over the requested h5 files, pulling out the velocity at
    # our raster points, then saving to xyz files for each output
    head, tail = os.path.split(input_file)
    func_name = mapping[tail]
    if verbose:
        PETSc.Sys.Print('Reading h5 file: ', input_file)
        PETSc.Sys.Print('Looking for: ', func_name)

    with DumbCheckpoint(input_file, mode=FILE_READ, comm=function.comm) as f:
        f.load(function,name=func_name)
        data_raster = function.at(raster_coords,dont_raise=True)

        if velocity:
            # if we have points outside mesh, we get a single none back
            # this messes up the stack, so swap them for a null array instead
            data_raster = np.array([ np.array([None,None]) if d is None else d for d in data_raster ])
            # uv raster is a list of 2-element np arrays
            # stck so we can slice and extract the u's and v's
            u_data_set = np.stack(data_raster,0)[:,0]
            v_data_set = np.stack(data_raster,0)[:,1]

            # write u to file
            u_filename = output_filexyz + "_u_" + ".xyz"
            with open(u_filename,"w") as f:
                for p,u in zip(raster_coords,u_data_set):
                    if u == None:
                        u = -9999
                    f.write(str(p[0]) + "\t" + str(p[1]) + "\t" + str(u) + "\n")
            # write v file
            v_filename = output_file + "_v_" + ".xyz"
            with open(v_filename,"w") as f:
                for p,v in zip(raster_coords,v_data_set):
                    if v == None:
                        v = -9999
                    f.write(str(p[0]) + "\t" + str(p[1]) + "\t" + str(v) + "\n")
        else:
            data_set = np.stack(data_raster,0)
            filename = output_file + ".xyz"
            with open(filename,"w") as f:
                for p,u in zip(raster_coords,data_set):
                    if u == None:
                        u = -9999
                    f.write(str(p[0]) + "\t" + str(p[1]) + "\t" + str(u) + "\n")


if __name__ == "__main__":
    main()
