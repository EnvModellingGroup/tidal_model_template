import shutil
from thetis import *
import os.path
import sys
import math
import numpy as np
import argparse

#create function
def main():

    parser = argparse.ArgumentParser(
         prog="create_raster_from_h5",
         description="""Will create set of velocity raster files (U and V) for each h5 file from a thetis run"""
    )
    # decision on resolution and extent of raster
    parser.add_argument(
            '--resolution',
            help="Set the resolution of the XY file. Default is 1000m",
            type=float,
            default=1000
            )
    # guessing I need the below as it was the same in create_raster_vel
    parser.add_argument(
            'input_directory',
            metavar='input_directory',
            help='The output directory of your run. Which will be our input for this script'
            )
    parser.add_argument(
            'output_filexyz',
            metavar='output_file',
            help='The output raster file stub. _u and _t will be added to this, e.g. ~/path/12ka_u_111.xyz'
            )
    parser.add_argument(
            'mesh',
            metavar='mesh',
            help='The mesh file.'
            )
    parser.add_argument(
            'min_max',
            metavar='min_max',
            type=float,
            nargs=4,
            help='The min and max of x and y in order: minx maxx miny maxy'
            )
    parser.add_argument(
            '-v',
            '--verbose',
            action='store_true',
            help="Verbose output: mainly progress reports.",
            default=False
            )
    parser.add_argument(
            '--tend',
            type=float,
            help="End time of simulation or when you want to stop. Default 3455100 seconds",
            default=3455100
            )
    parser.add_argument(
            '--tstart',
            type=float,
            help="start time of simulation or when you want to start. Default 1209600 seconds",
            default=1209600
            )
    parser.add_argument(
            '--texport',
            type=float,
            help="Export time of the simulation. Default 900 seconds",
            default=900
            )
    args = parser.parse_args()
    verbose = args.verbose
    resolution = args.resolution
    input_dir = args.input_directory
    meshfile = args.mesh
    min_max = args.min_max
    x_min = min_max[0]
    x_max = min_max[1]
    y_min = min_max[2]
    y_max = min_max[3]
    output_filexyz = args.output_filexyz
    t_end = args.tend
    t_start = args.tstart
    t_export = args.texport


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

    #load Mesh
    # for testing against the rectabgle_tsunami test (which has hard coded mesh)
    #mesh = RectangleMesh(500,250,100000.,50000.)
    mesh = Mesh(meshfile)

    # and --- create solver ---
    # Dummy: not actually used. Just use it to store the h5 data in
    solverObj = solver2d.FlowSolver2d(mesh, Constant(100))
    options = solverObj.options
    options.use_nonlinear_equations = True
    options.simulation_export_time = t_export
    options.simulation_end_time = t_end
    options.output_directory =  input_dir
    options.check_volume_conservation_2d = True
    options.fields_to_export = ['uv_2d', 'elev_2d']
    options.fields_to_export_hdf5 = ['uv_2d', 'elev_2d']

    options.manning_drag_coefficient = Constant(1.0)
    options.horizontal_viscosity = Constant(1.0)
    options.coriolis_frequency = Constant(1.0)
    options.timestep = 1.0
    options.use_grad_div_viscosity_term = True
    options.use_wetting_and_drying = False
    # THis is where the p1DG is set up. Here for both elev and uv
    options.element_family = "dg-dg"
    options.swe_timestepper_type = 'CrankNicolson'
    options.use_grad_div_viscosity_term = True
    options.use_grad_depth_viscosity_term = False

    # work out the export numbers, t_n and the times
    t_n = range(int(t_start/t_export),int((t_end/t_export) + 1))
    thetis_times = np.arange(t_start, t_end+(t_export/2), t_export)

    # create empty 2D arrays - same size as the raster. One for u, one for v
    # these arrays should match the size of the *output* raster file
    u_data_set = np.empty(len(raster_coords))
    print(np.shape(u_data_set))
    v_data_set = np.empty(len(raster_coords))

    # loop over the requested h5 files, pulling out the velocity at
    # our raster points, then saving to xyz files for each output
    for t in t_n:
        PETSc.Sys.Print('Reading h5 files. Time ',t,t*t_export)
        solverObj.load_state(t)
        uv, elev = solverObj.timestepper.solution.split()
        uv_raster = uv.at(raster_coords,dont_raise=True)
        # if we have points outside mesh, we get a single none back
        # this messes up the stack, so swap them for a null array instead
        uv_raster = np.array([ np.array([None,None]) if d is None else d for d in uv_raster ])
        print(np.shape(uv_raster))
        # uv raster is a list of 2-element np arrays
        # stck so we can slice and extract the u's and v's
        print(uv_raster)
        u_data_set = np.stack(uv_raster,0)[:,0]
        v_data_set = np.stack(uv_raster,0)[:,1]

        # write u to file
        u_filename = output_filexyz + "_u_" + f'{t:05}' + ".xyz"
        with open(u_filename,"w") as f:
            for p,u in zip(raster_coords,u_data_set):
                f.write(str(p[0]) + "\t" + str(p[1]) + "\t" + str(u) + "\n")

        # write v file
        v_filename = output_filexyz + "_v_" + f'{t:05}' + ".xyz"
        with open(v_filename,"w") as f:
            for p,v in zip(raster_coords,v_data_set):
                f.write(str(p[0]) + "\t" + str(p[1]) + "\t" + str(v) + "\n")

    # as a niceity, we also write out the time steps to another file
    t_filename = output_filexyz + "_times_thetis.txt"
    with open(t_filename,"w") as f:
        for t in thetis_times:
            f.write(str(t) + "\n")





if __name__ == "__main__":
    main()
