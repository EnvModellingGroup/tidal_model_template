import numpy as np
from helper_functions import *
import datetime
from thetis import *
import sys
import csv
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir, "sims")))
import params
from firedrake.petsc import PETSc
import gc
import argparse

petsc_options = PETSc.Options()
petsc_options["options_left"] = False

def main():

    parser = argparse.ArgumentParser(
         prog="extract gauges",
         description="""Extract elevation, velocity and BSS from a thetis run at a number of gauges"""
    )
    parser.add_argument(
        '-v', 
        '--verbose', 
        action='store_true', 
        help="Verbose output: mainly progress reports.",
        default=False
    )
    parser.add_argument(
        '-b', 
        '--bss', 
        action='store_true', 
        help="Also extract BSS? Assumues you have a bunch of bss.h5 files for each output step",
        default=False
    )
    parser.add_argument(
        '--velocity', 
        action='store_true', 
        help="Also extract velocity?",
        default=False
    )
    parser.add_argument(
        '-t',
        '--tide',
        help="The tide gauge file, including path. Default is '../data/tide_gauges.csv'",
        default="../data/tide_gauges.csv"
    )
    parser.add_argument(
        '-s',
        '--stub',
        help='short string to append onto output filenames, before the extension, e.g. using "7" would make the output "thetis_vs_obs_7.pdf"'
    )
    parser.add_argument(
        'model_dir',
        help='The model run directory, e.g. ../sims/base_case/'
    )

    args = parser.parse_args()
    verbose = args.verbose    
    tide_gauges = args.tide
    stub = args.stub
    thetis_dir = args.model_dir
    calc_bss = args.bss
    calc_velocity = args.velocity

    legacy_run = False # will not work with bss with this true...

    if verbose:
        print("Reading in tidal gauges")
    # tide gauge data comes back as two nested dicts
    #{location_name: {M2Amp: x, M2Phase:, y, etc, etc}
    tide_gauge_data = read_tide_gauge_data(tide_gauges)

    # coords to grab data
    gauge_locs = []
    for loc in tide_gauge_data:
        location = (float(tide_gauge_data[loc]['Y']),float(tide_gauge_data[loc]['X']))
        gauge_locs.append((location[1],location[0]))

    # How long does your simulations run for (s)
    t_end = params.end_time
    t_export = params.output_time
    t_start = params.spin_up
    t_n = int((t_end - t_start) / t_export) + 1
    thetis_times = t_start + t_export*np.arange(t_n)
    
    # this is a big assumption in terms of which h5 file exists. May need updating in future.
    chk = CheckpointFile(os.path.join(thetis_dir,params.output_dir,"hdf5/Elevation2d_00000.h5"),'r')

    thetis_mesh = chk.load_mesh()
    P1DG = FunctionSpace(thetis_mesh, "DG", 1)
    t_n = int((t_end - t_start) / t_export) + 1
    if calc_velocity:
        speed = Function(P1DG, name='speed')
        uv = Function(P1DG, name='vel_2d')  
        speed_data = np.empty((t_n,  len(gauge_locs)))
    if calc_bss:
        bss_data = np.empty((t_n,  len(gauge_locs)))

    elev = Function(P1DG, name='elev_2d')
    elev_data = np.empty((t_n,  len(gauge_locs)))

    count = int(t_start / t_export) #orig count = 0
    index = 0 
    for t in thetis_times:
        iexport = int(t/t_export)
        filename = '{0:s}_{1:05d}'.format("Elevation2d", iexport)
        if verbose:
            PETSc.Sys.Print("Elev",filename)
        with CheckpointFile(os.path.join(thetis_dir,params.output_dir,"hdf5",filename+".h5"), 'r') as afile:
            e = afile.load_function(thetis_mesh, "elev_2d")
            elev_data[index, :] = e.at(gauge_locs,dont_raise=True)
            PETSc.garbage_cleanup(comm=afile._comm)
            gc.collect()
            if calc_velocity:
                filename = '{0:s}_{1:05d}'.format("Velocity2d", iexport)
                if verbose:
                    PETSc.Sys.Print("Vel:",filename)
                with CheckpointFile(os.path.join(thetis_dir,params.output_dir,"hdf5",filename+".h5"), 'r') as afile:
                    uv = afile.load_function(thetis_mesh, "uv_2d")
                    u_data_set = uv.dat.data[:,0]
                    v_data_set = uv.dat.data[:,1]
                    speed.dat.data[:] = np.sqrt(u_data_set*u_data_set + v_data_set*v_data_set)
                    speed_data[index,:] = speed.at(gauge_locs,dont_raise=True)
                    PETSc.garbage_cleanup(comm=afile._comm)
                    gc.collect()
            if calc_bss:
                if verbose:
                    PETSc.Sys.Print('Reading BSS', index)
                with CheckpointFile(file_location + 'bss_{:05}.h5'.format(i), 'r') as chk:
                    bss = chk.load_function(thetis_mesh, "BSS")
                    bss_data[index, :] = bss.at(gauge_locs,dont_raise=True)
                    PETSc.garbage_cleanup(comm=chk._comm)
                    gc.collect()

            index = index + 1

    filename = "model_gauges_elev"
    if not stub is None:
        filename = filename + "_" + stub
    filename = filename + ".csv"
    elev_data = np.append(np.reshape(np.array(thetis_times),(-1,1)),elev_data,axis=1)
    header = ["Time"]
    header.extend(tide_gauge_data.keys())
    np.savetxt(os.path.join(thetis_dir,filename), elev_data, delimiter=",",header= ','.join(header),comments='')
    if calc_velocity:
        filename = "model_gauges_speed"
        if not stub is None:
            filename = filename + "_" + stub
        filename = filename + ".csv"
        speed_data = np.append(np.reshape(np.array(thetis_times),(-1,1)),speed_data,axis=1)
        np.savetxt(os.path.join(thetis_dir,filename), speed_data, delimiter=",",header= ','.join(header),comments='')
    if calc_bss:
        filename = "model_gauges_bss"
        if not stub is None:
            filename = filename + "_" + stub
        filename = filename + ".csv"
        bss_data = np.append(np.reshape(np.array(thetis_times),(-1,1)),bss_data,axis=1)
        np.savetxt(os.path.join(thetis_dir,filename), bss_data, delimiter=",",header= ','.join(header),comments='')

if __name__ == "__main__":
    main()
