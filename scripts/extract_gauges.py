import numpy as np
from matplotlib import pyplot as plt
import datetime
import uptide
from thetis import *
import sys
import csv
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir, "sims")))
import params
from firedrake.petsc import PETSc

# EDIT ME #
run_dir = "../sims/palaeoslip/"
tide_gauges = "../data/core_sites_sep2021.csv"
legacy_run = False # will not work with bss with this true...
extract_bss = False # only set true if you've done so in post-processing

# You *MAY* need to edit below this line
# Make sure below matches your main run file as much as possible
# *if* anything goes wrong with the analysis
#============================================================#

# now load in the tide gauge locations
# how long is the input? 
tide_gauge_data = {}
try:
    with open(tide_gauges, 'r') as csvfile:
        # need to read in a couple of lines, rather than a set of bytes
        # for sniffer to work properly
        temp_lines = csvfile.readline() + '\n' + csvfile.readline()
        dialect = csv.Sniffer().sniff(temp_lines, delimiters=",\t")
        csvfile.seek(0)
        # Expect file to be:
        # Name, X, Y, M2 amp, M2 phase, etc
        # Header should be as above, with capitalisation etc, but order is unimportant
        reader = csv.DictReader(csvfile,dialect=dialect)
        for row in reader:
            temp = dict(row) # copy
            temp.pop('Name') # remove name
            tide_gauge_data[row['Name']] = temp
except csv.Error:
    # it's not, so no idea what the heck has been thrown at us
    PETSc.Sys.Print("Sorry, I could not decipher your tide gauge data. Exiting.")
    PETSc.Sys.Print(tide_gauges[0])
    sys.exit(1)

# coords to grab data
gauge_locs = []
for loc in tide_gauge_data:
    location = (float(tide_gauge_data[loc]['Y']),float(tide_gauge_data[loc]['X']))
    gauge_locs.append((location[1],location[0]))


# How long does your simulations run for (s)
t_end = params.end_time
t_export = params.output_time
t_start = params.spin_up

# where are your thetis output files (note, do not include the hdf5 directory)
thetis_dir = params.output_dir

t_n = int((t_end - t_start + 1) / t_export)

# this is a big assumption in terms of which h5 file exists. May need updating in future.
chk = CheckpointFile(os.path.join(run_dir,"output/hdf5/Elevation2d_00000.h5"),'r')
mesh2d = chk.load_mesh()

# --- create dummy solver ---
solverObj = solver2d.FlowSolver2d(mesh2d, Constant(10.0))
options = solverObj.options
options.simulation_export_time = t_export
options.simulation_end_time = t_end
options.output_directory =  thetis_dir
options.manning_drag_coefficient = Constant(1.0)
options.horizontal_viscosity = Constant(1.0)
options.element_family = "dg-dg"
options.output_directory = os.path.join(run_dir,params.output_dir)


P1DG = FunctionSpace(mesh2d, "DG", 1)
speed = Function(P1DG, name='speed')

thetis_times = params.spin_up + params.output_time*np.arange(t_n)
elev_data = np.empty((t_n,  len(gauge_locs)))
speed_data = np.empty((t_n,  len(gauge_locs)))
bss_data = np.empty((t_n,  len(gauge_locs)))

count = int(t_start / t_export)
index = 0
for i in thetis_times:
    PETSc.Sys.Print('Reading h5 files. Time ',count,i)
    solverObj.load_state(count, legacy_mode=legacy_run)
    elev_data[index, :] = solverObj.fields.elev_2d.at(gauge_locs,dont_raise=True)
    u_data_set = solverObj.fields.uv_2d.dat.data[:,0]
    v_data_set = solverObj.fields.uv_2d.dat.data[:,1]
    speed_dat = np.sqrt(u_data_set[:]*u_data_set[:] + v_data_set[:]*v_data_set[:])
    speed.dat.data[:] = speed_dat
    speed_data[index, :] = speed.at(gauge_locs,dont_raise=True)
    count = count + 1
    index = index + 1

elev_data = np.append(np.reshape(np.array(thetis_times),(-1,1)),elev_data,axis=1)
speed_data = np.append(np.reshape(np.array(thetis_times),(-1,1)),speed_data,axis=1)

#add thetis times to file
np.savetxt(os.path.join(run_dir,"model_gauges_elev.csv"), elev_data, delimiter=",")
np.savetxt(os.path.join(run_dir,"model_gauges_speed.csv"), speed_data, delimiter=",")

count = 0
if (not legacy_run and calc_bss):
    # only if we're using the new checkpoint
    file_location = os.path.join(run_dir,'analysis/') #location of the BSS output files
    count = int(t_start / t_export)
    index = 0
    for i in thetis_times:
        PETSc.Sys.Print('Reading BSS files. Time ',count,i)
        with CheckpointFile(file_location + 'bss_{:05}.h5'.format(i), 'r') as chk:
            bss = chk.load_function(mesh2d, "BSS")
        bss_data[count, :] = bss.at(gauge_locs,dont_raise=True)
        count = count + 1
        index = index + 1

    bss_data = np.append(np.reshape(np.array(thetis_times),(-1,1)),bss_data,axis=1)
    np.savetxt(os.path.join(run_dir,"model_gauges_bss.csv"), bss_data, delimiter=",")

