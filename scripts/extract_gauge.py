import numpy as np
from matplotlib import pyplot as plt
import datetime
import uptide
from firedrake import *
from thetis import *
import utm
import sys
import csv

tide_gauges = "gbr_gauges_utm56.csv"
run_dir = "/home/jh1889/filestore/eem530/SLR110/"
mesh_file = "test.msh"
t_end = 2419200 # 28 days (so 14 day plot) 3455100 # 40 days-ish
start_file = 1344 # 14 days
t_export = 900

mesh = Mesh(os.path.join(run_dir,mesh_file))

# now load in the tide gauge locations
# how long is the input? 
tide_gauge_data = {}
try:
    with open(tide_gauges, 'r') as csvfile:
        # need to read in a couple of lines, rather thana set of bytes
        # for sniffer to work properly
        temp_lines = csvfile.readline() + '\n' + csvfile.readline()
        dialect = csv.Sniffer().sniff(temp_lines, delimiters=",\t")
        csvfile.seek(0)
        # Expect file to be:
        # Name, X, Y, M2 amp, M2 phase, etc
        # Header should be as above, with capitalisation etc, but order is unimportant
        reader = csv.DictReader(csvfile)
        for row in reader:
            temp = dict(row) # copy
            temp.pop('Name') # remove name
            tide_gauge_data[row['Name']] = temp
except csv.Error:
    # it's not, so no idea what the heck has been thrown at us
    print("Sorry, I could not decipher your tide gauge data. Exiting.")
    print(tide_gauges[0])
    sys.exit(1)

# coords to grab data
gauge_locs = []
for loc in tide_gauge_data:
    location = (float(tide_gauge_data[loc]['Y']),float(tide_gauge_data[loc]['X']))
    gauge_locs.append((location[1],location[0]))

file_location = os.path.join(run_dir,'hdf5/') #location of the Elevation2d output files
t_n = int(t_end/t_export + 1)
thetis_times = t_export*np.arange(t_n) + t_export
P1 = FunctionSpace(mesh, "CG", 1)
P2 = FunctionSpace(mesh, "DG", 1)
elev = Function(P2, name='elev_2d')
elev_data_set = np.empty((t_n-start_file, len(gauge_locs)))
for i in range(start_file,t_n):
    print('Reading h5 files. Time ',i,i*t_export)
    checkpoint_file = DumbCheckpoint(file_location + '/Elevation2d_{:05}'.format(i),    mode=FILE_READ)
    checkpoint_file.load(elev)
    checkpoint_file.close()
    elev_data_set[i-start_file, :] = elev.at(gauge_locs,dont_raise=True)


np.savetxt("model_gauges_slr110.csv", elev_data_set, delimiter=",")
