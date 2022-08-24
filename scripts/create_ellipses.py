import numpy as np
from matplotlib import pyplot as plt
import datetime
import uptide
from firedrake import *
from thetis import *
import utm
import sys
import csv
from utide._solve import solve
import pandas as pd


tide_gauges = "gbr_gauges_utm56.csv"
run_dir = "/home/jh1889/filestore/eem530/SLR43/"
mesh_file = "test.msh"
t_end = 2419200 # 28 days (so 14 day plot) 3455100 # 40 days-ish
start_file = 1344 # 14 days
t_export = 900

mesh = Mesh(os.path.join(run_dir,mesh_file))

def plot_ellipse(SEMA, ECC, INC, PHA):
    # from https://github.com/pwcazenave/PyFVCOM/blob/master/PyFVCOM/tidal_ellipse.py
    """
    Ellipse plot subfunction.
    Converted to Python by Pierre Cazenave, October 2012.
    """

    i = 1j

    SEMI = SEMA * ECC
    Wp = (1 + ECC) / 2 * SEMA
    Wm = (1 - ECC) / 2 * SEMA
    THETAp = INC - PHA
    THETAm = INC + PHA

    # Convert degrees into radians
    THETAp = THETAp / 180 * np.pi
    THETAm = THETAm / 180 * np.pi
    INC = INC / 180 * np.pi
    PHA = PHA / 180 * np.pi

    # Calculate wp and wm.
    wp = Wp * np.exp(i * THETAp)
    wm = Wm * np.exp(i * THETAm)

    dot = np.pi / 36
    ot = np.arange(0, 2 * np.pi, dot)
    a = wp * np.exp(i * ot)
    b = wm * np.exp(-i * ot)
    w = a + b

    wmax = SEMA * np.exp(i * INC)
    wmin = SEMI * np.exp(i * (INC + np.pi / 2))

    # tidy up - plot ellipse, but chnage colour depending on
    # cw or acw (is a > b, anti. b>a cw)

    plt.plot(np.real(w), np.imag(w))
    plt.axis('equal')
    plt.hold('on')
    plt.plot([0, np.real(wmax)], [0, np.imag(wmax)], 'm')
    plt.plot([0, np.real(wmin)], [0, np.imag(wmin)], 'm')
    plt.xlabel('u')
    plt.ylabel('v')
    plt.plot(np.real(a), np.imag(a), 'r')
    plt.plot(np.real(b), np.imag(b), 'g')
    plt.plot([0, np.real(a[0])], [0, np.imag(a[0])], 'ro')
    plt.plot([0, np.real(b[0])], [0, np.imag(b[0])], 'go')
    plt.plot([0, np.real(w[0])], [0, np.imag(w[0])], 'bo')
    plt.plot(np.real(a[0]), np.imag(a[0]), 'ro')
    plt.plot(np.real(b[0]), np.imag(b[0]), 'go')
    plt.plot(np.real(w[0]), np.imag(w[0]), 'bo')
    plt.plot(np.real([a[0], a[0]+b[0]]), np.imag([a[0], a[0]+b[0]]), linestyle='--', color='g')
    plt.plot(np.real([b[0], a[0]+b[0]]), np.imag([b[0], a[0]+b[0]]), linestyle='--', color='r')

    for n in range(len(ot)):
        plt.hold('on')
        plt.plot(np.real(a[n]), np.imag(a[n]), 'ro')
        plt.plot(np.real(b[n]), np.imag(b[n]), 'go')
        plt.plot(np.real(w[n]), np.imag(w[n]), 'bo')

    plt.hold('off')
    # save to file and asjust font, etc. 
    # we can then add ellipses to a map via the coords in the csv (easy, right?)
    plt.savefig('ellipse.png')

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

xvector = mesh.coordinates.dat.data
# sort this vector to make it easier later
X = []
Y = []
for n,xy in enumerate(xvector):
    X.append(xy[0])
    Y.append(xy[1])

file_location = os.path.join(run_dir,'hdf5') #location of the Velocity2d output files
t_n = int(t_end/t_export + 1)
thetis_times = t_export*np.arange(t_n) + t_export
P1 = VectorFunctionSpace(mesh, "CG", 1)
P2 = VectorFunctionSpace(mesh, "DG", 1)
vel = Function(P2, name='uv_2d')
u_data_set = np.empty((len(gauge_locs), t_n-start_file))
v_data_set = np.empty((len(gauge_locs), t_n-start_file))
for i in range(start_file,t_n):
    print('Reading h5 file '+file_location+'/Velocity2d_{:05}'.format(i)+'. Time ',i,i*t_export)
    checkpoint_file = DumbCheckpoint(file_location + '/Velocity2d_{:05}'.format(i), mode=FILE_READ)
    checkpoint_file.load(vel)
    checkpoint_file.close()
    vv = vel.at(gauge_locs,dont_raise=True)
    det = 0
    for v in vv:
        if (v is None):
            u_data_set[det,i-start_file] = None
            v_data_set[det,i-start_file] = None
        else:
            u_data_set[det,i-start_file] = v[0]
            v_data_set[det,i-start_file] = v[1]
        det = det + 1


np.savetxt("model_gauges_SLR43_u.csv", u_data_set, delimiter=",")
np.savetxt("model_gauges_SLR43_v.csv", v_data_set, delimiter=",")

