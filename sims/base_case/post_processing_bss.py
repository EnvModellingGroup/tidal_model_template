import numpy as np
import uptide
from thetis import *
import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
import params
from firedrake.petsc import PETSc

# where should the output of this analysis go
output_dir = 'analysis'
create_directory(output_dir)

# where is the output of your model?
thetis_dir = params.output_dir

# was this run created with the DumbCheckpoint code? If so, make this True
legacy_run = False

# You *MAY* need to edit below this line
# Make sure below matches your main run file as much as possible
# *if* anything goes wrong with the analysis
#============================================================#

# making an assumption here on what the hdf5 output is called
chk = CheckpointFile("output/hdf5/Elevation2d_00000.h5",'r')
thetis_mesh = chk.load_mesh()

chk = CheckpointFile('bathymetry.h5','r')
bathymetry2d = chk.load_function(thetis_mesh,'bathymetry')
chk.close()
chk = CheckpointFile('manning.h5','r')
manning = chk.load_function(thetis_mesh, 'manning')
chk.close()

# How long does your simulations run for (s)
t_end = params.end_time #params.end_time #40 days (i.e. 30 days of analysis)
# how often are exports produced in the main run?
t_export = params.output_time
# which is the start file?
t_start = params.spin_up 

# You shouldn't need to edit below here
#========================================
t_n = int((t_end - t_start) / t_export) + 1
thetis_times = t_start + t_export*np.arange(t_n)

P1 = FunctionSpace(thetis_mesh, "CG", 1)
# we need bathy and manning on the same mesh as the elev and vel
P1DG = FunctionSpace(thetis_mesh, "DG", 1)
manningdg = project(manning, P1DG)
bathydg = project(bathymetry2d, P1DG)

bathy = bathydg.dat.data[:].astype(np.single)
man = manningdg.dat.data[:].astype(np.single)
bss_data_set = np.empty((t_n,bathydg.dat.data.shape[0]),dtype="float32")
# we can now discard the functions
#del(manningdg)
print(bss_data_set.dtype, sys.getsizeof(bss_data_set))
#h.setrelheap()

count = 0
for t in thetis_times:
    iexport = int(t/t_export)
    filename = '{0:s}_{1:05d}'.format("Elevation2d", iexport)
    afile = CheckpointFile(os.path.join(thetis_dir,"hdf5",filename+".h5"), 'r')
    e = afile.load_function(thetis_mesh, "elev_2d")
    elev_data_set = e.dat.data[:].astype("float32")
    afile.close()
    # clean up after ourselves or we get a memory leak
    PETSc.garbage_cleanup(comm=afile._comm)
    filename = '{0:s}_{1:05d}'.format("Velocity2d", iexport)
    print("BSS: ",filename, end=" ", flush=True)
    afile = CheckpointFile(os.path.join(thetis_dir,"hdf5",filename+".h5"), 'r')
    uv = afile.load_function(thetis_mesh, "uv_2d")
    afile.close()
    PETSc.garbage_cleanup(comm=afile._comm)
    u_data_set = uv.dat.data[:,0].astype("float32")
    v_data_set = uv.dat.data[:,1].astype("float32")
    speed = np.sqrt(u_data_set*u_data_set + v_data_set*v_data_set)
    elev_bathy = elev_data_set + bathy
    elev_bathy[elev_bathy < 0.01] = 0.0
    tau_b = np.array(1024*9.81*man*man*speed*speed / (elev_bathy)**(1./3.))
    tau_b[ elev_bathy < 0.001] = 0.0 # we have < 1mm of water
    tau_b[ tau_b < 0.0 ] = 0.0 # we had no water (shouldn't happen due to above, but just in case)
    bss_data_set[count, :] = tau_b.astype("float32")
    elev_data_set = None
    speed = None
    u_data_set = None
    v_data_set = None
    elev_bathy = None
    tau_b = None
    e = None
    uv = None
    del(e,uv,elev_data_set,speed,u_data_set,v_data_set,elev_bathy,tau_b)
    count += 1

print(bss_data_set.dtype,sys.getsizeof(bss_data_set))
#print(h.heap())

print("Calculating averages", flush=True)
ave_bss = [] # ave of bss calc
max_bss = [] # max of bss calc
for i in range(bathydg.dat.data.shape[0]): # loop over nodes in the Function mesh
    tau_b = bss_data_set[:,i]
    ave_bss.append(np.mean(tau_b))
    max_bss.append(np.max(tau_b))

print("Saving checkpoints", flush=True)
# We then save all the scalar temporal stats in a single hdf5 file
with CheckpointFile(output_dir + '/temporal_stats_bss.h5', "w") as chk:
    chk.save_mesh(thetis_mesh)
    avebss = Function(P1DG, name="AveBSS")
    avebss.dat.data[:] = np.array(ave_bss)
    File( output_dir + '/average_bss.pvd').write(avebss)
    chk.save_function(avebss)
    maxbss = Function(P1DG, name="MaxBSS")
    maxbss.dat.data[:] = np.array(max_bss)
    File( output_dir + '/max_bss.pvd').write(maxbss)
    chk.save_function(maxbss, name='MaxBSS')

