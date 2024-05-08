import numpy as np
from thetis import *
import uptide.tidal_netcdf
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
t_end = params.end_time #40 days (i.e. 30 days of analysis)
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

bathy = bathydg.dat.data[:]
man = manningdg.dat.data[:]
# we can now discard the functions
del(manningdg)
del(bathydg)

uv = Function(P1DG, name='vel_2d')
u_data_set = np.empty((t_n, uv.dat.data.shape[0]),dtype=numpy.single)
v_data_set = np.empty((t_n, uv.dat.data.shape[0]),dtype=numpy.single)
count = 0
for t in thetis_times:
    iexport = int(t/t_export)
    filename = '{0:s}_{1:05d}'.format("Velocity2d", iexport)
    print("Vel:",filename, end=" ")
    with CheckpointFile(os.path.join(thetis_dir,"hdf5",filename+".h5"), 'r') as afile:
        uv = afile.load_function(thetis_mesh, "uv_2d")
        u_data_set[count, :] = uv.dat.data[:,0]
        v_data_set[count, :] = uv.dat.data[:,1]
        PETSc.garbage_cleanup(comm=afile._comm)
    count += 1


ave_speed = [] # average over speeds
max_speed = [] # max over speeds
ave_vel = [] # vector of ave u and ave v
max_vel = [] # vector of when max speed occurs

for i in range(uv.dat.data.shape[0]): # loop over nodes in the Function mesh
    u_vel = np.array(u_data_set[:,i]) # timeseries of u, v and elevation
    v_vel = np.array(v_data_set[:,i])
    speed = np.sqrt(u_vel*u_vel + v_vel*v_vel)
    ave_speed.append(np.mean(speed))
    max_speed.append(np.max(speed))
    ave_vel.append([np.mean(u_vel), np.mean(v_vel)])
    max_vel.append([u_vel[np.argmax(speed)],v_vel[np.argmax(speed)]])

# We then save all the scalar temporal stats in a single hdf5 file
with CheckpointFile(output_dir + '/temporal_stats_scal_vel.h5', "w") as chk:
    chk.save_mesh(thetis_mesh)
    avespeed = Function(P1DG, name="AveSpeed")
    avespeed.dat.data[:] = np.array(ave_speed)
    File( output_dir + '/ave_speed.pvd').write(avespeed)
    chk.save_function(avespeed)
    maxspeed = Function(P1DG, name="MaxSpeed")
    maxspeed.dat.data[:] = np.array(max_speed)
    File( output_dir + '/max_speed.pvd').write(maxspeed)
    chk.save_function(maxspeed)
   
# now the vectors
with CheckpointFile(output_dir + '/temporal_stats_velvec.h5', "w") as chk:
    chk.save_mesh(thetis_mesh)
    P1DG = VectorFunctionSpace(thetis_mesh, "DG", 1)
    avevel = Function(P1DG, name='ave_vel')
    avevel.dat.data[:] = np.array(ave_vel)
    File( output_dir + '/average_vel.pvd').write(avevel)
    chk.save_function(avevel, name='AveVel')
    maxvel = Function(P1DG, name='max_vel')
    maxvel.dat.data[:] = np.array(max_vel)
    File( output_dir + '/max_vel.pvd').write(maxvel)
    chk.save_function(maxvel, name='MaxVel')

