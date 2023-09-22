import numpy as np
import uptide
from thetis import *
import uptide.tidal_netcdf
import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
import params
from firedrake.petsc import PETSc
import gc

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

# Which tidal consituents to analyse?
constituent_order = ['M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1', 'M4', 'MS4', 'MN4'] #LEAVE
constituents = params.constituents
# alter the ones used based on the Rayleigh Criterion
constituents = uptide.select_constituents(constituents, t_end - params.spin_up)

# You shouldn't need to edit below here
#========================================
t_n = int((t_end - t_start) / t_export) + 1
thetis_times = t_start + t_export*np.arange(t_n)


P1 = FunctionSpace(thetis_mesh, "CG", 1)
# we need bathy and manning on the same mesh as the elev and vel
P1DG = FunctionSpace(thetis_mesh, "DG", 1)
manningdg = project(manning, P1DG)
bathydg = project(bathymetry2d, P1DG)

elev = Function(P1DG, name='elev_2d')
elev_data_set = np.empty((t_n, elev.dat.data.shape[0]),dtype=numpy.single)
bathy = bathydg.dat.data[:]
man = manningdg.dat.data[:]
# we can now discard the functions
del(manningdg)
del(bathydg)

count = 0
for t in thetis_times:
    iexport = int(t/t_export)
    filename = '{0:s}_{1:05d}'.format("Elevation2d", iexport)
    print("Elev",filename,end=" ")
    with CheckpointFile(os.path.join(thetis_dir,"hdf5",filename+".h5"), 'r') as afile:
        e = afile.load_function(thetis_mesh, "elev_2d")
        elev_data_set[count, :] = e.dat.data[:]
        PETSc.garbage_cleanup(comm=afile._comm)
        gc.collect()

    count += 1

print("Working out min/max")
max_fs = [] # maximum tide height
min_fs = [] # minimum tide height

for i in range(elev.dat.data.shape[0]): # loop over nodes in the Function mesh
    all_elev = np.array(elev_data_set[:,i])
    max_fs.append(np.max(all_elev))
    min_fs.append(np.min(all_elev))


print("about to do tidal amps")
# we now sort out the tidal components
detector_amplitudes = []
detector_phases = []
detector_maxfs = []
detector_minfs = []
detector_tidal_range = []

for i in range(elev.dat.data.shape[0]):
    thetis_elev = elev_data_set[:,i]
    tide = uptide.Tides(constituents)
    tide.set_initial_time(params.start_datetime)

    # Subtract mean
    thetis_elev = thetis_elev - thetis_elev.mean()
    thetis_amplitudes, thetis_phases = uptide.analysis.harmonic_analysis(tide, thetis_elev[:], thetis_times[:])
    
    detector_maxfs.append(max(thetis_elev[:]))
    detector_minfs.append(min(thetis_elev[:]))
    detector_tidal_range.append(max(thetis_elev) - min(thetis_elev))
    detector_amplitudes.append(thetis_amplitudes)
    detector_phases.append(thetis_phases)

print("tidal amps done")

# sort out the min, max and tidal range - save as h5 to rasterise
with CheckpointFile(output_dir + '/tidal_stats_scal.h5', "w") as chk:
    chk.save_mesh(thetis_mesh)
    tr = Function(P1DG, name="TidalRange")
    tr.dat.data[:] = np.array(detector_tidal_range)
    chk.save_function(tr)
    File( output_dir + '/tidal_range.pvd').write(tr)

    for i in constituents:
        amp = Function(P1DG, name= i +'_amp')
        phase = Function(P1DG, name= i +'_phase')
        phasepi = Function(P1DG, name = i+'_phasepi')
        amp.dat.data[:] = np.array(detector_amplitudes)[:,constituents.index(i)]
        phase.dat.data[:] = np.array(detector_phases)[:,constituents.index(i)]
        chk.save_function(amp)
        chk.save_function(phase)
        File( output_dir + '/' + i + '_amp.pvd').write(amp)
        File( output_dir + '/' + i + '_phase.pvd').write(phase)
        phasepi.dat.data[:] = np.arcsin(np.sin(phase.dat.data[:]))
        File( output_dir + '/' + i + '_phase_mod_pi.pvd').write(phasepi)
        chk.save_function(phasepi)

# We then save all the scalar temporal stats in a single hdf5 file
with CheckpointFile(output_dir + '/temporal_stats_elev.h5', "w") as chk:
    chk.save_mesh(thetis_mesh)
    maxfs = Function(P1DG, name="MaxFS")
    maxfs.dat.data[:] = np.array(max_fs)
    chk.save_function(maxfs, name='MaxFS')
    File( output_dir + '/max_fs.pvd').write(maxfs)
    minfs = Function(P1DG, name="MinFS")
    minfs.dat.data[:] = np.array(min_fs)
    chk.save_function(minfs, name='MinFS')
    File( output_dir + '/min_fs.pvd').write(minfs)

