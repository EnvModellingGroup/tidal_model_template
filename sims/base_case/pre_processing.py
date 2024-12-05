from thetis import *
import hrds
import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
import params

mesh2d = Mesh(os.path.join(os.path.pardir,os.path.pardir,params.mesh_file)) # mesh file

viscosity = 1.0 # viscosity, obvs. 1.0 is a decent value. 10 is high, 0.001 is very low. Lower means more 
# pretty eddies etc, but harder to solve, more likely to crash. Higher is more stable, but less realistic.

manning_drag = 0.025 # which value of drag?

# which boundary is the forced boundary. We increase the visc and manning there to aid stability
forcing_boundary = params.forcing_boundary

# what distance should be used for the boundary blending (m)
blend_dist = 50000

# first deal with bathymetry
with timed_stage('initialising bathymetry'):
    bathy = hrds.HRDS("../../data/gbr_400_utm56S.tif",rasters=['../../data/gbr_100_utm56S_cropped.tif','../../data/oti_bathy_utm56S_filled_cropped.tif'],distances=[500.0,10.0])
    bathy.set_bands()
    P1_2d = FunctionSpace(mesh2d, 'CG', 1)
    bathymetry2d = Function(P1_2d, name="bathymetry")
    xvector = mesh2d.coordinates.dat.data
    bvector = bathymetry2d.dat.data
    assert xvector.shape[0]==bvector.shape[0]
    for i, (xy) in enumerate(mesh2d.coordinates.dat.data):
        bvector[i] = -1.0 * bathy.get_val(xy)

chk = CheckpointFile('bathymetry.h5', 'w')
chk.save_mesh(mesh2d)
chk.save_function(bathymetry2d, name='bathymetry')
PETSc.Sys.Print("Done bathy")

# now create distance from boundary function
# typical length scale
L = Constant(1.0e3)
V = FunctionSpace(mesh2d, 'CG', 1)
v = TestFunction(V)
u = Function(V)
u.interpolate(Constant(0.0))
bcs = DirichletBC(V, Constant(0.0), forcing_boundary) #make sure this matches physicalID of open boundaries
solver_parameters = {
    'ksp_rtol': 1e-4,
}
# Before we solve the Eikonal equation, let's solve a Laplace equation to
# generate an initial guess
F = L**2*(inner(grad(u), grad(v))) * dx - v * dx
solve(F == 0, u, bcs)#, solver_parameters=solver_parameters)
solver_parameters = {
    'ksp_rtol': 1e-4,
}
# epss values set the accuracy (in meters) of the final "distance to boundary" function. To make
# more accurate add in extra iterations, eg, 500., 250., etc. This may result in the solver not
# converging.
epss = [100000., 50000., 10000., 7500., 5000., 2500.,1000.]
# solve Eikonal equations
for i, eps in enumerate(epss):
    PETSc.Sys.Print('Solving Eikonal with eps == ' + str(float(eps)))
    F = inner(sqrt(inner(grad(u), grad(u))), v) * dx - v * dx + eps*inner(grad(u), grad(v)) * dx
    solve(F == 0, u, bcs, solver_parameters=solver_parameters)

# for checking your blend distance
dist = Function(V, name='dist')
dist.interpolate(u)
File('dist.pvd').write(u)

# create a viscosity buffer
# you may need to tweak this a little (it should be as low as possible
# to get a stable model)
chk = CheckpointFile('viscosity.h5', 'w')
with timed_stage('initialising viscosity'):
    h_viscosity = Function(V, name='viscosity')
    h_viscosity.interpolate(max_value(viscosity, (viscosity*1000) * (1. - u / blend_dist)))
    chk.save_mesh(mesh2d)
    chk.save_function(h_viscosity, name='viscosity')

# create manning boundary of increased bottom friction. Note this is turned off
# but you may find the mode is more stable model with a small (10) increase
# however this may make your model less skillful
chk = CheckpointFile('manning.h5', 'w')
with timed_stage('initialising manning'):
    manning = Function(V, name='manning')
    manning.interpolate(max_value(manning_drag, (manning_drag*1) * (1. - u / blend_dist)))
    chk.save_mesh(mesh2d)
    chk.save_function(manning, name='manning')
