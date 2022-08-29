from thetis import *
# this imports our tidal forcing. If you want fes, comment out and uncomment the FES line
import tidal_forcing_tpxo as tidal_forcing
#import tidal_forcing_fes as tidal_forcing
import utm
import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
import params

mesh2d = Mesh(os.path.join(os.path.pardir,os.path.pardir,params.mesh_file))

#timestepping options
dt = 180 # reduce if solver does not converge
t_export = 900
t_end = params.end_time
output_dir = params.output_dir
utm_zone = params.utm_zone
utm_band=params.utm_band
P1 = FunctionSpace(mesh, "CG", 1)
cent_lat = params.cent_lat
cent_lon = params.cent_lon

# read bathymetry code
chk = DumbCheckpoint('bathymetry', mode=FILE_READ)
bathymetry2d = Function(P1)
chk.load(bathymetry2d,  name='bathymetry')
chk.close()

#read viscosity / manning boundaries code
chk = DumbCheckpoint('viscosity', mode=FILE_READ)
h_viscosity = Function(bathymetry2d.function_space(), name='viscosity')
chk.load(h_viscosity)
chk.close()
chk = DumbCheckpoint('manning', mode=FILE_READ)
manning = Function(bathymetry2d.function_space(), name='manning')
chk.load(manning)
chk.close()

# function to set up the Coriolis force
# Depends on a "central" lat/lon point in
# your mesh
def coriolis(mesh, lat, lon):
    R = 6371e3
    Omega = 7.292e-5
    lat_r = lat * pi / 180.
    f0 = 2 * Omega * sin(lat_r)
    beta = (1 / R) * 2 * Omega * cos(lat_r)
    x = SpatialCoordinate(mesh)
    x_0, y_0, utm_zone, zone_letter = utm.from_latlon(lat, lon)
    coriolis_2d = Function(FunctionSpace(mesh, 'CG', 1), name="coriolis_2d")
    coriolis_2d.interpolate(f0 + beta * (x[1] - y_0))

    return coriolis_2d

#account for Coriolis code - mesh, centre lat, centre lon
coriolis_2d = coriolis(mesh, cent_lat, cent_lon)

# --- create solver ---
solverObj = solver2d.FlowSolver2d(mesh, bathymetry2d)
options = solverObj.options
options.use_nonlinear_equations = True
options.simulation_export_time = t_export
options.simulation_end_time = t_end
options.output_directory =  output_dir
options.check_volume_conservation_2d = True
options.fields_to_export = ['uv_2d', 'elev_2d']
options.fields_to_export_hdf5 = ['uv_2d', 'elev_2d']

options.manning_drag_coefficient = manning #the manning function we created in initialisation & loaded above
options.horizontal_viscosity = h_viscosity #the viscosity 'cushion' we created in initialisation & loaded above
options.coriolis_frequency = coriolis_2d
options.timestep = dt
options.use_automatic_wetting_and_drying_alpha = True
options.wetting_and_drying_alpha_min = Constant(0.5)
options.wetting_and_drying_alpha_max = Constant(75.0)
options.use_wetting_and_drying = True
options.element_family = "dg-dg"
options.swe_timestepper_type = 'DIRK22'


options.swe_timestepper_options.solver_parameters = {
      'snes_type': 'newtonls',
      'snes_rtol': 1e-3,
      'ksp_rtol': 1e-3,
      'snes_monitor': None,
      'ksp_converged_reason': None,
      'ksp_type': 'gmres',
      'pc_type': 'fieldsplit',
  }

    
# set boundary/initial conditions code
tidal_elev = Function(bathymetry2d.function_space())
solverObj.bnd_functions['shallow_water'] = {
    params.forcing_boundary: {'elev': tidal_elev},
    #set open boundaries to tidal_elev function
    1000: {'un': 0.0},
    #set closed boundaries to zero velocity
}


solverObj.assign_initial_conditions(uv=Constant((1.0e-10,0.0)))

#work out our coords in lat/lon to save doing this every timestep.
mesh2d = tidal_elev.function_space().mesh()
xvector = mesh2d.coordinates.dat.data
llvector = []
for i,xy in enumerate(xvector):
    ll = utm.to_latlon(xy[0], xy[1], utm_zone, utm_band)
    llvector.append(ll)

def update_forcings(t):
  with timed_stage('update forcings'):
    print_output("Updating tidal field at t={}".format(t))
    tidal_forcing.set_tidal_field(tidal_elev, t, llvector)
    print_output("Done updating tidal field")


solverObj.iterate(update_forcings=update_forcings)
