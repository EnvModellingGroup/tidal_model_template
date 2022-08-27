import uptide
import uptide.tidal_netcdf
import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
import params

# Which constiuents do you want?
constituents = params.constituents
tide = uptide.Tides(constituents)
# set your start date and time
tide.set_initial_time(start_datetime)

# point me at your FES file
grid_file_name = "../../data/fes_2014.nc"
tnci = uptide.tidal_netcdf.AMCGTidalInterpolator(tide,  grid_file_name)
#tnci.set_mask_from_fill_value('m2amp', 1.844674e+19)

#linear increase
tt = 86400 # 24 hours

#No need to edit below here

def set_tidal_field(elev, t, llvector):
    tnci.set_time(t)
    evector = elev.dat.data
    # create a linear ramping over the tt seconds
    for i,xy in enumerate(llvector):
        lon = xy[1]
        lat = xy[0]


        multiplier = t / tt
        if t >= tt:
            multiplier = 1.0
        try:
            # FES is in cms, hence the 100
            evector[i] = multiplier * tnci.get_val((lat, lon)) / 100.
        except uptide.netcdf_reader.CoordinateError:
            evector[i] = 0.

    return evector

