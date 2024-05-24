#!/bin/bash
overwrite=true
directory="../sims/modern/output/hdf5/"
ncore=5
# these need to match the filelist below
# this is not a clever script, so easy to get this wrong
startTime=0
dt=900

resolution=100
projection=EPSG:32756
maskfile="../mesh/oti_smaller_extended_hires_mask.shp"
minmax="363000 430000 7380000 7430000"


function join { local IFS="$1"; shift; echo "$*"; }

function process_file {

    processing_file=${1}
    NUMBER=${processing_file%*.h5}  # retain the part before the h5 extension
    NUMBER=${NUMBER##*_}  # retain the part after the last _, which is the number
    echo $NUMBER

    # loop over variables
    file=${processing_file}
    # loop over variables with counter
    # create the raster ov the vtu
    python h5_2_raster.py --resolution ${resolution} --wd_mask ../sims/modern/bathymetry.h5 --velocity --min_max ${minmax} ${file} temp 
    gdalwarp -cutline ${maskfile} -s_srs ${projection} -of netCDF -r bilinear  -dstnodata -9999 -overwrite temp_u_${NUMBER}.xyz "${file}_coarse_u.nc"
    gdalwarp -cutline ${maskfile} -s_srs ${projection} -of netCDF -r bilinear  -dstnodata -9999 -overwrite temp_v_${NUMBER}.xyz "${file}_coarse_v.nc"
    rm temp_u_${NUMBER}.xyz
    rm temp_v_${NUMBER}.xyz

}

# loop through files twice; once to get the count and create the times
# then to process files, which is done in parallel.
FILES="${directory}/Velocity2d_0*.h5"
timeList=( )
count=0
for f in $FILES
do
    timeList+="$(($startTime + ($count * $dt))) "
    count=$[$count + 1]
done

for f in $FILES
do
    ((i=i%ncore)); ((i++==0)) && wait    
    # overwrite if the flag is set
    # only process if the output doesn't exist (and overwrite false)
    # we only check v, assume the equi u worked ok!
    if [[ ${overwrite} = true || ! -f "${f}_v.nc" ]] ; then
        echo ${f}        
        process_file ${f} &
    fi
done

# we can now merge all the .nc files created
# we can * as we've made the number 0ddd in Thetis so they should list in order
#ncecat ${directory}/*_500_u.nc -O "${directory}/Velocity2d_500_U.nc"
#ncecat ${directory}/*_500_v.nc -O "${directory}/Velocity2d_500_V.nc"

# then rename the record dimension as time
#ncrename -d record,time ${directory}/Velocity2d_500_U.nc
#ncrename -d record,time ${directory}/Velocity2d_500_V.nc

# adjust units of time to something sensible
#joined=$(join , ${timeList[@]})
#string="time[time]={${joined%,}}"
#ncap2 -s "${string}" ${directory}/Velocity2d_500_U.nc -O ${directory}/Velocity2d_500_U.nc
#ncap2 -s "${string}" ${directory}/Velocity2d_500_V.nc -O ${directory}/Velocity2d_500_V.nc
# we then have to alter the time to be a float
#ncap2 -s 'time=float(time)' --overwrite ${directory}/Velocity2d_500_U.nc ${directory}/Velocity2d_500_U.nc
#ncap2 -s 'time=float(time)' --overwrite ${directory}/Velocity2d_500_V.nc ${directory}/Velocity2d_500_V.nc


