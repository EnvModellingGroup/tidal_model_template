#!/bin/bash
# pull directory as argument
directory=$1
ncore=$2
mesh=$3

resolution=1000
projection=EPSG:32630
maskfile="../mesh/mask.shp"

# do this in two parts - the temporal stats, then tidal

# function names 
declare -a varname=(  "AveSpeed"
                      "MaxSpeed"
                      "AveBSS"
                      "MaxBSS"
                      "MinFS"
                      "MaxFS"
                     )

# The English equivalent of above - *same order*, include units
declare -a names=("Mean speed (m/s)"
                  "Max speed (m/s)"
                  "Mean BSS (kgm-1s-2)"
                  "Max BSS (kgm-1s-2)"
                  "Sim HAT (m)"
                  "Sim LAT (m)"
                 )

# loop over variables
for (( i=0; i<${#varname[@]}; i++));
do
    var=${varname[$i]}
    name=${names[$i]}
    file="${directory}/temporal_stats_scal.h5"
   	# loop over variables with counter
    echo "   Rasterising ${var}"
    # create the raster ov the vtu
    mpiexec -n ${ncore} python h5_2_raster.py --resolution ${resolution} ${file} ${mesh} temp --func ${var}
    # create a filename
    filename="${directory}/${var}".nc
    #mask it
    gdalwarp -cutline ${maskfile} -s_srs ${projection} -crop_to_cutline -of NetCDF -r bilinear  -dstnodata -9999 -overwrite temp.xyz "${filename}"
    # rename the variable to something sensible
    ncrename -v Band1,"${var}" "${filename}"		
    # change the long name variable to contain decent info
    ncatted -O -a long_name,"${var}",o,c,"${name}" "${filename}"

done

# function names
declare -a varname=("TidalRange" 
                      "M2_amp" 
                      "S2_amp"
                      "K1_amp"
                      "O1_amp"
                      "M2_phase"
                      "S2_phase"
                      "K1_phase"
                      "O1_phase"
                     )

# The English equivalent of above - *same order*, include units
declare -a names=("Tidal Range (m)"
                  "M2 amplitude (m)"
                  "S2 amplitude (m)"
                  "K1 amplitude (m)"
                  "O1 amplitude (m)"
                  "M2 phase (radians)"
                  "S2 phase (radians)"
                  "K1 phase (radians)"
                  "O1 phase (radians)"
                  "Mean speed (m/s)"
                  "Max speed (m/s)"
                  "Mean BSS (kgm-1s-2)"
                  "Max BSS (kgm-1s-2)"
                 )

# loop over variables
for (( i=0; i<${#varname[@]}; i++));
do
    var=${varname[$i]}
    name=${names[$i]}
    file="${directory}/tidal_stats_scal.h5"
   	# loop over variables with counter
    echo "   Rasterising ${var}"
    # create the raster ov the vtu
    mpiexec -n ${ncore} python h5_2_raster.py --resolution ${resolution} ${file} ${mesh} temp --func ${var}
    # create a filename
    filename="${directory}/${var}".nc
    #mask it
    gdalwarp -cutline ${maskfile} -s_srs ${projection} -crop_to_cutline -of NetCDF -r bilinear  -dstnodata -9999 -overwrite temp.xyz "${filename}"
    # rename the variable to something sensible
    ncrename -v Band1,"${var}" "${filename}"		
    # change the long name variable to contain decent info
    ncatted -O -a long_name,"${var}",o,c,"${name}" "${filename}"

done

