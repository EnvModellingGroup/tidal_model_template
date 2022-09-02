#!/bin/bash
# pull directory as argument
directory=$1
ncore=$2
mesh=$3

resolution=5000
projection=EPSG:32630

# filenames to rasterise, without .h5 extension
declare -a varname=("tidal_range" 
                      "M2_amp" 
                      "S2_amp"
                      "K1_amp"
                      "O1_amp"
                      "M2_phase"
                      "S2_phase"
                      "K1_phase"
                      "O1_phase"
                      "ave_speed"
                      "max_speed"
                      "ave_bss"
                      "max_bss"
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
for (( i=0; i<${#variables[@]}; i++));
do
    var=${variables[$i]}
    name=${names[$i]}
    echo ${var}
    file="${directory}/${var}"
    echo "Dealing with ${file}"
   	# loop over variables with counter
    echo "   Rasterising ${var}"
    # create the raster ov the vtu
    mpiexec -n ${ncore} python h5_2_raster.py --resolution ${resolution} ${file} ${mesh} temp
    # create a filename
    filename="${directory}/${var}".nc
    #mask it
    gdalwarp -cutline "${directory}/../mesh/mask.shp" -s_srs ${projection} -crop_to_cutline -of NetCDF -r bilinear  -dstnodata -9999 -overwrite temp.xyz "${filename}"
    # rename the variable to something sensible
    ncrename -v Band1,"${var}" "${filename}"		
    # change the long name variable to contain decent info
    ncatted -O -a long_name,"${var}",o,c,"${name}" "${filename}"

done

