#!/bin/bash
# pull directory as argument
directory=$1

# Variable names to pull out of the PVTUs
declare -a variables=("tidal_range" 
                      "M2_amp" 
                      "S2_amp"
                      "K1_amp"
                      "O1_amp"
                      "M2_phase"
                      "S2_phase"
                      "K1_phase"
                      "O1_phase"
                      "average_speed"
                      "max_speed"
                      "average_bss"
                      "max_bss"
                     )
# Variable names to pull out of the PVTUs
declare -a varname=("TidalRange" 
                      "M2_amp" 
                      "S2_amp"
                      "K1_amp"
                      "O1_amp"
                      "M2_phase"
                      "S2_phase"
                      "K1_phase"
                      "O1_phase"
                      "AveSpeed"
                      "MaxSpeed"
                      "AveBSS"
                      "MaxBSS"
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
    vname=${varname[$i]}
    echo ${var}
    file="${directory}/${var}_0.pvtu"
    echo "Dealing with ${file}"
   	# loop over variables with counter
    echo "   Rasterising ${var}"
    # create the raster ov the vtu
    python create_raster.py ${vname} ${file} temp.grd
    # create a filename
    filename="${directory}/${var}".nc
    #mask it
    gdalwarp -cutline "${directory}/../mesh/mask.shp" -s_srs EPSG:32630 -crop_to_cutline -of NetCDF -r bilinear  -dstnodata -9999 -overwrite temp.grd "${filename}"
    # rename the variable to something sensible
    ncrename -v Band1,"${var}" "${filename}"		
    # change the long name variable to contain decent info
    ncatted -O -a long_name,"${var}",o,c,"${name}" "${filename}"

done

