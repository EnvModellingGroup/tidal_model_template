#!/bin/bash
# pull directory as argument
directory=$1
ncore=$2
mesh=$3
velocity=true
bss=true 

resolution=8000
projection=EPSG:32630
maskfile="../../TM1mSLR/mesh/mask.shp"

function process_file {

    processing_file=${1}

    # loop over variables
    for (( i=0; i<${#varname[@]}; i++));
    do
        var=${varname[$i]}
        name=${names[$i]}
        file="${directory}/${processing_file}"
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

    rm temp.xyz

}

# do this in parts - the elev stats, then tidal. Then if user wants, velocity and BSS
# function names 
declare -a varname=("MinFS"
                    "MaxFS"
                     )
# The English equivalent of above - *same order*, include units
declare -a names=("Sim LAT (m)"
                  "Sim HAT (m)"
                     ); #moved this along to match others

# process this lot
process_file "temporal_stats_elev.h5" 

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
                 )

process_file "tidal_stats_scal.h5" #added =



if [ ${velocity} == true ]; then #added space between true and ]
    # function names 
    declare -a varname=("AveSpeed"
                    "MaxSpeed"
                     )
    # The English equivalent of above - *same order*, include units
    declare -a names=("Mean speed (m/s)"
                  "Max speed (m/s)"
                 )
    # process this lot
    process_file "temporal_stats_scal_vel.h5"
fi

if [ ${bss} == true ]; then #added space between true and ]
    # function names 
    declare -a varname=("AveBSS"
                    "MaxBSS"
                     )
    # The English equivalent of above - *same order*, include units
    declare -a names=("Mean BSS (kgm-1s-2)"
                  "Max BSS (kgm-1s-2)"
                 )
    # process this lot
    process_file "temporal_stats_bss.h5"
fi
