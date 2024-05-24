#!/bin/bash

# run in a firedrake venv
# more cores == faster processing (to limit of your machine and mesh)
# wraps the h5_2_pvtu.py script to sort elevation and velocity VTUs 
# for paraview visualisation from thetis h5 files.

directory="test/hdf5/"
ncore=6

function process_file {

    file=${1}
    output_file="${directory}/../${2}"
    # loop over variables with counter
    # create the raster ov the vtu
    mpiexec -n ${ncore} python h5_2_pvtu.py ${file} ${output_file}

}

FILES="${directory}/Velocity2d_*.h5"
for f in $FILES
do
    echo ${f}
    process_file ${f} "Velocity2d"
done

FILES="${directory}/Elevation2d_*.h5"
for f in $FILES
do
    echo ${f}
    process_file ${f} "Elevation2d"
done


