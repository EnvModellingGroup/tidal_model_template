import shutil
from thetis import *
import os.path
import sys
import math
import numpy as np
import argparse
import re


#create function
def main():

    parser = argparse.ArgumentParser(
         prog="create_vtu_from_h5",
         description="""Will create a vtu files for a h5 file"""
    )
    parser.add_argument(
            'input_file',
            metavar='input_file',
            help='The h5 file to convert to pvtu'
            )
    parser.add_argument(
            '-v',
            '--verbose',
            action='store_true',
            help="Verbose output: mainly progress reports.",
            default=False
            )
    parser.add_argument(
            'output_file',
            metavar='output_file',
            help='The output file PVD file. The vtu/pvtus will be treated magically.'
            )
    args = parser.parse_args()
    verbose = args.verbose
    input_file = args.input_file
    output_file = args.output_file
    if verbose:
        PETSc.Sys.Print("loading file")  

    # work out the number of the output
    head, tail = os.path.split(input_file)    
    p = re.compile(r'\d{5}')
    # need to post_process tail (the filename) to check if it's not part of a timestep
    m = p.search(tail)
    if m:
        tail = tail[0:m.start()]
        timestep = m.group()
    else:
        timestep = None
    
    func_name = ""
    if ("Elevation2d" in tail): 
        func_name = "elev_2d"
    elif ("Velocity2d" in tail):
        func_name = "uv_2d"
    else:
        PETSc.Sys.Print("I can't read that file")  

    # get output dir from the outputname
    outputdir, outputfile = os.path.split(output_file)    
    if (outputdir == ""):
        outputdir = "."
    if verbose:
        PETSc.Sys.Print("saving "+func_name+" to "+output_file+" ts: "+str(timestep))
    with CheckpointFile(input_file, "r") as f:
        mesh2d = f.load_mesh()
        function = f.load_function(mesh2d,func_name)
        visu_space = exporter.get_visu_space(function.function_space())
        # create our exporter
        e = exporter.VTKExporter(visu_space, func_name, outputdir, outputfile, next_export_ix=int(timestep))
        e.set_next_export_ix(int(timestep))
        e.export(function)
        



if __name__ == "__main__":
    main()
