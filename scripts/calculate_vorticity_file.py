from thetis import *
import os.path
import argparse
import re
import glob
import csv
from mpi4py import MPI

# Run like
#
# mpiexec -n 5 python calculate_vorticity_file.py ../sims/modern/output/hdf5/Velocity2d_01000.h5 vorticity
#
#

def main():

    parser = argparse.ArgumentParser(
         prog="calc_vorticity_from_h5",
         description="""Will calculate vorticity from a h5 velocity file"""
    )
    parser.add_argument(
            'input_file',
            metavar='input_file',
            help='The velocity (uv) h5 file'
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
            help='The output PVD without extension. The vtu/pvtus will be created magically.'
            )
    args = parser.parse_args()
    verbose = args.verbose
    input_file = args.input_file
    output_file = args.output_file
    
    if verbose:
        PETSc.Sys.Print("loading file: ", input_file)  

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
    if ("Velocity2d" in tail):
        func_name = "uv_2d"
    else:
        PETSc.Sys.Print("I can't read that file")  

    # get output dir from the outputname
    outputdir, outputfile = os.path.split(output_file)    
    if (outputdir == ""):
        outputdir = "."
    if verbose:
        PETSc.Sys.Print("    saving "+func_name+" to "+output_file+" ts: "+str(timestep))
    with CheckpointFile(input_file, "r") as f:
        mesh2d = f.load_mesh()
        function = f.load_function(mesh2d,func_name)
        # create a p1-cg space for vorticity
        vorticity = Function(FunctionSpace(mesh2d, 'CG', 1), name="vorticity")
        vorticity_calc = diagnostics.VorticityCalculator2D(function, vorticity)
        vorticity_calc.solve()
        L2_norm = sqrt(assemble(inner(vorticity,vorticity)*dx))
        PETSc.Sys.Print("L2 norm: ", L2_norm)
        visu_space = exporter.get_visu_space(vorticity.function_space())
        # create our exporter
        e = exporter.VTKExporter(visu_space, "vorticity", outputdir, outputfile, next_export_ix=int(timestep))
        e.set_next_export_ix(int(timestep))
        e.export(vorticity)
            


if __name__ == "__main__":
    main()
