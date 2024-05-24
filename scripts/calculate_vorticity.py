from thetis import *
import os.path
import argparse
import re
import glob
import csv
from mpi4py import MPI

# run like:
#
#mpiexec -n 5 python calculate_vorticity.py ../sims/modern/output/hdf5/ test.csv -v
#

def main():

    parser = argparse.ArgumentParser(
         prog="calc_vorticity_from_h5",
         description="""Will calculate L2 norm of vorticity from a h5 file"""
    )
    parser.add_argument(
            'input_dir',
            metavar='input_dir',
            help='The directory where the velocity (uv) h5 files are'
            )
    parser.add_argument(
            '-v',
            '--verbose',
            action='store_true',
            help="Verbose output: mainly progress reports.",
            default=False
            )
    parser.add_argument(
            '--output_file',
            metavar='output_file',
            help='The output PVD without extension. The vtu/pvtus will be created magically.'
            )
    parser.add_argument(
            "csv_output",
            metavar="csv_output",
            help="A csv file to store the time indeces and L2 norms"
            )
    args = parser.parse_args()
    verbose = args.verbose
    input_dir = args.input_dir
    output_file = args.output_file
    csv_output = args.csv_output

    rank = MPI.COMM_WORLD.Get_rank()
    

    output = False
    if output_file != None:
        output = True
    
    # get list of velocity h5 files
    h5_files = sorted(glob.glob(input_dir + '/Velocity2d*.h5'))
    l2_norms = []
    times_index = []
    
    for input_file in h5_files:
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

        if output:
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
            times_index.append(timestep)
            l2_norms.append(L2_norm)
            if output:
                visu_space = exporter.get_visu_space(vorticity.function_space())
                # create our exporter
                e = exporter.VTKExporter(visu_space, "vorticity", outputdir, outputfile, next_export_ix=int(timestep))
                e.set_next_export_ix(int(timestep))
                e.export(vorticity)
            
    if (rank == 0):
        with open(csv_output, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',',
                                quotechar='"', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(["Index","L2_norm"])
            for i, l in zip(times_index, l2_norms):
                writer.writerow([str(i), str(l)])
    

if __name__ == "__main__":
    main()
