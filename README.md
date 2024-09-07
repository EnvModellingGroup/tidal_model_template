# Tidal model template

This repo sets-up a series of standard tidal model as used by the group for a project.

Description of directories:
- data: put your rasters, gauge data, tidal forcing data, etc in here.
- mesh: put your shapefiles for meshing, rst files, .geo files and of course, your mesh here
- sims: put all of your model runs in here. 
- scripts: a bunch of scripts for analysing and visualising your tidal run

General workflow:
 - clone or fork this repo
 - Add it to your git or github
 - create a stagnant_water test in your sims folder
 - create directories in sims for each model you wish to run (e.g. with W&D, without W&D, mannings=0.03, mannings=0.02, etc)
 - commit these changes to your version
 - collect your data
 - create your mesh(es)
 - set up the various params.py
 - edit the job submission scripts to what your need
 - make sure everything is commited to your repo
 - if you're using viking, copy the directory structure there
 - once run, commit your final files. Use rsync to copy down from viking
 - carry out post-processing as required
 - run any scripts as required
 - commit any changes made

If you find any bugs in scripts, add an issue your fix to the original repo

