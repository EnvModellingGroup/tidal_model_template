#!/usr/bin/env python
from helper_functions import *
import numpy as np
import matplotlib.pyplot as plt
import datetime
import csv
import pandas as pd
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir, "sims")))
import params
import re
from ampang import *
import argparse
import uptide

def main():

    parser = argparse.ArgumentParser(
         prog="plot amp-phase",
         description="""Plot amp and phase against theoretical"""
    )
    parser.add_argument(
        '-v', 
        '--verbose', 
        action='store_true', 
        help="Verbose output: mainly progress reports.",
        default=False
    )
    parser.add_argument(
        '-m', 
        '--model', 
        help="the elevevation csv file (file only) from the extract_gague.py script. Defaults to model_gauges_elev.csv in the model dir",
        default="model_gauges_elev.csv"
    )
    parser.add_argument(
        '-t',
        '--tide',
        help="The tide gauge file, including path. Default is '../data/tide_gauges.csv'",
        default="../data/tide_gauges.csv"
    )
    parser.add_argument(
        '-s',
        '--stub',
        help='short string to append onto output filenames, before the extension, e.g. using "7" would make the output "thetis_vs_obs_7.pdf"'
    )
    parser.add_argument(
        'model_dir',
        help='The model run directory, e.g. ../sims/base_case/'
    )

    args = parser.parse_args()
    verbose = args.verbose    
    model_input = args.model
    tide_gauges = args.tide
    stub = args.stub
    model_dir = args.model_dir
    model_input = os.path.join(model_dir, model_input)

    #################################################
    # assumes you've run extract_guage.py and obtained the file for that

    # output dir is the run name, minus the csv file, which we discard
    output_dir, filename = os.path.split(model_input)
    # try make the output dir
    os.makedirs(os.path.join(output_dir,"Tidal_Validation"), exist_ok=True)

    constituents =  ['M2', 'S2', 'K1', 'O1']

    t_start = params.spin_up
    t_export = params.output_time
    t_end = params.end_time
    tide = uptide.Tides(constituents)  # select which constituents to use
    tide.set_initial_time(params.start_datetime + datetime.timedelta(hours=params.time_diff))

    if verbose:
        print("Reading in tidal gauges")
    # tide gauge data comes back as two nested dicts
    #{location_name: {M2Amp: x, M2Phase:, y, etc, etc}
    tide_gauge_data = read_tide_gauge_data(tide_gauges)

    model_data = {}
    df = pd.read_csv(model_input)
    thetis_times = df["Time"]

    for name in tide_gauge_data:
        # pull amplitude
        # Subtract mean
        thetis_elev = df[name]
        thetis_elev = thetis_elev - thetis_elev.mean()
        thetis_amplitudes, thetis_phases = uptide.analysis.harmonic_analysis(tide, thetis_elev, thetis_times)
        thetis_data = {}
        i = 0
        for c in constituents:
            thetis_data[c+"_amp"] = thetis_amplitudes[i]
            thetis_data[c+"_phase"] = thetis_phases[i]
            i += 1
        model_data[name] = thetis_data


    smin=0.5
    smax=1.5
    swdt=0.25
    tmin=-25.
    tmax=25
    twdt=15
    aaa='Ratio of Mag. (cal/obs)'
    bbb='Phase lag'
    refft=1.0
    refstd=1.+0.0j
        
    arng=np.array((smin-0.10,smax+0.1))
    prng=np.array((tmin-5.00,tmax+5.0))
    athc=np.array((smin,smax,swdt))
    pthc=np.array((tmin,tmax,twdt))

    # plot one per component
    for t in constituents:

        if verbose:
            print("Ploting: ",t)
        obs_amps = []
        obs_phases = []
        model_amps = []
        model_phases = []
        for l in tide_gauge_data:
            obs_amps.append(float(tide_gauge_data[l][t+" amp"]))
            obs_phases.append(np.radians(float(tide_gauge_data[l][t+" phase"])))
            model_amps.append(float(model_data[l][t+"_amp"]))
            model_phases.append(float(model_data[l][t+"_phase"]))

        obs_amps = np.array(obs_amps)
        obs_phases = np.array(obs_phases)
        model_amps = np.array(model_amps)
        model_phases = np.array(model_phases)
        index_to_remove = []
        for j in range(0,len(model_amps)):
            if np.isnan(model_amps[j]) or model_amps[j] > 10000 or model_amps[j] < 0.00001 or obs_amps[j] < 0.00001: #null value
                index_to_remove.append(j)

        model_amps = np.delete(model_amps, index_to_remove)
        obs_amps = np.delete(obs_amps, index_to_remove)
        model_phases = np.delete(model_phases, index_to_remove)
        obs_phases = np.delete(obs_phases, index_to_remove)

        # add check to make sure we have same number?

        
        # create comples array
        ratio_amps = model_amps/obs_amps
        phase_lag = model_phases-obs_phases
        ratios = ratio_amps*np.cos(phase_lag) -1j*ratio_amps*np.sin(phase_lag)
        fig = plt.figure(figsize=(8,8))
        dia = AmpPhsDiagram(reff=refft,fig=fig,amprange=arng,ampthck=athc,phsrange=prng,phsthck=pthc,amptitle=aaa,phstitle=bbb,rd_fmt="%.2f")
        contours=dia.add_contours(levs=[0.2,0.4,0.6],colors='0.5')
        plt.clabel(contours, inline=1, fontsize=10,fmt="%.1f")
        plt.title(t)
        plt.scatter(ratios.real,ratios.imag,marker='x', s=20, label='Mod/Obs')
        plt.scatter(1,0,marker='*',s=70, label='Reference')    
        plt.legend(loc="lower left",ncol=2,scatterpoints=1)
        name = t+"_ratio"
        if not stub is None:
            name = name +"_" + stub
        plt.savefig(os.path.join(output_dir,"Tidal_Validation",name+'.png'), dpi=180)

if __name__ == "__main__":
    main()
