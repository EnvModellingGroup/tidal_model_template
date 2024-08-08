#!/usr/bin/env python
from helper_functions import *
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib
from scipy.stats import linregress
import sys
import os
import uptide
import pandas as pd
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir, "sims")))
import params
import argparse

def main():

    parser = argparse.ArgumentParser(
         prog="analyse tides",
         description="""Analyse a thetis tidal run against known tidal gauge data"""
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
        '-n',
        '--names',
        action="store_true",
        help="Label points with gauge name",
        default=False
    )
    parser.add_argument(
        '-t',
        '--tide',
        help="The tide gauge file, including path. Default is '../data/tide_gauges.csv'",
        default="../data/tide_gauges.csv"
    )
    parser.add_argument(
        '-c',
        '--constituents',
        help="""Which constiuents to analyse. As a list, in order of importance. Default is M2, S2, K1, O1.
                Remember to add -- after this argument to stop counting""",
        nargs='+',
        default=["M2", "S2", "K1", "O1"]
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
    constituents_to_plot = args.constituents
    stub = args.stub
    model_dir = args.model_dir
    model_input = os.path.join(model_dir, model_input)
    add_names = args.names

    plt.rcParams.update({'font.size': 22})
    #################################################
    # assumes you've run extract_guage.py and obtained the file for that


    constituents = constituents_to_plot
    start_date = params.start_datetime
    tide = uptide.Tides(constituents)
    # TODO: and timezone
    tide.set_initial_time(start_date)

    # output dir is the run name, minus the csv file, which we discard
    output_dir, filename = os.path.split(model_input)
    # try make the output dir
    os.makedirs(os.path.join(output_dir,"Tidal_Validation"), exist_ok=True)

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

    fig = plt.figure(figsize=(15,15),dpi=180)
    n_plots = len(constituents_to_plot)
    # Now compare model to tide gagues, calculate error and create plot
    errors = []
    av_err = []
    average_amp = []
    r_vals = []
    p_vals = []
    std_errs = []
    constit = 0
    if verbose:
        print("Tidal analysis")

    for t in constituents_to_plot:
        if verbose:
            print("\tdealing with ", t)

        obs_amps = []
        obs_phases = []
        model_amps = []
        model_phases = []
        lons = []
        lats = []
        names = []
        for l in tide_gauge_data:
            obs_amps.append(float(tide_gauge_data[l][t+" amp"]))
            obs_phases.append(np.radians(float(tide_gauge_data[l][t+" phase"])))
            lats.append(float(tide_gauge_data[l]['Y']))
            lons.append(float(tide_gauge_data[l]['X']))
            names.append(l)
            model_amps.append(float(model_data[l][t+"_amp"]))
            model_phases.append(float(model_data[l][t+"_phase"]))
        
        obs_amps = np.array(obs_amps)
        obs_phases = np.array(obs_phases)
        model_amps = np.array(model_amps)
        model_phases = np.array(model_phases)
        lats = np.array(lats)
        lons = np.array(lons)
        names = np.array(names)
        
        index_to_remove = []
        for j in range(0,len(model_amps)):
            if np.isnan(model_amps[j]) or model_amps[j] > 10000 or model_amps[j] < 0.00001: #null value
                index_to_remove.append(j)

        model_amps = np.delete(model_amps, index_to_remove)
        model_phases = np.delete(model_phases, index_to_remove)
        obs_amps = np.delete(obs_amps, index_to_remove)
        obs_phases = np.delete(obs_phases, index_to_remove)
        lats = np.delete(lats, index_to_remove)
        lons = np.delete(lons, index_to_remove)
        names = np.delete(names, index_to_remove)
        average_amp.append(np.average(obs_amps))

        errors.append(uptide.analysis.error_analysis(model_amps, model_phases, obs_amps, obs_phases)[0])
        av_err.append(uptide.analysis.error_analysis(model_amps, model_phases, obs_amps, obs_phases)[1])
        ax = fig.add_subplot(math.ceil(n_plots/2),2,constit+1)

        # get StdDev for the observed amps
        std_dev = np.std(obs_amps)

        gradient, intercept, r_value, p_value, std_err = linregress(model_amps,obs_amps)
        ax.plot(model_amps,obs_amps,'bx')
        yLim = ax.get_ylim()
        xLim = ax.get_xlim()
        lim = max(np.amax(model_amps),np.amax(obs_amps))
        #for i, txt in enumerate(names):
        #    ax.annotate(txt, (model_amps[i], obs_amps[i]), fontsize=8)

        # 1:1 line and 1 sigma each way
        ax.plot((0.0,lim), (0,lim), 'k-')
        ax.plot((0.0,lim), (std_dev,lim+std_dev), 'k:')
        ax.plot((0.0,lim), (-std_dev,lim-std_dev), 'k:')

        # our line of best fit
        mn=np.min(model_amps)
        mx=np.max(model_amps)
        x1=np.linspace(mn,mx,50)
        y1=gradient*x1+intercept
        ax.plot(x1,y1,'-r')
        if add_names:
            for i, txt in enumerate(names):
                ax.annotate(txt, (model_amps[i], obs_amps[i]), fontsize=8)
        

        ax.set_title(t)
        # set equal, then adjust boundas
        plt.axis('equal')
        ax.set_xbound(lower=0,upper=lim)
        ax.set_ybound(lower=0,upper=lim)

        ax.set_xlabel("Model amp (m)")
        ax.set_ylabel("Gauges amp (m)")
        constit += 1 #  counter for model data

        r_vals.append(r_value)
        p_vals.append(p_value)
        std_errs.append(std_err)

    plt.subplots_adjust(wspace=0.4,hspace=0.4)
    filename = "thetis_vs_obs"
    if not stub is None:
        filename = filename + "_" + stub
    filename = filename + ".pdf"
    plt.savefig(os.path.join(output_dir,"Tidal_Validation",filename), dpi=180)

        
    average_amp = np.array(average_amp)
    errors = np.array(errors)
    if verbose:
        print("Error to Model:")
    # create a prettier table using Pandas (as we've loaded this already)
    dataframe = [errors/len(obs_amps), average_amp,  ((errors/len(obs_amps)) / average_amp), r_vals, p_vals, std_errs]
    skill = pd.DataFrame(dataframe, ["Error", "Gauge av. amp", "Relative err", "Pearson", "P-val", "Std. Err."], constituents_to_plot)
    pd.options.display.float_format = '{:,.4f}'.format
    if verbose:
        print(skill)
        print("\nNo. stations valid:", len(obs_amps), "\n")
    filename = "thetis_skill"
    if not stub is None:
        filename = filename + "_" + stub
    filename = filename + ".csv"
    with open(os.path.join(output_dir,"Tidal_Validation",filename), 'w') as f:
        skill.to_csv(f)

    model_tides = pd.DataFrame.from_dict(model_data, orient="index")
    filename = "thetis_tidal_data"
    if not stub is None:
        filename = filename + "_" + stub
    filename = filename + ".csv"
    with open(os.path.join(output_dir,"Tidal_Validation",filename), 'w') as f:
        model_tides.to_csv(f)

if __name__ == "__main__":
    main()
