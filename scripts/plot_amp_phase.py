#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import datetime
import csv
import pandas as pd
import os
import uptide
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir, "sims")))
import params
import re
from ampang import *

model_input = "../sims/base_case/model_gauges_elev.csv"
tide_gauges = "../data/uk_all_gagues_UTM30.csv"


constituents_to_plot = ["M2", "S2", "K1", "O1", "M4"]

#################################################
# assumes you've run extract_guage.py and obtained the file for that

# TODO: This is identical to analyse_tides. Perhaps someone can modulerise that?
constituents = params.constituents
start_date = params.start_datetime
tide = uptide.Tides(constituents)
tide.set_initial_time(start_date)

# output dir is the run name, minus the csv file, which we discard
output_dir, filename = os.path.split(model_input)
# try make the output dir
os.makedirs(os.path.join(output_dir,"Tidal_validation"), exist_ok=True)

# read in tide gauge data, so let's do easiest to hardest...
# how long is the input? 
tide_gauge_data = {}
model_data = {}
tg_order = [] # we need to do things in order, so this helps keep track
try:
    with open(tide_gauges, 'r') as csvfile:
        # need to read in a couple of lines, rather thana set of bytes
        # for sniffer to work properly
        temp_lines = csvfile.readline() + '\n' + csvfile.readline()
        dialect = csv.Sniffer().sniff(temp_lines, delimiters=",\t")
        csvfile.seek(0)
        # Expect file to be:
        # Name, X, Y, M2 amp, M2 phase, etc
        # Header should be as above, with capitalisation etc, but order is unimportant
        reader = csv.DictReader(csvfile)
        for row in reader:
            temp = dict(row) # copy
            temp.pop('Name') # remove name
            model_data[row['Name']] = {}            
            tide_gauge_data[row['Name']] = temp
            tg_order.append(row['Name'])
except csv.Error:
    # it's not, so no idea what the heck has been thrown at us
    print("Sorry, I could not decipher your tide gauge data. Exiting.")
    print(tide_gauges[0])
    sys.exit(1)

# tide gauge data comes back as two nested dicts
#{location_name: {M2Amp: x, M2Phase:, y, etc, etc}

ignore = []
model_data = {}
thetis_times = np.arange(params.spin_up, params.end_time, params.output_time)
df = pd.read_csv(model_input, header=None)

for name in tg_order:
    # pull amplitude
    idx = tg_order.index(name)
    # Subtract mean
    thetis_elev = df.iloc[:, idx]
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
tmin=-30.
tmax=30
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
for t in constituents_to_plot:

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
        if np.isnan(model_amps[j]) or model_amps[j] > 10000 or model_amps[j] < 0.00001: #null value
            index_to_remove.append(j)

    model_amps = np.delete(model_amps, index_to_remove)
    obs_amps = np.delete(obs_amps, index_to_remove)
    model_phases = np.delete(model_phases, index_to_remove)
    obs_phases = np.delete(obs_phases, index_to_remove)

    # add check to make sure we have same number?

    # create the complex array
    model = model_amps + 1j*model_phases
    obs = obs_amps + 1j*obs_phases
    ratio = model/obs

    fig = plt.figure(figsize=(8,8))
    dia = AmpPhsDiagram(reff=refft,fig=fig,amprange=arng,ampthck=athc,phsrange=prng,phsthck=pthc,amptitle=aaa,phstitle=bbb,rd_fmt="%.2f")
    contours=dia.add_contours(levs=[0.2,0.4,0.6],colors='0.5')
    plt.clabel(contours, inline=1, fontsize=10,fmt="%.1f")
    plt.title(t)
    plt.scatter(ratio.real,ratio.imag,marker='x', s=20, label='Mod/Obs')
    plt.scatter(1,0,marker='*',s=70, label='Reference')    
    plt.legend(loc="lower left",ncol=2,scatterpoints=1)
    plt.savefig(os.path.join(output_dir,"Tidal_validation",t+'_ratio.png'), dpi=180)


#plot all componenets on same?
