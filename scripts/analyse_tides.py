#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib
from scipy.stats import linregress
import csv
import sys
import os
import uptide
import pandas as pd
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir, "sims")))
import params

# EDIT ME #
plt.rcParams.update({'font.size': 22})

model_input = "../sims/base_case/model_gauges_elev.csv"
tide_gauges = "../data/uk_all_gagues_UTM30.csv"

constituents_to_plot = ["M2", "S2", "K1", "O1", "M4"]

#################################################
# assumes you've run extract_guage.py and obtained the file for that


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


# Now compare model to tide gagues
errors = []
av_err = []
i = 0
average_amp = []
for t in constituents:
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

    i+=1 
    
average_amp = np.array(average_amp)
errors = np.array(errors)
print("Error to Model:")
# create a prettier table using Pandas (as we've loaded this already)
dataframe = [errors/len(obs_amps), average_amp,  ((errors/len(obs_amps)) / average_amp)]
print(pd.DataFrame(dataframe, ["Error", "Gauge av. amp", "Relative err"], constituents))
print("\nNo. stations valid:", len(obs_amps), "\n")
with open(os.path.join(output_dir,"Tidal_validation","thetis_vs_stations.csv"), 'w') as f:
    writer = csv.writer(f)
    for lt,ln,name in zip(lats,lons,names):
        writer.writerow([name,ln,lt])


# fluidity to tide gauge plot
fig = plt.figure(figsize=(15,15),dpi=180)
n_plots = len(constituents_to_plot)
i = 0
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

    ax = fig.add_subplot(math.ceil(n_plots/2),2,i+1)

    # get StdDev for the observed amps
    std_dev = np.std(obs_amps)

    gradient, intercept, r_value, p_value, std_err = linregress(model_amps,obs_amps)
    ax.plot(model_amps,obs_amps,'bx')
    yLim = ax.get_ylim()
    xLim = ax.get_xlim()
    lim = max(np.amax(model_amps),np.amax(obs_amps))

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
    

    ax.set_title(t)
    # set equal, then adjust boundas
    plt.axis('equal')
    ax.set_xbound(lower=0,upper=lim)
    ax.set_ybound(lower=0,upper=lim)

    ax.set_xlabel("Model amp (m)")
    ax.set_ylabel("Gauges amp (m)")
    i += 1 #  counter for model data

    print(t, r_value, p_value, std_err)

plt.subplots_adjust(wspace=0.4,hspace=0.4)
plt.savefig(os.path.join(output_dir,"Tidal_validation","thetis_vs_obs.pdf"), dpi=180)
 

