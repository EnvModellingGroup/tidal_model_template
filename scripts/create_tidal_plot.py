import uptide
import numpy as np
import matplotlib.pyplot as plt
import datetime
import csv
import pandas as pd
import os

constituents = ['M2', 'S2', "K1", "O1", "P1", "N2", "K2", "Q1" ]
tide = uptide.Tides(constituents)  # select which constituents to use
tide.set_initial_time(datetime.datetime(2005,11,11,10,0,0)) 
model_input = "model_gauges_oti_hires.csv"
tide_gauges = "../data/tide_gauge_2.csv"
t_end =  1548900
t_start = 172800
t_export = 900
colours = ['#377eb8', '#ff7f00', '#4daf4a',
           '#a65628', '#984ea3', '#f780bf',
           '#999999', '#e41a1c', '#dede00']
params = {
  'legend.fontsize': 12,
  'xtick.labelsize': 12,
  'ytick.labelsize': 12,
  'axes.labelsize' : 14,
  'figure.subplot.left' : 0.18,
  'figure.subplot.top': 0.82,
  'figure.subplot.right': 0.95,
  'figure.subplot.bottom': 0.15,
  'text.usetex' : True
    }
plt.rcParams.update(params)


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



# list which gauges want plotting
model_times = np.arange(t_start,t_end,t_export)
df = pd.read_csv(model_input, header=None)
# cols are the tide gauges in same order and names above (tg_order)
print(df.shape)
# start and end indices
start_idx = int(t_start / t_export)
end_idx = int(t_end / t_export)

# now loop over tide gauges and plot them.
for name in tg_order:
    # pull amplitude
    obs_amps = []
    obs_phases = []
    idx = tg_order.index(name)
    for t in constituents:
        obs_amps.append(float(tide_gauge_data[name][t+" amp"]))
        obs_phases.append(np.radians(float(tide_gauge_data[name][t+" phase"])))

    t = np.arange(t_start, t_end, t_export)
    eta = tide.from_amplitude_phase(obs_amps, obs_phases, t)

    fig_summary=plt.figure(figsize=(12.0,6),dpi=360)
    ax=fig_summary.add_subplot(111)

    obs_ln = ax.plot(t / 86400., eta, color="grey", lw=2, label="Tide gauge", alpha=0.4)
    mod_ln = ax.plot(np.array(model_times) / 86400., df.iloc[start_idx:end_idx, idx], color=colours[0],lw=1, label="Model")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Water height (m)")
    lns = mod_ln + obs_ln
    labs = [l.get_label() for l in lns]
    leg = ax.legend(lns, labs, loc='lower right',ncol=1)
    leg.get_frame().set_edgecolor('k')
    ax.set_title(name.replace("_"," "))
    plt.savefig("tidal_plot_hi/tidal_plots_mod_all_"+name+".png", dpi=180)
    plt.close()
