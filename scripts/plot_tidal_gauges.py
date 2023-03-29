import uptide
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

model_input = "../sims/base_case/model_tide_gauges.csv"
tide_gauges = "../data/uk_all_gagues_UTM30.csv"

########################################################

def tex_escape(text):
    """
        :param text: a plain text message
        :return: the message escaped to appear correctly in LaTeX
    """
    conv = {
        '&': r'\&',
        '%': r'\%',
        '$': r'\$',
        '#': r'\#',
        '_': r'\_',
        '{': r'\{',
        '}': r'\}',
        '~': r'\textasciitilde{}',
        '^': r'\^{}',
        '\\': r'\textbackslash{}',
        '<': r'\textless{}',
        '>': r'\textgreater{}',
    }
    regex = re.compile('|'.join(re.escape(str(key)) for key in sorted(conv.keys(), key = lambda item: - len(item))))
    return regex.sub(lambda match: conv[match.group()], text)


constituents = params.constituents
tide = uptide.Tides(constituents)  # select which constituents to use
tide.set_initial_time(params.start_datetime) 

t_start = params.spin_up
t_export = params.output_time
t_end = params.end_time

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

# output dir is the run name, minus the csv file, which we discard
output_dir, filename = os.path.split(model_input)
# try make the output dir
os.makedirs(os.path.join(output_dir,"Tidal_Gauges"), exist_ok=True)

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

    # Prettier and fixes LaTeX issue
    nice_name = tex_escape(name.replace("_"," ").title())
    obs_ln = ax.plot(t / 86400., eta, color="grey", lw=2, label="Tide gauge", alpha=0.4)
    mod_ln = ax.plot(np.array(model_times)/86400., df.iloc[:, idx], color=colours[0], lw=1, label="Model")
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Water height (m)")
    lns = mod_ln + obs_ln
    labs = [l.get_label() for l in lns]
    leg = ax.legend(lns, labs, loc='lower right',ncol=1)
    leg.get_frame().set_edgecolor('k')
    ax.set_title(nice_name)
    plt.savefig(os.path.join(output_dir,"Tidal_Gauges","tidal_plot"+name+".png"), dpi=180)
    plt.close()
