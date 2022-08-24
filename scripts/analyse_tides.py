#!/usr/bin/env python
import datetime
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib
import argparse
from scipy.stats import linregress
import csv
import sys
import os
import uptide
import utm

constituents = ['M2','S2','O1','K1']
startdate = datetime.datetime(2000, 1, 1, 0, 0)
tide = uptide.Tides(constituents)
tide.set_initial_time(datetime.datetime(2000,1,1,0,0,0))
plt.rcParams.update({'font.size': 22})

def main():


    parser = argparse.ArgumentParser(
         prog="analyse tidal model run",
         description="""Analyse a thetis SWE model against gauges. Assumes you have netcdfs of your tidal constituents."""
    )
    parser.add_argument(
        '-v', 
        '--verbose', 
        action='store_true', 
        help="Verbose output: mainly progress reports.",
        default=False
    )
    parser.add_argument(
        '--tide_gauge',
        metavar='tide_gauge',
        help='Tide gauge data. A CSV with sensible headers.',
    )
    parser.add_argument(
        '--model_input',
        metavar='model_input',
        help='Directory where the netcdf files live. The script assumes the names...'
    )
    parser.add_argument(
        'output_stub',
        metavar='output_stub',
        help="The output file stub. Two files will be generated. A cross plot of amplitudes (Model against Obs, CSV file of error measures (D_n)."
    )

    
    args = parser.parse_args()
    verbose = args.verbose
    model_input = args.model_input
    output_stub = args.output_stub
    tide_gauges = args.tide_gauge
 
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
    except csv.Error:
        # it's not, so no idea what the heck has been thrown at us
        print("Sorry, I could not decipher your tide gauge data. Exiting.")
        print(tide_gauges[0])
        sys.exit(1)

    # tide gauge data comes back as two nested dicts
    #{location_name: {M2Amp: x, M2Phase:, y, etc, etc}

    ignore = []
    # create model data in same way
    for harmonic in constituents:
        nci_amp = uptide.netcdf_reader.NetCDFInterpolator(os.path.join(model_input,harmonic+'_amp.nc'), ('y','x'), ('y', 'x'))
        nci_amp.set_field(harmonic+"_amp")
        nci_phase = uptide.netcdf_reader.NetCDFInterpolator(os.path.join(model_input,harmonic+'_phase.nc'), ('y','x'), ('y', 'x'))
        nci_phase.set_field(harmonic+"_phase")
        for loc in tide_gauge_data:
            #location = (float(tide_gauge_data[loc]['Lat']),float(tide_gauge_data[loc]['Lon']))
            location = (float(tide_gauge_data[loc]['Y']),float(tide_gauge_data[loc]['X']))
            X = location[1]
            Y = location[0]
            # convert to UTM
            # NO LONGER NEEDED. DATA IN UTM55S
            #X,Y,zn,zl = utm.from_latlon(location[0], location[1])
            # grab data
            try:
                amp = nci_amp.get_val((Y,X),allow_extrapolation=False)
                phase = nci_phase.get_val((Y,X),allow_extrapolation=False)
                model_data[loc][harmonic+"_amp"] = amp
                model_data[loc][harmonic+"_phase"] = phase
            except uptide.netcdf_reader.CoordinateError:
                print("Can't find ",loc, X,Y)
                ignore.append(loc)

    # remove locations not in the model
    for k in ignore:
        tide_gauge_data.pop(k, None)

    #print(model_data)
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
    print("Error to Model")
    print("components:", constituents)
    print("error:", errors/len(obs_amps))
    print("Gauge av amp:", average_amp)
    print("Relative err:",  ((errors/len(obs_amps)) / average_amp))
    print("No. stations valid:", len(obs_amps))
    with open(output_stub+"_thetis_stations.csv", 'w') as f:
        writer = csv.writer(f)
        for lt,ln,name in zip(lats,lons,names):
            writer.writerow([name,ln,lt])


    # fluidity to tide gauge plot
    fig = plt.figure(figsize=(15,15),dpi=180)
    i = 0
    for t in constituents:
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

        ax = fig.add_subplot(2,2,i+1)

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
        

        ax.set_title(constituents[i])
        # set equal, then adjust boundas
        plt.axis('equal')
        ax.set_xbound(lower=0,upper=lim)
        ax.set_ybound(lower=0,upper=lim)

        ax.set_xlabel("Model amp (m)")
        ax.set_ylabel("Gauges amp (m)")
        i += 1 #  counter for model data

        print(r_value, p_value, std_err)

    plt.subplots_adjust(wspace=0.4,hspace=0.4)
    plt.savefig(output_stub+"_thetis_obs.pdf", dpi=180)
     

if __name__ == "__main__":
    main()


