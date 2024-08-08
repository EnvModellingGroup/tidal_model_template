import csv

def read_tide_gauge_data(tide_gauges):

    # read in tide gauge data, so let's do easiest to hardest...
    # how long is the input? 
    tide_gauge_data = {}
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
                tide_gauge_data[row['Name']] = temp
    except csv.Error:
        # it's not, so no idea what the heck has been thrown at us
        print("Sorry, I could not decipher your tide gauge data. Exiting.")
        print(tide_gauges[0])
        sys.exit(1)

    return(tide_gauge_data)
