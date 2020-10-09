#!/usr/bin python3

# correct vehicle instances
# if cmap instance file contains entries with
# ".H07.G08.E09.D10.B11.A12" or similar as vehicle
# run this to correct it
# also fixes common formatting errors in the data


import csv

path_to_cmap_instance_file = "path/to/file.csv"
path_to_fixed_output_file = "path/to/output/file.csv"

infile = open(path_to_cmap_instance_file, "r")
outfile = open("C:\\Users\\thoma\\Documents\\Uni\\bachelor\\cMap\\cmap_instances_corrected2.csv", "w", newline='')

inreader = csv.reader(infile, delimiter=',', quotechar='"')
outwriter = csv.writer(outfile, delimiter=',', quotechar='"')


for line in inreader:
    if line[9].startswith("."):
        vehicles = line[9][1:].split(".")
        prefix = line[8].split(".")[0]
        if prefix.startswith("\'"):
            prefix = prefix[1:]
            line[8] = line[8][1:]
        for i in vehicles:
            v = prefix + "." + i
            newline = line
            newline[9] = v
            outwriter.writerow(newline)
    else:
        # get rid of common formatting error
        if line[9].startswith("\'"):
            line[9] = line[9][1:]
        if line[8].startswith("\'"):
            line[8] = line[8][1:]
        outwriter.writerow(line)
