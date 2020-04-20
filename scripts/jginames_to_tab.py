#!/usr/bin/env python3
import csv,re,os


topdir="lib"
outfile="lib/jgi_names.tab"
with open(outfile,"w") as ofh:
    writer = csv.writer(ofh,delimiter="\t",lineterminator='\n')

    for infile in os.listdir(topdir):
        if not infile.endswith("_names.csv"):
            continue
        with open(os.path.join(topdir,infile),"r") as ifh:
            reader = csv.reader(ifh,delimiter=",")
            for line in reader:
                name=re.sub(" ","_",line[1])
                name=re.sub(";","",name)
                name=re.sub("\r","",name)
                writer.writerow([line[2],name])
