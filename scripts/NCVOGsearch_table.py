#!/usr/bin/env python3

import csv, os, sys, re

topdir="search"
outVOG="results/VOG_hits.tsv"
table = {}
orgs = {}
for infile in os.listdir(topdir):
    if not infile.endswith(".TFASTX.tab"):
        continue
    stem = infile.split('.')
    org = stem[0]
    orgs[org] = 1
    with open(os.path.join(topdir,infile),"r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            row = line.split("\t")
            VOG = row[0]
            VOG = re.sub(r'-consensus','',VOG)
            contig = row[1]
            evalue = row[10]
            if VOG not in table:
                table[VOG] = {}
            if org not in table[VOG]:
                table[VOG][org] = [contig,evalue] # store first hit
            elif table[VOG][org][1] > evalue:
                table[VOG][org] = [contig,evalue] # store best hit

orgnames = sorted(orgs.keys())
header = ['VOG']
header.extend(orgnames)
with open(outVOG,"w") as outVOGfh:
    outcsv = csv.writer(outVOGfh,delimiter="\t",lineterminator="\n")
    outcsv.writerow(header)
    for VOG in sorted(table.keys()):
        row = [VOG]
        for org in orgnames:
            if org in table[VOG]:
                row.append(table[VOG][org][1])
            else:
                row.append("")

        outcsv.writerow(row)
