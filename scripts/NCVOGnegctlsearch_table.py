#!/usr/bin/env python3

import csv, os, sys, re

jginamesfile ="lib/jgi_names.tab"
topdir="search/negctl_genome"
outVOG="results/NCVOG_negctl_hits_score.tsv"
outVOGctg="results/NCVOG_negctl_hits_ctgs.tsv"
table = {}
orgs = {}
jginames = {}
with open(jginamesfile,"r") as fh:
    tsvin = csv.reader(fh,delimiter="\t")
    for row in tsvin:
        jginames[row[0]] = row[1]

for infile in os.listdir(topdir):
    if not infile.endswith(".TFASTX.tab"):
        continue
    stem = infile.split('.')
    org = stem[0]
    if org in jginames:
        org = jginames[org]
    orgs[org] = 1
    with open(os.path.join(topdir,infile),"r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            row = line.split("\t")
            VOG = row[0]
            VOG = re.sub(r'-consensus','',VOG)
            contig = row[1]
            evalue = float(row[10])
            if VOG not in table:
                table[VOG] = {}
            if org not in table[VOG]:
                table[VOG][org] = [contig,evalue] # store first hit
            elif table[VOG][org][1] > evalue:
                table[VOG][org] = [contig,evalue] # store best hit

orgnames = sorted(orgs.keys())
header = ['VOG']
header.extend(orgnames)
with open(outVOG,"w") as outVOGfh, open(outVOGctg,"w") as outVOGctgfh:
    outcsv = csv.writer(outVOGfh,delimiter="\t",lineterminator="\n")
    outctgcsv = csv.writer(outVOGctgfh,delimiter="\t",lineterminator="\n")
    outcsv.writerow(header)
    outctgcsv.writerow(header)
    for VOG in sorted(table.keys()):
        row = [VOG]
        rowctg = [VOG]
        for org in orgnames:
            if org in table[VOG]:
                row.append(table[VOG][org][1])
                rowctg.append(table[VOG][org][0])
            else:
                row.append("")
                rowctg.append("")

        outcsv.writerow(row)
        outctgcsv.writerow(rowctg)
