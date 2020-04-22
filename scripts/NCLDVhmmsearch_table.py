#!/usr/bin/env python3

import csv, os, sys, re
evalue_cutoff = 1e-10
jginamesfile ="lib/jgi_names.tab"
topdir="search/protein_NCLDV"
outVOG="results/NCLDV_proteinhits_score.tsv"
outVOGctg="results/NCLDV_proteinhits_ctgs.tsv"
table = {}
orgs = {}
jginames = {}
with open(jginamesfile,"r") as fh:
    tsvin = csv.reader(fh,delimiter="\t")
    for row in tsvin:
        jginames[row[0]] = row[1]

for infile in os.listdir(topdir):
    if not infile.endswith(".domtbl"):
        continue
    org = re.sub('(\.LCG)?\.hmmscan\.domtbl$','',infile)
    if org in jginames:
        org = jginames[org]
    orgs[org] = 1
    with open(os.path.join(topdir,infile),"r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            line.strip("\n")
            row = line.split()
            VOG = row[0]
            hit = row[3]
            evalue = float(row[6])
            if evalue > evalue_cutoff:
                continue
            if VOG not in table:
                table[VOG] = {}
            if org not in table[VOG]:
                table[VOG][org] = [ [hit,evalue] ] # store first hit
            else:
                table[VOG][org].append( [hit,evalue] )# store next hit

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
                table[VOG][org].sort(key=lambda x: x[1])
                besthit = table[VOG][org][0]
                row.append("%s (%s hits)"%(besthit[1],len(table[VOG][org])))
#                rowctg.append(table[VOG][org][0])
            else:
                row.append("")
#                rowctg.append("")

        outcsv.writerow(row)
        outctgcsv.writerow(rowctg)
