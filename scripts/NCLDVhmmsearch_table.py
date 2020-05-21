#!/usr/bin/env python3

import csv, os, sys, re, gzip
evalue_cutoff = 1e-8
debug = 0
jginamesfile ="lib/jgi_names.tab"
topdir="search/protein_NCLDV"
gff="gff3"
outVOG="results/NCLDV_proteinhits_score.tsv"
outVOGctg="results/NCLDV_proteinhits_ctgs.tsv"
table = {}
orgs = {}
jginames = {}
with open(jginamesfile,"r") as fh:
    tsvin = csv.reader(fh,delimiter="\t")
    for row in tsvin:
        jginames[row[0]] = row[1]

gene2loc = {}
namematch = re.compile(r'ID=([^;]+)')
namepref  = re.compile(r'^[^\|]+\|\S+')
for infile in os.listdir(gff):
    if not infile.endswith(".gff3.gz"):
        continue
    i=0
    with gzip.open(os.path.join(gff,infile),'rb') as gff3:
        for line in gff3:
            gene = line.decode('ascii').strip()
            if gene.startswith('#'):
                continue
            row = gene.split("\t")
            m = namematch.match(row[-1])
            name = row[-1]
            if m:
                name = m.group(1)
            m= namepref.match(name)
            if not m:
                pref = name.split("_")[0]
                name = "%s|%s"%(pref,name)

            if name in gene2loc:
                print("found a gene name %s already stored - reading file %s"%(name,infile))
            gene2loc[name] = [row[0],int(row[3]),int(row[4]),row[6]]
            if debug and i % 10000 == 0:
                print("storing %s as %s (%s)"%(name,gene2loc[name],
                                                         infile))
            i += 1

for infile in os.listdir(topdir):
    if not infile.endswith(".domtbl"):
        continue
    org = re.sub('(\.LCG)?\.hmmscan\.domtbl$','',infile)
    if org in jginames:
        org = jginames[org]
    orgs[org] = 1
    with open(os.path.join(topdir,infile),"r") as fh:
        # parsing domtbl format from HMMER
        for line in fh:
            if line.startswith("#"):
                continue
            line.strip("\n")
            row = line.split()
            VOG = row[0]
            hit = row[3] # gene name
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
header = ['NCLDV_gene']
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
                name = re.sub(r'-T\d+$','',besthit[0])
                evalue  = besthit[1]
                besthitlocation = "NOLOC-%s"%(name)
                if name in gene2loc:
                    besthitlocation = "%s %s:%d..%d(%s)"%(name, gene2loc[name][0],
                                                       gene2loc[name][1],
                                                       gene2loc[name][2],
                                                       gene2loc[name][3])
                row.append("%s (%s hits)"%(besthit[1],len(table[VOG][org])))
                rowctg.append(besthitlocation)
            else:
                row.append("")
                rowctg.append("")

        outcsv.writerow(row)
        outctgcsv.writerow(rowctg)
