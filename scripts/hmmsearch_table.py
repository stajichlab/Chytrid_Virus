#!/usr/bin/env python3

import csv, os, sys, re, gzip, argparse

evalue_cutoff = 1e-8
debug = 0
top_search_dir = "search"
parser = argparse.ArgumentParser(description='summary stats to train on viral genomes')

parser.add_argument('-j','--jginames',default="lib/jgi_names.tab",
                    help='JGI prefix to species')
parser.add_argument('-u','--nofilter',action='store_true', help="Use Unfiltered SPAdes files")

parser.add_argument('-l','--library',default="NCLDV",choices=['NCLDV', 'NCVOG'], help="NCLDV or NCVOG HMM library results to parse")
parser.add_argument('-p','--prediction',default="prodigal",choices=['prodigal', 'funannotate'], help="prodigal or funannotate proteins")
#parser.add_argument('-d','--dir',default="search/prodigal_NCLDV",
#                    help='Search results dir')
parser.add_argument('-g','--genesgff',default="gff3",
                    help='Genes in GFF3 format')

parser.add_argument('--out',default="results",help='Out folder')

args = parser.parse_args()

topdir = os.path.join(top_search_dir)
if args.prediction is 'prodigal':
    if args.unfilt:
        topdir = os.path.join(top_search_dir,"%s_%s_nofilter"%(args.prediction,args.library))
    else:
        topdir = os.path.join(top_search_dir,"%s_%s"%(args.prediction,args.library))
else:
    topdir = os.path.join(top_search_dir,"%s_%s"%("protein",args.library))

outScore = os.path.join(args.out,"%s_protein_%s_score.tsv"%(args.library,args.prediction))
outCtg = os.path.join(args.out,"%s_protein_%s_ctg.tsv"%(args.library,args.prediction))

table = {}
orgs = {}
jginames = {}

with open(args.jginames,"r") as fh:
    tsvin = csv.reader(fh,delimiter="\t")
    for row in tsvin:
        jginames[row[0]] = row[1]

gene2loc = {}
namematch = re.compile(r'ID=([^;]+)')
namepref  = re.compile(r'^[^\|]+\|\S+')

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
            hit = row[0] # gene name
            q = row[3]
            evalue = float(row[6])
            if evalue > evalue_cutoff:
                continue
            if q not in table:
                table[q] = {}
            if org not in table[q]:
                table[q][org] = [ [hit,evalue] ] # store first hit
            else:
                table[q][org].append( [hit,evalue] )# store next hit

orgnames = sorted(orgs.keys())

for infile in os.listdir(args.genesgff):
    if not (infile.endswith(".gff3.gz") or
            infile.endswith(".gff.gz")):
        continue
    inspecies = re.sub(r'((\.prodigal\.gff\.gz)|\.genes\.gff3\.gz)$','',infile)
    inspecies = re.sub(r'\.((JGI)|(ncbi))','',inspecies)
    if inspecies in jginames:
        inspecies = jginames[inspecies]

    if inspecies not in orgs:
        print("skipping species %s not in the hmmsearch results"%(inspecies))
        continue

    print(inspecies)
    gene2loc[inspecies] = {}
    i=0
    with gzip.open(os.path.join(args.genesgff,infile),'rb') as gff3:
        ctg_orf_counter = {}
        for line in gff3:
            gene = line.decode('ascii').strip()
            if gene.startswith('#'):
                continue
            row = gene.split("\t")

            if row[1].startswith("Prodigal"):
                if row[0] not in ctg_orf_counter:
                    ctg_orf_counter[row[0]] = 0

                ctg_orf_counter[row[0]] += 1
                name = "%s_%d"%(row[0],ctg_orf_counter[row[0]])
            else:
                name = row[-1]
                m = namematch.match(name)
                if m:
                    name = m.group(1)

                m = namepref.match(name)
                if not m:
                    pref = name.split("_")[0]
                    name = "%s|%s"%(pref,name)

            if name in gene2loc[inspecies]:
                print("found a gene name %s already stored - reading file %s"%(name,infile))
            else:
                gene2loc[inspecies][name] = [row[0],int(row[3]),int(row[4]),row[6]]

            if debug and i % 10000 == 0:
                print("storing %s as %s (%s)"%(name,gene2loc[name],
                                                         infile))
            i += 1



for org in orgnames:
    if org not in gene2loc:
        print("expected to find organism %s in gene2loc, but did not"%(org))

for o in gene2loc:
    if o not in orgs:
        print("expected to find organism %s in orgs, but did not"%(o))

header = ["%s_gene"%(args.library)]
header.extend(orgnames)

with open(outScore,"w") as outScorefh, open(outCtg,"w") as outCtgfh:
    outcsv = csv.writer(outScorefh,delimiter="\t",lineterminator="\n")
    outctgcsv = csv.writer(outCtgfh,delimiter="\t",lineterminator="\n")
    outcsv.writerow(header)
    outctgcsv.writerow(header)
    for fam in sorted(table.keys()):
        row = [fam]
        rowctg = [fam]
        for org in orgnames:
            if org in table[fam]:
                table[fam][org].sort(key=lambda x: x[1])
                besthit = table[fam][org][0]
                name = re.sub(r'-T\d+$','',besthit[0])
                evalue  = besthit[1]
                besthitlocation = "NOLOC-%s"%(name)
                if name in gene2loc[org]:
                    besthitlocation = "%s %s:%d..%d(%s)"%(name, gene2loc[org][name][0],
                                                       gene2loc[org][name][1],
                                                       gene2loc[org][name][2],
                                                       gene2loc[org][name][3])
                row.append("%s (%s hits)"%(besthit[1],len(table[fam][org])))
                rowctg.append(besthitlocation)
            else:
                row.append("")
                rowctg.append("")

        outcsv.writerow(row)
        outctgcsv.writerow(rowctg)
