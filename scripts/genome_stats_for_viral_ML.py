#!/usr/bin/env python3

import sys, os, csv, argparse, re
from Bio import SeqIO
from Bio.SeqUtils import GC
import numpy as np
import numpy.ma as ma
from scipy.stats.mstats import mode, gmean, hmean

parser = argparse.ArgumentParser(description='summary stats to train on viral genomes')

parser.add_argument('-i','--infile', type=argparse.FileType('r'),nargs='?',default=sys.stdin,
                    help='Input file for reading genome DNA (FASTA)')
parser.add_argument('-o','--outbase',required=True,help="Output file base")
parser.add_argument('-p','--prodigal',type=argparse.FileType('r'),required=True,help="Prodigal GFF")

args = parser.parse_args()

contigs = {}
total_length = 0
for record in SeqIO.parse(args.infile, "fasta") :
        coverage = 1
        m = re.search(r'cov_(\d+\.\d+)',record.id)
        if m:
            coverage = float(m.group(1))

        contigs[record.id] = {
            'ID':  record.id,
            'LENGTH': len(record),
            'COVERAGE': coverage,
            'GC': "%.2f"%(GC(record.seq)),
            }
        total_length += len(record)

gff_parse = csv.reader(args.prodigal,delimiter="\t")
ORF_set = {}
ORF_count = 0
ORF_cols = ['length','strand','score','start','end']
for ORF in gff_parse:
    if ORF[0].startswith("#"):
        continue
    chrom = ORF[0]
    if chrom not in ORF_set:
        ORF_set[chrom] = []
    start = int(ORF[3])
    stop   = int(ORF[4])
    strand = 1
    if ORF[6] is '-':
        strand = -1

    ORF_count += 1
    ORF_set[chrom].append({'length': abs(stop - start) + 1,
                          'strand': strand, 'score': float(ORF[5]),
                          'start': start, 'stop': stop
                          })

for ctg in ORF_set:
    if ctg not in contigs:
        sys.stderr.write("cannot find chrom %s in FASTA file, but is in the Prodigal file, make sure these match"%(ctg))
        exit()
    ctgorfs = []
    intergenic = []
    contigs[ctg]['ORF_COUNT']        = 0
    minusgenes = plusgenes = 0
    last_stop = 1
    contigs[ctg]['FRACTION_OVERLAPPING_ORFS']  = 0
    for ORF in sorted(ORF_set[ctg], key=lambda orf: orf['start']):
        contigs[ctg]['ORF_COUNT'] += 1
        ctgorfs.append(ORF['length']) # keep track of all the lengths

        # plus vs minus strand ratio calculation
        if ORF['strand'] == 1:
            plusgenes += 1
        else:
            minusgenes += 1

        # examine intergenic space
        if ORF['start'] < last_stop:
            contigs[ctg]['FRACTION_OVERLAPPING_ORFS'] += 1
        else:
            intergenic.append( ORF['start'] - last_stop + 1)
        last_stop = ORF['stop']

    # calculate ratio of + to - strand genes on a contig: total number of plus (calculated)
    if plusgenes > 0 and minusgenes > 0:
        contigs[ctg]['ORF_PLUS_STRAND_RATIO'] = "%.3f"%(plusgenes / minusgenes)
    else:
        contigs[ctg]['ORF_PLUS_STRAND_RATIO'] = 'NA'
    ctgorfs_np = np.array(ctgorfs)
    contigs[ctg]['ORF_PER_KB'] = "%.2f"%(1000 * contigs[ctg]['ORF_COUNT'] / contigs[ctg]['LENGTH'])
    contigs[ctg]['ORF_MEAN'] = "%.2f"%(np.mean(ctgorfs_np))
    contigs[ctg]['ORF_MEDIAN'] = "%.2f"%(np.median(ctgorfs_np))
    intergenic_np = np.array(intergenic)
    contigs[ctg]['INTERGENIC_MEAN'] = "%.2f"%(np.median(intergenic))
    contigs[ctg]['INTERGENIC_MEDIAN'] = "%.2f"%(np.median(intergenic))
    if contigs[ctg]['ORF_COUNT'] > 0:
        contigs[ctg]['FRACTION_OVERLAPPING_ORFS'] = "%.2f"%(contigs[ctg]['FRACTION_OVERLAPPING_ORFS'] / contigs[ctg]['ORF_COUNT'])

sys.stderr.write("There are %d contigs total length = %d\n"%(len(contigs),total_length))
sys.stderr.write("There are %d ORFs\n"%(ORF_count))

with open(args.outbase + ".sum_stat.tsv","w",newline='') as sumstat:
    fieldnames = ['ID', 'LENGTH','GC', 'COVERAGE','ORF_COUNT','ORF_PER_KB','ORF_MEAN','ORF_MEDIAN','ORF_SD','ORF_PLUS_STRAND_RATIO','INTERGENIC_MEAN','INTERGENIC_MEDIAN','INTERGENIC_SD','FRACTION_OVERLAPPING_ORFS']
    outtsv = csv.DictWriter(sumstat,delimiter="\t",quoting=csv.QUOTE_MINIMAL,fieldnames=fieldnames)
    outtsv.writeheader()
    for ctgname in contigs:
        outtsv.writerow(contigs[ctgname])
