
#!/usr/bin/env python3

import sys, os, csv, argparse, re
from Bio import SeqIO
from Bio.SeqUtils import GC
import numpy as np
import numpy.ma as ma
from scipy.stats.mstats import mode, gmean, hmean

def parse_hmmer_dombl(filehandle):
    "Parse HMMER hmmsearch domtbl for VOG/DLV db search"
    table = {}
    for line in args.search_ncldv:
        if line[0].startswith("#"):
            continue
        line = line.strip("\n")
        row  = line.split()
        for row in search:
            hit = row[0]
            evalue = float(row[6])
            q = row[1]
            if hit not in table:
                table[hit] = {}
            if (q not in table[hit] or
                table[hit][q]['evalue'] > evalue):
                score = 0.0
                if row[7] is not '.':
                    score = float(row[7])
                table[row[1]][row[0]] = {'evalue': evalue,
                                            'score': score,
                                            'qalnlen': abs(int(row[20]) - int(row[19]))+1,
                                            }
        return table


ribosomal_evalue_cutoff = 1e-10
mitochondria_evalue_cutoff = 1e-5
parser = argparse.ArgumentParser(description='summary stats to train on viral genomes')

parser.add_argument('-i','--infile', type=argparse.FileType('r'),nargs='?',default=sys.stdin,
                    help='Input file for reading genome DNA (FASTA)')
parser.add_argument('-o','--outbase',required=True,help="Output file base")
parser.add_argument('-p','--prodigal',type=argparse.FileType('r'),required=True,help="Prodigal GFF")
parser.add_argument('-sd','--search_ncldv',type=argparse.FileType('r'),required=False,help="NCLDV HMMER domtbl (from prodigal peps)")
parser.add_argument('-sv','--search_ncvog',type=argparse.FileType('r'),required=False,help="NCVOG HMMER domtbl (from prodigal peps)")
parser.add_argument('-sv','--pfam',type=argparse.FileType('r'),required=False,help="Pfam HMMER domtbl (from prodigal peps)")
parser.add_argument('-rrna','--ribosomal',type=argparse.FileType('r'),required=False,help="rRNA BLASTN TAB")
parser.add_argument('-mt','--mitochondria',type=argparse.FileType('r'),required=False,help="mitochondria TFASTX TAB")

args = parser.parse_args()

MT = {}
if args.mitochondria:
    search = csv.reader(args.mitochondria,delimiter="\t")
    for row in search:
        evalue = float(row[10])
        if evalue < mitochondria_evalue_cutoff:
            if row[1] not in MT or MT[row[1]][2] < row[2]:
                MT[row[1]] = row

rRNA = {}
if args.ribosomal:
    search = csv.reader(args.ribosomal,delimiter="\t")
    for row in search:
        evalue = float(row[10])
        if evalue < ribosomal_evalue_cutoff:
            if row[0] not in rRNA or rRNA[row[0]][2] < row[2]:
                rRNA[row[0]] = row

ncldvhits = parse_hmmer_dombl(args.search_ncldv)
ncvoghits = parse_hmmer_dombl(args.search_ncvog)



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
            'ORF_COUNT': 0,
            'ORF_PER_KB': 0,
            'ORF_MEAN':   0,
            'ORF_MEDIAN': 0,
            'ORF_PLUS_STRAND_RATIO': 1,
            'INTERGENIC_MEAN': 0,
            'INTERGENIC_MEDIAN': 0,
            'FRACTION_OVERLAPPING_ORFS': 0,
            'NCVOG_HITS': 0,
            'NCLDV_HITS': 0,
            'CAPSID': 0,
            'RIBOSOMAL': 0,
            'MT': 0,
            }
        if record.id in ncvoghits:
            contigs[record.id]['NCVOG_HITS'] = len(ncvoghits[record.id])
        if record.id in ncldvhits:
            contigs[record.id]['NCLDV_HITS'] = len(ncldvhits[record.id])
            if 'capsid' in ncldvhits[record.id]:
                contigs[record.id]['CAPSID'] = 1
        if record.id in rRNA:
            contigs[record.id]['RIBOSOMAL'] = 1
        if record.id in MT:
            contigs[record.id]['MT'] = 1
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
    score = 0.0
    if ORF[5] is not '.':
        score = float(ORF[5])
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
    contigs[ctg]['ORF_COUNT']            = 0
    contigs[ctg]['NUM_ORF_NO_PFAM_HITS'] = 'NA'
    minusgenes = plusgenes = 0
    last_stop = 1
    contigs[ctg]['FRACTION_OVERLAPPING_ORFS']  = 0
    for ORF in sorted(ORF_set[ctg], key=lambda orf: orf['start']):
        contigs[ctg]['ORF_COUNT'] += 1
        ctgorfs.append(ORF['length'] / 1000) # keep track of all the lengths in kb

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
    contigs[ctg]['ORF_MEAN'] = "%.3f"%(np.mean(ctgorfs_np))
    contigs[ctg]['ORF_MEDIAN'] = "%.3f"%(np.median(ctgorfs_np))
    intergenic_np = np.array(intergenic)
    contigs[ctg]['INTERGENIC_MEAN'] = "%.2f"%(np.median(intergenic))
    contigs[ctg]['INTERGENIC_MEDIAN'] = "%.2f"%(np.median(intergenic))
    if contigs[ctg]['ORF_COUNT'] > 0:
        contigs[ctg]['FRACTION_OVERLAPPING_ORFS'] = "%.2f"%(contigs[ctg]['FRACTION_OVERLAPPING_ORFS'] / contigs[ctg]['ORF_COUNT'])

sys.stderr.write("There are %d contigs total length = %d\n"%(len(contigs),total_length))
sys.stderr.write("There are %d ORFs\n"%(ORF_count))

with open(args.outbase + ".sum_stat.tsv","w",newline='') as sumstat:
    fieldnames = ['ID', 'LENGTH','GC', 'COVERAGE','ORF_COUNT','ORF_PER_KB','ORF_MEAN','ORF_MEDIAN','ORF_PLUS_STRAND_RATIO',
    'INTERGENIC_MEAN','INTERGENIC_MEDIAN','FRACTION_OVERLAPPING_ORFS','NCLDV_HITS','NCVOG_HITS','CAPSID','RIBOSOMAL','MT',
    'NUM_ORF_NO_PFAM_HITS']
    outtsv = csv.DictWriter(sumstat,delimiter="\t",quoting=csv.QUOTE_MINIMAL,fieldnames=fieldnames)
    outtsv.writeheader()
    for ctgname in contigs:
        outtsv.writerow(contigs[ctgname])
