
#!/usr/bin/env python3

import sys, os, csv, argparse, re
from Bio import SeqIO
from Bio.SeqUtils import GC
import numpy as np
import numpy.ma as ma
from scipy.stats.mstats import mode, gmean, hmean

def parse_hmmersearch_domtbl(filehandle,cutoff):
    "Parse HMMER hmmsearch domtbl for VOG/DLV db search"
    table = {}
    if not filehandle:
        print("no open filehandle")
        return table

    for line in filehandle:
        if line[0].startswith("#"):
            continue
        row = line.strip("\n").split()
        hit = row[0]
        evalue = float(row[6])
        q = row[3]
        if evalue > cutoff:
            continue
        score = 0.0
        if row[7] is not '.':
            score = float(row[7])

        if hit not in table:
            table[hit] = {}

        if (q not in table[hit] or
            table[hit][q]['evalue'] > evalue):
            aln_len = abs(int(row[20]) - int(row[19])) + 1
            table[hit][q] = {'evalue': evalue,
                            'score': score,
                            'qalnlen': aln_len,
                            }
    return table

def parse_blast (filehandle,cutoff):
    "Parse BLAST Or FASTA 12 column output"
    result_table = {}
    if not filehandle:
        print("no valid input filehandle for parse_blast")
        return result_table
    search = csv.reader(filehandle,delimiter="\t")

    for row in search:
        evalue = float(row[10])
        if evalue > cutoff:
            continue
        q = row[0]
        h = row[1]
        pid = row[2]
        if q not in result_table:
            result_table[q] = row
    return result_table


def gff_parse(filehandle):
    ctg_orf_counter = {}
    ORF_set = {}
#    ORF_cols = ['length','strand','score','start','end']
    if not filehandle:
        print("No GFF filehandle provided")
        return ORF_set
    gff_parse = csv.reader(filehandle,delimiter="\t")
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

#        ORF_count += 1
        score = 0.0
        if ORF[5] is not '.':
            score = float(ORF[5])

        name = ORF[8]

        if ORF[1].startswith("Prodigal"):
            if chrom not in ctg_orf_counter:
                ctg_orf_counter[chrom] = 0
            ctg_orf_counter[chrom] += 1
            name = "%s_%d"%(chrom,ctg_orf_counter[chrom])
        else:
            print("Expected Prodigal format GFF using %s"%(name))

#        if name in ORF_names:
#            print("warning already saw orf %s in set"%(name))

#        ORF_names[name] = {'chrom': chrom,
#                           'name': name,
#                           'length': abs(stop - start) + 1,
#                           'strand': strand, 'score': float(ORF[5]),
#                           'start': start, 'stop': stop
#                           }
        ORF_set[chrom].append({'name': name,
                               'length': abs(stop - start) + 1,
                               'strand': strand, 'score': float(ORF[5]),
                               'start': start, 'stop': stop
                              })
    return ORF_set

ribosomal_evalue_cutoff = 1e-10
mitochondria_evalue_cutoff = 1e-5
viral_domain_cutoff = 1e-5
pfam_cutoff = 0.01
parser = argparse.ArgumentParser(description='summary stats to train on viral genomes')

parser.add_argument('-i','--infile', type=argparse.FileType('r'),nargs='?',default=sys.stdin,
                    help='Input file for reading genome DNA (FASTA)')
parser.add_argument('-o','--outbase',required=True,help="Output file base")
parser.add_argument('-p','--prodigal',type=argparse.FileType('r'),required=True,help="Prodigal GFF")
parser.add_argument('-sd','--search_ncldv',type=argparse.FileType('r'),required=False,help="NCLDV HMMER domtbl (from prodigal peps)")
parser.add_argument('-sv','--search_ncvog',type=argparse.FileType('r'),required=False,help="NCVOG HMMER domtbl (from prodigal peps)")
parser.add_argument('--pfam',type=argparse.FileType('r'),required=False,help="Pfam HMMER domtbl (from prodigal peps)")
parser.add_argument('-rrna','--ribosomal',type=argparse.FileType('r'),required=False,help="rRNA BLASTN TAB")
parser.add_argument('-mt','--mito',type=argparse.FileType('r'),required=False,help="mitochondria TFASTX TAB")


args = parser.parse_args()

MT = parse_blast(args.mito,mitochondria_evalue_cutoff)
rRNA = parse_blast(args.ribosomal,ribosomal_evalue_cutoff)

ncldvhits = parse_hmmersearch_domtbl(args.search_ncldv,viral_domain_cutoff)
ncvoghits = parse_hmmersearch_domtbl(args.search_ncvog,viral_domain_cutoff)

print("%d ncldv hits hits parsed"%(len(ncldvhits.keys())))
print("%d ncvog hits parsed"%(len(ncvoghits.keys())))
pfamhits = parse_hmmersearch_domtbl(args.pfam,pfam_cutoff)

print("%d pfam hits parsed"%(len(pfamhits.keys())))

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
            'ORF_COUNT':  0,
            'ORF_PLUS_STRAND_RATIO': 1,
            'INTERGENIC_MEAN': 0,
            'INTERGENIC_MEDIAN': 0,
            'FRACTION_OVERLAPPING_ORFS': 0,
            'NCVOG_HITS': 0,
            'NCLDV_HITS': 0,
            'NCVOG_PERCENT': 0,
            'NCLDV_PERCENT': 0,
            'CAPSID': 0,
            'RIBOSOMAL': 0,
            'MT': 0,
            'PFAM_HITS': 0,
            'PFAM_PERCENT': 0
            }

        if record.id in rRNA:
            contigs[record.id]['RIBOSOMAL'] = 1

        if record.id in MT:
            contigs[record.id]['MT'] = 1

        total_length += len(record)


ORF_set   = gff_parse(args.prodigal)

for ctg in ORF_set:
    if ctg not in contigs:
        sys.stderr.write("cannot find chrom %s in FASTA file, but is in the Prodigal file, make sure these match"%(ctg))
        exit()
    ctgorfs = []
    intergenic = []
    minusgenes = plusgenes = 0
    last_stop = 1
    last_strand = 1
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
            if ORF['strand'] == last_strand:
                contigs[ctg]['FRACTION_OVERLAPPING_ORFS'] += 1
        else:
            intergenic.append( ORF['start'] - last_stop + 1)
        last_stop = ORF['stop']
        last_strand = ORF['strand']

        if ORF['name'] in ncvoghits:
            contigs[ctg]['NCVOG_HITS'] += 1
        if ORF['name'] in ncldvhits:
            contigs[ctg]['NCLDV_HITS'] += 1
            #print(ncldvhits[ORF['name']].keys())
            if 'capsid' in ncldvhits[ORF['name']]:
                contigs[ctg]['CAPSID'] += 1
        if ORF['name'] in pfamhits:
            contigs[ctg]['PFAM_HITS'] += 1

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
    contigs[ctg]['INTERGENIC_MEAN'] = "%.2f"%(np.mean(intergenic))
    contigs[ctg]['INTERGENIC_MEDIAN'] = "%.2f"%(np.median(intergenic))
    if contigs[ctg]['ORF_COUNT'] > 0:
        contigs[ctg]['FRACTION_OVERLAPPING_ORFS'] = "%.2f"%(contigs[ctg]['FRACTION_OVERLAPPING_ORFS'] / contigs[ctg]['ORF_COUNT'])

    if contigs[ctg]['ORF_COUNT'] > 0:
        # how many ORFs (fraction) are VOGs or LDV  hits
        contigs[ctg]['NCVOG_PERCENT'] = "%.4f"%(100.0 * contigs[ctg]['NCVOG_HITS']/ contigs[ctg]['ORF_COUNT'])
        contigs[ctg]['NCLDV_PERCENT'] = "%.4f"%(100.0 * contigs[ctg]['NCLDV_HITS']/ contigs[ctg]['ORF_COUNT'])
        # how many ORFs (fraction) have PFAM hits
        contigs[ctg]['PFAM_PERCENT'] = "%.4f"%(100.0 * contigs[ctg]['PFAM_HITS'] / contigs[ctg]['ORF_COUNT'])


sys.stderr.write("There are %d contigs total length = %d\n"%(len(contigs),total_length))
#sys.stderr.write("There are %d ORFs\n"%(ORF_count))

with open(args.outbase + ".sum_stat.tsv","w",newline='') as sumstat:
    fieldnames = ['ID', 'LENGTH','GC', 'COVERAGE','ORF_COUNT','ORF_PER_KB','ORF_MEAN','ORF_MEDIAN','ORF_PLUS_STRAND_RATIO',
    'INTERGENIC_MEAN','INTERGENIC_MEDIAN','FRACTION_OVERLAPPING_ORFS',
    'NCLDV_HITS','NCLDV_PERCENT','NCVOG_HITS','NCVOG_PERCENT','CAPSID',
    'RIBOSOMAL','MT',
    'PFAM_HITS','PFAM_PERCENT']
    outtsv = csv.DictWriter(sumstat,delimiter="\t",quoting=csv.QUOTE_MINIMAL,fieldnames=fieldnames)
    outtsv.writeheader()
    for ctgname in contigs:
        outtsv.writerow(contigs[ctgname])
