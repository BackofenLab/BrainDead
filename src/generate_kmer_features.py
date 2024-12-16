#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO
import subprocess
from subprocess import PIPE
import re
import argparse
import os.path
import sys
BINDIR = os.path.dirname(os.path.realpath(__file__)) 
def find_kmer_hits(sequence, kmer):
    return [m.start() for m in re.finditer('(?='+kmer+')', sequence)] # re with look-ahead for overlaps

def call_command (cmd):
    p = subprocess.Popen(cmd,shell=True,stdin=None, stdout=PIPE)
    (result, error) = p.communicate()
    if error:
        raise RuntimeError("Error in calling cmd or perl script\ncmd:{}\nstdout:{}\nstderr:{}".format(cmd, result, error))
    cs_array = result.strip().decode()
    return cs_array

def get_subopt_intarna_strs(sequence, minE_subopt, minE_intarna):
    cmd_subopt = "echo \"{}\" | RNAsubopt -e 100 | perl {}/RNAsubopt-minEperPos.pl".format(sequence,BINDIR)
    csstr_subopt = call_command(cmd_subopt)
    assert (1 == len([l for l in csstr_subopt.split('\n')])) # 1 line expected
    arr_subopt = csstr_subopt.split(';')[1:] # First entry is always empty!
    assert len(sequence) == len(arr_subopt)

    cmd_intarna = "IntaRNA -q \"{}\" -t \"{}\" -n 1000 --noSeed -m M --outMaxE={} --outMode=C --outOverlap=B --outCsvCols=id1,id2,E,hybridDBfull |  grep -v \"^#\" |perl {}/IntaRNA-minEperPos.pl".format(sequence, sequence,str(minE_intarna+0.1), BINDIR)
    csstr_intarna = call_command(cmd_intarna)
    
    if len([l for l in csstr_intarna.split('\n')]) == 3:
        arr_intarna = csstr_intarna.split('\n')[1].split(';')[1:] # Second line is the target line and desired, First entry is "target" and discarded
    elif len([l for l in csstr_intarna.split('\n')]) == 1:
        print("csstr_intarna:", csstr_intarna)
        arr_intarna = [0 for i in range(len(sequence)) ]
        print("arr_intarna:", arr_intarna)
    else:
        print(cmd_intarna)
        print("csstr_intarna:", csstr_intarna)
        raise RuntimeError ("Unexpected IntaRNA output number of lines")

    assert len(sequence) == len(arr_intarna)


    str_subopt_filtE = ''.join([ sequence[i] if float(e) >= minE_subopt else 'N' for  i, e in enumerate(arr_subopt)])


    str_intarna_filtE = ''.join([ sequence[i] if float(e) >= minE_intarna else 'N' for  i, e in enumerate(arr_intarna)])
    
    return str_subopt_filtE, str_intarna_filtE

def find_hits_all(sequence, seq_subopt, seq_intarna, kmer):
    hits_seq = find_kmer_hits(sequence, kmer)
    hits_subopt = find_kmer_hits(seq_subopt, kmer)
    hits_intarna = find_kmer_hits(seq_intarna, kmer)
    hits_intarna_subopt = sorted(list(set(hits_subopt).intersection(hits_intarna)))
    return hits_seq, hits_subopt, hits_intarna, hits_intarna_subopt

def is_valid_file(file_name):
    if os.path.isfile(file_name):
        return os.path.abspath(file_name)
    else:
        raise FileNotFoundError(os.path.abspath(file_name))
        
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Generate kmer-based count features based on sequence plus RNAsubopt and IntaRNA position-wise energies..'
        '\nSample calls: \"python generate_kmer_features.py --kmers "AGA,GC,GGG" --fasta test.fa --out-csv counts.csv\"'\
        '\"python generate_kmer_features.py --kmers "AGA,GC,GGG" --fasta test.fa --out-csv "stdout" --report-counts --minE-subopt -5 --minE-intarna -2\"')

    parser.add_argument('--kmers', required=True, type=str, help='List of kmers as a comma separated string e.g. \"AGG,GA,GG\"')
    parser.add_argument('--fasta', required=True, type=is_valid_file, help='Sequences to extract features from as a FASTA file')
    parser.add_argument('--report-counts', action='store_true', help='Whether to report counts as integer, default is binary nohit(0)-hit(1)'),
    parser.add_argument('--out-csv', type=str, default='stdout', help='CSV File name to write counts, pass "stdout" for stdout ')
    parser.add_argument('--minE-subopt', default=-3, type=int, help='Minimum free energy of the position on RNAsubopt result')
    parser.add_argument('--minE-intarna', default=-3, type=int, help='Minimum free energy of the position on IntaRNA result')

    args = parser.parse_args()
    print(args)
    out_csv_str = "id"
    kmers_list = args.kmers.split(',')
    for kmer in kmers_list:
        out_csv_str += ",{}_any,{}_intra,{}_dimer,{}_free".format(kmer,kmer,kmer,kmer)
    out_csv_str += '\n'
    for r in SeqIO.parse(args.fasta, format='fasta'):
        print(r.id)
        out_csv_str += r.id
        seq_subopt, seq_intarna = get_subopt_intarna_strs(str(r.seq), minE_subopt=args.minE_subopt, minE_intarna=args.minE_intarna)
        for kmer in kmers_list:
            hseq, hsubopt, hintarna, hsubopt_intarna = find_hits_all(str(r.seq),seq_subopt,seq_intarna, kmer)
            print(kmer, hseq, hsubopt, hintarna, hsubopt_intarna)
            cseq, csubopt, cintarna, csubopt_intarna  = len(hseq), len(hsubopt), len(hintarna), len(hsubopt_intarna)
            if args.report_counts is True:
                out_csv_str += ',{},{},{},{}'.format(cseq, csubopt, cintarna, csubopt_intarna)
            else:
                binary_hits = ['0' if c==0 else '1' for c in [cseq, csubopt, cintarna, csubopt_intarna]] 
                out_csv_str += ","+','.join(binary_hits)
        out_csv_str += '\n'

    if args.out_csv == "stdout":
        print(out_csv_str)
    else:
        with open(args.out_csv, 'w') as outfile:
            outfile.write(out_csv_str)
