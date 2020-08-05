#!/bin/python

import sys
import random

# input fasta file
inFasta = sys.argv[1]


# parse dna sequence
def fasta_parser(fasta):
    dna = ""
    with open(fasta) as f:
        for line in f:
            dna += line.strip("\n").lower()
    dna = dna[6:]
    return dna


# generate random coordinates
def random_seq_generator():
    dna_seqs = []
    begin = 9411250
    end = 48020000
    while begin < end:
        rand_begin = random.randint(begin, end)
        rand_end = rand_begin + random.randint(30, 2000)
        dna_seqs.append([rand_begin, rand_end])
        begin = rand_end
    return dna_seqs


# create bed file format
def seq_to_bed(seqs):
    bed_file = ""
    for each_seq in seqs:
        bed_file += "chr21\t" + str(each_seq[0]) + "\t" + str(each_seq[1]) + "\n"
    bed_file = bed_file[:-1]
    return bed_file


dna = fasta_parser(inFasta)
bed = seq_to_bed(random_seq_generator())


# write to bedfile
bed_file = open("bedfile.bed", "w+")
bed_file.write(bed)
bed_file.close()
