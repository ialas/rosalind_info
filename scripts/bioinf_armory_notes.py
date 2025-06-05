# -*- coding: utf-8 -*-
"""
Created on Fri May 16 13:47:35 2025

@author: imraa
"""

#%% Conda Environment
# I used python-class

#%% Libraries
import numpy as np
import pandas as pd
import math
import os

from Bio import Entrez
from Bio import SeqIO

#%% Functions Created:

def textParser(filename):
    # filename is a string
    f = open(filename,'r');
    textStorage = [];
    for line in f:
        textStorage.append(line.strip()) # remove newline immediately
    return textStorage

def dna_codon_table(filename, invert=False):
    # Automap codon table
    f = open(filename, 'r')
    lineStorage = [];
    for line in f:
        lineStorage.append(line);
    # remove \n
    lineStorage = [i.split('\n') for i in lineStorage]
    lineStorage = [i[0] for i in lineStorage]
    # Separate
    lineStorage = [i.split(' ') for i in lineStorage]
    lineStorage2 = [];
    for j in lineStorage:
        lineStorage2.append(list(filter(None, j)))
    # Create codon dictionary
    dna_codon_dict = {};
    for i in lineStorage2:
        count = 0;
        while count <= 7:
            dna_codon_dict[i[count]] = i[count+1];
            count += 2
    # dna_codon_dict['']
    if invert == False:
        return dna_codon_dict
    elif invert == True:
        flip_dna_codon = {};
        for i in lineStorage2:
            count = 0;
            while count <= 7:
                if i[count+1] not in flip_rna_codon:
                    flip_dna_codon[i[count+1]] = [i[count]]
                else:
                    flip_dna_codon[i[count+1]].append(i[count])
                count += 2;
        return flip_dna_codon

def dna_reverse_complement(dna_string):
    complement = {};
    complement['A'] = 'T';
    complement['G'] = 'C';
    complement['T'] = 'A';
    complement['C'] = 'G';
    complement_seq = ''.join([complement[i] for i in dna_string][::-1])
    return complement_seq

def fasta_seq_cleaner(filename):
    # Grab information, remove newlines.
    textStorage = textParser(filename)
    # Find the headers.
    count=0;
    header_count = [];
    for i in textStorage:
        if i[0] == '>':
            header_count.append(count)
        count += 1;
    # Combine everything that is between headers.
    count = 0;
    textStorage2 = [];
    for i in header_count:
        # Everything except the last sequence (avoid index issues)
        if i != header_count[-1]:
            beginning = i+1;
            end = header_count[count+1];
            textStorage2.append(textStorage[i]) # Make a rosalind one.
            textStorage2.append([''.join(textStorage[beginning:end])][0]);
        if i == header_count[-1]: # for the last one (parser was failing to get the last strain)
            beginning=i+1
            end = len(textStorage)
            textStorage2.append(textStorage[i])
            textStorage2.append([''.join(textStorage[beginning:end])][0]);    
        count += 1
    # Generate dictionary of dna.
    dna_dict = {};
    count=0;
    while count < len(textStorage2):
        dna_dict[textStorage2[count]] = textStorage2[count+1];
        count += 2;
    # return the dictionary of DNA.
    return dna_dict

#%% 001: Intro to Bioinformatics Armory

# Given: DNA sequence
# Find: Nucleotide counts

# Program: Sequence Manipulation Suite
    # Does:
        # Format conversion
            # Fasta, EMBL, GenBank, etc
        # Sequence analysis
        # Sequence figures
        # Random sequences
# Used program to solve question. 
# Below is a code-based implementation
        
# Rapid Introduction to Molecular Biology
f = open('rosalind_info/data/rosalind_ini.txt', 'r')
count=0
for line in f:
    string1 = line
# string1 = 'AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC'
countA = 0;
countT = 0;
countC = 0;
countG = 0;
for i in string1:
    if i == 'A':
        countA +=1
    elif i == 'T':
        countT +=1
    elif i == 'C':
        countC +=1
    elif i == 'G':
        countG +=1
sol001 = [countA,countC,countG,countT]

with open('rosalind_info/outputs/rosalind_ini_answer_05.txt','w') as o:
    for line in sol001:
        o.write(f"{line}\n")
        
#%% 002: GenBank Introduction

from Bio import Entrez # Entrez is a data retrieval system offered by the NCBI

# Entrez.esearch() # searchs any of the NCBI databases
    # Arguments:
        # db # This can be "nucleotide" for GenBank or "pubmed" for Pubmed
        # term # This is the query field

# Example search:
Entrez.email = 'Email Address'; # you have to put your email in here

handle = Entrez.esearch(db='nucleotide', term=' "Rhytisma"[Organism] AND ("2002/06/09"[PDAT] : "2010/07/17"[PDAT])')
record = Entrez.read(handle)
fsolution = record["Count"]
with open('rosalind_info/outputs/rosalind_gbk_answer_004.txt','w') as o:
    o.write(fsolution)
    
#%% 003: Data Formats

# Popular file formats (from 2013)
# fasta, nexus, phylip
# genbank is its own file format

# Given: <10 GenBank entry IDs
# Return: The shortest of strings associated with the IDs in fasta format

# The Sequence Manipulation Suite (SMS) had a GenBank -> Fasta converter
# Manual Option: For each GenBank entry ID, download the relevant genbank file. Then convert them all to fasta. Then find the shortest.
    # Steps 2 and 3 can be swapped, to avoid the time-consuming step.
    
# Coding Option:
    
f = open('rosalind_info/data/rosalind_frmt(2).txt', 'r')
count=0
for line in f:
    string1 = line
string_list = string1.split(' '); # separate
# remove newline
last_entr = string_list[-1]
last_entr = last_entr.split('\n')[0];
# replace newline entry
string_list[-1] = last_entr
gbk_id_list = string_list;

from Bio import Entrez
from Bio import SeqIO

Entrez.email = 'example email'
handle = Entrez.efetch(db='nucleotide', id = gbk_id_list, rettype = 'fasta');
records = list(SeqIO.parse(handle, 'fasta')) # get a list of SeqIO objets in Fasta Format
seq_list = [len(i.seq) for i in records]
min_seq = min(seq_list)
min_seq_index = seq_list.index(min_seq)
min_seq_name = records[min_seq_index].description
min_seq_full = list(records[min_seq_index].seq)

# issue: fasta seq are formatted with max 70 characters per line
# soln: separate fasta seq by every 70 characters
split_min_seq_full = [min_seq_full[i:i+70] for i in range(0,len(min_seq_full),70)]

with open('rosalind_info/outputs/rosalind_frmt_answer_003.txt','w') as o:
    o.write(">" + min_seq_name + '\n')
    for line in split_min_seq_full:
        line = ''.join(line) + '\n' # add newline to separate
        o.write(line)
        
#%% 004: New Motif Discovery

# ProSite database: helds find motifs in proteins that have already been discovered in other sequences
# MEME:
    # input: collection of DNA/protein sequences
    # output: motifs exceeding a user-specified statistical confidence threshold
    
    # does not take into account gaps. motifs that use indels are not observed. GLAM2 can do that.
    
# Run using MEME website.
# Command line argument: meme rosalind_meme.txt -protein -oc . -nostatus -time 14400 -mod zoops -nmotifs 3 -minw 6 -maxw 50 -objfun classic -markov_order 0

# Issue I noted:
    # The text version of the output from MEME. The sequence it found, it doesn't include the ambiguous areas. So if there was a part of the seequence that was different, lets say W on one and K on the others, it would choose 1. I had to manually put that into Rosalind
    
#%% 005: Pairwise Global Alignment

# EMBOSS's Needle tool aligns DNA to DNA or RNA to RNA (protein) sequences.
    # Normally, Scoring is defined as:
        # +1 for matching symbols
        # -1 for all mismatches
        # -5 for a gap opening (each new gap)
        # -1 for each gap extension (ex: a gap with 3 extensions is -3)
    # Needle has:
        # -10 for gap opening
        # -1 for each gap extension
        
# Installed to conda environment (seen above)

# Coding implementation:
    # Given: 2 genbank IDs
    # Output: maximum global alignment score between the two DNA strings

f = open('rosalind_info/data/rosalind_need.txt', 'r')
count=0
for line in f:
    string1 = line
string_list = string1.split(' '); # separate
# remove newline
last_entr = string_list[-1]
last_entr = last_entr.split('\n')[0];
# replace newline entry
string_list[-1] = last_entr
gbk_id_list = string_list;

Entrez.email = 'imraanalas@gmail.com'
handle = Entrez.efetch(db='nucleotide', id = gbk_id_list, rettype="fasta");
records = list(SeqIO.parse(handle, "fasta")) # get a list of SeqIO objets in Fasta Format
for i in records:
    print(i.seq)
    
# Use Emboss's Needle. 
# Settings:
    # output: pair
    # matrix: DNAfull
    # gap open, end gap open = 10, 10
    # gap extend, end gap extend = 1.0, 1.0
    # end gap = TRUE
    
#%% 006: FastQ format introduction

# Sequencing quality is defined using the Phred quality score.
    # Q = -10*log_10_(P)
        # P is the probability that the corresponding base call is incorrect.
        # It's per base. You can graph the phred score per base as a scatter plot.
    # Q = 10 -> 0.1 probability (90% accuracy)
    # Q = 20 -> 0.01 probability (99% accuracy)
    # Q = 30 -> 0.001 probability (99.9% accuracy)

# In a fastQ, the quality score (Q) is scored from 0-40.
# Each quality score is represented by the ASCII characters 33-73.
    # For instance, ! is ascii code 33, and represents a Q of 0.
    # + is ascii code 43, and represents a Q of 10
    
# Converting FastQ to FASTA
# Given: FastQ file
# Output: FASTA file

# Declare file locations (existing and to exist)
input_file = 'rosalind_info/data/rosalind_tfsq.txt'
output_file = 'rosalind_info/outputs/rosalind_tsfq_fasta.fasta'
if not os.path.exists(output_file): # if the file doesn't exist
    SeqIO.convert(input_file,'fastq', output_file,'fasta') # convert file to fasta and put in output location
    
#%% 007: Read Quality Distribution

# Phred quality score is classic, but we would want to look at it in aggregate instead of per base.
# Per-read quality score distribution:
    # The distribution of the average quality of each read.
    
# FastQC is what is typically used for quick quality assessment.

# Given: Quality Threshold & FastQ entries for multiple reads (of the same strain)
# Return: # of reads with an average quality BELOW the threshold.
    # They want # of reads BELOW the threshold

# Implementation

input_file = 'rosalind_info/data/rosalind_phre.txt'
f = open(input_file, 'r')
count=0
string_all = [i.strip('\n') for i in f];

# quality score
qual_score = int(string_all.pop(0)) # the first entry is the quality, assume integer
f.close() # had to close to make changes to the file to account for below

# Issue: I don't think SeqIO likes when there's a random integer as the first line in the file.
# I could code it such that I take the integer, store it as qual_score, then rewrite the file back without it.
# But I don't plan on it.

failed_qc = 0;
for record in SeqIO.parse(input_file,"fastq"): # for each fastq record in the fastq file
    phred = record.letter_annotations['phred_quality']; # get the list of phred scores, per base in the sequence
    # print(phred)
    ave_phred = sum(phred)/len(phred); # get the average phred score for one sequence
    # print(ave_phred)
    if ave_phred < qual_score: # BELOW the threshold
        failed_qc += 1; # only on low threshold
        
#%% 008: Protein Translation

# SMS suite has the Translate Tool to translate nucleotides into amino acids
# BioPython possesses the translate() method for converting DNA or RNA strings to protein strings

# Given: DNA string. Protein string.
# Output: Determine which genetic code variant was used to translate the DNA into the protein string

# Import
from Bio.Seq import translate

input_file = 'rosalind_info/data/rosalind_ptra(1).txt'
f = open(input_file, 'r')
count=0
string_all = [i.strip('\n') for i in f];

# dna
coding_dna = string_all[0]
# coding_dna = 'ATGGCCATGGCGCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA'
# protein 
protein_str = string_all[1]
# protein_str = 'MAMAPRTEINSTRING'
# genetic code variant list
gcv_list = [1,2,3,4,5,6,9,10,11,12,13,14,15]
prot_1 = [];
prot_1_yes = [];
prot_1_conf = [];
prot_2 = [];
prot_2_yes = [];
prot_2_conf = [];

for gcv in gcv_list:
    prot_trans = translate(coding_dna, table=gcv,to_stop=True) # translate DNA using 1 of 15 gcv tables. continue until stop codon found
    prot_trans_no_stop = translate(coding_dna, table=gcv,to_stop=False) # just try it
    prot_1.append(prot_trans);
    prot_2.append(prot_trans_no_stop);
    if prot_trans == protein_str:
        prot_1_yes.append(gcv)
    if prot_trans_no_stop == protein_str:
        prot_2_yes.append(gcv)
    # if len(prot_trans) == len(protein_str):
    #     prot_1_conf.append(gcv)
    # if len(prot_trans_no_stop) == len(protein_str):
    #     prot_2_conf.append(gcv)
        
if len(prot_1_yes) > 0:
    print('With Stop:')
    print(prot_1_yes)
if len(prot_2_yes) > 0:
    print('Without Stop:')
    print(prot_2_yes)
    
# Note: My first attempt failed to find any exact matches. My second attempt did.
# Theory: First attempt did not have any 100% identity matches, but still had 100% coverage.
# The second attempt had 100% identity coverage (every base aligned perfectly).
    
#%% 009: Read Filtration by Quality

# Filter out poor-quality reads using FastQ Quality Filter from fASTX (toolkit).
    # Actually, keep and report the good-quality reads.

# I think I can do it programmatically

# Given: FastQ File that contains 1 additional row on top.
    # First row has two integers. The first is the quality threshold. The second is the % of bases.
    # Basically: If X is the % of bases value. Then X% of bases must be above the quality threshold for their phred score.
    # Try a list comprehension to get True/Fase for greater than or less than teh quality score, then sum, then get a % relative to total length, then compare to % value X.
    
# Implementation

# Manually changed file extension to fastq
input_file = 'rosalind_info/data/rosalind_filt.fastq' 
f = open(input_file, 'r')
count=0
string_all = [i.strip('\n') for i in f];

# get quality threshold & percentage bases
metric_vals = string_all[0];
qual_thresh = metric_vals.split(' ')[0];
qual_thresh = int(qual_thresh)
perc_base = metric_vals.split(' ')[1];
perc_base = int(perc_base)

# example
qual_thresh = 19;
perc_base = 57;

# manually remove the first row of the filt file

# solve 
filter_count = 0;
phred_eval = [];
phred_test = []
phred_fail = []
print(perc_base)
for record in SeqIO.parse(input_file,"fastq"): # for each fastq record in the fastq file
    phred = record.letter_annotations['phred_quality']; # get the list of phred scores, per base in the sequence
    filter_score = sum([(i >= qual_thresh) for i in phred])/len(phred); # get % of above threshold bases
    phred_eval.append(filter_score)
    # phred_test.append(sum([(i > qual_thresh) for i in phred])/len(phred))
    if filter_score >= (perc_base/100):
        # print(sum([(i >= qual_thresh) for i in phred])/len(phred))
        # phred_fail.append(sum([(i >= qual_thresh) for i in phred])/len(phred))
        
        # [(i > qual_thresh) for i in phred] returns a list of T/F, where T is if the phred count for that base was above the quality threshold
        # sum[above] counts up all the phred scores that were above the threshold
        # sum[above]/len(phred) tells us what % of phred scores were above the threshold
        # we then compare that to the perc_base value (divided by 100)
        # if the % of phred scores above the threshold is greater than the perc base, we would keep that 
        # so we increment the filter_count
        filter_count += 1;

succ_phred = [b for b,i in enumerate(phred_eval) if i >= perc_base/100]
succ_phred_2 = [b for b,i in enumerate(phred_eval) if i > perc_base/100]

print('FASTQ Filtered Based On: Quality Threshold of ' + str(qual_thresh) + ' and Percentage of Bases (' + str(perc_base) + '%) is: ' + str(filter_count))

# my first attempt assumed that they were looking for the reads that would be filtered as bad. i will now try the opposite
    # Yup. needed to get the opposite value. Ambiguous question.

#%% 010: Complimenting a Strand of DNA

# Imports

# Implementation
input_file = 'rosalind_info/data/rosalind_rvco.txt'
textStorage = textParser(input_file)

# get the names
names = [i for i in textStorage if i[0] == '>']
# get the sequences
sequences = [i for i in textStorage if i[0] != '>']
# each seuqnece is on 2 rows, combine
sequence_new = [];
counter = 0;
while counter < len(sequences):
    sequence_new.append([sequences[counter] + sequences[counter+1]])
    counter += 2;

# dictionary comprehension to create a dictionary of names:sequences
text_dict = {}
# text_dict = {k:v for (k,v) in zip(names, sequences)}
text_dict = {k:v for (k,v) in zip(names, sequence_new)}

# Generate a dictionary of complements (back when I thought I had to compare parts of strings against complements of all the strings)
# comp_dict = {}
# comp_seq = [];
# for i in text_dict:
#     reverse_seq = dna_reverse_complement(text_dict[i]);
#     comp_seq.append(reverse_seq)
# comp_dict = {k:v for (k,v) in zip(names, comp_seq)}

# oh. all I have to do is check to see if the complement of the string is identical to the string.

comp_dict = dict.fromkeys(list(text_dict.keys()))

counter=0;
for i in text_dict:
    comp_dict[i] = [dna_reverse_complement(text_dict[i][0])]
    if text_dict[i][0] == dna_reverse_complement(text_dict[i][0]): # [0] is because dna_reverse_complement takes a string and also returns a string
        counter += 1;
        
print(counter)
    

for i in text_dict:
    print(text_dict[i][0])
    print(comp_dict[i][0])
    
#%% 011: Subotimal Local Alignment

# Transposons are relatively short intervals of DNA that can, during the process of recombination, allow for duplication/deletion of regions 
    
# Typically, a transposon don't necessarily align spatially in a consistent manner, so we must search across the the DNA string to find matches that don't classically align well.
# To do this, we use LAlign, to find suboptimal alignments that would allow for transposons to insert/delete during recombination.

# Given:
    # Two DNA strings in fasta format that share an inexact repeat of about 32-40 bp (no more than 3 changes/indels).
# Return:
    # Total number of occurrences of the repeat, in the first string, then in the second string.
    
# So we need to find the substring of about 32-40 bp, then search for inexact matches of it in the two DNA strings.

# Don't implement this here

# Issues:
    # The tool they linked from EBI, i think the link they had was out-of-date. There was no guidance on parameters to set.
    # The FASTA package manual was unable to be downloaded, as the bioch.virginia.edu link doesn't seem to work anymore.
# Solution:
    # Khang Tsung Fei from the Questions section of the board recommended the usage of DotLet, which visualizes the alignment of two sequences against each other.
    # By looking for the short (~32-40 bp) regions that shared alignment, I was able to find that there was only 7 and 5 of them.
    
#%% 012 Base Quality Distribution

# Quality of reads often degrade over the course of a sequence.
# Per Base Sequence Quality module of FastQC program genereates a quality distribution.

# Given:
    # Fastq file, quality threshold q
# Return:
    # Number of positions where the mean base quality falls below the threshold q (BELOW)
    
# Clarification: Is the threshold q different from the "accepted base quality threshold"?
    # As in, is this simply just find the # of positions (bases) that hav ea phred_quality below q?
    # or do I have to determine what positions have a phred quality score that is below the Q3 (25th percentile)?
# Secondary Clarification:
    # I think, that they:
        # 1) Find the phred_quality per base
        # 2) Average the phred_quality per base per sequence?
            # So, if the phred_quality for base 1 is [10] for #1, and then [20] for #2, and so on. Average those, determine if the average falls below the quality threshold.
            

# Implementation

# Imports
import scipy
from scipy import stats
import numpy as np

# 
input_file = 'rosalind_info/data/rosalind_bphr.txt'
f = open(input_file, 'r')
count=0
string_all = [i.strip('\n') for i in f];

# quality score
qual_score = int(string_all.pop(0)) # the first entry is the quality, assume integer
f.close() # had to close to make changes to the file to account for below

# Issue: I don't think SeqIO likes when there's a random integer as the first line in the file.
# I could code it such that I take the integer, store it as qual_score, then rewrite the file back without it.
# But I don't plan on it.

# Manually take out the QC score from the original file before proceeding.

failed_qc = qual_score;
# failed_qc = 23;
phred_scores = [];
num_records = 0;
# Get all the phred_scores
for record in SeqIO.parse(input_file,"fastq"): # for each fastq record in the fastq file
    phred = record.letter_annotations['phred_quality']; # get the list of phred scores, per base in the sequence
    # print(phred)
    phred_scores.append(phred);
    num_records += 1; # Store the amount of sequences you have.

phred_scores = np.array(phred_scores) # Convert the list of phred_scores into an array.
phred_ave = (phred_scores.sum(axis=0)/num_records) # For each "base", as in each part of the sequence, average the phred scores across the sequences. So all of the 1st bp phred scores will be summed and then averaged to get the mean phred_score for the 1st bp.

# Evaluate, for each position (bp), which positions (bp) failed the QC threhsold.
failed_position_qc = 0;
for ave_phred_score in phred_ave:
    if ave_phred_score < failed_qc:
        failed_position_qc += 1

print('Number of positions where mean base quality falls below given threshold: ' + str(failed_position_qc))

#%% 013: Base Filtration Quality

# Instead of deleting poor quality bp in the reads, trim them.

# Use the FastQ Quality Trimmer tool on Galaxy: https://usegalaxy.org/root?tool_id=fastq_quality_trimmer
    # Set window_size to 1
