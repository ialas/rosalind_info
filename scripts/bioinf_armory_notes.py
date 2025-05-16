# -*- coding: utf-8 -*-
"""
Created on Fri May 16 13:47:35 2025

@author: imraa
"""
#%% Libraries
import numpy as np
import pandas as pd
import math

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
    
