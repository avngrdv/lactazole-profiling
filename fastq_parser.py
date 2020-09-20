from __future__ import division

# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 01:37:10 2020
@author: Vinogradov Alex
"""

import numpy as np
import os, tables

def fetch_DNA(filename):
    '''
    Fetch DNA as a list of sequences from a .fastq file.
    .fastq files are base call .fastqs 
    from single-read run Illumina data.
    
    Parameters
    ----------
    filename : full path to a .fastq file

    Returns
    -------
    DNA : a list of sequences each 
          containing a single read

    '''
    
    with open(filename, 'r') as f:
        content = f.readlines()
        
        DNA = content[1::4]
        DNA = [x.rstrip('\n') for x in DNA]
                
        f.close()
        
    return DNA

def translate(seq):
    '''
    In silico translation function which converts a DNA sequence 
    into a peptide according to the genetic code as specified in tables.py. 
    
    This is NOT a multipurpose translation function that finds
    a start codon in every frame, etc. According to our experimental design,
    sequencing reads start in the middle of the underlying ORF,
    in-frame with it, so a straightforward reading of nucleotide
    triples suffices here.
    
    
    Parameters
    ----------
    seq : DNA sequence (dtype=str)

    Returns
    -------
    translated peptide (dtype=str)

    '''

    #ORFs containing ambigious symbols or no stop codon
    #are labelled with a '+' sign.     
    def find_stop(peptide):
        ind = peptide.find('_')
        if ind == -1:
            return peptide + '+'
        else:
            return peptide[:ind]
        
    pep = ''
    for i in range(0, len(seq), 3):
        try:
            pep += tables.translation_table[seq[i:i+3]]
        except:
            pep += '+'
    
    pep = find_stop(pep)
    return pep

def transform_peps(peps):
    '''
    P matrix holds a peptide dataset as a 2D numpy array of shape
    (number of peptides x length of the longest sequence in the dataset).
    Shorter sequences are right-padded to the longest peptide in the
    dataset.
    
    Transformation to a numpy array is done to speed-up and 
    simplify downstream analysis. 

    Parameters
    ----------
    peps : a list of peptide sequences represented as strings

    Returns
    -------
    P : P matrix, a 2D-numpy array

    '''

    #find the longest peptide in the dataset
    L = max([len(x) for x in peps])

    #create an empty array of approapriate dimensions     
    P = np.zeros((len(peps), L), dtype='<U1')
    
    #iteratively fill it
    for i,pep in enumerate(peps):
        P[i,:len(pep)] = list(pep)
        
    return P

def hamming_distance(P, pep):
    '''
    A custom Hamming distance (HD) calculator. Calculates the HD
    between peptide pep and each sequence in the P-array.
    
    ----------
    P : P-matrix as a 2D-numpy array
    
    pep : Peptide sequence to compare P against. dtype=str

    Returns
    -------
    
    1D-numpy array of Hamming distances between elements of P and pep
    '''
    return np.sum(P != pep, axis=1)


def parse_sample(fname, design):
    '''
    Parser for a single .fastq file.
    Fetches a .fastq file, reads DNA from it, and fills a P matrix. 
    Peptides in P are filtered to get rid of junk sequences, etc.
    P as a 2D-numpy array is written to a file at the end of the operation.

    Parameters
    ----------
    fname : path to a .fastq file to be analyzed

    design :  an instance of PepDesign, specifying the sequence
              of the analyzed library
        
    Returns
    -------
    None.

    '''
    
    print('Parsing', fname, '. . .')
    
    #Load DNA, translate it and represent the peptide as a numpy array (P)
    DNA = fetch_DNA(fname)           
    peps = [translate(x) for x in DNA]
    P = transform_peps(peps)
    
    #Fetch peptides which are as long as specified by the library design
    P = P[np.sum(P != '', axis=1) == design.L]
    
    #Make sure peptides do not contain ambiguous symbols, etc
    P = P[np.all(P != '+', axis=1)]

    #Fetch variable regions
    P = P[:,np.nonzero(design.vr_mask)[0]]

    #Fetch sequences containing less than two mutations from the parent
    #i.e. fetch the parent and its single-point mutants
    HD = hamming_distance(P, design.wt)
    P = P[HD < 2]
    
    #Write to file as an .npy
    #Infer sample name from the .fastq file name
    sample_name = fname.split('.')[0]
    fname = sample_name + '_as_P'
    np.save(fname, P.astype(np.string_))
    return


def parse_cwd(design):
    '''
    Parse all .fastq files in the current working directory

    Parameters
    ----------
    design :  an instance of PepDesign, specifying the sequence
              of the analyzed library
        
    
    Returns
    -------
    None.

    '''
    
    cwd = os.getcwd()
    fnames = [os.path.join(v) for v in os.listdir(cwd) if v.endswith(".fastq")]
    
    for fname in fnames:
        parse_sample(fname, design)

    return