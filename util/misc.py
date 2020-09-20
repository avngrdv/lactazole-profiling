# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 23:56:58 2020
@author: Vinogradov Alex
"""

import numpy as np
import tables

def _compute_F(P):
    '''
    Compute the frequency for each single-point mutant in the dataset P.

    Parameters
    ----------
    P : 2D-numpy array holding peptides as rows

    Returns
    -------
    F : Frequencies for each single-point mutant as 2D-numpy array with
        the shape (number of amino acids in the alphabet, 
                   number of randomized positions in the library).
        
        Note that cells corresponding to the wild-type amino acids
        in the F matrix do not represent the frequency for the wild-type
        peptide. This needs to be calculated separately.
    '''
    
    C = np.zeros((len(tables.aas),P.shape[1]))
    for i,aa in enumerate(tables.aas):
        C[i] = np.sum(P == aa, axis=0)
    
    F = np.divide(C, P.shape[0])
    return F

def _compute_wt_freq(P, wt):
    '''
    Compute the frequency of the parent peptide in the sample.

    Parameters
    ----------
    P : P-array of peptides (2D-numpy array)
    wt : amino acid sequence of the parent peptide (only variable region)
         1D-numpy array

    Returns
    -------
    Frequency of the wild-type peptide in the sample
    
    '''
    
    #c: count of parent peptide reads in the dataset
    c = np.sum(np.all(P == wt, axis=1))
    return np.divide(c, P.shape[0])

def _wt_score(P_pos, P_neg, wt):
    '''
    Calculate the Y score for the parent peptide.
    
    Parameters
    ----------
    P_pos : P-array for positive ("reactive") peptides (2D-numpy array)
    P_neg : P-array for negative ("unreactive") peptides (2D-numpy array)
    
    wt : amino acid sequence of the parent peptide (only variable region)
         1D-numpy array

    Returns
    -------
    Y score for the parent peptide

    '''
    
    wt_p = _compute_wt_freq(P_pos, wt)
    wt_n = _compute_wt_freq(P_neg, wt)
    return np.divide(wt_p, wt_n)

def compute_Y(P_pos, P_neg, wt):
    '''
    Calculate the Y score matrix, represented as a 2D numpy array 
    of shape  (number of amino acids in the alphabet, 
               number of randomized positions in the library)
    
    containing Y scores for individual mutants. Cells corresponding
    to wild type sequence are filled in "manually".

    Parameters
    ----------
    P_pos : P-array for positive ("reactive") peptides (2D-numpy array)
    P_neg : P-array for negative ("unreactive") peptides (2D-numpy array)
    
    wt : amino acid sequence of the parent peptide (only variable region)
         1D-numpy array

    Returns
    -------
    Y : Y scores

    '''
    
    #Compute mutant frequency matrices
    F_pos = _compute_F(P_pos)
    F_neg = _compute_F(P_neg)
    
    #Get the Y scores. At this stage the cells corresponding
    #to the parent peptides contain nonsense values
    Y = np.divide(F_pos, F_neg)
    
    #Correct the parent peptide cells to correct values
    wt_Y = _wt_score(P_pos, P_neg, wt)
    for i,x in enumerate(wt):
        Y[tables.aas.index(x), i] = wt_Y
        
    return Y
    
    

def _Agresti_Coull_interval(C, n, z=1.96):
    '''
    Agresti-Coull confidence interval calculator.

    Parameters
    ----------
    C : 2D-numpy array with
        the shape (number of amino acids in the alphabet, 
                   number of randomized positions in the library),
        which holds read counts for each library member.
        
    n : total number of reads in the sample
    z : quantile of a standard normal distribution

    Returns
    -------
    CI : One-sided confidence interval at z

    '''
    
    #Fill in a p-hat matrix
    p_hat = np.divide(C + 0.5*z, n + z**2)        

    #Return the one-sided confidence interval value    
    return z * np.sqrt(np.divide(p_hat * (1 - p_hat), n + z**2))  

def Y_sampling_error(P_pos, P_neg, z=1.96):
    '''
    Calculate 95% confidence intervals for Y scores computed
    from P_pos and P_neg datasets. 

    Parameters
    ----------
    P_pos : P-array for positive ("reactive") peptides (2D-numpy array)
    P_neg : P-array for negative ("unreactive") peptides (2D-numpy array)
    z: quantile of a standard normal distribution

    Returns
    -------
    I : A 2D-numpy array, same dimension as Y, with one-sided confidence
        intervals for each library mutant.

    '''
    

    #Compute mutant frequency matrices
    F_pos = _compute_F(P_pos)
    F_neg = _compute_F(P_neg)
    
    #Number of trials
    n_pos = P_pos.shape[0]
    n_neg = P_neg.shape[0]
    
    #Calculate sampling confidence intervals both for
    #positive and negative samples
    CI_pos = _Agresti_Coull_interval(F_pos*n_pos, n_pos, z)
    CI_neg = _Agresti_Coull_interval(F_neg*n_neg, n_neg, z)
    
    #Get the confidence interval for Y scores by propagating the error
    CI_pos_rel = np.divide(CI_pos, F_pos)
    CI_neg_rel = np.divide(CI_neg, F_neg)  
    I = np.sqrt(np.square(CI_pos_rel) + np.square(CI_neg_rel))
    I = np.divide(I, np.log(2))
    
    return I
