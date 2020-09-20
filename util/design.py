# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 20:59:03 2020
@author: Vinogradov Alex
"""

import numpy as np

class PepDesign(object):
    
    def __init__(self, pep_sequence, wt):
        '''
        The PepDesign constructor holds the information about the library design,
        i.e. the parent peptide sequence, amino acid positions subject to 
        randomization, constant regions, etc. Simplifies the parser.
    
        At this stage, the library design constructor only works for libraries
        of fixed length (i.e. every peptide is the same length). 
        
        Parameters
        ----------
        pep_sequence : generic library peptide sequence as dtype=str.
                       Amino acids subject to randomization should be
                       indicated as numerals.
                       e.g. the string 'LPENGA1111111111111111YPYDVPDYAGELARP'
                       indicates that the sequence 'LPENGA' is fixed for all
                       library members, whereas the following 16 amino acids
                       are subject to randomization (mutagenesis, etc).
                       
        wt : amino acid sequence of the parent peptide. dtype=str
            
        Returns
        -------
        None.

        '''
        
        self.design = np.array(list(pep_sequence))
        self.wt_full = np.array(list(wt))
        self._setup_params()
        
    def _setup_params(self):
        
        #L is expected length of the library peptide
        self.L = self.design.size
        
        #vr_mask is masking array for variable region positions
        #cr_mask is masking array for constant region positions
        self.vr_mask = [x.isdigit() for x in self.design]
        self.cr_mask = [not x.isdigit() for x in self.design]

        #vr: variable region
        #cr: constant region      
        self.vr = self.design[self.vr_mask]
        self.cr = self.design[self.cr_mask]
        
        #self.wt is the variable region of the parent peptide
        self.wt = self.wt_full[self.vr_mask]
        
        return    