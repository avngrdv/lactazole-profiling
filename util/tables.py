# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 01:37:10 2020
@author: Vinogradov Alex
"""

translation_table = {
                    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
                    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
                    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
                    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
                    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
                    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
                    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
                    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
                    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
                    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
                    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
                    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
                    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
                    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
                    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
                    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
                    }
         
aas = ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
      'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')
      
bases = ('T', 'C', 'A', 'G')

codons = ('TTT', 'TTG', 'TCT', 'TCG', 'TAT', 'TAG', 'TGT', 'TGG', 'CTT',
          'CTG', 'CCT', 'CCG', 'CAT', 'CAG', 'CGT', 'CGG', 'ATT', 'ATG',
          'ACT', 'ACG', 'AAT', 'AAG', 'AGT', 'AGG', 'GTT', 'GTG', 'GCT',
          'GCG', 'GAT', 'GAG', 'GGT', 'GGG', 'TTC', 'TTA', 'TCC', 'TCA',
          'TAC', 'TAA', 'TGC', 'TGA', 'CTC', 'CTA', 'CCC', 'CCA', 'CAC',
          'CAA', 'CGC', 'CGA', 'ATC', 'ATA', 'ACC', 'ACA', 'AAC', 'AAA',
          'AGC', 'AGA', 'GTC', 'GTA', 'GCC', 'GCA', 'GAC', 'GAA', 'GGC',
          'GGA')