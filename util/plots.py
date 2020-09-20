# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 17:17:00 2019
@author: Vinogradov Alex
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import tables

def Y_score(Y, lib, basename):
    '''
    Plot the Y score matrix
    '''

    fig, ax = plt.subplots(1, 1, figsize=(10, 10.5), dpi=300)

    scale_min = 1.02 * np.nanmin(Y)
    scale_max = 1.02 * np.nanmax(Y)
    norm = mpl.colors.Normalize(vmin=scale_min,vmax=scale_max)
    
    c = ax.pcolor(Y, cmap=plt.cm.magma, norm=norm,
                  edgecolors='w', linewidths=4)
    cbar = fig.colorbar(c, ax=ax)

    cbar.ax.set_ylabel("Y score", rotation=-90, va="bottom", fontsize=22)
    cbar.ax.tick_params(labelsize=20)    

    #set ticks
    ax.set_xticks(np.arange(Y.shape[1])+0.5)
    ax.set_yticks(np.arange(Y.shape[0])+0.5)   
    ax.set_xticklabels(lib.wt)
    ax.set_yticklabels(tables.aas)

    #set labels
    ax.set_xlabel('Wild type amino acid', fontsize=25)
    ax.set_ylabel('Mutated to', fontsize=25)
    ax.tick_params(axis='both', which='major', labelsize=21)
    ax.set_title('Y-score map', fontsize=27)
    
    #save png and svg, and close the file
    svg = basename + '.svg'
    png = basename + '.png'
    fig.savefig(svg, bbox_inches = 'tight')
    fig.savefig(png, bbox_inches = 'tight')
    return
