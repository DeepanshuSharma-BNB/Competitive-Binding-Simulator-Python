#!/usr/bin/env python

# competitive_binding_sim.py

# simulate binding curves for competitive binding events (fluorescence polarization) 

# reference: 
# Michael H. A. Roehrl, Julia Y. Wang, and Gerhard Wagner, Biochemistry 43, 16056 (2004) 

# SNF Winter 2014/2015

import sys, os, argparse
import numpy as np
import matplotlib.pyplot as plt
from itertools import repeat

parser = argparse.ArgumentParser(description="Plot simulated competitive binding curves for the input values")

parser.add_argument("-o", "--output", help="Output filename for plot", default="binding_sim.png") 
parser.add_argument("-kd", "--kd", help="Binding constant for the labeled ligand (micromolar)", default=1, type=float)
parser.add_argument("-ki", "--ki", help="Binding constant for the competitor (micromolar)", default=0.1, type=float) 
parser.add_argument("-rt", "--receptor_total", help="Concentration of receptor (normally protein; micromolar)", default=1, type=float)
parser.add_argument("-lt", "--ligand_total", help="Concentration of labeled ligand; micromolar", default=0.01, type=float)

args = parser.parse_args()

print "------------------------------------"
print "|    competitive_binding_sim.py    |"
print "|  plot simulated binding curves   |"
print "|      run with -h for help        |"
print "|           snf 1/2015             |"
print "------------------------------------\n\n"

print "Arguments:"
print args

# test for reasonability of assumptions?? e.g. r_t >> l_st 

xvals = np.logspace(-9,-3,100)

def um_to_molar(val): 
    return val/1000000.

def molar_to_um(val): 
    return val*1000000.

k_d = um_to_molar(args.kd)
#k_i = nm_to_molar(args.ki )
k_is = [um_to_molar(args.ki / i) for i in (0.01, 0.1, 1, 10, 100)]
r_t = um_to_molar(args.receptor_total)
l_st = um_to_molar(args.ligand_total )

# equation 17 in the above reference - for "complete competition"; i.e. ligand and competitor are mutex  
def eqn17(l_t, k_i):
    d = k_d + k_i + l_st + l_t - r_t
    e = (l_t - r_t) * k_d + (l_st - r_t) * k_i + k_d * k_i
    f = -1 * k_d * k_i * r_t
    theta = np.arccos( (-2 * d**3 + 9 * d * e - 27 * f) / (2 * np.sqrt( (d**2 - 3*e)**3 ) ) )
    fsb = ( 2 * np.sqrt( d**2 - 3*e ) * np.cos(theta/3.) - d ) / ( 3 * k_d + 2 * np.sqrt( d**2 - 3*e ) * np.cos(theta/3.) - d ) 
    return fsb


plots = [ map(eqn17, xvals, repeat(k_i, len(xvals))) for k_i in k_is ]

fig, ax = plt.subplots()
for result in plots: 
    plt.plot(xvals, result)
plt.title("Kd (uM): %3.2f Ki (uM): %3.2f Recep (uM): %3.2f Labeled_lig (uM): %3.2f" % (args.kd, args.ki, args.receptor_total, args.ligand_total), fontsize=10)
plt.axis([min(xvals), max(xvals), 0, 1]) 
ax.set_xscale("log")
ax.set_xticks([um_to_molar(.001), um_to_molar(.01), um_to_molar(.1), um_to_molar(1), um_to_molar(10), um_to_molar(100), um_to_molar(1000)])
ax.set_xticklabels(["1nM", "10nM", "100nM", "1uM", "10uM", "100uM", "1mM"], fontsize=8)
ax.set_yticklabels(["0", "0.2", "0.4", "0.6", "0.8", "1"], fontsize=8)
ax.set_ylabel("Fraction bound") 
ax.set_xlabel("[Competitor]")
plt.legend(["Ki (uM): %3.2f" % molar_to_um(k_i) for k_i in k_is], fontsize=8)
plt.savefig(args.output)
plt.close(fig) 
