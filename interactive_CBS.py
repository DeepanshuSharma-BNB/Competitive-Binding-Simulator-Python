#!/usr/bin/env python3

# competitive_binding_sim.py

# Simulate binding curves for competitive binding events (fluorescence polarization)

# Reference:
# Michael H. A. Roehrl, Julia Y. Wang, and Gerhard Wagner, Biochemistry 43, 16056 (2004)

# Updated for Python 3 compatibility and interactive visualization

import argparse
import numpy as np
import matplotlib.pyplot as plt
from itertools import repeat
from matplotlib.widgets import Slider

parser = argparse.ArgumentParser(description="Plot simulated competitive binding curves for the input values")

parser.add_argument("-o", "--output", help="Output filename for plot", default="binding_sim.png") 
parser.add_argument("-kd", "--kd", help="Binding constant for the labeled ligand (micromolar)", default=1, type=float)
parser.add_argument("-ki", "--ki", help="Binding constant for the competitor (micromolar)", default=0.1, type=float) 
parser.add_argument("-rt", "--receptor_total", help="Concentration of receptor (normally protein; micromolar)", default=1, type=float)
parser.add_argument("-lt", "--ligand_total", help="Concentration of labeled ligand; micromolar", default=0.01, type=float)

args = parser.parse_args()

def um_to_molar(val): 
    return val / 1_000_000.

def molar_to_um(val): 
    return val * 1_000_000.

def eqn17(l_t, k_i, k_d, l_st, r_t):
    d = k_d + k_i + l_st + l_t - r_t
    e = (l_t - r_t) * k_d + (l_st - r_t) * k_i + k_d * k_i
    f = -1 * k_d * k_i * r_t
    theta = np.arccos((-2 * d**3 + 9 * d * e - 27 * f) / (2 * np.sqrt((d**2 - 3*e)**3)))
    fsb = (2 * np.sqrt(d**2 - 3*e) * np.cos(theta / 3.) - d) / (3 * k_d + 2 * np.sqrt(d**2 - 3*e) * np.cos(theta / 3.) - d)
    return fsb

def update(val):
    k_d = um_to_molar(kd_slider.val)
    k_is = [um_to_molar(ki_slider.val / i) for i in (0.01, 0.1, 1, 10, 100)]
    r_t = um_to_molar(rt_slider.val)
    l_st = um_to_molar(lt_slider.val)
    plots = [[eqn17(x, k_i, k_d, l_st, r_t) for x in xvals] for k_i in k_is]
    ax.clear()
    for result in plots:
        ax.plot(xvals, result)
    ax.set_xscale("log")
    ax.set_xticks([um_to_molar(i) for i in [0.001, 0.01, 0.1, 1, 10, 100, 1000]])
    ax.set_xticklabels(["1nM", "10nM", "100nM", "1uM", "10uM", "100uM", "1mM"], fontsize=8)
    ax.set_ylabel("Fraction bound") 
    ax.set_xlabel("[Competitor]")
    ax.legend([f"Ki (uM): {molar_to_um(k_i):.2f}" for k_i in k_is], fontsize=8)
    fig.canvas.draw_idle()

xvals = np.logspace(-9, -3, 100)
k_d = um_to_molar(args.kd)
k_is = [um_to_molar(args.ki / i) for i in (0.01, 0.1, 1, 10, 100)]
r_t = um_to_molar(args.receptor_total)
l_st = um_to_molar(args.ligand_total)

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.1, bottom=0.3)
for k_i in k_is:
    ax.plot(xvals, [eqn17(x, k_i, k_d, l_st, r_t) for x in xvals])
ax.set_xscale("log")
ax.set_ylabel("Fraction bound")
ax.set_xlabel("[Competitor]")

axcolor = 'lightgoldenrodyellow'
ax_kd = plt.axes([0.1, 0.15, 0.65, 0.03], facecolor=axcolor)
ax_ki = plt.axes([0.1, 0.1, 0.65, 0.03], facecolor=axcolor)
ax_rt = plt.axes([0.1, 0.05, 0.65, 0.03], facecolor=axcolor)
ax_lt = plt.axes([0.1, 0.0, 0.65, 0.03], facecolor=axcolor)

kd_slider = Slider(ax_kd, 'Kd (uM)', 0.01, 10.0, valinit=args.kd)
ki_slider = Slider(ax_ki, 'Ki (uM)', 0.01, 10.0, valinit=args.ki)
rt_slider = Slider(ax_rt, 'Receptor (uM)', 0.1, 10.0, valinit=args.receptor_total)
lt_slider = Slider(ax_lt, 'Ligand (uM)', 0.001, 1.0, valinit=args.ligand_total)

kd_slider.on_changed(update)
ki_slider.on_changed(update)
rt_slider.on_changed(update)
lt_slider.on_changed(update)

plt.show()