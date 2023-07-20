import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description='Visualize vartigs in a vartig-coverage plot. ')

parser.add_argument("vartig_file", help="the vartig file to visualize.", type=str)
parser.add_argument("-q", "--min-hapq", help = "minimum HAPQ threshold (default = 1)", type = int, default = 1) 
parser.add_argument("-c", "--min-cov", help = "minimum COV threshold (default = 1.5)", type = int, default = 1.5) 
parser.add_argument("-l", "--min-len", help = "minimum vartig length threshold (default = 1000)", type = int, default = 1000) 
args = parser.parse_args()

hapq_cut = int(args.min_hapq)
len_cutoff = int(args.min_len)
cov_cutoff =  float(args.min_cov)

def normal (x):
    return np.log2(x+np.array(1))


import numpy as np
from matplotlib import collections  as mc

import os
import re
from sys import argv
import subprocess
import shlex
p = re.compile('COV:(\d*\.?\d+)')
snp_p = re.compile('BASERANGE:(\d+)-(\d+)')
hapq_p = re.compile('HAPQ:(\d+)')

threshold = [0.00,1.00]
mult = 1.0
plt.rcParams.update({'font.size': 7})
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams.update({'font.family':'arial'})
cm = 1/2.54  # centimeters in inches

#cmap = plt.cm.cool
#cmap = plt.cm.viridis
cmap = plt.cm.jet
#cmap = plt.cm.plasma
cmap_hq = plt.cm.binary
#cmap_hq = plt.cm.plasma
hap = args.vartig_file
hap_cov = []
lines = []
line_colors = []
line_colors_hapq = []
widths = []
max_br = 0
cont = False
hapq_fail = False
for line in open(hap,'r'):
    if line[0] == '>':
        name = line.split()[0][1:]
        ar = p.findall(line)
        br = snp_p.findall(line)
        hapr = hapq_p.findall(line)
        if len(ar) == 0:
            cont = True
            #print(line)
            continue
        curr_cov = float(ar[0])
        start = int(br[0][0])
        end = int(br[0][1])
        hapq = int(hapr[0])
        if hapq > 60:
            hapq = 60
        #print(hapq)
        baserange = (start,end)
        if end - start < len_cutoff:
            cont = True
            continue
        if hapq < hapq_cut:
            hapq_fail = True
        if end > max_br:
            max_br = end
    else:
        if cont:
            cont = False
            continue
        total_alt = 0
        al_c = 0
        num_zero = 0
        for x in line:
            if x == '0' or x == '1' or x == '2':
                al_c += 1
            if x == '1' or x == '2':
                total_alt += 1 
            else:
                num_zero += 1
        if curr_cov >= cov_cutoff:
            if al_c == 0:
                #print(line, curr_cov)
                continue
            hap_cov.append(curr_cov)
            lines.append([(baserange[0], (curr_cov+1)), (baserange[1], (curr_cov+1))])
            line_colors.append(total_alt / al_c)
            line_colors_hapq.append(hapq)
            if hapq_fail:
                hapq_fail = False
                widths.append(0)
            else:
                widths.append(4)
        else:
            hapq_fail = False

lc = mc.LineCollection(lines, linewidths=widths, array = np.array(line_colors), cmap = cmap, alpha = 1.0)
lc_hq = mc.LineCollection(lines, linewidths=widths, array = np.array(line_colors_hapq), cmap = cmap_hq, alpha = 1.0)
lc_hq.set_clim(vmin=0, vmax=60)
lc.set_clim(vmin=0, vmax=1)
#print(lc)

fig, ax = plt.subplots(2)
ax[0].add_collection(lc)
ax[0].set_xlim(0, max_br)
ax[0].set_ylim(np.min(hap_cov)-5, np.max(hap_cov) + 5)
fig.colorbar(lc, ax=ax[0], orientation='vertical')
fig.colorbar(lc_hq, ax=ax[1], orientation='vertical')



ax[1].add_collection(lc_hq)
ax[1].set_xlim(0, max_br)
ax[1].set_ylim(np.min(hap_cov)-5, np.max(hap_cov) + 5)
ax[0].set_title("Vartigs colored by alternate allele ratio")
ax[1].set_title("Vartigs colored by HAPQ")
ax[1].set_xlabel("Base position")
ax[0].set_ylabel("Coverage")
ax[1].set_ylabel("Coverage")
plt.show()

