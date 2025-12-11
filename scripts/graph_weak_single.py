#!/usr/bin/env python3

# Usage: 
#   module load scipy-stack/2023b
#   PYTHONNOUSERSITE=1 python3 graph_weak_single.py <data_file>
# Or use teachsetup which includes scipy-stack:
#   source scripts/teachsetup
#   PYTHONNOUSERSITE=1 python3 graph_weak_single.py <data_file>
import matplotlib.pyplot as plt
import sys

if len(sys.argv) != 2:
    print("Please give path to results.txt")
    sys.exit(1)

file=sys.argv[1]

def parse_row(ls):
    l=ls.split('|')
    c=int(l[1])
    # Check if second column is eps or width
    try:
        e=float(l[2])
        t=float(l[3])
    except (ValueError, IndexError):
        # If parsing fails, try alternative format
        t=float(l[2])
        e=None
    return (c,t,e)

core=[]
time=[]
workload=[]

i=0
with open(file) as f:
    for line in f:
        if line[0] != '|':
            continue
        if i==0:
            i+=1
            continue
    
        (c,t,e) = parse_row(line)
    
        if i==1:
            i+=1
    
        core.append(c)
        time.append(t)
        if e is not None:
            workload.append(e)

bloo='#1f77b4'
gren='#2ca02c'

# Determine workload label based on header or data
workload_label = "# of Pixels"
if 'width' in file.lower() or (len(workload) > 0 and workload[0] > 1000):
    workload_label = "Image Width"

fig, ax1 = plt.subplots(figsize=(6,5))
ax2 = ax1.twinx()
ax1.plot(core,time,'o-',label='time',color=bloo)
ax1.set_xlabel('Cores', fontsize=14)
ax1.set_ylabel('Time (seconds)', fontsize=14,color=bloo)
if len(workload) > 0:
    ax2.plot(core,workload,'-',label=workload_label,color=gren)
    ax2.set_ylabel(workload_label, fontsize=12,color=gren)
ax1.legend()

plt.savefig(f'{file}.png', bbox_inches='tight', dpi=400)
