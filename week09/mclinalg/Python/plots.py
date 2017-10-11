#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

s = 20
seq_length=2**16
max_scrambled_digits=30
iflag=0

rand = Faur(s,seq_length,max_scrambled_digits,iflag)


fig = plt.figure()
ax1 = fig.add_subplot(121)

## the data
N=1000
x = np.random.randn(N)
y = np.random.randn(N)

## left panel
ax1.scatter(x,y,color='blue',s=5,edgecolor='none')
ax1.set_aspect(1./ax1.get_data_ratio()) # make axes square


## right panel
ax2 = fig.add_subplot(122)
props = dict(alpha=0.5, edgecolors='none' )

handles = []
colors = ['blue', 'green', 'magenta', 'cyan']
for color in colors:
    x = np.random.randn(N)
    y = np.random.randn(N)
    s = np.random.randint(50,200)
    handles.append(ax2.scatter(x, y, c=color, s=s, **props))

ax2.set_ylim([-5,11])
ax2.set_xlim([-5,11])

ax2.legend(handles, colors)
ax2.grid(True)
ax2.set_aspect(1./ax2.get_data_ratio())
plt.show()