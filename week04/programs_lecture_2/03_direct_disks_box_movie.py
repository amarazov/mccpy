import random

import numpy as np
from matplotlib import animation
from matplotlib import pyplot as plt


def direct_disks_box(N, sigma):
    condition = False
    while condition == False:
        L = [(random.uniform(sigma, 1.0 - sigma), random.uniform(sigma, 1.0 - sigma))]
        for k in range(1, N):
            a = (random.uniform(sigma, 1.0 - sigma), random.uniform(sigma, 1.0 - sigma))
            min_dist = min(np.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) for b in L)
            if min_dist < 2.0 * sigma:
                condition = False
                break
            else:
                L.append(a)
                condition = True
    return L


fig = plt.figure()

ax = plt.axes(aspect='equal', autoscale_on=True)

circles = []


def init():
    plt.axis([0, 1, 0, 1])
    plt.setp(plt.gca(), xticks=[0, 1], yticks=[0, 1])
    pos = direct_disks_box(N, sigma)
    for (x, y), c in zip(pos, colors):
        circle = plt.Circle((x, y), radius=sigma, fc=c)
        circles.append(circle)
        ax.add_patch(circle)
    ax.set_title('t = 0')
    return circles


def animate(i):
    pos = direct_disks_box(N, sigma)
    for circle, c in zip(circles, pos):
        circle.center = c
    ax.set_title('t = ' + str(i))
    return circles


N = 4
colors = ['r', 'b', 'g', 'orange']
sigma = 0.2
n_runs = 10

anim = animation.FuncAnimation(fig, animate,
                               init_func=init,
                               frames=n_runs,
                               interval=1000,
                               repeat=False)

plt.show()
