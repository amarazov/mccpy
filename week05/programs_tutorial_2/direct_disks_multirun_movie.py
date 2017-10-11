import random

import numpy as np
from matplotlib import animation
from matplotlib import pyplot as plt

fig = plt.figure()
border_color = 'k'
plt.setp(plt.gca(), xticks=[0, 1], yticks=[0, 1], aspect='equal')
ax = plt.axes(aspect='equal', autoscale_on=True)


def dist(x, y):
    d_x = abs(x[0] - y[0]) % 1.0
    d_x = min(d_x, 1.0 - d_x)
    d_y = abs(x[1] - y[1]) % 1.0
    d_y = min(d_y, 1.0 - d_y)
    return np.sqrt(d_x ** 2 + d_y ** 2)


def direct_disks(N, sigma):
    n_iter = 0
    condition = False
    while condition == False:
        n_iter += 1
        L = [(random.random(), random.random())]
        for k in range(1, N):
            a = (random.random(), random.random())
            min_dist = min(dist(a, b) for b in L)
            if min_dist < 2.0 * sigma:
                condition = False
                break
            else:
                L.append(a)
                condition = True
    return n_iter, L


circles = []


def init(arrow_scale=.2):
    plt.axis([0, 1, 0, 1])
    plt.setp(plt.gca(), xticks=[0, 1], yticks=[0, 1])
    iterations, config = direct_disks(N, sigma)
    print('run', -1)
    print(iterations - 1, 'tabula rasa wipe-outs before producing the following configuration')
    print(config)
    config_per = periodicize(config)
    for (x, y), c in zip(config_per, colors):
        circle = plt.Circle((x, y), radius=sigma, fc=c)
        circles.append(circle)
        ax.add_patch(circle)
    ax.set_title('t = -1')
    return circles


def periodicize(config):
    images = [-1.0, 0.0, 1.0]
    return [(x + dx, y + dy) for (x, y) in config for dx in images for dy in images]


N = 16
eta = 0.28
sigma = np.sqrt(eta / N / np.pi)
n_runs = 8
colors = ['r' for i in range(8 * N)]


def animate(run):
    iterations, config = direct_disks(N, sigma)
    print('run', run)
    print(iterations - 1, 'tabula rasa wipe-outs before producing the following configuration')
    print(config)
    config_per = periodicize(config)
    for circle, c, in zip(circles, config_per):
        circle.center = c
    ax.set_title('t = ' + str(round(run, 4)))
    return circles


anim = animation.FuncAnimation(fig, animate,
                               init_func=init,
                               frames=n_runs,
                               interval=3000,
                               repeat=False)

plt.show()
