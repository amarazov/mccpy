import numpy as np
from matplotlib import animation
from matplotlib import pyplot as plt

plt.style.use('ggplot')

xvec = {1: 3, 2: 2, 3: 1, 4: 3, 5: 2, 6: 1, 7: 3, 8: 2, 9: 1}
yvec = {1: 1, 2: 1, 3: 1, 4: 2, 5: 2, 6: 2, 7: 3, 8: 3, 9: 3}

neighbor = {1: [2, 4, 1, 1], 2: [3, 5, 1, 2], 3: [3, 6, 2, 3],
            4: [5, 7, 4, 1], 5: [6, 8, 4, 2], 6: [6, 9, 5, 3],
            7: [8, 7, 7, 4], 8: [9, 8, 7, 5], 9: [9, 9, 8, 6]}

N_runs = 20

fig = plt.figure()

ax = plt.axes(aspect='equal', autoscale_on=True)
ax.set_xticks([])
ax.set_yticks([])
H, xedges, yedges = np.histogram2d([], [], bins=(3, 3),
                                   range=[[1, 3], [1, 3]], normed=True)
extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
histo = plt.imshow(H, extent=extent, interpolation='nearest', vmin=0, vmax=1.00)
histo.set_cmap('hot')
plt.colorbar()


def animate(i):
    list_vec = []
    if i < 10:
        run_str = '0' + str(i)
    else:
        run_str = str(i)
    for n_runs in range(100000):
        pos = 9
        for iter in range(i):
            pos = neighbor[pos][np.random.randint(0, 4)]
        list_vec.append(pos)

    x = [xvec[k] for k in list_vec]
    y = [yvec[k] for k in list_vec]

    plt.xticks([])
    plt.yticks([])
    H, xedges, yedges = np.histogram2d(x, y, bins=(3, 3),
                                       range=[[1, 3], [1, 3]], normed=True)
    extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
    histo = plt.imshow(H, extent=extent, interpolation='nearest', vmin=0, vmax=1.00)
    histo.set_cmap('hot')
    ax.set_title('t = ' + run_str)


anim = animation.FuncAnimation(fig, animate,
                               frames=N_runs,
                               interval=2000,
                               repeat=False,
                               )
plt.show()
