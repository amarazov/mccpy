import numpy as np
from matplotlib import animation
from matplotlib import pyplot as plt

plt.style.use('ggplot')

sigma = 0.4
epsilon = 0.1
s_map = [(1.0, 1.0), (2.0, 1.0)]
neighbor = [[1], [0]]
pos = 0
tmax = 20

fig = plt.figure()

ax = plt.axes(aspect='equal', autoscale_on=True)
ax.set_xlim(0.5, 2.5), ax.set_xticks([])
ax.set_ylim(0.5, 1.5), ax.set_yticks([])
patch = plt.Circle(s_map[pos], radius=sigma, fc='r')
plt.plot([1.5, 1.5], [0.5, 1.5], 'b')


def init():
    patch.center = s_map[pos]
    ax.add_patch(patch)
    ax.set_title('t = 0')
    return patch,


def animate(i):
    global pos
    number_string = str(i).zfill(len(str(tmax)))
    newpos = neighbor[pos][0]
    if i % 2 == 1:
        patch.set_color('y')
    else:
        patch.set_color('r')
    if np.random.random() < epsilon:
        newpos = pos
    pos = newpos
    patch.center = s_map[pos]

    ax.set_title('t = ' + number_string)
    return patch


anim = animation.FuncAnimation(fig, animate,
                               init_func=init,
                               frames=tmax,
                               interval=1000,
                               repeat=False)
plt.show()
