import numpy as np
from matplotlib import animation
from matplotlib import pyplot as plt

plt.style.use('ggplot')

sigma = 0.4  # sigma and s_map are needed for the graphical output
s_map = [(1.0, 1.0), (2.0, 1.0), (3.0, 1.0),
         (1.0, 2.0), (2.0, 2.0), (3.0, 2.0),
         (1.0, 3.0), (2.0, 3.0), (3.0, 3.0)]
neighbor = [[1, 3, 0, 0], [2, 4, 0, 1], [2, 5, 1, 2],
            [4, 6, 3, 0], [5, 7, 3, 1], [5, 8, 4, 2],
            [7, 6, 6, 3], [8, 7, 6, 4], [8, 8, 7, 5]]
site = [8, ]
N_runs = 100

fig = plt.figure()

ax = plt.axes(aspect='equal', autoscale_on=True)
ax.set_xlim(0.0, 4.0), ax.set_xticks([])
ax.set_ylim(0.0, 4.0), ax.set_yticks([])
patch = plt.Circle(s_map[site[0]], radius=sigma, fc='r')
plt.plot([0.5, 3.5], [1.5, 1.5], 'b')
plt.plot([0.5, 3.5], [2.5, 2.5], 'b')
plt.plot([1.5, 1.5], [0.5, 3.5], 'b')
plt.plot([2.5, 2.5], [0.5, 3.5], 'b')


def init():
    patch.center = s_map[site[0]]
    ax.add_patch(patch)
    ax.set_title('t = 0')
    return patch,


def animate(i):
    site[0] = neighbor[site[0]][np.random.randint(0, 4)]
    patch.center = s_map[site[0]]
    if i % 2 == 1:
        patch.set_color('y')
    else:
        patch.set_color('r')
    ax.set_title('t = ' + str(i))
    print(i, patch.center)
    return patch


anim = animation.FuncAnimation(fig, animate,
                               init_func=init,
                               frames=N_runs,
                               interval=1000,
                               repeat=False)

plt.show()
