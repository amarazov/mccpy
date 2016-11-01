import numpy as np
from matplotlib import animation
from matplotlib import pyplot as plt

plt.style.use('ggplot')
np.random.seed(1234)

sigma = 0.4
epsilon = 0.0  # probability to switch from red to blue pebble, and vice versa

fig = plt.figure()

ax = plt.axes(aspect='equal', autoscale_on=True)
ax.set_xlim(0.5, 6.5), ax.set_xticks([])
ax.set_ylim(-2.5, 3.5), ax.set_yticks([])

s_map_red = [(1.0, 1.0), (2.0, 1.0), (3.0, 1.0),
             (1.0, 2.0), (2.0, 2.0), (3.0, 2.0),
             (1.0, 3.0), (2.0, 3.0), (3.0, 3.0)]
offset = 3.0
s_map_blue = [(x + offset, y - offset) for (x, y) in s_map_red]
neighbor = [[1, 3, 0, 0], [2, 4, 0, 1], [2, 5, 1, 2],
            [4, 6, 3, 0], [5, 7, 3, 1], [5, 8, 4, 2],
            [7, 6, 6, 3], [8, 7, 6, 4], [8, 8, 7, 5]]
color = ['blue', ]  # chose 'red' or 'blue'
site = [8, ]
tmax = 240
patch = plt.Circle(s_map_red[site[0]], radius=sigma, fc='r')

plt.plot([0.5, 3.5], [0.5, 0.5], 'r')
plt.plot([0.5, 3.5], [1.5, 1.5], 'r')
plt.plot([0.5, 3.5], [2.5, 2.5], 'r')
plt.plot([1.5, 1.5], [0.5, 3.5], 'r')
plt.plot([2.5, 2.5], [0.5, 3.5], 'r')
plt.plot([3.5, 3.5], [0.5, 3.5], 'r')
plt.plot([0.5 + offset, 3.5 + offset], [1.5 - offset, 1.5 - offset], 'b')
plt.plot([0.5 + offset, 3.5 + offset], [2.5 - offset, 2.5 - offset], 'b')
plt.plot([0.5 + offset, 3.5 + offset], [3.5 - offset, 3.5 - offset], 'b')
plt.plot([0.5 + offset, 0.5 + offset], [0.5 - offset, 3.5 - offset], 'b')
plt.plot([1.5 + offset, 1.5 + offset], [0.5 - offset, 3.5 - offset], 'b')
plt.plot([2.5 + offset, 2.5 + offset], [0.5 - offset, 3.5 - offset], 'b')


def init():
    patch.center = s_map_red[site[0]]
    ax.add_patch(patch)
    ax.set_title('t = 0')
    return patch,


def animate(i):
    maxlength = len(str(tmax - 1))
    number_string = str(i).zfill(maxlength)
    if color[0] == 'red':
        patch.center = s_map_red[site[0]]
    if color[0] == 'blue':
        patch.center = s_map_blue[site[0]]
    ax.set_title('t = ' + number_string)
    # End of graphical output
    newsite = neighbor[site[0]][np.random.randint(0, 4)]
    newcolor = color[0]
    if i % 2 == 1:
        patch.set_color('g')
    else:
        patch.set_color('r')
    if (color[0] == 'red') and (site[0] == 2) and (newsite == 2):
        if np.random.random() < epsilon:
            newcolor = 'blue'
            newsite = 6
            print("transition red->blue at time = ", i)
    if (color[0] == 'blue') and (site[0] == 6) and (newsite == 6):
        if np.random.random() < epsilon:
            newcolor = 'red'
            newsite = 2
            print("transition blue->red at time = ", i)
    site[0] = newsite
    color[0] = newcolor
    return patch


anim = animation.FuncAnimation(fig, animate,
                               init_func=init,
                               frames=tmax,
                               interval=1000,
                               repeat=False)

plt.show()
