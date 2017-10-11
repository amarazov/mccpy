import random

from matplotlib import animation
from matplotlib import pyplot as plt

fig = plt.figure()

ax = plt.axes(aspect='equal', autoscale_on=True)

circles = []


def init(arrow_scale=.2):
    plt.axis([0, 1, 0, 1])
    plt.setp(plt.gca(), xticks=[0, 1], yticks=[0, 1])
    for (x, y), c in zip(L, colors):
        circle = plt.Circle((x, y), radius=sigma, fc=c)
        circles.append(circle)
        ax.add_patch(circle)
    ax.set_title('t = 0')
    return circles


L = [[0.25, 0.25], [0.75, 0.25], [0.25, 0.75], [0.75, 0.75]]
sigma = 0.15
sigma_sq = sigma ** 2
delta = 0.1
colors = ['r', 'b', 'g', 'orange']
n_steps = 100


def animate(i):
    a = random.choice(L)
    b = [a[0] + random.uniform(-delta, delta), a[1] + random.uniform(-delta, delta)]
    min_dist = min((b[0] - c[0]) ** 2 + (b[1] - c[1]) ** 2 for c in L if c != a)
    box_cond = min(b[0], b[1]) < sigma or max(b[0], b[1]) > 1.0 - sigma
    if not (box_cond or min_dist < 4.0 * sigma ** 2):
        a[:] = b
    for circle, c, in zip(circles, L):
        circle.center = c
    ax.set_title('t = ' + str(i))


anim = animation.FuncAnimation(fig, animate,
                               init_func=init,
                               frames=n_steps,
                               interval=300,
                               repeat=False)

plt.show()
